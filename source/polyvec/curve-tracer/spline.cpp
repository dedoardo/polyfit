// Header
#include <polyvec/curve-tracer/spline.hpp>

// Polyvec
#include <polyvec/utils/system.hpp>
#include <polyvec/core/options.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/utils/potrace.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/curve-tracer/curve_line.hpp>
#include <polyvec/curve-tracer/find-fitting-points.hpp>
#include <polyvec/utils/directions.hpp>
#include <polyvec/curve-tracer/measure_accuracy_polygon.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/curve-tracer/regularity_handling.hpp>
#include <polyvec/curve-tracer/fit_options.hpp>

// c++ stl
#include <cstdlib> // abs
#include <algorithm> // find_first_of
#include <sstream>
#include <iostream>

using namespace Eigen;
using namespace std;
using namespace polyfit;

NAMESPACE_BEGIN(polyvec) 

std::string fitting_suffix() {
#if !USE_CURVATURE_VARIATION_MINIMIZING_INITIALIZATION && !USE_CURVATURE_VARIATION_MINIMIZATION && !REDUCED_DOF
	return "-old";
#else
	std::stringstream ss;
#if USE_CURVATURE_VARIATION_MINIMIZING_INITIALIZATION
	ss << "-curveVarInit";
#endif
#if REDUCED_DOF
	ss << "-reduceddof";
#endif
#if USE_CURVATURE_VARIATION_MINIMIZATION
	ss << "-curvvar";
#endif

#if REGULARIZE_BEZIER
	ss << "-bezierreg";
#endif
	return ss.str();
#endif
}

const char* tangent_fit_type_to_string(const TangentFitType v) {
	switch (v) {
	case TANGENT_FIT_CONSTANT:
		return "TANGENT_FIT_CONSTANT";
	case TANGENT_FIT_LERP:
		return "TANGENT_FIT_LERP";
	case TANGENT_FIT_LERP_SYM:
		return "TANGENT_FIT_LERP_SYM";
	default:
		return "Unknown";
	}
}

// -------------------------------------------------------------------------
Index poly_linidx(const Eigen::Matrix2Xd& poly, const Index idx) {
	return (idx + poly.cols() * 2) % poly.cols();
}

// -------------------------------------------------------------------------

bool polyvec::CurveSequenceFitter::edge_has_important_tangent(int i_edge) const
{
	const Eigen::Vector2d d = CircularAt(polygon, i_edge + 1) - polygon.col(i_edge);

	// axis-aligned
	if (std::abs(d.x()) < PF_EPS || std::abs(d.y()) < PF_EPS)
		return true;

	// 45
	if (std::abs(std::abs(d.x()) - std::abs(d.y())) < PF_EPS)
		return true;

	return false;
}

double cross2(const Eigen::Vector2d& v1, const Eigen::Vector2d& v2) {
	return v1.x() * v2.y() - v1.y() * v2.x();
}

//Returns a Bezier curve that is minimizes curvature variation and interpolates the two endpoints and endtangent directions.
//t1 and t2 must be unit vectors.
BezierCurve* init_curvature_variation_minimization_bezier(const Eigen::Vector2d& p1, const Eigen::Vector2d& t1, const Eigen::Vector2d& p2,
	const Eigen::Vector2d& t2) {
	//Jaklic et al., Curvature variation minimizing cubic Hermite interpolants
	Matrix2Xd control_points(2, 4);

	control_points.col(0) = p1;
	control_points.col(3) = p2;

	Eigen::Vector2d dp = p2 - p1;
	double c01 = cross2(t1, t2);
	double c0 = cross2(t1, dp);
	double c1 = cross2(t2, dp);

	double d1 = -2 * c1 / c01 / 3;
	double d2 = 2 * c0 / c01 / 3;

	if (std::abs(c01) <= 0.01) { //tangents are parallel
		d1 = d2 = dp.norm() / 3;
	}

	assert_break(d1 >= 0 && d2 >= 0); //this should hold for polygon corners

	control_points.col(1) = p1 + d1 * t1;
	control_points.col(2) = p2 - d2 * t2;

	auto curve = new BezierCurve();
	curve->set_control_points(control_points);
	return curve;
}

template<TangentFitType FIT>
void polyvec::CurveSequenceFitter::init_all_fits() {
	init_fit<TangentFitType(FIT - 1) >();
	//recursively init all previous types
	init_all_fits<TangentFitType(FIT - 1) >();
}

template<>
void polyvec::CurveSequenceFitter::init_all_fits<TangentFitType(0) >()
{ }

template<TangentFitType FIT>
void CurveSequenceFitter::init_fit() {
	initial_fits[FIT].resize(polygon.cols());

	for (int i = 0; i < polygon.cols(); ++i) {
		init_corner_fit<FIT>(i);
	}
}

template<>
void polyvec::CurveSequenceFitter::init_corner_fit<TANGENT_FIT_LERP>(int cornerId) {
	auto corner_prev_id = poly_linidx(polygon, cornerId - 1);
	auto corner_next_id = poly_linidx(polygon, cornerId + 1);

	const Vector2d& corner_prev = polygon.col(corner_prev_id);
	const Vector2d& corner = polygon.col(cornerId);
	const Vector2d& corner_next = polygon.col(corner_next_id);

	const Vector2d mid_prev = misc::lerp(corner, corner_prev, .5);
	const Vector2d mid_next = misc::lerp(corner, corner_next, .5);

#if USE_CURVATURE_VARIATION_MINIMIZING_INITIALIZATION
	auto curve = init_curvature_variation_minimization_bezier(mid_prev, (corner - corner_prev).normalized(), mid_next, (corner_next - corner).normalized());
#else
	const double len_prev = (corner - corner_prev).norm();
	const double len_next = (corner_next - corner).norm();

	Matrix2Xd control_points(2, 4);
	control_points.col(0) = mid_prev;
	control_points.col(3) = mid_next;

	if (len_prev < len_next) {
		control_points.col(1) = misc::lerp(corner, mid_prev, len_prev / (len_prev + len_next));
		control_points.col(2) = misc::lerp(mid_next, corner, len_next / (len_prev + len_next));
	}
	else {
		control_points.col(1) = misc::lerp(mid_prev, corner, len_prev / (len_prev + len_next));
		control_points.col(2) = misc::lerp(corner, mid_next, len_next / (len_prev + len_next));
	}

	auto curve = new BezierCurve();
	curve->set_control_points(control_points);
#endif

	initial_fits[TANGENT_FIT_LERP][cornerId].emplace_back();

    bool invert_front_coordinate_system = invert_edge_coordinate_system[Circular(invert_edge_coordinate_system, (Index)cornerId - 1)];
    bool invert_back_coordinate_system = invert_edge_coordinate_system[cornerId];
	auto param = new GlobFitBezierAngleBasedParametrization(std::shared_ptr<BezierCurve>(curve), invert_front_coordinate_system, invert_back_coordinate_system);
	initial_fits[TANGENT_FIT_LERP][cornerId].back().curve = param;

	auto cornerT = curve->project(corner);
	initial_fits[TANGENT_FIT_LERP][cornerId].back().fitting_info.add_tangent_samples(param, cornerT);
	initial_fits[TANGENT_FIT_LERP][cornerId].back().fitting_info.add_midpoint_samples(points, polygonV, corner_prev_id, 0.5, cornerId, 0.5, circular);

    // Computing the edge ids
    initial_fits[TANGENT_FIT_LERP][cornerId].back().fitting_info.edge_src = BoundaryGraph::make_edge_id(corner_prev_id, cornerId);
    initial_fits[TANGENT_FIT_LERP][cornerId].back().fitting_info.edge_dst = BoundaryGraph::make_edge_id(cornerId, corner_next_id);
}

template<>
void polyvec::CurveSequenceFitter::init_corner_fit<TANGENT_FIT_LERP_SYM>(int cornerId) {
	auto corner_prev_id = poly_linidx(polygon, (Index)cornerId - 1);
	auto corner_next_id = poly_linidx(polygon, (Index)cornerId + 1);

	const Vector2d& corner_prev = polygon.col(corner_prev_id);
	const Vector2d& corner = polygon.col(cornerId);
	const Vector2d& corner_next = polygon.col(corner_next_id);

	const double len_prev = (corner - corner_prev).norm();
	const double len_next = (corner_next - corner).norm();

	Vector2d mid_prev, mid_next;
	double prevT, nextT;

	if (len_prev < len_next) {
		prevT = .5;
		nextT = .5 * len_prev / len_next;
		if ((0.5 - nextT) * len_next < 0.1)
			nextT = 0.5; // don't create extremely short line segments
		mid_prev = misc::lerp(corner_prev, corner, prevT);
		mid_next = misc::lerp(corner, corner_next, nextT);
	}
	else {
		prevT = 1 - .5 * len_next / len_prev;
		nextT = .5;
		if ((prevT - 0.5) * len_prev < 0.1)
			prevT = 0.5; // don't create extremely short line segments
		mid_prev = misc::lerp(corner_prev, corner, prevT);
		mid_next = misc::lerp(corner, corner_next, nextT);
	}

	//center of circular arc
	auto tangent_prev = (corner - mid_prev).normalized();
	auto tangent_next = (mid_next - corner).normalized();
	double det = cross2(tangent_next, tangent_prev);
	auto dp = mid_next - mid_prev;
	double tangent_length;

	if (std::abs(det) < 0.001) {
		//degenerate corner
		tangent_length = std::min(len_prev, len_next) / 3;
	}
	else {
		double circle_radius = 1.0 / det * (tangent_next.x() * dp.x() + tangent_next.y() * dp.y());
		Eigen::Vector2d center(mid_prev.x() + circle_radius * tangent_prev.y(), mid_prev.y() - circle_radius * tangent_prev.x());
		double sector_angle = std::acos((mid_prev - center).dot(mid_next - center) / (circle_radius * circle_radius));
		tangent_length = 4.0 / 3.0 * std::tan(sector_angle / 4.0) * std::abs(circle_radius);
	}

	Matrix2Xd control_points(2, 4);
	control_points.col(0) = mid_prev;
	control_points.col(3) = mid_next;
	control_points.col(1) = mid_prev + tangent_length * tangent_prev;
	control_points.col(2) = mid_next - tangent_length * tangent_next;

		auto curve = new BezierCurve();
		curve->set_control_points(control_points);

	//Add a line before the bezier
	if (prevT > .5) {
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].emplace_back();

		GlobFitCurve_Line* line = new GlobFitCurve_Line();
		line->set_points(misc::lerp(corner_prev, corner, 0.5), mid_prev);

        const bool invert_coordinate_system = invert_edge_coordinate_system[Circular(invert_edge_coordinate_system, (Index)cornerId - 1)];
		auto param = new GlobFitLineParametrization(std::shared_ptr<GlobFitCurve_Line>(line), invert_coordinate_system, invert_coordinate_system);

		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().curve = param;
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.add_tangent_samples(line);
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.add_midpoint_samples(points, polygonV, corner_prev_id, 0.5, corner_prev_id, prevT,
			circular);
        initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.edge_src = BoundaryGraph::make_edge_id(corner_prev_id, cornerId);
        initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.edge_dst = BoundaryGraph::make_edge_id(corner_prev_id, cornerId);
	}

	//Add Bezier
	{
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].emplace_back();

        const bool invert_front_coordinate_system = invert_edge_coordinate_system[Circular(invert_edge_coordinate_system, (Index)cornerId - 1)];
        const bool invert_back_coordinate_system = invert_edge_coordinate_system[cornerId];
		auto param = new GlobFitBezierAngleBasedParametrization(std::shared_ptr<BezierCurve>(curve), invert_front_coordinate_system, invert_back_coordinate_system);
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().curve = param;

		auto cornerT = 0.5;
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.add_tangent_samples(param, cornerT);
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.add_midpoint_samples(points, polygonV, corner_prev_id, prevT, cornerId, nextT, circular);

		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.edge_src = BoundaryGraph::make_edge_id(corner_prev_id, cornerId);
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.edge_dst = BoundaryGraph::make_edge_id(cornerId, corner_next_id);
	}

	//Add a line after the bezier
	if (nextT < .5) {
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].emplace_back();

		GlobFitCurve_Line* line = new GlobFitCurve_Line();
		line->set_points(mid_next, misc::lerp(corner, corner_next, 0.5));

        const bool invert_coordinate_system = invert_edge_coordinate_system[cornerId];
		auto param = new GlobFitLineParametrization(std::shared_ptr<GlobFitCurve_Line>(line), invert_coordinate_system, invert_coordinate_system);

		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().curve = param;
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.add_tangent_samples(line);
		initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.add_midpoint_samples(points, polygonV, cornerId, nextT, cornerId, 0.5, circular);

        initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.edge_src = BoundaryGraph::make_edge_id(cornerId, corner_next_id);
        initial_fits[TANGENT_FIT_LERP_SYM][cornerId].back().fitting_info.edge_dst = BoundaryGraph::make_edge_id(cornerId, corner_next_id);
	}
}

template<>
void polyvec::CurveSequenceFitter::init_corner_fit<TANGENT_FIT_CONSTANT>(int cornerId) {
	auto corner_prev_id = poly_linidx(polygon, cornerId - 1);
	auto corner_next_id = poly_linidx(polygon, cornerId + 1);

	const Vector2d& corner_prev = polygon.col(corner_prev_id);
	const Vector2d& corner = polygon.col(cornerId);
	const Vector2d& corner_next = polygon.col(corner_next_id);

	const Vector2d mid_prev = misc::lerp(corner, corner_prev, .5);
	const Vector2d mid_next = misc::lerp(corner, corner_next, .5);

	GlobFitCurve_Line* line1 = new GlobFitCurve_Line();
	line1->set_points(mid_prev, corner);

	GlobFitCurve_Line* line2 = new GlobFitCurve_Line();
	line2->set_points(corner, mid_next);

    const bool invert_edge_before_coordinate_system = invert_edge_coordinate_system[Circular(invert_edge_coordinate_system, (Index)cornerId - 1)];
    const bool invert_edge_after_coordinate_system = invert_edge_coordinate_system[cornerId];
	auto param1 = new GlobFitLineParametrization(std::shared_ptr<GlobFitCurve_Line>(line1), invert_edge_before_coordinate_system, invert_edge_before_coordinate_system);
	auto param2 = new GlobFitLineParametrization(std::shared_ptr<GlobFitCurve_Line>(line2), invert_edge_after_coordinate_system, invert_edge_after_coordinate_system);
	param1->set_back_primary_axis(corner_normals[cornerId], true);
	param2->set_front_primary_axis(corner_normals[cornerId], true);	

#if !ALLOW_CORNERS_TO_MOVE
	param1->fix_parameter(2, 0.0);
	param1->fix_parameter(3, 0.0);
	param2->fix_parameter(0, 0.0);
	param2->fix_parameter(1, 0.0);
#endif

	initial_fits[TANGENT_FIT_CONSTANT][cornerId].emplace_back();
	initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().curve = param1;
	initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.add_midpoint_samples(points, polygonV, corner_prev_id, 0.5, corner_prev_id, 1.0, circular);
	initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.add_tangent_samples(line1);
    initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.edge_src = BoundaryGraph::make_edge_id(corner_prev_id, cornerId);
    initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.edge_dst = BoundaryGraph::make_edge_id(corner_prev_id, cornerId);
	if (edge_has_important_tangent(corner_prev_id))
		initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().hold_front_tangent = true;	

	initial_fits[TANGENT_FIT_CONSTANT][cornerId].emplace_back();
	initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().curve = param2;
	initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.add_midpoint_samples(points, polygonV, cornerId, 0.0, cornerId, 0.5, circular);
	initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.add_tangent_samples(line2);
    initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.edge_src = BoundaryGraph::make_edge_id(cornerId, corner_next_id);
    initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().fitting_info.edge_dst = BoundaryGraph::make_edge_id(cornerId, corner_next_id);
	if (edge_has_important_tangent(cornerId))
		initial_fits[TANGENT_FIT_CONSTANT][cornerId].back().hold_end_tangent = true;
}

// -------------------------------------------------------------------------
void CurveSequenceFitter::add_polygon_normals() {
	assert_break(circular);

	corner_normals.resize(polygon.cols());

	for (Index i = 0; i < polygon.cols(); ++i) {
		const Vector2d corner_prev_pt = polygon.col(poly_linidx(polygon, i - 1));
		const Vector2d corner_pt = polygon.col(i);
		const Vector2d corner_next_pt = polygon.col(poly_linidx(polygon, i + 1));

		auto dir_prev = (corner_pt - corner_prev_pt).normalized();
		auto dir_next = (corner_next_pt - corner_pt).normalized();

		// corner normals must allow fixing tangents at important edges (they provide the x-axis for lines).
		// At non-important edges, they can be arbitrary

		auto prev_is_important = edge_is_important[Circular(edge_is_important, i - 1)] || edge_has_important_tangent(Circular(edge_is_important, i - 1));
		auto next_is_important = edge_is_important[i] || edge_has_important_tangent(i);

		if ((prev_is_important && next_is_important))
		{
			auto n = polyvec::util::normal_dir(dir_prev + dir_next);
			
			corner_normals[i] = -n / n.dot(polyvec::util::normal_dir(dir_prev)); // allows to fix both tangents
		}
		else if (prev_is_important)
			corner_normals[i] = -polyvec::util::normal_dir(dir_prev);
		else if (next_is_important)
			corner_normals[i] = -polyvec::util::normal_dir(dir_next);
		else
			corner_normals[i] = -polyvec::util::normal_dir(dir_prev); //arbitrary
		
	}
}

void CurveSequenceFitter::add_corner_fit(CurvePrimitiveSequence& seq, const Eigen::Index corner, TangentFitType fitType) {
	static int next_primitive_id = 0;
	auto& primitiveConfiguration = initial_fits[fitType][corner];
	auto firstPrimitive = seq.primitives.size();

	for (auto& initialP : primitiveConfiguration) {
		CurvePrimitive primitive;
		primitive.corner = (int)corner;
		primitive.corner_pt = polygon.col(corner);
		primitive.curve.reset(initialP.curve->clone());
		primitive.fitting_info = initialP.fitting_info;

		primitive.hold_front_tangent = initialP.hold_front_tangent;
		primitive.hold_end_tangent = initialP.hold_end_tangent;		

		//debug
		primitive.endpoint_src = primitive.curve->get_curve()->pos(0.);
		primitive.endpoint_dst = primitive.curve->get_curve()->pos(1.);

		seq.primitives.emplace_back(std::move(primitive));
	}
}

template <typename T>
double calculate_objective_energy(std::vector<T>& objectives) {
	double sum = 0;

	for (auto& o : objectives) {
		Eigen::VectorXd obj;
		Eigen::MatrixXd deriv;

		o.compute_objective_and_jacobian(obj, deriv);

		for (int i = 0; i < obj.rows(); ++i) {
			sum += o.get_weight() * obj.row(i) * obj.row(i);
		}
	}

	return sum;
}

// -------------------------------------------------------------------------
CurveSequenceFitter::CurveSequenceFitter(
    const Eigen::Matrix2Xd& points, 
    const Eigen::Matrix2Xd& polygon, 
    const std::vector<Eigen::Index>& polygon_vertices,
	const polyfit::Regularity::RegularityInformation& regularity, 
    const int polygon_id,
    const std::vector<bool>& invert_edge_coordinate_system,
    const bool circular
) : 
    points(points), 
    polygon(polygon), 
    polygon_vertices(polygon_vertices), 
    circular(circular), 
    regularity(regularity), 
    polygon_id(polygon_id),
    invert_edge_coordinate_system(invert_edge_coordinate_system) {
	assert_break(polygon.cols() == polygon_vertices.size());

    if (this->invert_edge_coordinate_system.empty()) {
        this->invert_edge_coordinate_system.resize(polygon_vertices.size(), false);
    }

	polygonV.resize(polygon_vertices.size());	

	for (int i = 0; i < polygon_vertices.size(); ++i) {
		polygonV(i) = polygon_vertices[i];
        PF_VERBOSE_F("Index %d - %f %f - %f %f", polygonV(i), polygon(0, i), polygon(1, i), points(0, polygonV(i)), points(1, polygonV(i)));
	}

	geom::aabb raster_aabb;
	for (int i = 0; i < points.cols(); ++i)
		raster_aabb.add(points.col(i));
	raster_aabb_diagonal = (raster_aabb.max - raster_aabb.min).norm();

	edge_is_important = get_important_or_axis_aligned_edges(polygon, regularity);	

	// Measure the polygon
	polygon_measures.resize(polygon.cols());
	for (int i = 0; i < polygon.cols(); ++i) {
		polyfit::CurveTracer::measure_accuracy_signed_extended_polygon_fit_asymmetric(points, polygonV, i, polygon_measures[i], circular);
	}

	add_polygon_normals();
	init_all_fits();

	edge_angles = find_edge_angles_from_parallel(polygon, regularity);

	PF_VERBOSE_F("curve-sequence-fitter: points %lld", points.cols());
	PF_VERBOSE_F("curve-sequence-fitter: polygon corners %lld ", polygon.cols());
}

bool can_merge_consecutive_lines(
    const Eigen::Vector2d& tangent_0,
    const Eigen::Vector2d& tangent_1
) {
    return (tangent_0.normalized() - tangent_1.normalized()).squaredNorm() < 0.001;
}

void merge_consecutive_parallel_lines(CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, int& primitiveTo, const std::vector<bool>& prevent_merge) {
	for (size_t i = primitiveFrom; i < primitiveTo; ++i) {
		if (circular || i < primitiveTo - 1) {
			int curve_id = (int)i;
			int next_curve_id = (int)((i + 1) % seq.primitives.size());

			//merge consecutive parallel lines
			auto current_line = dynamic_cast<GlobFitLineParametrization*>(seq.primitives[curve_id].curve.get());
			auto next_line = dynamic_cast<GlobFitLineParametrization*>(seq.primitives[next_curve_id].curve.get());
			if (current_line && next_line)
			{
                if (!prevent_merge.empty() &&
                    (prevent_merge[BoundaryGraph::unpack_edge_id(seq.primitives[curve_id].fitting_info.edge_dst)(1)] ||
                    prevent_merge[BoundaryGraph::unpack_edge_id(seq.primitives[curve_id].fitting_info.edge_src)(0)])) {
                    PF_VERBOSE_F("Skipping merge for primitive %d", i);
                    continue;
                }

				auto dir1 = current_line->get_curve()->dposdt(0.0);
				auto dir2 = next_line->get_curve()->dposdt(0.0);
				if (can_merge_consecutive_lines(dir1, dir2))
				{
					// lines are parallel
					double old_end_t = dir1.norm() / (dir1.norm() + dir2.norm());
					auto& fit1 = seq.primitives[curve_id].fitting_info;
					auto& fit2 = seq.primitives[next_curve_id].fitting_info;					
					current_line->merge_with(next_line);
					seq.primitives[curve_id].fits_next_corner = true;
					seq.primitives[curve_id].hold_front_tangent |= seq.primitives[next_curve_id].hold_front_tangent;
					seq.primitives[curve_id].hold_end_tangent |= seq.primitives[next_curve_id].hold_end_tangent;
					fit1.fit_midpoints.insert(fit1.fit_midpoints.end(), fit2.fit_midpoints.begin(), fit2.fit_midpoints.end());
					fit1.fit_midpoint_normals.insert(fit1.fit_midpoint_normals.end(), fit2.fit_midpoint_normals.begin(), fit2.fit_midpoint_normals.end());

					auto merge_tangents = [&](FittingInfo::Tangents& source, FittingInfo::Tangents& add)
					{
						source.fit_tangents.insert(source.fit_tangents.end(), add.fit_tangents.begin(), add.fit_tangents.end());
						for (auto& t : source.fit_tangent_ts)
							t *= old_end_t;
						for (auto& t : add.fit_tangent_ts)
							t = old_end_t + t * (1 - old_end_t);
						source.fit_tangent_ts.insert(source.fit_tangent_ts.end(), add.fit_tangent_ts.begin(), add.fit_tangent_ts.end());
					};					

					merge_tangents(fit1.dense_tangents, fit2.dense_tangents);
					merge_tangents(fit1.sparse_tangents, fit2.sparse_tangents);

					seq.primitives.erase(seq.primitives.begin() + next_curve_id);
					--primitiveTo;					
				}
			}
		}
	}
}

// -------------------------------------------------------------------------
void CurveSequenceFitter::solve(CurvePrimitiveSequence& seq, bool circular, bool allow_parallel_handling)
{
	solve(seq, circular, 0, allow_parallel_handling);
}
void CurveSequenceFitter::solve(CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, bool allow_parallel_handling)
{
	int primitiveTo = seq.primitives.size();
	solve(seq, circular, primitiveFrom, primitiveTo, allow_parallel_handling);
}

std::vector<bool> get_important_or_axis_aligned_edges(const Eigen::Matrix2Xd& polygon, const polyfit::Regularity::RegularityInformation& regularity)
{
	std::vector<bool> edge_is_important(polygon.cols(), false);
	for (auto& r : regularity.important_edges())
	{
		edge_is_important[r.v0] = true;
	}
	//also add axis-aligned edges
	for (int i = 0; i < polygon.cols(); ++i)
	{
		const Eigen::Vector2d d = CircularAt(polygon, i + 1) - polygon.col(i);
		if (d.cwiseAbs().minCoeff() <= PF_EPS)
			edge_is_important[i] = true;
	}
	return edge_is_important;
}

std::vector<double> find_edge_angles_from_parallel(const Eigen::Matrix2Xd& polygon, const polyfit::Regularity::RegularityInformation& regularity)
{
	struct EntryData
	{
		vec2 dir;
		int merge_count = 0;
	};

	UnionFind<EntryData> uf(polygon.cols());

	// initialize union-find with original edge directions
	for (int i = 0; i < polygon.cols(); ++i)
		uf[i].dir = CircularAt(polygon, i + 1) - polygon.col(i);

	// find edge angles by averaging connected edges (weighted somewhat by edge length)
	for (auto& p : regularity.parallels())
	{
		auto edge0 = p.v00;
		auto edge1 = p.v10;

		auto dir0 = uf[edge0].dir;
		auto dir1 = uf[edge1].dir;

		uf.merge(edge0, edge1);
		auto merged = uf.getRepresentative(edge0);

		uf[merged].merge_count++;
		if (dir0.dot(dir1) > 0)
			uf[merged].dir = dir0 + dir1;
		else
			uf[merged].dir = dir0 - dir1;
	}

	// find the final edge angles
	std::vector<double> edge_angles(polygon.cols(), std::numeric_limits<double>::quiet_NaN());
	for (int i = 0; i < polygon.cols(); ++i)
	{
		auto rep = uf.getRepresentative(i);

		if (uf[rep].merge_count == 0)
			continue;

		Vector2d avg_dir = uf[rep].dir;
		const Vector2d orig_dir = CircularAt(polygon, i + 1) - polygon.col(i);

		// find the correct orientation
		if (orig_dir.dot(avg_dir) < 0)
			avg_dir *= -1;

		// calculate the corresponding angle
		edge_angles[i] = std::atan2(avg_dir.y(), avg_dir.x());
	}

	return edge_angles;
}

void addCurveSequenceToFitter
	(polyvec::CurveFitter& curve_fitter, CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, int primitiveTo, int polygon_corners, 
		const std::vector<bool>& is_edge_important_or_axisaligned, const std::vector<double>& edge_angles, bool allow_parallel_handling)
{	
	for (int i = primitiveFrom; i < primitiveTo; ++i)
		curve_fitter.add_curve(&seq.primitives[i]);

	if (!circular)
	{
		seq.primitives[primitiveFrom].curve->reduce_degrees_of_freedom(DofOptions::FIX_FRONT_TANGENT | DofOptions::KEEP_FRONT_ON_BISECTOR);
		seq.primitives[primitiveTo - 1].curve->reduce_degrees_of_freedom(DofOptions::FIX_BACK_TANGENT | DofOptions::KEEP_BACK_ON_BISECTOR);
	}	

	for (size_t i = primitiveFrom; i < primitiveTo; ++i) {
		int curve_id = (int)i;
		auto& prim = seq.primitives[curve_id];
		GlobFitCurveParametrization* curve = prim.curve.get();		

		if (prim.fitting_info.edge_src == prim.fitting_info.edge_dst)
		{
			auto src_edge = BoundaryGraph::unpack_edge_id(prim.fitting_info.edge_src);
			if (is_edge_important_or_axisaligned[src_edge(0)])
			{
				// This is a line segment on an important edge: fix the tangent
				curve->reduce_degrees_of_freedom(DofOptions::FIX_FRONT_TANGENT);				
			}

			if (curve_fitter.get_fit_options().consider_parallel && allow_parallel_handling)
			{
				if (!std::isnan(edge_angles[src_edge(0)]))
				{
					auto lineParam = static_cast<GlobFitLineParametrization*>(curve);
					auto line = std::static_pointer_cast<GlobFitCurve_Line>(lineParam->get_curve());
					auto d = line->get_points().col(1) - line->get_points().col(0);
					auto current_angle = std::atan2(d.y(), d.x());
					// how much must the endpoints move to get to the prescribe angle
					auto distance = std::abs(current_angle - edge_angles[src_edge(0)]) * d.norm() / 2;

					if (distance < 0.75)
						// we want to prescribe the angle of this line
						curve_fitter.prescribe_angle(lineParam, edge_angles[src_edge(0)]);
				}
			}
		}

		if (circular || i < primitiveTo - 1) {
			int next_curve_id = (int)((i + 1) % seq.primitives.size());
			auto& prim_next = seq.primitives[next_curve_id];
			GlobFitCurveParametrization* curve_next = prim_next.curve.get();

			curve_fitter.make_g0(curve, curve_next);

			//C0-corners will not be affected by make_g1(), so we can just call it
			curve_fitter.make_g1(curve, curve_next);

			curve_fitter.make_g2(curve, curve_next);


			if (curve_fitter.get_fit_options().keep_endpoints_on_bisectors)
			{
				if (prim.fits_same_corner(prim_next, polygon_corners) ||
					(dynamic_cast<GlobFitBezierAngleBasedParametrization*>(curve)
						&& dynamic_cast<GlobFitBezierAngleBasedParametrization*>(curve_next)))
				{
					curve->reduce_degrees_of_freedom(DofOptions::KEEP_BACK_ON_BISECTOR);
					curve_next->reduce_degrees_of_freedom(DofOptions::KEEP_FRONT_ON_BISECTOR);
				}
				else
				{
					// restrict how far the curves can move - not yet used

					auto line_cur = dynamic_cast<GlobFitLineParametrization*>(curve);
					auto line_next = dynamic_cast<GlobFitLineParametrization*>(curve_next);

					// tangent is secondary axis
					if (line_cur)
					{
						auto points = static_pointer_cast<GlobFitCurve_Line>(line_cur->get_curve())->get_points();
						auto& coordinate_system = line_cur->get_back_coordinate_system();
						auto hard_limit = (points.col(0) - points.col(1)).dot(coordinate_system.secondary());
						auto soft_limit = 0.1 * hard_limit;
					}

					if (line_next)
					{
						auto points = static_pointer_cast<GlobFitCurve_Line>(line_next->get_curve())->get_points();
						auto& coordinate_system = line_next->get_front_coordinate_system();
						auto hard_limit = (points.col(1) - points.col(0)).dot(coordinate_system.secondary());
						auto soft_limit = 0.1 * hard_limit;
					}

				}
			}
		}
	}
}

void CurveSequenceFitter::solve ( CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, int& primitiveTo, bool allow_parallel_handling) {
    if ( circular ) {
        assert_break ( primitiveFrom == 0 && primitiveTo == seq.primitives.size() );
    }

	merge_consecutive_parallel_lines(seq, circular, primitiveFrom, primitiveTo);	

	CurveFitter curve_fitter(raster_aabb_diagonal, circular);

	addCurveSequenceToFitter(curve_fitter, seq, circular, primitiveFrom, primitiveTo, polygon.cols(), edge_is_important, edge_angles, allow_parallel_handling);
    
    curve_fitter.solve_with_stiffening();
}

const auto isBezier = [] ( const CurvePrimitive& p ) {
    return dynamic_cast<GlobFitBezierAngleBasedParametrization*> ( p.curve.get() ) != nullptr;
};
CurvePrimitive* CurveSequenceFitter::first_bezier(CurvePrimitiveSequence& seq, int firstPrimitive) const
{
	//find the first Bezier
	auto firstBezier = std::find_if(seq.primitives.begin() + firstPrimitive, seq.primitives.end(), isBezier);

	if (firstBezier != seq.primitives.end()) {
		return &(*firstBezier);
	}
	else
		return nullptr;
};
CurvePrimitive* CurveSequenceFitter::last_bezier(CurvePrimitiveSequence& seq, int firstPrimitive) const
{
	auto primitivesInSequence = seq.primitives.size() - firstPrimitive;
	auto end = seq.primitives.rbegin() + primitivesInSequence;
	auto lastBezier = std::find_if(seq.primitives.rbegin(), end, isBezier);

	if (lastBezier != end) {
		return &(*lastBezier);
	}
	else
		return nullptr;
};

void CurveSequenceFitter::fix_front ( CurvePrimitiveSequence& seq, int firstPrimitive ) {
    //find the first Bezier
	auto firstBezier = first_bezier(seq, firstPrimitive);
	if(firstBezier)
        firstBezier->curve->fix_parameter ( 2, 0.0 );
}

void CurveSequenceFitter::fix_back ( CurvePrimitiveSequence& seq, int firstPrimitive ) {
    //find the last Bezier
	auto lastBezier = last_bezier(seq, firstPrimitive);

    if ( lastBezier ) {
        lastBezier->curve->fix_parameter ( 4, 0.0 );
    }
}

CurvePrimitiveSequence CurveSequenceFitter::all_with_type ( TangentFitType type, bool optimize, bool fixFront, bool fixBack ) {
    CurvePrimitiveSequence seq;

    for ( int i = 0; i < polygon.cols(); ++i ) {
        const int first_prim = (int) seq.primitives.size();
        add_corner_fit ( seq, i, type );

        if ( fixFront ) {
            fix_front ( seq, first_prim );
        }

        if ( fixBack ) {
            fix_back ( seq, first_prim );
        }		

        if ( optimize ) {
            solve ( seq, false, first_prim, false );
        }
    }

    return std::move ( seq );
}

CurvePrimitiveSequence CurveSequenceFitter::fit_individual_corners(
	std::vector<std::pair<int, TangentFitType>> corners,
	bool circular,
	bool _fix_front,
	bool _fix_back,
	bool optimize,
	bool allow_parallel_handling
) {
	CurvePrimitiveSequence seq;

	for (auto& c : corners) {
		add_corner_fit(seq, c.first, c.second);
	}
	/*add_corner_fit(seq, 0, TANGENT_FIT_CONSTANT);
	add_corner_fit(seq, 1, TANGENT_FIT_LERP_SYM);
	circular = false;*/

	if (_fix_front) {
		fix_front(seq, 0);
	}

	if (_fix_back) {
		fix_back(seq, 0);
	}

	if (optimize) {
		solve(seq, circular, allow_parallel_handling);
	}

#if ADAPTIVE_POINT_WEIGHT
	if (circular)
	{
		bool changed = true;
		while (changed)
		{
			changed = false;
			for (auto& p : seq.primitives)
			{
				if (p.error.accuracy.max_error() > ADAPTIVE_ACCURACY_CHECK_THRESHOLD && p.point_weight_multiplier <= 128)
				{								
					p.point_weight_multiplier *= 2;
					changed = true;
				}
			}
			if (changed)
				solve(seq, circular, allow_parallel_handling);
		}
	}
#endif

    return std::move ( seq );
}

// --------  Evolutionary fitting  -----------------------
//Represents a consecutive sequence of corners on the polygon using the
//indices of the first (inclusive) and last (inclusive) corner. The sequence
//is circular, i.e. the first index can be greater than the last index.
CurveSequenceFitter::CornerSequence::CornerSequence()
    : _first_incl ( -1 ), _last_incl ( -1 )
{ }

CurveSequenceFitter::CornerSequence::CornerSequence ( int first_incl, int last_incl )
    : _first_incl ( first_incl ), _last_incl ( last_incl )
{ }

int CurveSequenceFitter::CornerSequence::first_incl() const {
    return _first_incl;
}
int CurveSequenceFitter::CornerSequence::last_incl() const {
    return _last_incl;
}

//TFunc: void(int)
template <typename TFunc>
void CurveSequenceFitter::CornerSequence::for_each_corner ( int polygon_corners, TFunc&& callback ) {
    int i = first_incl();

    while ( true ) {
        if ( i >= polygon_corners ) {
            break;    //circular
        }

        std::forward<TFunc> ( callback ) ( i );

        if ( i == last_incl() || ( last_incl() == polygon_corners && i + 1 == polygon_corners ) ) {
            break;
        }

        i = ( i + 1 ) % polygon_corners;
    }
}

int CurveSequenceFitter::CornerSequence::included_corners ( int polygon_corners ) const {
    int corners = _last_incl + 1 - _first_incl;

    if ( corners <= 0 ) {
        corners += polygon_corners;
    }

    return corners;
}

CurveSequenceFitter::EvolutionaryFittingState::EvolutionaryFittingState ( const Eigen::Matrix2Xd& rasterPoints, const Eigen::VectorXi& polygonPoints,
		bool circular )
    : cornerBelief ( polygonPoints.size() ) {
    for ( auto& e : cornerBelief ) {
        e.set();    //make all options viable
    }

    polygon_corner_angles.resize ( polygonPoints.size() );
    polygon_raster_corner_sizes.resize ( polygonPoints.size() );
    corner_good_probability.resize ( polygonPoints.size() );

    for ( int i = 0; i < polygonPoints.size(); ++i ) {
        auto corner = rasterPoints.col ( polygonPoints ( i ) );
        auto prev   = rasterPoints.col ( polygonPoints ( ( i - 1 + polygonPoints.size() ) % polygonPoints.size() ) );
        auto next   = rasterPoints.col ( polygonPoints ( ( i + 1 ) % polygonPoints.size() ) );
        auto to_prev = ( prev - corner ).normalized();
        auto to_next = ( next - corner ).normalized();
        polygon_corner_angles[i] = std::acos ( std::min ( 1.0, std::max ( -1.0, to_prev.dot ( to_next ) ) ) );
        polygon_raster_corner_sizes[i] = polyfit::GeomRaster::raster_corner_size ( rasterPoints, polygonPoints ( i ) );
    }
}

bool CurveSequenceFitter::EvolutionaryFittingState::is_corner_certain ( int corner ) const {
    return cornerBelief[corner].count() == 1;
};
TangentFitType CurveSequenceFitter::EvolutionaryFittingState::get_highest_priority_possible_type ( int corner ) const {
    for ( int i = 0; i < TANGENT_FIT_SAMPLES_COUNT; ++i )
        if ( cornerBelief[corner][i] ) {
            return TangentFitType ( i );
        }

    assert_break ( false );
}

std::vector<polyfit::CurveTracer::AccuracyMeasurement> CurveSequenceFitter::get_polygon_measures() const {
    std::vector<polyfit::CurveTracer::AccuracyMeasurement> polygonMeasures ( polygon.cols() );

    for ( int i = 0; i < polygon.cols(); ++i ) {
        polyfit::CurveTracer::measure_accuracy_signed_extended_polygon_fit_asymmetric ( points, polygonV, i, polygonMeasures[i], circular );
    }

    return polygonMeasures;
}

std::vector<CurveSequenceFitter::CornerSequence> CurveSequenceFitter::EvolutionaryFittingState::get_connected_sequences() const {
    std::vector<CornerSequence> connectedSequences;
    int firstBreakPoint = -1; //break before this corner
    int currentSequenceStart = -1;
    int i = 0;

    while ( ! ( firstBreakPoint == -1 && i >= cornerBelief.size() ) ) {
        auto prev = ( i - 1 + cornerBelief.size() ) % cornerBelief.size();

        if ( !is_corner_certain ( prev % cornerBelief.size() ) && !is_corner_certain ( i % cornerBelief.size() ) ) {
            //we need to break here
            if ( currentSequenceStart != -1 ) {
                connectedSequences.push_back ( CornerSequence(currentSequenceStart, ( i - 1 ) % cornerBelief.size()) );
            }

            currentSequenceStart = i;

            if ( firstBreakPoint == -1 ) {
                firstBreakPoint = i;
            }
        }

        if ( i == firstBreakPoint + cornerBelief.size() ) {
            break;
        }

        ++i;
    }

    if ( connectedSequences.empty() ) {
        connectedSequences.emplace_back ( 0, cornerBelief.size() );    //cornerBelief.size() on purpose to denote no break
    }

    return connectedSequences;
}    

CurvePrimitiveSequence CurveSequenceFitter::fit_evolutionary_simple ( std::vector<TangentFitType>& outCornerTypes,
        const std::function<void ( const EvolutionaryFittingState& state, const FittingAttempt* fits ) >& callback ) {
        
    FitClassifierRandomForest classifier;
    classifier.load_from_default_file();
	return fit_evolutionary_simple(classifier, outCornerTypes, callback);
}

CurvePrimitiveSequence CurveSequenceFitter::fit_evolutionary_simple(FitClassifier& classifier, std::vector<TangentFitType>& outCornerTypes, const std::function<void(const EvolutionaryFittingState& state, const FittingAttempt* fits)>& callback) {
    classify_evolutionary_simple(classifier, outCornerTypes, callback);

    std::vector<std::pair<int, TangentFitType>> corner_fits(polygon.cols());
    for (int i = 0; i < corner_fits.size(); ++i) {
        corner_fits[i] = std::make_pair(i, outCornerTypes[i]);
    }

	return fit_individual_corners(corner_fits, true, false, false, true, true);
}

void CurveSequenceFitter::classify_evolutionary_simple(
    FitClassifier& classifier,
    std::vector<TangentFitType>& outCornerTypes,
    const std::function<void(const EvolutionaryFittingState& state, const FittingAttempt* fits)>& callback
) {
    EvolutionaryFittingState state(points, polygonV, circular);

    // probability for every type and every corner
    std::vector<std::vector<double>> fit_probabilities(TANGENT_FIT_SAMPLES_COUNT);

	std::vector<FitClassifier::PolygonCornerInfo> polygon_corner_info(polygon.cols());
	for (int corner = 0; corner < polygon.cols(); ++corner) {
		polygon_corner_info[corner] = FitClassifier::PolygonCornerInfo(
			polygon_measures[corner],
			(polygon.col(corner) - polygon.col((corner - 1 + polygon.cols()) % polygon.cols())).norm(),
			(polygon.col(corner) - polygon.col((corner + 1) % polygon.cols())).norm(),
			state.polygon_corner_angles[corner], state.polygon_raster_corner_sizes[corner],
			polyfit::AngleUtils::number_of_neighbors_with_different_convexity(polygon, corner));
	}

#if POLYVEC_INCLUDE_CONSTANT_TANGENTS_TO_FIT_ATTEMPTS
    for (int type = TANGENT_FIT_LERP; type <= TANGENT_FIT_CONSTANT; ++type) {
#else
    for (int type = TANGENT_FIT_LERP; type < TANGENT_FIT_CONSTANT; ++type) {
#endif
        FittingAttempt fits[3];

        for (int j = 0; j < 3; ++j) {
            bool fixFront = j & 1 == 1;
            bool fixEnd = (j >> 1) & 1 == 1;

            fits[j].primitive_seq = all_with_type(TangentFitType(type), true, fixFront, fixEnd);
			fits[j].measure_sequence(polygon);
        }

        fit_probabilities[type].resize(polygon.cols());

        for (int corner = 0; corner < polygon.cols(); ++corner) {            

            FitClassifier::FitInfo fit_info[4];

            for (int j = 0; j < 3; ++j) {
                int source_idx = j;

                fit_info[j].error = fits[source_idx].corner_errors[corner];
                fit_info[j].distance_to_corner = fits[source_idx].distances_to_corner[corner];
            }

            state.corner_good_probability[corner] = classifier.evaluate_fit_probability(fit_info, polygon_corner_info[corner], polygon, corner);
            fit_probabilities[type][corner] = state.corner_good_probability[corner];
        }

		if(callback)
			callback(state, fits);
    }

    //make a decision on the final fit
    outCornerTypes.resize(polygon.cols());

    for (int corner = 0; corner < polygon.cols(); ++corner) {
        if (fit_probabilities[TANGENT_FIT_LERP][corner] > 0.75) {
            outCornerTypes[corner] = TANGENT_FIT_LERP;
        }
        else if (fit_probabilities[TANGENT_FIT_LERP_SYM][corner] > 0.75) {
            outCornerTypes[corner] = TANGENT_FIT_LERP_SYM;
        }
        else if (fit_probabilities[TANGENT_FIT_LERP][corner] < 0.25 && fit_probabilities[TANGENT_FIT_LERP_SYM][corner] < 0.25) {
            outCornerTypes[corner] = TANGENT_FIT_CONSTANT;
        }
        else {
            PF_VERBOSE_F("Corner %d is unclear", (int)corner);
            bool found_fit = false;

            //unclear case
            for (int type = 0; type < TANGENT_FIT_CONSTANT; ++type)
                if (fit_probabilities[type][corner] > 0.5) {
                    outCornerTypes[corner] = TangentFitType(type);
                    found_fit = true;
                    break;
                }

            if (!found_fit) {
                outCornerTypes[corner] = TANGENT_FIT_CONSTANT;
            }
        }
    }

    outCornerTypes = get_regularity_actions<HandledRegularityNewClassifications>(polygon, regularity, outCornerTypes, polygon_id);

#if POST_ACCURACY_CHECK
	{
		bool changed = true;
		while (changed)
		{
			changed = false;
			// Check if the final fit is accurate enough. If not, downgrade asymmetric fits.
			std::vector<std::pair<int, TangentFitType>> corner_fits(polygon.cols());
			for (int i = 0; i < corner_fits.size(); ++i) {
				corner_fits[i] = std::make_pair(i, outCornerTypes[i]);
			}
			FittingAttempt fit;
			fit.primitive_seq = fit_individual_corners(corner_fits, true, false, false);
			fit.measure_sequence(polygon);


			for (int i = 0; i < polygon.cols(); ++i)
			{
				if (outCornerTypes[i] > TANGENT_FIT_LERP_SYM)
					continue;

				FitClassifier::FitInfo fit_info[3];

				for (int j = 0; j < 3; ++j) {
					fit_info[j].error = fit.corner_errors[i];
					fit_info[j].distance_to_corner = fit.distances_to_corner[i];
				}

				auto good_probability = classifier.evaluate_fit_probability(fit_info, polygon_corner_info[i], polygon, i);

				if (good_probability < 0.5)
				{
					outCornerTypes[i] = (TangentFitType)(outCornerTypes[i] + 1);
					changed = true;
				}
			}
			if(changed)
				outCornerTypes = get_regularity_actions<HandledRegularityNewClassifications>(polygon, regularity, outCornerTypes, polygon_id);
		}
	}
#endif
}

CurvePrimitiveSequence CurveSequenceFitter::fit_initial_guess(
    const std::vector<TangentFitType>& tangentFits
) {
    std::vector<std::pair<int, TangentFitType>> corner_fits(polygon.cols());
    for (int i = 0; i < corner_fits.size(); ++i) {
        corner_fits[i] = std::make_pair(i, tangentFits[i]);
    }

    CurvePrimitiveSequence seq = fit_individual_corners(corner_fits, true, false, false, false);
    int primitiveTo = seq.primitives.size() - 1;

    // This makes things more painful than necessary
    // merge_consecutive_parallel_lines(seq, true, 0, primitiveTo);

    return move(seq);
}

void CurveSequenceFitter::FittingAttempt::measure_sequence(const Eigen::Matrix2Xd& polygon)
{
	resize(polygon.cols());
	for (int i = 0; i < primitive_seq.primitives.size(); ++i) {
		auto& p = primitive_seq.primitives[i];
		corner_errors[p.corner].combine(p.error);

		auto corner_t = p.curve->get_curve()->project(polygon.col(p.corner));
		distances_to_corner[p.corner] = std::min(distances_to_corner[p.corner],
			(polygon.col(p.corner) - p.curve->get_curve()->pos(corner_t)).norm());
	}
}


NAMESPACE_END(polyvec)