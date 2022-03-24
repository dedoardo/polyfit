// polyvec
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/regularity/check_opposite_points_winding_number.hpp>
#include <polyvec/regularity/are_points_facing_inside.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curve_line.hpp>
#include <polyvec/curve-tracer/measure_hausdorff.hpp>
#include <polyvec/curve-tracer/curve_objectives.hpp>
#include <polyvec/curve-tracer/measure_curvature.hpp>
#include <polyvec/options.hpp>
#include <polyvec/curve-tracer/curvature_variation_optimize_coordinate_descent.hpp>
#include <polyvec/regularity/continuations.hpp>
#include <polyvec/utils/directions.hpp>
#include <polyvec/core/constants.hpp>

// libc++
#include <algorithm>
#include <unordered_set>

using namespace std;
using namespace Eigen;
using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

std::string Parallel::to_string() const
{
	return StringUtils::fmt("%Parallel (%d %d) (%d %d)", v00, v01, v10, v11);
}

std::string Continuation::to_string() const
{
	return StringUtils::fmt("%Continuation (%d -> %d -> %d -> %d)", v0, v1);
}

std::string Symmetry::to_string() const
{
	return StringUtils::fmt("%Symmetry (%d %d), relaxed: %d", v0, v1, is_relaxed);
}

std::string ImportantEdge::to_string() const
{
	return StringUtils::fmt("%Important edge (%d)", v0);
}

void order_symmetry_vertices_to_canonical(
	const mat2x& B,
	std::vector<polyfit::Symmetry::SubPath>& S,
	bool circular
) {
	for (size_t i = 0; i < S.size(); ++i) {
		auto& s = S[i];

		// reordering the symmetry points consistently with the orientation of the boundary
		bool joint0 = polyfit::Symmetry::is_subpath_joint0(B, s, circular);
		bool joint1 = polyfit::Symmetry::is_subpath_joint1(B, s, circular);

		// making v00 and v10 the consecutive vertices
		if (joint1 && !joint0) {
			swap(s.v00, s.v01);
			swap(s.v10, s.v11);
		}

		// making the direction +1 from v11 -> v01
		if (CircularDist(B, s.v01, s.v00) > CircularDist(B, s.v00, s.v01)) {
			swap(s.v01, s.v11);
			swap(s.v00, s.v10);
		}
	}
}


// Given two edges, it tests whether the opposite points are aligned with respect
// to the edge interpolating the two parallel edges
// ----------------------------------------------------------------------------
void test_and_set_opposite_points_alignment(
	const mat2x& B,
	const std::vector<int>& C,
	const vecXi& V,
	Parallel& r
) {
	const vec2 p00 = B.col(V(r.v00));
	const vec2 p01 = B.col(V(r.v01));
	const vec2 p10 = B.col(V(r.v10));
	const vec2 p11 = B.col(V(r.v11));

	const vec2 d_axis_0 = (p00 - p01).normalized();
	const vec2 d_axis_1 = (p11 - p10).normalized();
	const double a_axis_0 = atan2(d_axis_0.y(), d_axis_0.x());
	const double a_axis_1 = atan2(d_axis_1.y(), d_axis_1.x());
	const double a_axis_lerp = .5 * (a_axis_0 + a_axis_1);
	vec2 d_axis_lerp(cos(a_axis_lerp), sin(a_axis_lerp));
	d_axis_lerp.normalize(); // ehm...

	const double t_project_00 = abs(LineUtils::project_t(p00, p11, p11 + d_axis_lerp));
	const double t_project_11 = abs(LineUtils::project_t(p11, p00, p00 + d_axis_lerp));
	const double t_project_01 = abs(LineUtils::project_t(p01, p10, p10 + d_axis_lerp));
	const double t_project_10 = abs(LineUtils::project_t(p10, p01, p01 + d_axis_lerp));

	double t_project_scaled = PV_REGULARITY_MAX_ALIGNMENT_DISTANCE;
	t_project_scaled *= GeomRaster::get_resolution_scaling_factor_32x32(B);

	r.aligned_00_11 = t_project_00 < t_project_scaled && t_project_11 < t_project_scaled;
	r.aligned_01_10 = t_project_01 < t_project_scaled && t_project_10 < t_project_scaled;

	r.aligned_00_11 = r.aligned_00_11 || GeomRaster::are_points_axis_aligned(p00, p11);
	r.aligned_01_10 = r.aligned_01_10 || GeomRaster::are_points_axis_aligned(p01, p10);

	// The need to be both convex
	r.aligned_00_11 = r.aligned_00_11 && C[r.v00] == -1 && C[r.v11] == -1;
	r.aligned_01_10 = r.aligned_01_10 && C[r.v01] == -1 && C[r.v10] == -1;
}

// Represents a uniform grid in angle space
class AngularGrid
{
public:
	AngularGrid(int bucket_count)
		: buckets(bucket_count), bucket_angle(2 * M_PI / bucket_count)
	{

	}

	void add_entry(Index identifier, double angle)
	{
		angle = normalize_angle(angle);
		auto bucket = get_bucket(angle);
		buckets[bucket].emplace_back(identifier, angle);
	}

	// TFunc: void(int identifier, double angle)
	template <typename TFunc>
	void for_each_entry(TFunc&& callback) const
	{
		for (auto& bucket : buckets)
			for(auto& entry : bucket)
				std::forward<TFunc>(callback)(entry.first, entry.second);
	}

	// TFunc: void(int identifier, double angle)
	template <typename TFunc>
	void for_each_entry_in_angle_range(double lower_angle, double upper_angle, TFunc&& callback) const
	{
		lower_angle = normalize_angle(lower_angle);
		upper_angle = normalize_angle(upper_angle);

		int lower_bucket = get_bucket(lower_angle);
		int upper_bucket = get_bucket(upper_angle);

		int bucket = lower_bucket;
		while (true)
		{
			for (auto& entry : buckets[bucket])
			{
				if (is_angle_in_range(entry.second, lower_angle, upper_angle))
					std::forward<TFunc>(callback)(entry.first, entry.second);
			}

			if (bucket == upper_bucket)
				break;

			bucket = Circular(buckets, bucket + 1);
		}
	}

private:

	double normalize_angle(double angle) const
	{
		if (angle < 0)
			angle += 2 * M_PI;
		if (angle >= 2 * M_PI)
			angle -= 2 * M_PI;
		return angle;
	}

	bool is_angle_in_range(double angle, double lower_angle, double upper_angle) const
	{
		if (lower_angle < upper_angle)
			return angle >= lower_angle && angle <= upper_angle;
		else
			return angle >= lower_angle || angle <= upper_angle;
	}

	int get_bucket(double normalized_angle) const
	{
		return (int)(normalized_angle / bucket_angle);
	}

	const double bucket_angle;
	std::vector<std::vector<std::pair<Index, double>>> buckets; // entries are pairs of identifier and angle
};

// ----------------------------------------------------------------------------
void find_all_parallel_edges(
	const mat2x& B,
	const vecXi& V,
	const std::vector<int>& C,
	RegularityInformation& reg,
	const bool circular
) {
	mat2x P;
	BoundaryGraph::trace_to_points(B, V, P);

	// first, we insert all edges into angle buckets, so we can find parallel edges faster
	// make sure that we have an even number of buckets and each bucket fills an angle of at least PV_REGULARITY_PARALLEL_MAX_DEVIATION_ANGLE
	int bucket_count = 2 * (int)std::floor(M_PI / PV_REGULARITY_PARALLEL_MAX_DEVIATION_ANGLE);
	AngularGrid angular_grid(bucket_count);
	
	bool is_ccw;
	WindingNumber::compute_orientation(P, is_ccw);

	// pre-calculate edge directions, skip flat vertices
	std::vector<bool> is_flat_corner(V.size(), false);
	for (Index i = 0; i < V.size(); ++i)
	{
		if (C[i] == CircularAt(C, i + 1) && C[i] == CircularAt(C, i - 1))
			continue; //do not skip flat angles when they are part of a gradual direction change
		auto angle = AngleUtils::spanned_shortest(B.col(CircularAt(V, i - 1)), B.col(V(i)), B.col(CircularAt(V, i + 1)));
		is_flat_corner[i] = angle > PF_RAD(170);
	}

	for (Index i = 0; i < V.size(); ++i)
	{
		Index vi0 = i;
		Index vi1 = Circular(V, i + 1);
		//skip flat vertices
		if (is_flat_corner[vi0])
			vi0 = Circular(V, vi0 - 1); 
		if (is_flat_corner[vi1])
			vi1 = Circular(V, vi1 + 1);

		const Vector2d p0 = B.col(V(vi0));
		const Vector2d p1 = B.col(V(vi1));
		const Vector2d d = p1 - p0;

		// skip tiny edges
		if (d.squaredNorm() < 2. + PF_EPS)
			continue;

		// calculate the angle of the edge
		const double angle = std::atan2(d.y(), d.x());
		angular_grid.add_entry(i, angle);	
	}

	angular_grid.for_each_entry([&](int i, double angle_i)
	{
		Index vi0 = i;
		Index vi1 = Circular(V, i + 1);

		//skip flat vertices
		if (is_flat_corner[vi0])
			vi0 = Circular(V, vi0 - 1);
		if (is_flat_corner[vi1])
			vi1 = Circular(V, vi1 + 1);

		const vec2 pi0 = B.col(V(vi0));
		const vec2 pi1 = B.col(V(vi1));
		const vec2 di = (pi1 - pi0);
		auto length_i = di.norm();

		const Vector2d normal = polyvec::util::normal_dir(di) * (is_ccw ? 1.0 : -1.0);
		
		// find all edges that fulfill the angle criterion for parallelity
		angular_grid.for_each_entry_in_angle_range(angle_i + M_PI - PV_REGULARITY_PARALLEL_MAX_DEVIATION_ANGLE, angle_i + M_PI + PV_REGULARITY_PARALLEL_MAX_DEVIATION_ANGLE,
			[&](int j, double angle_j)
		{
			if (i >= j)
				return;

            Parallel r;
            r.v00 = i;
            r.v01 = Circular(V, i + 1);
            r.v10 = j;
            r.v11 = Circular(V, j + 1);
            const int dist_11_00 = min(CircularDist(V, r.v11, r.v00), CircularDist(V, r.v00, r.v11));
            const int dist_01_10 = min(CircularDist(V, r.v01, r.v10), CircularDist(V, r.v10, r.v01));

            if (dist_11_00 == 0 || dist_01_10 == 0) {
                return;
            }

            r.connected_00_11 = dist_11_00 <= 1;
            r.connected_01_10 = dist_01_10 <= 1;

			Index vj0 = j;
			Index vj1 = Circular(V, j + 1);

			//skip flat vertices
			if (is_flat_corner[vj0])
				vj0 = Circular(V, vj0 - 1);
			if (is_flat_corner[vj1])
				vj1 = Circular(V, vj1 + 1);

			const vec2 pj0 = B.col(V(vj0));
			const vec2 pj1 = B.col(V(vj1));
			const vec2 dj = (pj1 - pj0);
			auto length_j = dj.norm();

			// Connecting through outside
			if ((pj0 - pi0).dot(normal) < 0)
				return;

			auto di1j0 = (pj0 - pi1); //the vector connecting the parallel vertex pair	
			auto di0j1 = (pj1 - pi0);

			if (!are_points_facing_inside(P, vi0, vj0) || !are_points_facing_inside(P, vi0, vj1) ||
				!are_points_facing_inside(P, vi1, vj0) || !are_points_facing_inside(P, vi1, vj1)) {
				// the points don't see each other, there is something in between them
				return;
			}			

			if ((di.cwiseAbs().minCoeff() < PF_EPS && dj.cwiseAbs().minCoeff() >= PF_EPS) || (dj.cwiseAbs().minCoeff() < PF_EPS && di.cwiseAbs().minCoeff() >= PF_EPS))
				return; //if one of them is axis-aligned and the other is not			

			// Test overlap
			const double t00 = LineUtils::project_t(pi0, pj0, pj1);
			const double t01 = LineUtils::project_t(pi1, pj0, pj1);
			const double t10 = LineUtils::project_t(pj0, pi0, pi1);
			const double t11 = LineUtils::project_t(pj1, pi0, pi1);

			double overlap0, overlap1;
			if (!Num::test_and_calculate_interval_overlap(0., 1., min(t00, t01), max(t00, t01), overlap0) ||
				!Num::test_and_calculate_interval_overlap(0., 1., min(t10, t11), max(t10, t11), overlap1)) {
				return;
			}

			const double overlap0_relative = overlap0 * length_j / min(length_j, length_i);
			const double overlap1_relative = overlap1 * length_i / min(length_j, length_i);
			const double min_overlap_thr = .75;
			if (overlap0_relative < min_overlap_thr || overlap1_relative < min_overlap_thr) {
				return;
			}

			// checking if the match is visible in the local neighborhood
			double distance = std::numeric_limits<double>::infinity();
			if (t00 >= -PF_EPS && t00 <= 1 + PF_EPS)
				distance = std::min(distance, ((1 - t00) * pj0 + t00 * pj1 - pi0).norm());
			if (t01 >= -PF_EPS && t01 <= 1 + PF_EPS)
				distance = std::min(distance, ((1 - t01) * pj0 + t01 * pj1 - pi1).norm());

			if (t10 >= -PF_EPS && t10 <= 1 + PF_EPS)
				distance = std::min(distance, ((1 - t10) * pi0 + t10 * pi1 - pj0).norm());
			if (t11 >= -PF_EPS && t11 <= 1 + PF_EPS)
				distance = std::min(distance, ((1 - t11) * pi0 + t11 * pi1 - pj1).norm());

			if (distance > max(length_i, length_j) * .75) {
				return;
			}

			// check if the raster supports this parallelity
			// find the primary direction of the two edges (the axis with maximum extent)
			auto get_primary_direction = [](const vec2& d) {
				if (std::abs(std::abs(d.x()) - std::abs(d.y())) < PF_EPS)
					return -1; //45�
				else if (std::abs(d.x()) > std::abs(d.y()))
					return 0;
				else
					return 1;
			};

			int primary_i = get_primary_direction(di);
			int primary_j = get_primary_direction(dj);
			if (primary_i == -1)
				primary_i = primary_j;
			if (primary_j == -1)
				primary_j = primary_i;			
			int primary = primary_i == -1 ? 0 : primary_i;

			// now walk along the raster in the primary direction and count the 
			// steps

			// counts the maximum and minimum step in the primary direction
			auto count_max_step = [&](int v_from, int v_to)
			{
				vec2 maxSteps(0, 0);
				int last_step = v_from;
				int v = v_from;
				vec2 minmax(0, 0);
				while (v != v_to)
				{
					auto last_v = v;
					v = Circular(B, v + 1);
					if (std::abs((B.col(last_v) - B.col(v))(1 - primary)) > PF_EPS || v == v_to)
					{
						double current_steps = (B.col(v) - B.col(last_step))(primary);
						last_step = v;
						if (current_steps < minmax(0))
							minmax(0) = current_steps;
						if (current_steps > minmax(1))
							minmax(1) = current_steps;
					}
				}
				return minmax;
			};

			auto minmaxi = count_max_step(V(vi0), V(vi1));
			auto minmaxj = count_max_step(V(vj0), V(vj1));

			// The primary directions can be different if we are at a 45� line (1-1-1-1-1-....)
			if (primary_i != primary_j && (minmaxi.cwiseAbs().maxCoeff() > 1 || minmaxj.cwiseAbs().maxCoeff() > 1))
				return; //edges lie in different octants

			// one of the lines is in the opposite direction
			minmaxj = -minmaxj;
			std::swap(minmaxj(0), minmaxj(1));

			if (di(1 - primary) != 0 && dj(1 - primary) != 0) //if they are axis-aligned, they are definitely parallel
			{

				bool matching_steps = ((minmaxi - minmaxj).cwiseAbs().maxCoeff() <= 1. + PF_EPS);
				auto maxi = minmaxi.cwiseAbs().maxCoeff();
				auto maxj = minmaxj.cwiseAbs().maxCoeff();
				//if ((maxi == 1 && maxj > 1) || (maxj == 1 && maxi > 1))
				//	matching_steps = false;
				if(!matching_steps)
					return;
			}

			// Checking which of the pair of corners are aligned
			if (r.connected_00_11)
			{
				// Connecting line must be (almost) axis-aligned or at a 45� angle
				if (di0j1.cwiseAbs().minCoeff() > 1. + PF_EPS && std::abs(std::abs(di0j1.x()) - std::abs(di0j1.y())) > 1.0 + PF_EPS) {
					return;
				}
			}

			if (r.connected_01_10)
			{
				// Connecting line must be (almost) axis-aligned or at a 45� angle
				if (di1j0.cwiseAbs().minCoeff() > 1. + PF_EPS && std::abs(std::abs(di1j0.x()) - std::abs(di1j0.y())) > 1.0 + PF_EPS) {
					return;
				}
			}

			// sets r.aligned_00_11 and r.aligned_01_10
			test_and_set_opposite_points_alignment(B, C, V, r);

			reg.add(r);
		});
	});
}

// ----------------------------------------------------------------------------
void find_all_symmetric_points(
	const mat2x& B,
	const vecXi& V,
	const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries,
    const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries_local,
	RegularityInformation& reg,
	const bool circular
) {
	const double same_angle_threshold = PF_EPS;

	// find a mapping from boundary indices to polygon indices
	std::vector<int> boundary_to_polygon(B.cols(), -1);
	for (int i = 0; i < V.size(); ++i)
		boundary_to_polygon[V(i)] = i;

	for (size_t k = 0; k < raster_symmetries.size(); ++k) {
		const polyfit::Symmetry::SubPath& s = raster_symmetries[k];

		const Eigen::Index size = CircularDist(B, s.v01, s.v00);

		// walk along the two subregions in opposite directions
		Eigen::Index bi = s.v01;
		Eigen::Index bj = s.v11;

		int last_vi = -1;
		int last_vj = -1;

		while (true) {
			const int vi = boundary_to_polygon[bi];
			const int vj = boundary_to_polygon[bj];

			if (vi != -1 && vj != -1)
			{
				// there is a polygon vertex at both boundary positions

				Symmetry e;
				e.v0 = vi;
				e.v1 = vj;
				e.axis = s.axis;
				e.size = size;
				e.region = k;

				// check if this is a strong or weak symmetry
				const double angle_i = AngleUtils::spanned_shortest(B.col(CircularAt(V, vi - 1)), B.col(bi), B.col(CircularAt(V, vi + 1)));
				const double angle_j = AngleUtils::spanned_shortest(B.col(CircularAt(V, vj - 1)), B.col(bj), B.col(CircularAt(V, vj + 1)));

				const bool same_angle = abs(angle_i - angle_j) < same_angle_threshold;

				e.is_relaxed = !same_angle;	
				reg.add(e);

				// check edge symmetry
				if (last_vi != -1 && last_vj != -1) {
					Symmetry edge_symmetry = e;
					edge_symmetry.v0 = last_vi;
					edge_symmetry.v1 = vj;
					reg.add_edge_symmetry(edge_symmetry);
				}

				last_vi = vi;
				last_vj = vj;
			}

			// we are at the end of the region
			if (bi == s.v00)
				break;

			// advance in the regions
			bi = Circular(B, bi + 1);
			bj = Circular(B, bj - 1);
		}
	}

    for (size_t i = 0; i < raster_symmetries_local.size(); ++i) {
        const polyfit::Symmetry::SubPath& s = raster_symmetries_local[i];
        PF_VERBOSE_F("Local symmetry %d %d %d %d", s.v01, s.v00, s.v10, s.v11);
        const int poly_v0 = boundary_to_polygon[s.v00];
        const int poly_v1 = boundary_to_polygon[s.v10];

        if (poly_v0 == -1 || poly_v1 == -1) {
            continue;
        }

        Symmetry sym;
        sym.v0 = poly_v0;
        sym.v1 = poly_v1;
        sym.axis = s.axis;
        sym.size = CircularDist(B, s.v01, s.v11);
        sym.region = i;
		sym.is_relaxed = false;
        reg.add(sym);
    }
}

void find_all_important_edges(const mat2x & P, RegularityInformation& reg)
{
	auto is_cornerish = [&](int i)
	{
		auto next = polyfit::PathUtils::find_next_non_degenerate_corner<polyfit::PathUtils::DegenerateClassifier2Px>(P, i, polyfit::Circular(P, i + 1), true);
		auto prev = polyfit::PathUtils::find_next_non_degenerate_corner<polyfit::PathUtils::DegenerateClassifier2Px>(P, i, polyfit::Circular(P, i - 1), true);

		if (polyfit::AngleUtils::spanned_shortest(P.col(prev), P.col(i), P.col(next)) < PF_RAD(120))
			return true;

		return false;
	};

	// Calculate the bounding box of the polygon
	geom::aabb bounding_box;
	for (int i = 0; i < P.cols(); ++i) {
		bounding_box.add(P.col(i));
	}

	for (int i = 0; i < P.cols(); ++i) {
		auto& p = P.col(i);
		const auto& next = polyfit::CircularAt(P, i + 1);
		auto nextnext = polyfit::CircularAt(P, i + 2);
		auto prev = polyfit::CircularAt(P, i - 1);

		auto d = next - p;

		bool is_long_axis_aligned = false;
		is_long_axis_aligned |= (std::abs(d.x()) < PF_EPS && std::abs(d.y()) >= 0.5 * bounding_box.height());
		is_long_axis_aligned |= (std::abs(d.y()) < PF_EPS && std::abs(d.x()) >= 0.5 * bounding_box.width());
		is_long_axis_aligned |= std::abs(std::abs(d.x()) - std::abs(d.y())) <= 1 + PF_EPS && d.norm() >= 0.5 * std::min(bounding_box.width(), bounding_box.height());

		auto d_next = nextnext - next;
		auto d_prev = p - prev;

		bool is_short_diagonal =
			std::abs(std::abs(d.x()) - std::abs(d.y())) < PF_EPS
			&& d.cwiseAbs().sum() <= 4. + PF_EPS
			&& (d_next.cwiseAbs().minCoeff() < PF_EPS)
			&& (d_prev.cwiseAbs().minCoeff() < PF_EPS)
			&& (std::abs(d_next.x()) < PF_EPS) != (std::abs(d_prev.x()) < PF_EPS);

		if (is_long_axis_aligned)
			reg.add_important_edge(i);

		if (is_short_diagonal) {

			if (d_prev.cwiseAbs().maxCoeff() > d.cwiseAbs().sum() * 3 && is_cornerish(polyfit::Circular(P, i - 1)))
				reg.add_important_edge((int)polyfit::Circular(P, i - 1));

			if (d_next.cwiseAbs().maxCoeff() > d.cwiseAbs().sum() * 3 && is_cornerish(polyfit::Circular(P, i + 2)))
				reg.add_important_edge((int)polyfit::Circular(P, i + 1));
		}
	}
}

// ----------------------------------------------------------------------------
void RegularityInformation::update (
	const mat2x& B,        // Fitting boundary
	const vecXi& V,        // Polygonal fit on B
	const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries,	
    const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries_local, 
	const bool circular
) {
#if !POLYVEC_REGULARIZE
	return;
#endif

	PF_VERBOSE_S("construct regularity graph");

	mat2x P;
	BoundaryGraph::trace_to_points(B, V, P);	

    // Find parallel edges
	if (parallels_dirty)
	{
		std::vector<int> C;
		PathUtils::compute_convexities(P, C);

		clear_parallels();
		find_all_parallel_edges(B, V, C, *this, circular);
		parallels_dirty = false;
	}
	
    // Find symmetries	
	if (symmetries_dirty)
	{
		clear_symmetries();
		per_vertex_symmetries.resize(V.size());
		find_all_symmetric_points(B, V, raster_symmetries, raster_symmetries_local, *this, circular);
		symmetries_dirty = false;
	}

	if (important_edges_dirty)
	{
		clear_important_edges();
		find_all_important_edges(P, *this);
		important_edges_dirty = false;
	}

	if (per_edge_continuations_forward.size() != V.size())
		per_edge_continuations_forward.resize(V.size(), -1);
	if (per_edge_continuations_backward.size() != V.size())
		per_edge_continuations_backward.resize(V.size(), -1);

	polygon_vertices = V.size();
}

void RegularityInformation::add(const Parallel& p)
{
	_parallels.push_back(p);
}

void RegularityInformation::add(const Symmetry& s)
{
	_vertex_symmetries.push_back(s);
	per_vertex_symmetries[s.v0].push_back(_vertex_symmetries.size() - 1);
	if(s.v0 != s.v1)
		per_vertex_symmetries[s.v1].push_back(_vertex_symmetries.size() - 1);
}

void RegularityInformation::add(const Continuation& s)
{
	_continuations.push_back(s);

	if (s.v0 == Circular<size_t>(polygon_vertices, s.v0_prev + 1))
	{
		per_edge_continuations_forward[s.v0] = _continuations.size() - 1;
		per_edge_continuations_backward[s.v1] = _continuations.size() - 1;
	}
	else
	{
		per_edge_continuations_backward[s.v0] = _continuations.size() - 1;
		per_edge_continuations_forward[s.v1] = _continuations.size() - 1;
	}
}

const Continuation* RegularityInformation::get_edge_continuation(Vertex edge_src, Vertex edge_dst) const
{
	auto get_continuation = [this](int i) -> const Continuation*
	{
		if (i == -1)
			return nullptr;
		else
			return &_continuations[i];
	};

	if (edge_dst == Circular<size_t>(polygon_vertices, edge_src + 1))
		return get_continuation(per_edge_continuations_forward[edge_dst]);
	else if(edge_dst == Circular<size_t>(polygon_vertices, edge_src - 1))
		return get_continuation(per_edge_continuations_backward[edge_dst]);
	else throw std::runtime_error("The specified edge does not exist.");
}

void RegularityInformation::add_important_edge(int v0)
{
	ImportantEdge e;
	e.v0 = v0;
	_important_edges.push_back(e);
}

void RegularityInformation::add_edge_symmetry(const Symmetry & s)
{
	_edge_symmetries.push_back(s);
}

nse::util::IteratorRange<RegularityInformation::SymmetryIterator> RegularityInformation::vertex_symmetries(Vertex v) const
{
	return nse::util::IteratorRange<SymmetryIterator>(
		SymmetryIterator(*this, per_vertex_symmetries[v].begin()), 
		SymmetryIterator(*this, per_vertex_symmetries[v].end()));
}

void RegularityInformation::reset()
{
	clear_parallels();
	clear_continuations();
	clear_symmetries();
	clear_important_edges();	
}

void RegularityInformation::add_continuations(const mat2x & boundary, vecXi & polygon, std::vector<polyfit::BoundaryGraph::Edge>& E, bool circular, bool allow_polygon_modification)
{
	_continuations.clear();
	auto candidates = Regularity::find_continuation_candidates(boundary, polygon, *this, circular, allow_polygon_modification);
	Regularity::pick_continuations(boundary, polygon, candidates, *this, circular, E);
}

void RegularityInformation::add_continuations(const mat2x & boundary, vecXi & polygon, bool circular)
{
	std::vector<BoundaryGraph::Edge> unused;
	add_continuations(boundary, polygon, unused, circular, false);
}

void RegularityInformation::reindex_after_vertex_deletion(const std::vector<int>& deleted_vertices_ordered, int old_vertex_count)
{
	// i : Old index
	// allow_deletes_direction: 0 if re-indexing deleted vertices should throw an exception
	//                          otherwise, +1 or -1 to determine the direction of the newly corresponding vertex after deletion
	std::function<int(int, int)> new_index = [&](int i, int allow_deletes_direction = 0) -> int
	{
		// count how many deleted vertices are before i
		auto it = std::lower_bound(deleted_vertices_ordered.begin(), deleted_vertices_ordered.end(), i);
		if (it != deleted_vertices_ordered.end() && *it == i)
		{
			if(allow_deletes_direction == 0)				
				throw std::runtime_error("Trying to re-index a vertex that is deleted.");
			else
				return new_index(Circular(old_vertex_count, i + allow_deletes_direction), allow_deletes_direction);
		}
		return i - std::distance(deleted_vertices_ordered.begin(), it);
	};
	
	for (auto& r : _continuations)
		r.reindex(new_index);

	clear_parallels();
	clear_symmetries();
	clear_important_edges();
}


NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)