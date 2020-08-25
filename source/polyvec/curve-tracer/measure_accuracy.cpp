// polyvec
#include <polyvec/curve-tracer/measure_accuracy.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/core/macros.hpp>
#include <polyvec/utils/directions.hpp>

#include <algorithm>

using namespace polyvec;
using namespace std;

// this shouldn't be changed unless a pixel doesn't have a length 1 edge
#define FALLBACK_PIXEL_CENTER (.5 + PF_EPS_LARGE)

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

void AccuracyMeasurement::combine(const AccuracyMeasurement& m) {
	e_pos = max(e_pos, m.e_pos);
	e_neg = max(e_neg, m.e_neg);
	e_l1 += m.e_l1;
	e_l2_sq += m.e_l2_sq;
	count += m.count;
}

void AccuracyMeasurement::reset() {
	new (this) AccuracyMeasurement;
	e_neg = -INFINITY;
	e_pos = -INFINITY;
}

void AccuracyMeasurement::add_pos(const double d) {
	PF_ASSERT(d >= 0.);
	e_pos = std::max(e_pos, d);
	e_l1 += d;
	e_l2_sq += d * d;
	++count;
}

void AccuracyMeasurement::add_neg(const double d) {
	//PF_ASSERT(d >= 0.);
	e_neg = std::max(e_neg, d);
	e_l1 += d;
	e_l2_sq += d * d;
	++count;
}

void AccuracyMeasurement::add(const double d) {
	PF_ASSERT(d >= 0.);
	e_neg = std::max(e_neg, d);
	e_pos = std::max(e_pos, d);
	e_l1 += d;
	e_l2_sq += d * d;
	++count;
}

void AccuracyMeasurement::finish() {
	if (e_neg == -INFINITY) {
		e_neg = 0.;
	}

	if (e_pos == -INFINITY) {
		e_pos = 0.;
	}
}

double AccuracyMeasurement::max_error() const {
	return std::max(e_neg, e_pos);
}

double AccuracyMeasurement::min_error() const {
	return std::min(e_neg, e_pos);
}

//Measures the accuracy of a single point in a given direction (keeps the orientation). Can return nan
double measure_accuracy(polyvec::GlobFitCurve & curve, const Eigen::Vector2d& p, const Eigen::Vector2d& n, const double normal_sign)
{
	auto n_norm = n.normalized();
	auto pixel_center = p + 0.5 * n_norm;

	//classify pixel center with respect to curve
	auto t = curve.project(pixel_center);

	if (t == 0 || t == 1)
	{
		//pixel center projects onto endpoint - check if this is stable
		double epsilon = 0.05;
		auto tPlus = curve.project(pixel_center + epsilon * n);
		auto tMinus = curve.project(pixel_center - epsilon * n);
		if (tPlus == t && tMinus == t)
			return std::numeric_limits<double>::quiet_NaN(); //we cannot classify this point
	}

	auto p_proj = curve.pos(t);
	auto t_proj = curve.dposdt(t);
	Eigen::Vector2d curve_normal = normal_sign * polyvec::util::normal_dir(t_proj);

	double classify_center_vs_curve = (pixel_center - p_proj).dot(curve_normal);

	if (classify_center_vs_curve == 0.)
		return 0.5; //pixel center lies on curve

	if (classify_center_vs_curve < 0)
		return (p_proj - pixel_center).norm() + 0.5; //pixel center lies between curve and raster

	//pixel center lies on the same side w.r.t. raster and w.r.t. curve
	//shoot ray
	double ray_t, curve_t;
	int intersections = curve.approximate_intersect_ray_closest(p, n_norm, ray_t, curve_t);
	if(intersections ==  0)
		//move a bit to the side
		intersections = curve.approximate_intersect_ray_closest(p + 0.05 * polyvec::util::normal_dir(n_norm), n_norm, ray_t, curve_t);
	if (intersections == 0)
		//move a bit to the other side
		intersections = curve.approximate_intersect_ray_closest(p - 0.05 * polyvec::util::normal_dir(n_norm), n_norm, ray_t, curve_t);
	if (intersections == 0 || ray_t > 0.5 || ray_t < 0)
		return std::numeric_limits<double>::quiet_NaN();
	else
		return ray_t;
}

AccuracyMeasurement measure_accuracy(polyvec::GlobFitCurve & curve, const mat2x & P, const mat2x & N)
{
	AccuracyMeasurement result;
	result.reset();

	PF_ASSERT(P.cols() == N.cols());

	for (size_t i = 0; i < P.cols(); ++i)
	{
		auto measure_pos = measure_accuracy(curve, Eigen::Vector2d(P.col(i)), Eigen::Vector2d(N.col(i)), 1.0);
		auto measure_neg = measure_accuracy(curve, Eigen::Vector2d(P.col(i)), Eigen::Vector2d(-N.col(i)), -1.0);

		if (std::isnan(measure_pos) && std::isnan(measure_neg))
		{
			//both measures are invalid, use direct distance
			auto proj_dir = curve.pos(curve.project(P.col(i))) - P.col(i);
			if (proj_dir.dot(N.col(i)) >= 0)
				result.add_pos(proj_dir.norm());
			else
				result.add_neg(proj_dir.norm());
		}
		else if (!std::isnan(measure_pos) && !std::isnan(measure_neg))
		{
			//both measures are valid, take the minimum
			if (measure_pos < measure_neg)
				result.add_pos(measure_pos);
			else
				result.add_neg(measure_neg);
		}
		else if (!std::isnan(measure_pos))
			result.add_pos(measure_pos);
		else
			result.add_neg(measure_neg);
	}

	result.finish();
	return result;
}

void measure_accuracy_signed_extended(
	polyvec::GlobFitCurve& curve,  // primitive tested
	const mat2x& P,                      // fitting points
	const mat2x& N,                      // fitting point normals
	AccuracyMeasurement& m               // results
) {
	PF_ASSERT(P.cols() == N.cols());

	// Constructing two lines that extend the current curve endpoints along the same 
	// tangent direction. This is to make the test more robust to midpoints exactly at 
	// the endpoints, but outside an epsilon range.
	// A valid alternative would be to move the first and last midpoint closer to the 
	// previous, but it can happen to have single lines which fit a single midpoint and
	// for which this is not necessarily well defined
	GlobFitCurve_Line line_beg_ext, line_end_ext;
	const double line_ext_len = 10.;
	mat2 line_pts;
	
	line_pts.col(0) = curve.pos(0.);
	line_pts.col(1) = line_pts.col(0) - curve.dposdt(0.) * line_ext_len;
	line_beg_ext.set_points(line_pts);

	line_pts.col(0) = curve.pos(1.);
	line_pts.col(1) = line_pts.col(0) + curve.dposdt(1.) * line_ext_len;
	line_end_ext.set_points(line_pts);

	m.reset();
	for (size_t i = 0; i < P.cols(); ++i) {
		const vec2& ray_o = P.col(i);
		const vec2& ray_d = N.col(i).normalized(); // todo assume... 
	
		// falling back to l2 from offset midpoint
		double ray_at, curve_at;
		if (curve.approximate_intersect_ray_closest(ray_o, ray_d, ray_at, curve_at) == 0 ||
			abs(ray_at) > FALLBACK_PIXEL_CENTER) {
			
			// if this is the first or last point, let's attempt to intersect the respective line
			// we need to expli
			if (i == 0 && 
				line_beg_ext.approximate_intersect_ray_closest(ray_o, ray_d, ray_at, curve_at) &&
				abs(ray_at) < FALLBACK_PIXEL_CENTER) {

				if (ray_at > 0.) {
					m.add_pos(ray_at);
				} else {
					m.add_neg(abs(ray_at));
				}

				continue;
			}

			if (i == P.cols() - 1 &&
				line_end_ext.approximate_intersect_ray_closest(ray_o, ray_d, ray_at, curve_at) && 
				abs(ray_at) < FALLBACK_PIXEL_CENTER) {

				if (ray_at > 0.) {
					m.add_pos(ray_at);
				} else {
					m.add_neg(abs(ray_at));
				}

				continue;
			}

			// otherwise attempting the test from the respective pixel centers
			const vec2 p_mid_off_out = ray_o + ray_d * .5;
			const vec2 p_mid_off_in = ray_o - ray_d * .5;
			const double e_pos = .5 + (curve.pos(curve.project(p_mid_off_out)) - p_mid_off_out).norm();
			const double e_neg = .5 + (curve.pos(curve.project(p_mid_off_in)) - p_mid_off_in).norm();
			const double e = std::min(e_pos, e_neg);

			if (e_pos < e_neg) {
				m.add_pos(e);
			} else {
				m.add_neg(e);
			}
		}
		else {
			if (ray_at > 0.) {
				m.add_pos(ray_at);
			}
			else {
				m.add_neg(abs(ray_at));
			}
		}
	}

	m.finish();
	if (m.e_neg == -INFINITY) {
		m.e_neg = 0.;
	}

	if (m.e_pos == -INFINITY) {
		m.e_pos = 0.;
	}
}

void measure_accuracy_unsigned(const polyvec::CurvePrimitive* prim, AccuracyMeasurement& m) {
	PF_ASSERT(!prim->fitting_info.fit_midpoints.empty());

	new (&m) AccuracyMeasurement;
	m.e_neg = -INFINITY;
	m.e_pos = -INFINITY;

	const auto& curve = prim->curve.get();

	for (Vertex i = 0; i < (Vertex)prim->fitting_info.fit_midpoints.size(); ++i) {
		const vec2& p_mid = prim->fitting_info.fit_midpoints.at(i);

		const double e = (curve->get_curve()->pos(curve->get_curve()->project(p_mid)) - p_mid).norm();
		m.e_pos = std::max(m.e_pos, e);
		m.e_l1 += e;
		m.e_l2_sq += e * e;
		++m.count;
	}
	
	m.e_neg = m.e_pos;
}

void measure_accuracy_unsigned_extended(const polyvec::CurvePrimitive* prim, AccuracyMeasurement& m) {
	 PF_ASSERT(prim && prim->curve)
	const auto& curve = prim->curve->get_curve().get();

	new (&m) AccuracyMeasurement;
	m.e_neg = -INFINITY;
	m.e_pos = -INFINITY;

	const double d_thr = .5;

	for (int i = 0; i < prim->fitting_info.fit_midpoints.size(); ++i) {
		// attempt 1. distance from midpoint to closest point on the curve
		const vec2& pmid = prim->fitting_info.fit_midpoints.at(i);
		double d = (curve->pos(curve->project(pmid)) - pmid).norm();

		if (d > d_thr) {
			const vec2& pmid_off_pos = pmid + prim->fitting_info.fit_midpoint_normals.at(i) * .5;
			const vec2& pmid_off_neg = pmid - prim->fitting_info.fit_midpoint_normals.at(i) * .5;

			const double d_pos = (curve->pos(curve->project(pmid_off_pos)) - pmid_off_pos).norm();
			const double d_neg = (curve->pos(curve->project(pmid_off_neg)) - pmid_off_neg).norm();
			d = d_thr + min(d_pos, d_neg);
		}

		m.e_pos = max(m.e_pos, d);
		m.e_l1 += d;
		m.e_l2_sq += d * d;
		++m.count;
	}

	m.e_neg = m.e_pos;

	if (m.e_neg == -INFINITY) {
		m.e_neg = 0.;
	}

	if (m.e_pos == -INFINITY) {
		m.e_pos = 0.;
	}
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)