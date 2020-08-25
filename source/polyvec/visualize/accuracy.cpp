#include <polyvec/visualize/accuracy.hpp>
#include <polyvec/curve-tracer/measure_accuracy.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/io/draw.hpp>
#include <polyvec/curve-tracer/curve_line.hpp>

#define ACCURACY_MEASURE_STYLE Style::outline(colors::red, 2.)
#define LABEL_FONT draw::font_pdf / 2.5, Style::text()
#define CURVE_COLOR_ALT0 colors::talking_orange
#define CURVE_COLOR_ALT1 colors::forest_green
#define CURVE_STROKE_WIDTH 5.

using namespace std;
using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Visualize)

#if 0

double measure_l1_signed_midpoints_with_slack(
	const vec2& p,
	const vec2& n,
	GlobFitCurve* c,
	bool draw = true) {

	if ((c->pos(c->project(p)) - p).norm() < 0.025) {
		return 0.;
	}

	double ray_at, curve_at;
	if (c->approximate_intersect_ray_closest(p, n, ray_at, curve_at) == 0 ||
		abs(ray_at) > .5) {

		const vec2 p_mid_off_out = p + n * .5;
		const vec2 p_mid_off_in = p - n * .5;
		const vec2 p_curve_out = c->pos(c->project(p_mid_off_out));
		const vec2 p_curve_in = c->pos(c->project(p_mid_off_in));
		const double e_out = (p_mid_off_out - p_curve_out).norm();
		const double e_in = (p_mid_off_in - p_curve_in).norm();

		if (e_out < e_in) {
			if (draw) {
				draw::line(p, p_mid_off_out, ACCURACY_MEASURE_STYLE);
				draw::line(p_mid_off_out, p_curve_out, ACCURACY_MEASURE_STYLE);
			}
			return .5 + e_out;
		}
		else {
			if (draw) {
				draw::line(p, p_mid_off_in, ACCURACY_MEASURE_STYLE);
				draw::line(p_mid_off_in, p_curve_in, ACCURACY_MEASURE_STYLE);
			}
			return -(.5 + e_in);
		}
	}
	else {
		if (draw) {
			draw::line(p, p + ray_at * n, ACCURACY_MEASURE_STYLE);
		}
		return ray_at;
	}
}

// raster -> curve l1 measurements
void accuracy_measure_l1_signed_midpoints_with_slack(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	accuracy_measure_unsigned_midpoints_with_slack(prims);

	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			continue;
		}

		CurveTracer::AccuracyMeasurement m;
		CurveTracer::measure_accuracy_signed_extended(&prim, m);
		const string label = StringUtils::fmt("%i, in: %f out: %.3f l1: %.3f l1-avg: %.3f l2: %.3f l2-avg: %.3f",
			i, m.e_in, m.e_out, m.e_l1, m.e_l1_normalized, m.e_l2, m.e_l2_normalized);
		draw::text(curve->get_curve()->pos(.35), label, LABEL_FONT);
	}
}

// raster -> curve l2 measurements
void accuracy_measure_l2_unsigned_midpoints(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			continue;
		}

		for (size_t j = 0; j < prim.fit_midpoints->size(); ++j) {
			const vec2 p_mid = prim.fit_midpoints->at(j);
			const vec2 p_curve = curve->get_curve()->pos(curve->get_curve()->project(p_mid));
			draw::line(p_mid, p_curve, ACCURACY_MEASURE_STYLE);
		}

		CurveTracer::AccuracyMeasurement m;
		CurveTracer::measure_accuracy_unsigned(&prim, m);
		const string label = StringUtils::fmt("max: %.3f l1: %.3f l1-avg: %.3f l2: %.3f l2-avg: %.3f",
			m.e_out, m.e_l1, m.e_l1_normalized, m.e_l2, m.e_l2_normalized);
		draw::text(curve->get_curve()->pos(.35), label, LABEL_FONT);
	}
}

void accuracy_measure_unsigned_midpoints_with_slack(
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve->get_curve().get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			continue;
		}

		for (size_t j = 0; j < prim.fit_midpoints->size(); ++j) {
			vec2 pmid = prim.fit_midpoints->at(j);
			vec2 pmid_n = prim.fit_midpoint_normals->at(j);
			vec2 pcurve = curve->pos(curve->project(pmid));

			const double d = (pcurve- pmid).norm();
			if (d >= .5) {
				const vec2 pmid_off_pos = pmid + .5 * pmid_n;
				const vec2 pmid_off_neg = pmid - .5 * pmid_n;
				const double d_pos = (pmid_off_pos - curve->pos(curve->project(pmid_off_pos))).norm();
				const double d_neg = (pmid_off_neg - curve->pos(curve->project(pmid_off_neg))).norm();

				if (d_pos < d_neg) {
					draw::line(pmid, pmid_off_pos, ACCURACY_MEASURE_STYLE);
					pmid = pmid_off_pos;
					pcurve = curve->pos(curve->project(pmid_off_pos));
				} else {
					draw::line(pmid, pmid_off_neg, ACCURACY_MEASURE_STYLE);
					pmid = pmid_off_neg;
					pcurve = curve->pos(curve->project(pmid_off_neg));
				}
			}

			draw::line(pmid, pcurve, ACCURACY_MEASURE_STYLE);
		}
	}
}

// raster -> polygon(curve) l1 measurements
void accuracy_measure_l1_signed_midpoints_with_slack_polygon(
	const mat2x& P, // polygon
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			continue;
		}

		if (prim.corner == -1 || prim.fit_midpoints == nullptr || prim.fit_midpoint_normals == nullptr) {
			continue;
		}

		const vec2 p0 = prim.endpoint_src;
		const vec2 p1 = P.col(prim.corner);
		const vec2 p2 = prim.endpoint_dst;

		mat2 L0, L1;
		L0.col(0) = p0;
		L0.col(1) = p1;
		L1.col(0) = p1;
		L1.col(1) = p2;

		GlobFitCurve_Line* l0 = new GlobFitCurve_Line;
		l0->set_points(L0);
		GlobFitCurve_Line* l1 = new GlobFitCurve_Line;
		l1->set_points(L1);

		double e_in_max = -INFINITY;
		double e_out_max = -INFINITY;
		double e_in_acc = 0.;
		double e_out_acc = 0.;
		double e_l1 = 0.;
		double e_l2 = 0.;

		for (size_t j = 0; j < prim.fit_midpoints->size(); ++j) {
			const vec2 p_mid = prim.fit_midpoints->at(j);
			const vec2 n_mid = prim.fit_midpoint_normals->at(j).normalized();

			const double d0 = (l0->pos(l0->project(p_mid)) - p_mid).norm();
			const double d1 = (l1->pos(l1->project(p_mid)) - p_mid).norm();

			double e;
			if (d0 < d1) {
				e = measure_l1_signed_midpoints_with_slack(p_mid, n_mid, l0);
			} else {
				e = measure_l1_signed_midpoints_with_slack(p_mid, n_mid, l1);
			}

			if (e < 0) {
				e_in_max = max(e_in_max, abs(e));
				e_in_acc += abs(e);
				e_l1 += abs(e);
				e_l2 += abs(e) * abs(e);
			}
			else {
				e_out_max = max(e_out_max, abs(e));
				e_out_acc += abs(e);
				e_l1 += abs(e);
				e_l2 += abs(e) * abs(e);
			}
		}

		if (e_in_max == -INFINITY) {
			e_in_max = 0.;
		}
		if (e_out_max == -INFINITY) {
			e_out_max = 0.;
		}

		double e_l1_avg = e_l1 / prim.fit_midpoints->size();
		double e_l2_avg = e_l2 / prim.fit_midpoints->size();
		draw::line(p0, p1, Style::outline(i % 2 ? CURVE_COLOR_ALT0 : CURVE_COLOR_ALT1, 2.));
		draw::line(p1, p2, Style::outline(i % 2 ? CURVE_COLOR_ALT0 : CURVE_COLOR_ALT1, 2.));

		string label = StringUtils::fmt("in-max %.3f out-max %.3f l1: %.3f l1-avg: %.3f l2: %.3f l2-avg: %.3f",
			e_in_max, e_out_max, e_l1, e_l1_avg, e_l2, e_l2_avg);
		draw::text(curve->get_curve()->pos(.3), label, LABEL_FONT);
	}
}

// raster -> polygon(curve) l2 measurements
void accuracy_measure_l2_unsigned_midpoints_polygon(
	const mat2x& P, // polygon
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			draw::line(prim.endpoint_src, prim.endpoint_dst, Style::outline(colors::black, CURVE_STROKE_WIDTH));
			continue;
		}

		if (prim.corner == -1 || prim.fit_midpoints == nullptr || prim.fit_midpoint_normals == nullptr) {
			continue;
		}

		const vec2 p0 = prim.endpoint_src;
		const vec2 p1 = P.col(prim.corner);
		const vec2 p2 = prim.endpoint_dst;

		mat2 L0, L1;
		L0.col(0) = p0;
		L0.col(1) = p1;
		L1.col(0) = p1;
		L1.col(1) = p2;

		GlobFitCurve_Line* l0 = new GlobFitCurve_Line;
		l0->set_points(L0);
		GlobFitCurve_Line* l1 = new GlobFitCurve_Line;
		l1->set_points(L1);

		double e_l1 = 0.;
		double e_l2 = 0.;

		double e_max = -INFINITY;

		for (size_t j = 0; j < prim.fit_midpoints->size(); ++j) {
			const vec2 p_mid = prim.fit_midpoints->at(j);
			const vec2 p_curve0 = l0->pos(l0->project(p_mid));
			const vec2 p_curve1 = l1->pos(l1->project(p_mid));

			double d0 = (p_mid - p_curve0).norm();
			double d1 = (p_mid - p_curve1).norm();
			if (d0 < d1) {
				e_max = max(e_max, d0);
				e_l1 += d0;
				e_l2 += d0 * d0;
				draw::line(p_mid, p_curve0, ACCURACY_MEASURE_STYLE);
			} else {
				e_max = max(e_max, d1);
				e_l1 += d1;
				e_l2 += d1 * d1;
				draw::line(p_mid, p_curve1, ACCURACY_MEASURE_STYLE);
			}
		}

		draw::line(p0, p1, Style::outline(i % 2 ? CURVE_COLOR_ALT0 : CURVE_COLOR_ALT1, 2.));
		draw::line(p1, p2, Style::outline(i % 2 ? CURVE_COLOR_ALT0 : CURVE_COLOR_ALT1, 2.));

		const double e_l1_avg = e_l1 / prim.fit_midpoints->size();
		const double e_l2_avg = e_l2 / prim.fit_midpoints->size();
		const string label = StringUtils::fmt("max: %.3f l1: %.3f l1-avg: %.3f l2: %.3f l2-avg: %.3f",
			e_max, e_l1, e_l1_avg, e_l2, e_l2_avg);
		draw::text(curve->get_curve()->pos(.3), label, LABEL_FONT);
	}
}

void accuracy_error_polygon_greater_l1(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			draw::line(prim.endpoint_src, prim.endpoint_dst, Style::outline(colors::black, CURVE_STROKE_WIDTH));
			continue;
		}

		if (prim.corner == -1 || prim.fit_midpoints == nullptr || prim.fit_midpoint_normals == nullptr) {
			continue;
		}

		const vec2 p0 = prim.endpoint_src;
		const vec2 p1 = P.col(prim.corner);
		const vec2 p2 = prim.endpoint_dst;

		mat2 L0, L1;
		L0.col(0) = p0;
		L0.col(1) = p1;
		L1.col(0) = p1;
		L1.col(1) = p2;

		GlobFitCurve_Line* l0 = new GlobFitCurve_Line;
		l0->set_points(L0);
		GlobFitCurve_Line* l1 = new GlobFitCurve_Line;
		l1->set_points(L1);

		double e_in_max = -INFINITY;
		double e_out_max = -INFINITY;
		double e_in_acc = 0.;
		double e_out_acc = 0.;

		for (size_t j = 0; j < prim.fit_midpoints->size(); ++j) {
			const vec2 p_mid = prim.fit_midpoints->at(j);
			const vec2 n_mid = prim.fit_midpoint_normals->at(j).normalized();

			const double d0 = (l0->pos(l0->project(p_mid)) - p_mid).norm();
			const double d1 = (l1->pos(l1->project(p_mid)) - p_mid).norm();

			double e;
			if (d0 < d1) {
				e = measure_l1_signed_midpoints_with_slack(p_mid, n_mid, l0, false);
			}
			else {
				e = measure_l1_signed_midpoints_with_slack(p_mid, n_mid, l1, false);
			}

			if (e < 0) {
				e_in_max = max(e_in_max, abs(e));
				e_in_acc += abs(e);
			}
			else {
				e_out_max = max(e_out_max, abs(e));
				e_out_acc += abs(e);
			}
		}

		if (e_in_max == -INFINITY) {
			e_in_max = 0.;
		}
		if (e_out_max == -INFINITY) {
			e_out_max = 0.;
		}

		double e_in_acc_norm = e_in_acc / prim.fit_midpoints->size();
		double e_out_acc_norm = e_out_acc / prim.fit_midpoints->size();

		CurveTracer::AccuracyMeasurement m;
		CurveTracer::measure_accuracy_signed_extended(&prim, m);

		double e_poly_max = max(e_out_max, e_in_max);
		double e_curve_max = max(m.e_out, m.e_in);

		vec3 color;
		if (e_curve_max - 0.025 > e_poly_max) {
			color = colors::red;
		}
		else {
			color = colors::calm_blue;
		}

		draw::curve(curve->get_curve().get(), CURVE_STROKE_WIDTH, color);
	}
}

void accuracy_error_polygon_greater_l2(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
) {
	for (size_t i = 0; i < prims.size(); ++i) {
		const auto& prim = prims[i];
		const auto& curve = prim.curve.get();

		if (!prim.fit_midpoints || !prim.fit_midpoint_normals) {
			draw::line(prim.endpoint_src, prim.endpoint_dst, Style::outline(colors::black, CURVE_STROKE_WIDTH));
			continue;
		}

		if (prim.corner == -1 || prim.fit_midpoints == nullptr || prim.fit_midpoint_normals == nullptr) {
			continue;
		}

		const vec2 p0 = prim.endpoint_src;
		const vec2 p1 = P.col(prim.corner);
		const vec2 p2 = prim.endpoint_dst;

		mat2 L0, L1;
		L0.col(0) = p0;
		L0.col(1) = p1;
		L1.col(0) = p1;
		L1.col(1) = p2;

		GlobFitCurve_Line* l0 = new GlobFitCurve_Line;
		l0->set_points(L0);
		GlobFitCurve_Line* l1 = new GlobFitCurve_Line;
		l1->set_points(L1);

		double e_acc = 0.;
		double e_max = -INFINITY;

		for (size_t j = 0; j < prim.fit_midpoints->size(); ++j) {
			const vec2 p_mid = prim.fit_midpoints->at(j);
			const vec2 p_curve0 = l0->pos(l0->project(p_mid));
			const vec2 p_curve1 = l1->pos(l1->project(p_mid));

			double d0 = (p_mid - p_curve0).norm();
			double d1 = (p_mid - p_curve1).norm();
			if (d0 < d1) {
				e_max = max(e_max, d0);
				e_acc += d0;
			}
			else {
				e_max = max(e_max, d1);
				e_acc += d1;
			}
		}

		CurveTracer::AccuracyMeasurement m;
		CurveTracer::measure_accuracy_unsigned(&prim, m);

		vec3 color;
		if (max(m.e_in, m.e_out) - .025 > e_max) {
			color = colors::red;
		}
		else {
			color = colors::calm_blue;
		}

		draw::curve(curve->get_curve().get(), CURVE_STROKE_WIDTH, color);
	}
}

void accuracy_error_polygon_greater(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
) {

}

#endif

NAMESPACE_END(Visualize)
NAMESPACE_END(polyfit)