#include <polyvec/polygon-tracer/error-metrics.hpp>

// polyfit
#include <polyvec/api.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/core/options.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/core/options.hpp>

// c++ stl
#include <algorithm>

// Eigen
#include <Eigen/Geometry>

#define TEST_PARAM_DEFAULT(v) (v)
#define TEST_PARAM_20_UP(v) ((v) + (.1 * (v)))
#define TEST_PARAM_20_DOWN(v) ((v) - (.1 * (v)))

// 1.
#define SMOOTHNESS_LIMIT TEST_PARAM_DEFAULT(VectorOptions::get()->error_smoothness_limit)

// 2.
#define CONTINUITY_LIMIT TEST_PARAM_DEFAULT(VectorOptions::get()->error_continuity_limit)

// 3.
#define INFLECTION_LIMIT TEST_PARAM_DEFAULT(VectorOptions::get()->error_inflection_limit)

// 4.
#define INFLECTION_PENALTY TEST_PARAM_DEFAULT(VectorOptions::get()->error_inflection_penalty)

// unit weight
#define SMOOTHNESS_WEIGHT VectorOptions::get()->error_smoothness_weight
#define ACCURACY_WEIGHT VectorOptions::get()->error_accuracy_weight

// 5.
#define CONTINUITY_WEIGHT TEST_PARAM_DEFAULT(VectorOptions::get()->error_continuity_weight)

using namespace polyvec;
using namespace std;

namespace polyfit {
	namespace ErrorMetrics {
		double smoothness(const double rad0, const double rad1, const int flags) {
			const bool is_first = flags & EDGE_FIRST_IN_PATH;
			const bool is_last = flags & EDGE_LAST_IN_PATH;

			double e = 0.;
			e += .5 * (!(is_first) ? smoothness(rad0) : 0.);
			e += .5 * (!(is_last) ?  smoothness(rad1) : 0.);
			return e;
		}

		double smoothness(const double rad0) {
			return SMOOTHNESS_WEIGHT * min((M_PI - rad0) / SMOOTHNESS_LIMIT, 1.0);
		}

		double accuracy(const double d_min, const double d_max, const int flags) {
			const double overflow = max(d_max - 1., 0.) + max(d_min - .5, 0.);
			return ACCURACY_WEIGHT * d_min + max(d_max - .5, 0.) + overflow;
		}

		double accuracy_implicit(const double d_min, const double d_max, const int flags) {
			return ACCURACY_WEIGHT * 2 * max(0., d_max - .5);
		}

		double continuity(const double rad0, const double rad1, const int flags) {
			if ((flags & EDGE_LAST_IN_PATH) || (flags & EDGE_FIRST_IN_PATH)) {
				return 0.;
			}

			if (flags & EDGE_HAS_INFLECTION) {
                return CONTINUITY_WEIGHT * min(1., INFLECTION_PENALTY + (M_PI - min(rad0, rad1)) / INFLECTION_LIMIT  * (1. - INFLECTION_PENALTY));
                //const double a = min(2 * M_PI - rad0 - rad1, CONTINUITY_LIMIT);
                //return CONTINUITY_WEIGHT * a / CONTINUITY_LIMIT;
			}
			else {
				return CONTINUITY_WEIGHT * min(abs(rad0 - rad1) / CONTINUITY_LIMIT, 1.0);
			}
		}

		bool accuracy_within_bounds(const double d_min, const double d_max, const vec2& edge_dir) {
			// This is only checking the distance measures, not the error

			const double phi = atan2(edge_dir.cwiseAbs().minCoeff(), edge_dir.cwiseAbs().maxCoeff());
			const double d_eps = .1;
			// adaptive accuracy
			const double max_threshold = std::min(1.0, 0.5 + (.5 + d_eps) * cos(phi) - .5 * sin(phi));

			return (d_min <= .5 + PF_EPS && d_max <= max_threshold - PF_EPS);

			//return (d_min <= (.5 + PF_EPS) && d_max <= (1. - PF_EPS));
		}

		bool accuracy_within_bounds_relaxed(const double d_min, const double d_max) {
			return d_min < 1. - PF_EPS && d_max < 1. - PF_EPS;
		}

		vec3 calculate_at_edge(const mat2x& points, const mat24& p, const vec4i& v, int flags) {
#if 0
			assert_break(p.cols() == 4);

			const double rad0 = PathUtils::shortest_angle_spanned(p.col(0), p.col(1), p.col(2));
			const double rad1 = PathUtils::shortest_angle_spanned(p.col(1), p.col(2), p.col(3));
			flags |= PathUtils::are_consecutive_angles_opposite(p.col(0), p.col(1), p.col(2), p.col(3)) ? EDGE_HAS_INFLECTION : 0x0;

#if 1
			dbg::info(FMT("error-metric - angles:(%.3f %.3f) inflection: %d", geom::degrees(rad0), geom::degrees(rad1), flags & EDGE_HAS_INFLECTION));
#endif

			const vec2 d_err = PathUtils::distance_bounds_from_points(points, p.block(0, 1, 2, 2), v.segment(1, 2), );

			vec3 error;
			error(0) = smoothness(rad0, rad1, flags);
			error(1) = accuracy(d_err.x(), d_err.y(), flags);
			error(2) = continuity(rad0, rad1, flags);
			return error;
#endif
			return vec3(0., 0., 0.);
		}

		vec3 calculate_inner_sep(double d_min, double d_max, const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3) {
			const int flags = AngleUtils::have_opposite_convexity(p0, p1, p2, p3) ? EDGE_HAS_INFLECTION : 0x0;
			const double rad0 = AngleUtils::spanned_shortest(p0, p1, p2);
			const double rad1 = AngleUtils::spanned_shortest(p1, p2, p3);
			vec3 e;
			e(SMOOTHNESS) = ErrorMetrics::smoothness(rad0, rad1, flags);
			e(CONTINUITY) = ErrorMetrics::continuity(rad0, rad1, flags);

			if (VectorOptions::get()->accuracy_metric == AccuracyMetric::OneSided) {
				e(ACCURACY) = ErrorMetrics::accuracy(d_min, d_max, flags);
			}
			else if (VectorOptions::get()->accuracy_metric == AccuracyMetric::Implicit) {
				e(ACCURACY) = ErrorMetrics::accuracy_implicit(d_min, d_max, flags);
			}
			else {
				PF_ABORT;
			}
			return e;
		}

		vec3 calculate_first_sep(double d_min, double d_max, const vec2 p1, const vec2 p2, const vec2 p3) {
			const int flags = EDGE_FIRST_IN_PATH;
			const double rad1 = AngleUtils::spanned_shortest(p1, p2, p3);
			vec3 e;
			e(SMOOTHNESS) = ErrorMetrics::smoothness(M_PI, rad1, flags);
			e(CONTINUITY) = ErrorMetrics::continuity(M_PI, rad1, flags);
			if (VectorOptions::get()->accuracy_metric == AccuracyMetric::OneSided) {
				e(ACCURACY) = ErrorMetrics::accuracy(d_min, d_max, flags);
			}
			else if (VectorOptions::get()->accuracy_metric == AccuracyMetric::Implicit) {
				e(ACCURACY) = ErrorMetrics::accuracy_implicit(d_min, d_max, flags);
			}
			else {
				PF_ABORT;
			}
			return e;
		}

		vec3 calculate_last_sep(double d_min, double d_max, const vec2 p1, const vec2 p2, const vec2 p3) {
			const int flags = EDGE_FIRST_IN_PATH;
			const double rad1 = AngleUtils::spanned_shortest(p1, p2, p3);
			vec3 e;
			e(SMOOTHNESS) = ErrorMetrics::smoothness(M_PI, rad1, flags);
			e(CONTINUITY) = ErrorMetrics::continuity(M_PI, rad1, flags);
			if (VectorOptions::get()->accuracy_metric == AccuracyMetric::OneSided) {
				e(ACCURACY) = ErrorMetrics::accuracy(d_min, d_max, flags);
			}
			else if (VectorOptions::get()->accuracy_metric == AccuracyMetric::Implicit) {
				e(ACCURACY) = ErrorMetrics::accuracy_implicit(d_min, d_max, flags);
			}
			else {
				PF_ABORT;
			}
			return e;
		}

		double calculate_inner(double d_min, double d_max, const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3) {
			const vec3 e = calculate_inner_sep(d_min, d_max, p0, p1, p2, p3);
			return e.sum();
		}

		double calculate_first(double d_min, double d_max, const vec2 p1, const vec2 p2, const vec2 p3) {
			const vec3 e = calculate_first_sep(d_min, d_max, p1, p2, p3);
			return e.sum();
		}

		double calculate_last(double d_min, double d_max, const vec2 p1, const vec2 p2, const vec2 p3) {
			const vec3 e = calculate_last_sep(d_min, d_max, p1, p2, p3);
			return e.sum();
		}
	}
}