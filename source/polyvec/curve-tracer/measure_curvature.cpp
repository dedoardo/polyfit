// polyvec
#include <polyvec/curve-tracer/measure_curvature.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/curve-tracer/spline.hpp>

// libc++
#include <cstdlib> // min, max

using namespace polyvec;
using namespace std;

// How densely should the curvature be evaluated with respect to the parametric t
#define K_EVAL_STEP_T .05

// If true, the curvature is only evaluated between [K_EVAL_RELAXED_SLACK_T 1-K_EVAL_RELAXED_SLACK_T]
#define K_EVAL_RELAXED_BOUNDS 1
#define K_EVAL_RELAXED_SLACK_T .125 // too big?

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

void measure_curvature_radius(
	polyvec::CurvePrimitive& prim,
	CurvatureMeasurement& m
) {
    measure_curvature_radius_2(*prim.curve->get_curve(), m);
}

void measure_curvature_radius_2(
    polyvec::GlobFitCurve& curve,
    CurvatureMeasurement& m
) {
    double t_eval_k = K_EVAL_STEP_T;

    while (t_eval_k < 1. + PF_EPS) {
#if K_EVAL_RELAXED_BOUNDS
        if (t_eval_k <= K_EVAL_RELAXED_SLACK_T || t_eval_k >= (1. - K_EVAL_RELAXED_SLACK_T)) {
            t_eval_k += K_EVAL_STEP_T;
            continue;
        }
#endif

        double k;
        SmoothCurveUtil::curvature_eval(curve.dposdt(t_eval_k), curve.dposdtdt(t_eval_k), k);
        const double r = 1. / abs(k);

        m.r_min = min(m.r_min, r);
        m.r_max = max(m.r_max, r);

        t_eval_k += K_EVAL_STEP_T;
    }

}

void CurvatureMeasurement::combine(const CurvatureMeasurement& other)
{
	r_min = std::min(r_min, other.r_min);
	r_max = std::max(r_max, other.r_max);
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)
