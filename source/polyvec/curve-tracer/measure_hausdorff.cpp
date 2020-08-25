// Polyvec
#include <polyvec/curve-tracer/measure_hausdorff.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/core/constants.hpp>

// libc++
#include <algorithm>

using namespace polyfit;
using namespace std;

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(CurveTracer)

void HausdorffDistance::combine(const HausdorffDistance& other) {
    l2_sq += other.l2_sq;
    l2_sq_min = min(l2_sq_min, other.l2_sq_min);
    l2_sq_max = max(l2_sq_max, other.l2_sq_max);
    count += other.count;
}

HausdorffDistance measure_hausdorff(
    GlobFitCurve* c0,
    GlobFitCurve* c1,
    const double t0_min, const double t0_max,
    const double t1_min, const double t1_max
) {
    PF_ASSERT(c0);
    PF_ASSERT(c1);

    HausdorffDistance measure;

    double t = t0_min;

    double l2_sq;
    vec2   p_t;

    while (t < t0_max + PF_EPS) {
        p_t = c0->pos(max(min(t, 1.), 0.));
        l2_sq = (p_t - c1->pos(c1->project(p_t))).squaredNorm();
        measure.l2_sq += l2_sq;
        measure.l2_sq_max = max(measure.l2_sq_max, l2_sq);
        measure.l2_sq_min = min(measure.l2_sq_min, l2_sq);

        ++measure.count;
        t += PV_MEASURE_HAUSDORFF_SAMPLING_RATE;
    }

    t = t1_min;

    while (t < t1_max + PF_EPS) {
        p_t = c1->pos(max(min(t, 1.), 0.));
        l2_sq = (p_t - c0->pos(c0->project(p_t))).squaredNorm();
        measure.l2_sq += l2_sq;
        measure.l2_sq_max = max(measure.l2_sq_max, l2_sq);
        measure.l2_sq_min = min(measure.l2_sq_min, l2_sq);

        ++measure.count;
        t += PV_MEASURE_HAUSDORFF_SAMPLING_RATE;
    }

    return measure;
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvec)