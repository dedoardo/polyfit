/*
    Measures the Hausdorff distance metric between two curves
*/
#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>

// math
#include <cmath>

NAMESPACE_BEGIN(polyvec)
class GlobFitCurve;

NAMESPACE_BEGIN(CurveTracer)

struct HausdorffDistance {
    double l2_sq = 0.;            // sum of squared distances
    double l2_sq_max = -INFINITY; // maximum squared distance
    double l2_sq_min = INFINITY;  // minimum squared distance
    int    count = 0;             // number of samples taken


    void combine( const HausdorffDistance& other );
};

HausdorffDistance measure_hausdorff (
    GlobFitCurve* c0,
    GlobFitCurve* c1, 
    const double t0_min = 0., const double t0_max = 1.,
    const double t1_min = 0., const double t1_max = 1.
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvec)