#pragma once

#include <polyvec/core/macros.hpp>

NAMESPACE_BEGIN(polyvec)

// Forward declaration
class BezierCurve;

NAMESPACE_BEGIN(CurveTracer)

void minimize_curvature_variation_coordinate_descent (
    const double p1x, const double p1y,
    const double t1x, const double t1y,
    const double p2x, const double p2y,
    const double t2x, const double t2y,
    double& t1_scale, double& t2_scale
);

void minimize_curvature_variation_coordinate_descent(
    polyvec::BezierCurve& bezier
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvec)