#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

namespace polyvec {
	struct CurvePrimitive;
    class GlobFitCurve;
}

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

struct CurvatureMeasurement {
	// r = 1 / |k|
	double r_min = INFINITY;
	double r_max = INFINITY;

	void combine(const CurvatureMeasurement& other);
};

void measure_curvature_radius(
	polyvec::CurvePrimitive& prim,
	CurvatureMeasurement& m
);

void measure_curvature_radius_2 (
    polyvec::GlobFitCurve& curve,
    CurvatureMeasurement& m
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)