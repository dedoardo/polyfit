#pragma once

#include <polyvec/api.hpp>
#include <polyvec/curve-tracer/spline_types.hpp>

NAMESPACE_BEGIN(polyvec)


namespace FitValidity
{
	using Type = int;

	const Type FAILED_ALL = 0;
	const Type PASSED_ACCURACY = 1;
	const Type PASSED_CURVATURE = 2;
	const Type PASSED_ALL = PASSED_ACCURACY | PASSED_CURVATURE;

	Type evaluate_validity(const CurveFitError& error, const polyfit::CurveTracer::AccuracyMeasurement& polygonMeasure);
}

// Consolidates the per-primitive fitting errors into per-polygon corner fitting errors.
std::vector<CurveFitError> get_corner_errors(size_t corners, const std::vector<CurvePrimitive>& primitives);

NAMESPACE_END(polyvec)
