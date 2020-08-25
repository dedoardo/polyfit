#include <polyvec/curve-tracer/spline_validation.hpp>

polyvec::FitValidity::Type polyvec::FitValidity::evaluate_validity(const CurveFitError& error, const polyfit::CurveTracer::AccuracyMeasurement& polygonMeasure) {
	Type result = 0;

	if (error.accuracy.max_error() <= 0.51 || error.accuracy.max_error() <= polygonMeasure.max_error() + 0.05) {
		result |= PASSED_ACCURACY;
	}

	if (error.curvature.r_min >= 1.5) {
		result |= PASSED_CURVATURE;
	}

	return result;
};

std::vector<polyvec::CurveFitError> polyvec::get_corner_errors(size_t corners, const std::vector<CurvePrimitive>& primitives) {
	std::vector<CurveFitError> errors(corners);

	for (auto& p : primitives) {
		errors[p.corner].combine(p.error);
	}

	return errors;
}