#pragma once

#include <polyvec/curve-tracer/measure_accuracy.hpp>
#include <polyvec/curve-tracer/measure_curvature.hpp>

namespace polyvec {
	// CurveFitter
	// -------------------------------------------------------------------------
	struct CurveFitError {
		polyfit::CurveTracer::AccuracyMeasurement accuracy;
		polyfit::CurveTracer::CurvatureMeasurement curvature;

		Eigen::Vector2d failure;

		void combine(const CurveFitError& other);
		void reset();
	};
}