/*
	This shouldn't probably exist, I have no idea where to put it right now.
*/
#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/measure_accuracy.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

void measure_accuracy_signed_extended_polygon_fit_symmetric(
	const mat2x& B,       // raster boundary
	const vecXi& P,        // polygon
	const Vertex corner,      // corner
	AccuracyMeasurement& m,// result
	const bool circular
);

void measure_accuracy_signed_extended_polygon_fit_asymmetric(
	const mat2x& B,       // raster boundary
	const vecXi& P,        // polygon
	const Vertex corner,      // corner
	AccuracyMeasurement& m,// result
	const bool circular
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)
