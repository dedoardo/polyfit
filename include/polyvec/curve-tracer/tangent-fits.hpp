/*
	Enumerates the possible tangent interpolation types and contains
	the basic routines for initializing a curve given a raster, a polygon
	corner and it's ring-1 neighborhood
*/
#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

enum class TangentFit {
	ASYMMETRIC,
	SYMMETRIC,
	DISCONTINUOUS
};

void construct_asymmetric_bezier3(
	const mat2x& B,		  // raster boundary
	const vec3i& corners, // polygon corners
	mat24& C              // control points of the curve bezier
);

void construct_symmetric_bezier3(
	const mat2x& B,		  // raster boundary
	const vec3i& corners, // polygon corners
	mat24& C              // control points of the curve bezier
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)