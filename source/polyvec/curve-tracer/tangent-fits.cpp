// polyvec
#include <polyvec/curve-tracer/tangent-fits.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

void construct_asymmetric_bezier3(
	const mat2x& B,		  // raster boundary
	const vec3i& corners, // polygon corners
	mat24& C              // control points of the curve bezier
	) {

}

void construct_symmetric_bezier3(
	const mat2x& B,		  // raster boundary
	const vec3i& corners, // polygon corners
	mat24& C              // control points of the curve bezier
	) {

}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)