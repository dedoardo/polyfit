#pragma once

#include <polyvec/api.hpp>
#include <polyvec/curve-tracer/spline.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Visualize)

void primitives_indices(
	const std::vector<polyvec::CurvePrimitive>& prims
);

void primitives_alternated(
	const std::vector<polyvec::CurvePrimitive>& prims
);

// polygon
void primitives_alternated(
	const mat2x& P,
	bool circular 
);

void primitives_midpoints(
	const std::vector<polyvec::CurvePrimitive>& prims
);

void primitives_midpoints_normals(
	const std::vector<polyvec::CurvePrimitive>& prims
);

NAMESPACE_END(Visualize)
NAMESPACE_END(polyfit)