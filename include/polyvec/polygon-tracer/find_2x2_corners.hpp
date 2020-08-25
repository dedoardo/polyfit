/*
	Given a list of raster points, it finds two consecutive edges
	which are at least 2 pixels long
*/
#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void find_2x2_corners(
	const mat2x& B, // grid aligned points
	vecXi&       C, // corners
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)