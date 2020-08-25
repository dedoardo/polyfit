#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Given a raster boundary B (including midpoints), it finds constant
	monotone sequences which correspond to perfectly rasterized lines and
	stores the indices to all points which are *NOT* plausible candidates to
	be the original sources/destination.

	For each point the distances and convexities of the two neighboring points
	in each direction are recorded. The points which fall in one of this two 
	cases are considered inner points of the line sequence:

	Only non-flat neighboring points are considered

	i) convexities: -1 +1 -1 +1 -1
	   distances:    N  1  0  N  1

	ii) convexities: +1 -1 +1 -1 +1
	    distances:    1  N  0  1  N

	.... (check the implementation)

	The number of regions is increased every time a new point is added which 
	is not the immediate on-flat neighbor of the previous one in the sequence.
	todo: Currently the count is off by one as circularity is not handled, but
	this measure is only used indicatively.

	Returns the number of lines (not points) which have been detected.
*/
int detect_perfect_lines_inner_points(
	const mat2x& B,
	vecXi&		 L,
	const bool   circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)