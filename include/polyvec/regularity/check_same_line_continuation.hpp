#pragma once

#include <polyvec/core/types.hpp>
#include <polyvec/core/macros.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

/*
	Checks if two corners are part of the same line, the edges which are 
	would form this line are are P(i-1) P(i) P(j) P(j+1).
	All the vertices are assumed to exist (no need to check for circularity)

	The routine checks the distance from P(i) and P(j) to the segment 
	connecting P(i-1) P(j+1).

	The angles are P(i) and P(j) should have opposite convexity for the match
	to be considered valid.

	Returns true if the corners are a valid candidate and a weight in [0 1]
	range.
*/
bool check_same_line_continuation(
	const mat2x& B,
	const vecXi& P,
	const Vertex i,
	const Vertex j,
	double& weight
);

bool check_line_continuity(
	const vec2 pi0,
	const vec2 pi1,
	const vec2 pj0,
	const vec2 pj1,
	double& weight
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)