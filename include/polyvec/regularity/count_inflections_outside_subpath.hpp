#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

/*
	Given a polypath and a list of vertices, it counts the 
	number of edges connecting consecutive vertices with angles of
	opposite convexity.
	Only the inflections which are outside the specified closed interval
	are counted.

	in:
		P is a polygonal path
		V is the list of vertices defining the polygon edges
		B pair of indices in P which determines the closed interval
	
*/
int count_inflections_outside_subpath(
	const mat2x& P,
	const vecXi& V,
	const vec2i& B,
	bool circular
	);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)
