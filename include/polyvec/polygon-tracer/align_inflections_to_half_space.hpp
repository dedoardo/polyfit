#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Finds all inflections longer than two pixels and attempts to make
	the edge after/before axis-aligned if it's nearly axis aligned.
	This is only attempted if the outgoing direction lies in the opposite 
	half space defined by the normal with respect to the previous edge.
	(a) is a good candidate, while (b) is not

	     ____  ____
	    |          |
	    |n         |n
	____|      ____|
	(a)        (b)

	The operation is accepted if not new inflections are introduced.

	Depending on where the next point is place with respect to the candidate
	axis aligned edge, the operation can be a relocation or an insertion.
	Relocations are essentially done to prevent inserting degenerate edges

	Returns the number of edges which have been aligned.

	After all the edges have been aligned, the graph edges are pruned in
	order for successive shortest_paths() to not undo the changes.

	Ideally, the polygon should be recompute after calling this method, even
	if the polygon is guaranteed to not contain degenerate edges it will
	likely not be optimal with respect to the energy (especially smoothness)
*/
int align_inflections_to_half_space(
	const mat2x& B,                      // raster boundary (in)
	std::vector<BoundaryGraph::Edge>& E, // graph edges (in/out)
	vecXi& P,							 // polygon vertices (in/out)
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)