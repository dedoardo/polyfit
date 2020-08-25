#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Removes the edges from the graph in order for the shortest path to contain them.
	Pixel points along the boundary which form a 90 degree angle and have length >= 2
	edges are visually perceived as clear corners, though not necessarily sharp.

	The input points don't have to be grid aligned, but should only contain 90 degree
	turns, otherwise the behavior is undefined.

	Returns the number of edges which have been removed
*/
int remove_edges_crossing_implicit_corners(
	const mat2x& B,                      // Raster points
	std::vector<BoundaryGraph::Edge>& E, // Edges to be pruned
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)