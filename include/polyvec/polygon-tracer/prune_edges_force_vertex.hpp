#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Prunes the list of edges in such a way that the successive shortest path will go through the corner
	the passability of the corner is preserved.

	All the edges which contain the vertex in question without overlapping it are marked for deletion.

	Returns the number of edges which have been removed
*/
int prune_edges_force_vertex(
	const mat2x& B,						 // raster boundary
	std::vector<BoundaryGraph::Edge>& E, // graph to be pruned
	const Eigen::Index v,				 // vertex the path will go through
	const bool circular				 
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)