#pragma once

#include <polyvec/api.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Finds the cycle of minimum cost for the boundary B, without any pruning

	in:
		R is a clockwise oriented raster polyline (grid spacing)

	out:
		V are the vertices of the polygonal approximation (indexes B)
		B is the fitting polyline to which the polygon refers to
		E are the list of edges which connecting other plausible polygon edges

	Returns 0 on failure
*/
int minimum(const mat2x& R, vecXi& P, mat2x& B, std::vector<BoundaryGraph::Edge>& E, const bool circular);

int minimum(
    ShortestPath::State& s, 
    const mat2x& B, 
    const std::vector<Vertex>& M, 
    const std::vector<BoundaryGraph::Edge>& E, 
    vecXi& P, 
    const bool circular,
    const std::vector<bool>& corners_c0 = std::vector<bool>(),
    const std::vector<bool>& B_flat = std::vector<bool>()
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
