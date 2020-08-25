#pragma once

// polyvec
#include <polyvec/api.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

#include <polyvec/regularity/regularity_action.hpp>

// libc++


NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Finds the shortest cycle for the boundary B on its fully connected graph
	Detects axis-aligned edges and prunes the graph to force them in future iterations
	For each locally symmetric symmetric region, it regularizes the path only if
	the shortest path on the updated graph is not degenerate.

	in:
		R is a clockwise oriented polyline where points are grid-aligned

	out:
		P is written to and will contain the polygon trace
		B is the list of points to which P refers to and it can be a superset of R
		E is the list of edges connecting points in B which are also plausible candidates to be polygon edges

	Returns the number of iterations
*/
//int iterative_local(const mat2x& R, 
//	vecXi& P, 
//	mat2x& B, 
//	std::vector<BoundaryGraph::Edge>& E, 
//	const bool circular,
//	IterationCallback on_iter = nullptr);

//int iterative_local(
//	ShortestPath::State& s, 
//	const mat2x& B, 
//	const std::vector<Vertex>& M, 
//	std::vector<BoundaryGraph::Edge>& E, 
//	vecXi& P, 
//	const bool circular,
//	IterationCallback on_iter = nullptr,
//    const std::vector<bool>& corners_C0 = std::vector<bool>()
//);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
