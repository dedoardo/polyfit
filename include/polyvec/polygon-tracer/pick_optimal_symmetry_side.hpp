/*
	Given two symmetric raster regions which are approximated by two polygonal segments,
	it checks if they are already symmetric and in case which one of the two sides 
	is best to copy onto the other.

	The optimality of the side first uses the length of the two paths and if 
	that fails then it uses the same error metrics used during the shortest path.

	It also returns the ordered list of points for each of the symmetric regions.
*/
#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

// Returns
//   	0 -> already symmetric
//   	-1 -> symmetry 0 is optimal
//   	+1 -> symmetry 1 is optimal
int pick_optimal_symmetry_side(
	const mat2x& B,          // raster boundary
	const vecXi& P,          // polygonal approximation
	const std::vector<BoundaryGraph::Edge>& E, // edges
	const std::vector<Vertex> PE, // polygon edges
	const vec2i& sym0,       // (first, last) points in the raster for symmetry 0
	const vec2i& sym1,       // (first, last) points in the raster for symmetry 1
	std::vector<Vertex>& S0, // points of P contained in sym0
	std::vector<Vertex>& S1, // points of P contained in sym1
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)