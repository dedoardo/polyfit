#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/polygon-tracer/iterative-local.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	A locally symmetric fit is obtained with polygon-tracer/iterative-local
	Longer and disjoint symmetries are detected in the raster and
	they are applied iteratively as long as no inflections are introducted
	and no other global symmetries are disrupted

	in:
		R is a clockwise oriented polyline where points are grid-aligned

	out:
		P is the polygon fit whose vertices reference points in B
		B is the fitting path constructed from R
		E is the list of all possible edges in R

	returns the number of iterations of the local and global fitting
*/
//int iterative_global(
//	const mat2x& R,
//	vecXi& P,
//	mat2x& B,
//	std::vector<BoundaryGraph::Edge>& E,
//	const bool circular,
//	IterationCallback on_iter = nullptr
//);

int iterative_global(
	ShortestPath::State& s,
	const mat2x& B,
	const std::vector<Vertex>& M,
	const std::vector<BoundaryGraph::Edge>& E,
	vecXi& P,
	std::vector<Symmetry::SubPath>& raster_symmetries,
	const bool circular,
	IterationCallback on_iter = nullptr,
    const std::vector<bool>& corners_C0 = std::vector<bool>()
);

/*
	Fills in the list of edge indices from their values.
	The edge array is scanned linearly and the output array is cleared
	before inserting any other edge.

	Returns 0 if any of the edges in V could not be found in V, regardless
	I will contain all the valid edges
*/
int find_edge_indices(
	const std::vector<BoundaryGraph::Edge>& E, // List of edges where the index is computed from
	const std::vector<BoundaryGraph::Edge>& V, // Values for which the respective indices need to be computed
	std::vector<size_t>& I                     // Where the indices are written
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
