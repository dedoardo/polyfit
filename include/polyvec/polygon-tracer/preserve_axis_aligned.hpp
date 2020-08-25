#pragma once

#include <polyvec/api.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void preserve_axis_aligned_edges(ShortestPath::State& G_state,
	const mat2x& B,
	const std::vector<Vertex>& M,
	std::vector<BoundaryGraph::Edge>& E,
	vecXi& P,
	const bool circular
);


NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)