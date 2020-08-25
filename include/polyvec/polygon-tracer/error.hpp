#pragma once

#include <polyvec/api.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

// Calculates the error of each edge for each component without recomputing 
// the error measures
vec3 error_per_metric(const mat2x& B, 
	const std::vector<BoundaryGraph::Edge>& E,
	const vecXi& P,
	const bool circular);

// returns the sum of each individual cost. it's cost_individual.sum(), but
// it's better to encapsulate this.
double error_total(const mat2x& B, 
	const std::vector<BoundaryGraph::Edge>& E,
	const vecXi& P, 
	const bool circular);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)