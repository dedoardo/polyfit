#pragma once

#include <polyvec/api.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

#include <polyvec/regularity/regularity_action.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void add_parallel_edges(
	const mat2x& raster,
	const bool circular,
	std::vector<RegularityAction*>& regularity_actions);


NAMESPACE_END(polyfit)
NAMESPACE_END(PolygonTracer)