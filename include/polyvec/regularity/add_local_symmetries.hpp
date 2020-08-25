#pragma once

// polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/regularity/regularity_action.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

// Adds all local symmetries to the list of regularity actions
void add_local_symmetries(
	const mat2x& raster,
    std::vector<Symmetry::SubPath>& raster_symmetries_local,
    const std::vector<bool>& raster_flat, // make this a set pls
	const bool circular,
	std::vector<RegularityAction*>& regularity_actions);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
