#pragma once

// polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

// PLEASE!
// Pass the accuracy map, we don't want to loop through all the edges
// Also re-indexes the regularities to account for removed nodes
void replace_flat_nodes_with_edges (
    const mat2x& raster,
    vecXi& polygon,
    const std::vector<bool>& raster_flat,
	polyfit::Regularity::RegularityInformation& regularities,
    std::vector<polyfit::Symmetry::SubPath>& raster_symmetries,
    std::vector<polyfit::Symmetry::SubPath>& raster_symmetries_local
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)