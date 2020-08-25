#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/polygon-tracer/iterative-global.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/regularity/symmetry.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int regularized(
    ShortestPath::State& s,
    const mat2x& B,
    const std::vector<Vertex>& M,
    std::vector<BoundaryGraph::Edge>& E,
    polyfit::Regularity::RegularityInformation& RE,
    vecXi& P,
    std::vector<Symmetry::SubPath>& raster_symmetries,
    std::vector<Symmetry::SubPath>& raster_symmetries_local,
    const bool circular,
    IterationCallback on_iter = nullptr,
    const std::vector<bool>& corners_C0 = std::vector<bool>(),
    const std::vector<bool>& B_flat = std::vector<bool>(),
    const std::unordered_set<int>& kopf_junctions = std::unordered_set<int>()
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)