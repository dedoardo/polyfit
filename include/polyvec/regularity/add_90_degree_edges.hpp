#pragma once

// polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/regularity/regularity_action.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void add_90_degree_edges(
    const Eigen::Matrix2Xd& raster,
    const Eigen::VectorXi& polygon,
    const bool circular,
    std::vector<RegularityAction*>& regularity_actions
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
