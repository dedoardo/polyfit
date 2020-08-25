#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyvec)

std::vector<polyfit::Regularity::Continuation> find_continuation_candidates_2018(
    const Eigen::Matrix2Xd& B, 
    const Eigen::VectorXi& V,
    polyfit::Regularity::RegularityInformation& reg,
    const bool circular,
    const bool allow_move
);

NAMESPACE_END(polyvec)