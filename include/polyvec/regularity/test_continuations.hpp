#pragma once

#include <polyvec/regularity/construct_regularity_graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

// Returns true if the two segments are a continuation
// The vertices in the edge are supposed to be initialized.
bool test_continuation(
    const Eigen::Matrix2Xd& PP_0, // First polygon points
    const Eigen::Matrix2Xd& PP_1, // Second polygon edge
    Continuation& r     // Output data for the continuation
);

void compute_continuation_accuracy_metrics(
    const Eigen::Matrix2Xd& PP_0, const Eigen::Matrix2Xd& PP_1,
    Continuation& r
);

void compute_continuation_curvature_metrics(
    const Eigen::Matrix2Xd& PP_0, const Eigen::Matrix2Xd& PP_1,
    Continuation& r
);

void filter_same_side_continuations_greedy(
    const std::vector<const Eigen::Matrix2Xd*>& polygons,
    std::vector<Continuation>& r,
    const bool circular
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)