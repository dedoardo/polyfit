#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

void test_and_set_parallel_alignment(
    const Eigen::Matrix2Xd& points,
    const std::vector<int>& C,
    const Eigen::VectorXi& V,
    Parallel& r
);

bool test_parallel(
    const Eigen::Matrix2Xd& B,
    const Eigen::VectorXi& V,
    const int v00,
    const int v01,
    const int v10,
    const int v11,
    Parallel&  r
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)