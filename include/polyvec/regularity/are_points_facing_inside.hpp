#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyvec)

/*
    Returns true if the line connecting the two points does not intersect any other line in the polygon
*/
bool are_points_facing_inside (
    const Eigen::Matrix2Xd& P,
    const int v0,
    const int v1
);

/*
    Returns true if the line connecting the two points does not intersect any other line in any other polygon
*/
bool are_points_facing_inside(
    const std::vector<const Eigen::Matrix2Xd*>& P,
    const int polygon_0,
    const int v0,
    const int polygon_1,
    const int v1
);

NAMESPACE_END(polyvec)