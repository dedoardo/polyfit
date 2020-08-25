#pragma once

#include <polyvec/api.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// NOTE: not being used, can be thrown away.
//
// Sort points  to get a canonical representation
// this way we can quickly check if two pixel polygons are identical (duplicates)
Eigen::Matrix2Xi sort_points(const Eigen::Matrix2Xi& points);


NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

