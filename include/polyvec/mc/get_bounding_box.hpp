#pragma once

#include <polyvec/api.hpp>

#include <Eigen/Core>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// Get the bounding box of a bunch of paths
Eigen::AlignedBox<int, 2> get_bounding_box(const std::vector<Eigen::Matrix2Xi>& regions);

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)
