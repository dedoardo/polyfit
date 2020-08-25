#pragma once

#include <polyvec/api.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// checks if two closed pixel paths are idenctical (takes into consideration that they can be shifted)
// is not really used in the current pipeline.
bool are_pixel_polygons_identical(const Eigen::Matrix2Xi& polygon1, const Eigen::Matrix2Xi& polygon2);

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

