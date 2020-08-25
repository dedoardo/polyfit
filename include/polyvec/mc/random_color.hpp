#pragma once

#include <polyvec/api.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// returns a random color in RGBA (values in [0 1])
// Currently, alpha channel is always one.
Eigen::Vector4d random_color();


NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

