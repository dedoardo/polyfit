#ifndef POLYFIT_COLOR_H_
#define POLYFIT_COLOR_H_

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Color)

vec3 error(double t_min, double t_max, double t);

NAMESPACE_END(Color)
NAMESPACE_END(polyfit)

#endif // POLYFIT_COLOR_H_