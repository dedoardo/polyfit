#pragma once

#include <polyvec/core/types.hpp>
#include <polyvec/core/macros.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

bool check_opposite_points_winding_number (
	const mat2x& P,
	const Vertex i,
	const Vertex j
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)