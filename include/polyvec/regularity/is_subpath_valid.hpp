#ifndef POLYFIT_SYMMETRY_IS_SUBPATH_VALID_H_
#define POLYFIT_SYMMETRY_IS_SUBPATH_VALID_H_

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Symmetry)

/*
	Given polygonal path and two symmetry subpaths, it checks if the
	vertices for a symmetric polygons within that region

	S0 and S1 are expected to be ordered as follows
	S0(0) < S0(1) <= S1(0) <= S1(1)

	in:
		P  polygonal path
		S0 vertices of P defining the bounds of the first symmetric region
		S1 vertices of P defining the bounds of the second symmetric region
		V vertices of P defining the polygon to be checked for symmetry within S0 and S1
*/
bool is_subpath_valid(
	const mat2x& P,
	const vec2i S0,
	const vec2i S1,
	const vecXi V,
	const bool circular
);

NAMESPACE_END(Symmetry)
NAMESPACE_END(polyfit)

#endif // POLYFIT_SYMMETRY_IS_SUBPATH_VALID_H_