#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

/*
	Given a polygonal approximation V to a raster boundary B, it returns all the 
	boundary midpoints which lie below the specified subpath. The boundary is not
	expected to contain only grid aligned points, but only those will be 
	considered as the midpoints are meant to be pixel corners

	The midpoints interval spanned by the subpath is closed and the first
	and last midpoints are also considered if overlapping.

	V_fit contains the pair of vertices which form the midpoint, while P_fit
	optionally stores the position. The vertices are not guaranteed to be consecutive
	as described above.

	This methods effectively guesses which ones are the midpoints that the subpath 
	is meant to approximate as we only have an exact mapping at the vertices V.
	
	It tests at which parametric value the closest point on the polygon lies
	w.r.t. the input bounds.

	The number of midpoints found is returned
*/
int find_pixel_centers_for_subpath(
	const mat2x& B,     // raster boundary
	const vecXi& V,     // vertices of the polygonal approximation
	const size_t v_beg, // index of the first edge (source endpoint)
	const double t_beg, // t along the first edge where the subpath begins,
	const size_t v_end, // index of the last edge (source endpoint)
	const double t_end, // t along the first edge where the subpath begin
	mat2xi& V_fit,		// output vertex previous to midpoint
	mat2x* P_fit, 		// output midpoints 
	bool circular
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)