#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>

#include <Eigen/Geometry>

#define PF_RAD(deg) ((double)deg * M_PI / 180.)
#define PF_DEG(rad) ((double)rad * 180. / M_PI)
#define PF_ISNAN(v) (!((v) == (v)))

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(AngleUtils)

/*
	Returns true if both points are in the same quadrant. The test is inclusive
	if any of the two points lie exactly along one of the axis.
*/
PV_INLINE bool lie_in_the_same_quadrant(
	const vec2& p0,
	const vec2& p1
);

PV_INLINE double deviation_from_horizontal(const vec2& d);
PV_INLINE double deviation_from_horizontal(double angle);

PV_INLINE double deviation_from_vertical(const vec2& d);
PV_INLINE double deviation_from_vertical(double angle);

PV_INLINE double spanned_shortest(const vec2 & p0, const vec2 & p1, const vec2 & p2);

PV_INLINE double spanned_shortest_between(const vec2 & d0, const vec2 & d1);

PV_INLINE double spanned_clockwise_between(const vec2& d0, const vec2& d1);

PV_INLINE bool have_opposite_convexity(const vec2 & p0, const vec2 & p1, const vec2 & p2, const vec2 & p3);

PV_INLINE bool have_opposite_convexity_with_tol(const vec2 & p0, const vec2 & p1, const vec2 & p2, const vec2 & p3);

PV_INLINE int count_inflections(const mat2x & points, bool circular = true);

PV_INLINE bool is_flat(const vec2 pp, const vec2 p, const vec2 pn);

PV_INLINE int convexity(const vec2 pp, const vec2 p, const vec2 pn, int orientation);

PV_INLINE int number_of_neighbors_with_different_convexity(const Eigen::Matrix2Xd& polygon, int corner);

// convexity should be writable for 1 + 2 * neighborhood_size elements
// the polygon should contain enough elements before and after the first/last neighboring corners to be able to compute the convexity
PV_INLINE void neighborhood_convexity(const Eigen::Matrix2Xd& P, const int corner, const int neighborhood_size, int* convexity);

/*
	Returns the number of edges in P which overlap 1-pixel axis-aligned steps in B.

	B is assumed to be the raster boundary which P is approximating, though there are 
	no requirements on the points except that they should have the same ordering.
*/
PV_INLINE int count_visually_degenerate_edges(
	const mat2x& B, 
	const vecXi& P, 
	const bool circular = false
);

PV_INLINE int count_visually_inconsistent_edges(const mat2x & P, bool circular = false);

#include "angle.inl"

NAMESPACE_END(AngleUtils)
NAMESPACE_END(polyfit)