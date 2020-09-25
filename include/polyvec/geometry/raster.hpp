#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/options.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/line.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(GeomRaster)

// Returns true if p is within an epsilon of a intersection point between
// lines in a unit grid.
PV_INLINE bool is_grid_aligned(const vec2& p, const double tol = 1e-2);

// If P(v) is grid aligned, then it returns v. Otherwise it returns the next vertex
// along P which is grid aligned.
// If P contains no axis aligned points from v to itself (or to the last vertex if the
// path is not circular) then -1 is returned.
PV_INLINE Vertex find_next_grid_aligned(const mat2x& P, const Vertex v, const bool circular);

// Returns true if the two points are within an epsilon distance. L1 norm
PV_INLINE bool are_overlapping(const vec2& p0, const vec2& p1, const double tol = 1e-2);

// Returns true if the edge is length 1 and axis-aligned (assumes that points can only be placed
// at pixel centers or midpoints)
PV_INLINE bool is_edge_degenerate(const vec2& src, const vec2& dst);

// Returns the maximum inner and maximum outer distance from the line segment to the pixel centers
// of the approximated boundary.
// If the distance along the axis is above 0.5, the distance to the closest point on the line from
// the pixel center offset along the normal is considered.
PV_INLINE vec2 distance_bounds_from_points_with_slack (
    const mat2x& B,	  // raster boundary
    const vec2 p_src, // segment source point
    const vec2 p_dst, // segment destination point
    const int  v_src, // segment source vertex
    const int  v_dst  // segment destination vertex
);

// Given a value defined for 32x32 images it returns the scaling factor to be applied for an image with the given bounding box
PV_INLINE double get_resolution_scaling_factor_32x32(
    const mat2x& B
);

// Returns the size of the raster corner at the given raster vertex. Returns 0 if the vertex is on
// a straight raster segment. Otherwise, the size of the raster corner is the minimum length the
// raster continues in a constant direction to either side of the vertex.
PV_INLINE double raster_corner_size(const polyfit::mat2x& raster_boundary, Vertex v);

// Returns true if the points are axis-aligned
PV_INLINE bool are_points_axis_aligned(const vec2& p0, const vec2& p1, const double tol = 1e-4);

NAMESPACE_END(GeomRaster)
NAMESPACE_END(polyfit)

#include <polyvec/api.hpp>
#include "raster.inl"
