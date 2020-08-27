#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PathUtils)

// if the points have counter-clockwise ordering it reverses the path
void order_clockwise(mat2x& P); 


struct DegenerateClassifierAxisAligned1Px { static bool should_skip_edge(const vec2&, const vec2&); };
struct DegenerateClassifier2Px { static bool should_skip_edge(const vec2&, const vec2&); };

// Given a polygonal path and a corner, it finds the corner in the direction which is not degenerate
// and not nearly flat.
// The direction is determined by the pair of consecutive corners rather than a -1/+1
template <typename TDegenerateClassifier = DegenerateClassifierAxisAligned1Px>
int find_next_non_degenerate_corner (
    const Eigen::Matrix2Xd& P,
    const int vertex,          // Current vertex
    const int vertex_next,     // Next expected vertex
    const bool circular
);

// todo: use the WindingNumber::* function
int compute_orientation(const Eigen::Matrix2Xd& points);

// if direction is +1 it computes the convex hull of the points
// if direction is -1 it computes the concave hull of the the points
// the overlaps refers to the fact that this kind of configurations |_| 
// will result in duplicate vertices which won't be removed
void compute_hull_with_overlaps(const Eigen::Matrix2Xd& points, Eigen::Matrix2Xd& boundary, const int direction);

// Returns the length of the diagonal of the bounding box of a set of points
double bounding_box_diagonal(const Eigen::Matrix2Xd& points);

// fills the array with 0 if the points is a flat vertex, +1 if the point is convex and -1 if concave
// it assumes that the points are ordered clockwise
void compute_convexities(const Eigen::Matrix2Xd& points, std::vector<int>& convexities);

// Returns true if the edge causes a change of concavity between v0-v1 and v2-v3
bool is_edge_inflected(const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& v3);

// finds the next two consecutive non-flat points which have the same convexity and returns the most distant one
// if such edge does not exists it returns the last point regardless
// v0 can be flat
Vertex end_of_convex_region(const std::vector<int>& C, const Vertex v0, const int d, bool circular);

// returns true if v is in the closed interval [v0 v1]
bool contains_closed(const int size, Vertex v0, Vertex v1, Vertex v, bool circular);

// returns true if v is in the closed interval [v0 v1]
bool contains_closed(const mat2x& P, Vertex v0, Vertex v1, Vertex v, bool circular);

// return true if v is in the open interval (v0 v1)
bool contains_open(const mat2x& P, Vertex v0, Vertex v1, Vertex v, bool circular);

// Returns true if the distance between v0 and v equals v and v1, returning false
// if the path is not circular and a gap between any of the points exists
bool has_equal_distance(const mat2x& P, Vertex v0, Vertex v, Vertex v1, bool circular);

// It returns the next vertex in the array which is not flat (either +1 or -1)
Vertex next_transition_vertex(const std::vector<int>& convexities, const Eigen::Index from, const int direction, const bool circular = true);

// Similarly to next_transition_vertex it returns the distance between from and the next non-flat vertex
Vertex distance_to_next_corner(const std::vector<int>& convexities, const Eigen::Index from, const int direction);

// counts the number of vertices which have convexity == 0 from 'from' in direction 'direction', it doesn't
// include the last vertex
int count_flat_vertices_before_transition(const std::vector<int>& convexities, const Eigen::Index from, const Eigen::Index direction);

// Counts the non-flat vertices between, but excluding the source and destination
int count_transitions_between_vertices(const std::vector<int>& convexities, const Eigen::Index src, const Eigen::Index dst);

// checks if either the x/y coordinates are within epsilon distance
// move this somewhere else
bool are_axis_aligned(const Eigen::Vector2d& pt0, const Eigen::Vector2d& pt1);

// returns the angle spanned by the two vectors in range [0 PI]
double shortest_angle_spanned(const vec2& p0, const vec2& p1, const vec2& p2);

// checks if the edge p1p2 is an inflection. It checks if the sign of the cross product 
// to rotate the previous/next edge has a different sign
bool are_consecutive_angles_opposite(const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3);
	
// evaluates the distance between the midpoints of P in range vp(0) vp(1) to the edge. The edge doesn't
// necessarily have to start at specific points of P, hence the separated indices
vec2 distance_bounds_from_points (const mat2x& P, const mat2& E, vec2i vp);

// More compact version of the function above which returns true if the edge connecting the two vertices
// is within the relaxed accuracy threshold
bool edge_is_within_relaxed_distance_bounds(const mat2x& P, const int vsrc, const int vdst);

// returns true if exp is the first vertex encountered traversing the polygon P from i in direction d
bool points_in_same_partition(const mat2x& P, Vertex i, Vertex exp, Vertex fail, int d);

// returns the vertex previous to vsym in the opposite direction of v
Vertex opposite(const Eigen::Index& size, const Eigen::Index& v, const Eigen::Index& vsym);

// given an unordered list of vertices, it orders them in the same direction of the boundary
// placing the first point after vstart. vstart is not required to be in the list of vertices.
// The list *cannot* contain duplicate elements.
void order(const mat2x& P,   // points
	std::vector<Vertex>& V,  // vertices to be ordered
	Vertex vstart,           // suggested starting point
	bool circular = true);

// A wrapper for AngleUtils::spanned_shortest that returns the angle at a polygon corner
// attempting to handle circularity consistently.
// Returns the angle in radians or 0 if the angle does not exist.
double angle_at_corner(
	const mat2x& B,
	const vecXi& P,
	const int v,
	const bool circular
);

// Returns the number of degenerate edges (length 1 inflections aligned with the raster)
// which occurr outside the specified intervals.
int count_degenerate_edges_exclusive(
	const mat2x& P, // polygon vertices
	const mat2xi& E, // open intervals to exclude
	const bool circular
);

// Given a list of points and a list of vertices with the same orientation, it finds the first 
// vertex of P which lies before the queried point, excluding the point itself.
// Returns -1 if the path is empty or only contains the corner in question
Eigen::Index find_closest_corner_before_strict(
	const mat2x& B,      // list of points
	const vecXi& P,      // polygon vertices
	const Eigen::Index q // corner limit (index to B)
);

// Given a list of points and a list of vertices with the same orientation, it finds the first 
// vertex of P which lies after the queried point, excluding the point itself.
// Returns -1 if the path is empty or only contains the corner in question
Eigen::Index find_closest_corner_after_strict(
	const mat2x& B,      // list of points
	const vecXi& P,      // polygon vertices
	const Eigen::Index q // corner limit (index to B)
);

NAMESPACE_END(PathUtils)
NAMESPACE_END(polyfit)