#pragma once

#include <polyvec/core/types.hpp>
#include <polyvec/core/macros.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Symmetry)

enum Axis {
	AXIS_HORIZONTAL = 0,
	AXIS_VERTICAL,
	AXIS_SEC_QUAD02,
	AXIS_SEC_QUAD13,
	AXIS_COUNT
};

// Defines two raster subsequences which are symmetric with respect to 'axis'
// The two sequences are not necessarily sequential, but the vertices are
// consistenly ordered. 
// if the points lay on a flat edge with an even number of vertices, the vertices
// won't match and is_joint0, is_joint1 should be called to properly check whether
// the subpath is sequential (if the raster is not circular this could cause segfaults).
// find_longest
// v00 - v10 are symmetric w.r.t. axis
// v01 - v11 are symmetric w.r.t. axis
// find shortest
// the vertices are symmetric as above, but for ease of use they point to the convex vertices 
// and not the flat ones.
// in practice the points can be reconstructed by traversing the boundary in the direction +1 between:
// v01 -> v00  v10 -> v11
struct SubPath {
	int axis = AXIS_COUNT;
	Vertex v00 = -1;
	Vertex v01 = -1;
	Vertex v10 = -1;
	Vertex v11 = -1;

	operator bool() const {
		return v00 != -1 && v01 != -1 && v10 != -1 && v11 != -1 && axis >= 0 && axis < AXIS_COUNT && v00 != v01 && v10 != v11;
	}

	bool operator == (const SubPath& other) const {
		bool same_axis = axis == other.axis;
		bool sameflip = v00 == other.v01 && v01 == other.v00 && v11 == other.v10 && v10 == other.v11;
		bool same0 = v00 == other.v00 && v01 == other.v01 && v10 == other.v10 && v11 == other.v11;
		bool same1 = v00 == other.v10 && v01 == other.v11 && v10 == other.v00 && v11 == other.v01;
		bool sameinv = v00 == other.v11 && v01 == other.v10 && v10 == other.v01 && v11 == other.v00;
		return same_axis && (same0 || same1 || sameflip || sameinv);
	}
};

// extracts the longest symmetric subpaths found in the list of raster points
std::vector<SubPath> find_longest(const mat2x& P, bool circular = true);

// extracts the smallest symmetric subpaths found int the list of points (3 edges)
std::vector<SubPath> find_shortest(const mat2x& P, bool circular = true);

// Reorders the vertices v* of the sequential symmetry s guaranteeing that v00/v10 are the points closest to each other
// and v10/v11 are the ones farthest away from each other. 
void order_vertices(std::vector<SubPath>& S);

// s.v00 matches/overlaps s.v10
bool is_subpath_joint0(const mat2x& P, const SubPath& s, bool circular = true);

// s.v01 matches/overlaps s.v11
bool is_subpath_joint1(const mat2x& P, const SubPath& s, bool circular = true);

// all symmetry endpoints differ
bool is_subpath_disjoint(const mat2x& P, const SubPath& s, bool circular = true);

// hides the logic behind the ordering of the symmetries returning true if s is
// a sequential reflective symmmetry and v is the first point on the reflection axis.
bool is_point_symmetric(const SubPath& s, const Vertex v, bool circular = true);

// todo: address this in the filtering
// todo: can remove, check references first
vec2i subpath_bounds0(const mat2x& P, const SubPath& s);
vec2i subpath_bounds1(const mat2x& P, const SubPath& s);

// todo: this should be integrated directly in the code for the shortest symmetries
// returns the number of points that a sequential symmetric region spans from start to end
Vertex calculate_length(const mat2x& P, const std::vector<int>& turns, const SubPath& s, bool circular);

NAMESPACE_END(Symmetry)
NAMESPACE_END(polyfit)
