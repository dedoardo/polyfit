// polyvec
#include <polyvec/polygon-tracer/align_inflections_to_half_space.hpp>
#include <polyvec/polygon-tracer/prune_edges_force_vertex.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/path.hpp>

// Minimum length for an inflected edge to to be a valid candidate
// This number should be 2
#define MIN_CANDIDATE_LENGTH 2.
#define MIN_CANDIDATE_LENGTH_SQ (MIN_CANDIDATE_LENGTH * MIN_CANDIDATE_LENGTH - PF_EPS)

// Maximum angle spanned between two consecutive segments that want
// to be regularized to 90 degrees.
// This angle should be 120
#define MAX_CANDIDATE_ANGLE (M_PI_2 + PF_RAD(30) + PF_EPS)

// Log the polygon vertices after each alignment operation
#define LOG_POLYGON_ON_EDIT 1

using namespace Eigen;
using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

bool edge_is_candidate(
	const mat2x& B,			  // raster points
	const std::vector<int> C, // raster point convexities
	const int c, // shared corner
	const int i, // other inflection endpoint
	const int j  // other candidate edge endpoint
) {
	const vec2 pc = B.col(c);
	const vec2 pi = B.col(i);
	const vec2 pj = B.col(j);

	// check axis aligned
	if (!PathUtils::are_axis_aligned(pc, pi)) {
		return false;
	}

	// check reasonable length
	const double ci_norm_sq = (pc - pi).squaredNorm();
	const double cj_norm_sq = (pc - pj).squaredNorm();

	if (ci_norm_sq < MIN_CANDIDATE_LENGTH_SQ || cj_norm_sq < MIN_CANDIDATE_LENGTH_SQ) {
		return false;
	}

	// check nearly axis-aligned
	const double a_spanned = AngleUtils::spanned_shortest_between(pi - pc, pj - pc);
	if (a_spanned > MAX_CANDIDATE_ANGLE) {
		return false;
	}

	int d = CircularDist(B, c, j) > CircularDist(B, j, c) ? -1 : +1;
	const int len = PathUtils::count_flat_vertices_before_transition(C, c, d);
	
	if (len == 1) {
		return false;
	}

	return true;
}

// e0/e1 are not required to be ordered
// Returns 0 if not aligned 
// If a valid alignment is possible v_align is written with the index (of B) of the
// vertex which has been inserted
int align_edge_to_halfspace(
	const mat2x& B,            // raster boundary
	const std::vector<int>& C, // per pixel convexities |C| = |B|
	vecXi& P,                  // polygon vertices
	const Index iP_edit,	   // index in P where the new corner would ideally be placed
	const Index v_prev,        // point preceeding either e0 or e1
	const Index e0,			   // index of B of one of the two candidate edge endpoints
	const Index e1,			   // index of B of one of the two candidate edge endpoints
	Index&      v_align		   // new vertex inserted in P
) {
	PF_LOGF("e0 %d e1 %d edit %d", e0, e1, iP_edit);

	// ordering is not assumed, but since the points are consecutive the circular distance
	// between the two is a good enough discriminant
	const Index d01 = CircularDist(C, e0, e1);
	const Index d10 = CircularDist(C, e1, e0);

	// finding the next non-flat point
	const int e_dir = d01 > d10 ? -1 : +1;
	int v = Circular(B, e0 + e_dir);
	while (v != e1 && C[v] == 0) {
		v = Circular(B, v + e_dir);
	}

	const Index v_next = Circular(B, v + e_dir);

	if (v == e1) {
		PF_VERBOSE_F("candidate %d already aligned", v);
		return 0;
	}

	// Testing if the outgoing edge lies in the opposite halfspace of the previous edge
	const Vertex dprev0 = ShortestDist(C, v_prev, e0);
	const Vertex dprev1 = ShortestDist(C, v_prev, e1);
    PF_VERBOSE_F("dprev0 %d dprev1 %d", dprev0, dprev1);
	const Vertex e_src = dprev0 < dprev1 ? e0 : e1;
	const Vertex e_dst = e_src == e0 ? e1 : e0;
	const vec2 dprev = (B.col(v_prev) - B.col(e_src));
	const vec2 dnext = (B.col(e_dst) - B.col(e_src));
	
	if (dprev.dot(dnext) > 0) {
		PF_VERBOSE_F("candidate %d lies in the wrong half-space (%d %d) dot (%d %d)", v, e_src, v_prev, e_src, e_dst);
		return 0;
	}

	// Testing if the new edge creates an inflection
	const Vertex e_next = CircularAt(P, iP_edit + e_dir + e_dir);
	
    PF_VERBOSE_F("test inflection (%d %d %d %d)", e_src, v, e_dst, e_next);

	if (AngleUtils::have_opposite_convexity(
		B.col(e_src), B.col(v), B.col(e_dst), B.col(e_next)
	)) {
        PF_VERBOSE_F("candidate %d created inflection (%d %d %d %d)", v, e_src, v, e_dst, e_next);
		return 0;
	}

	// If the candidate is the pixel corner right before the next vertex
	// then we relocate, as otherwise we would introduce a degenerate edge
	// corresponding to the grid step. 
	// This is guaranteed to preserve the convexity
	if (v_next == e1) {
        PF_VERBOSE_F("relocate candidate %d", v);
		P(iP_edit) = v;
	}
	// Otherwise we insert a vertex right before the new one and check if
	// this violates the convexity. (can it?)
	else {
        PF_VERBOSE_F("insert candidate %d", v);
		InsertAt(P, iP_edit, v);
	}

	v_align = v;
	return 1;
}

int align_inflections_to_half_space(
	const mat2x& B,
	std::vector<BoundaryGraph::Edge>& E,
	vecXi& P,
	const bool circular
) {
	PF_LOG_DIVIDER;
	PF_LOGS("attempting to align inflections to half space");

	// corners for which the shortest path should be forced to go through
	std::vector<Index> corners_force;

	// per pixel convexities
	std::vector<int> C;
	PathUtils::compute_convexities(B, C);

	int aligned_tot = 0;
	int aligned;
	do {
		aligned = 0;

#if LOG_POLYGON_ON_EDIT
		PF_LOG_DIVIDER;
		PF_LOGF("polygon changed %d", P.size());
		for (Index i = 0; i < P.size(); ++i) {
			PF_LOGF("\t[%d] %d", i, P(i));
		}
#endif

		for (Index i = 0; i < P.size() && !aligned; ++i) {
			Index corner, inflection, candidate, i_edit; // candidate information
			Index v_align; // new vertex

			if (circular || i > 0) {
				corner = P(i);
				inflection = CircularAt(P, i + 1);
				candidate = CircularAt(P, i - 1);
				i_edit = i;

				if (edge_is_candidate(B, C, corner, inflection, candidate)) {
					if (align_edge_to_halfspace(B, C, P, i_edit, inflection, corner, candidate, v_align)) {
						corners_force.emplace_back(v_align);
						++aligned;
						continue;
					}
				}
			}

			if (circular || i < P.cols() - 2) {
				corner = CircularAt(P, i + 1);
				inflection = P(i);
				candidate = CircularAt(P, i + 2);
				i_edit = Circular(P, i + 1);

				if (edge_is_candidate(B, C, corner, inflection, candidate)) {
					if (align_edge_to_halfspace(B, C, P, i_edit, inflection, corner, candidate, v_align)) {
						corners_force.emplace_back(v_align);
						++aligned;
						continue;
					}
				}
			}
		}

		aligned_tot += aligned;
	} while (aligned);

	PF_LOGF("forcing %d corners in the graph", corners_force.size());
	for (int i = 0; i < corners_force.size(); ++i) {
		int pruned = prune_edges_force_vertex(B, E, corners_force[i], circular);
		PF_LOGF("pruned %d edges", pruned);
	}
	
	return aligned_tot;
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)