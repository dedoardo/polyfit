#include <polyvec/regularity/symmetry.hpp>

#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/debug.hpp>

#define SHRINK_SYMMETRY_TO_CORNER_SAME_CONVEXITY 1
#define DEBUG_TRAVERSE 0
#define DEBUG_CRITERIA 0
#define INVERTED_Y 1 // (0, 0) is at the top-left of the image

#if DEBUG_TRAVERSE
	#define TRAVERSE_LOG(fmt, ...) PF_LOGF(fmt, __VA_ARGS__)
#else
	#define TRAVERSE_LOG(...)
#endif

#if DEBUG_CRITERIA
#define CRITERIA_LOG(fmt, ...) PF_LOGF(fmt, __VA_ARGS__)
#else
#define CRITERIA_LOG(...)
#endif

using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Symmetry)

enum Directions {
	DIR_LEFT = 0,
	DIR_RIGHT,
	DIR_DOWN,
	DIR_UP,
	DIR_COUNT
};

const int axis_horz_flip[DIR_COUNT] = {
	DIR_LEFT, DIR_RIGHT, DIR_UP, DIR_DOWN
};

const int axis_vert_flip[DIR_COUNT] = {
	DIR_RIGHT, DIR_LEFT, DIR_DOWN, DIR_UP
};

const int axis_quad02_flip[DIR_COUNT] = {
	DIR_DOWN, DIR_UP, 
#if INVERTED_Y
	DIR_LEFT, DIR_RIGHT
#else
	DIR_RIGHT, DIR_LEFT
#endif
};

const int axis_quad13_flip[DIR_COUNT] = {
	DIR_UP, DIR_DOWN,
#if INVERTED_Y
	DIR_RIGHT, DIR_LEFT
#else
	DIR_LEFT, DIR_RIGHT

#endif
};

const int* sym_axis_flip[AXIS_COUNT] = {
	axis_horz_flip,
	axis_vert_flip,
	axis_quad02_flip,
	axis_quad13_flip
};

inline int direction(const vec2& p0, const vec2& p1) {
	const vec2 p01 = p1 - p0;
	PF_ASSERT(p01.cwiseAbs().minCoeff() < PF_EPS); // axis-aligned only
	PF_ASSERT(p01.cwiseAbs().maxCoeff() > PF_EPS); // separate points
			
	if (p01.cwiseAbs().x() > p01.cwiseAbs().y()) {
		return p01.x() > 0 ? DIR_RIGHT : DIR_LEFT;
	} else {
		return p01.y() > 0 ? DIR_DOWN : DIR_UP;
	}
}

// grows the region starting from v0 and v1 in both direction as long as the turns match.
// this assumes that v0 and v1 are part of the same boundary (not necessarily circular), that is
// turns are checked in opposite directions.
// This can be extended to multiple boundaries by having arbitrary directions for the two start points
void find_symmetric_subpath(const mat2x& P, const Vertex v0, const Vertex v1, SubPath& S, int sym_axis, bool circular) {
	S.axis = sym_axis;
	S.v00 = S.v01 = v0;
	S.v10 = S.v11 = v1;

	TRAVERSE_LOG("find symmetric subpath v0 %lld v1 %lld axis %d", v0, v1, sym_axis);
	const int max_turns = (int)(P.cols() / 2);

	// walking along the the positive direction from v0 and negative from v1
	Vertex v0n = v0, v1p = v1;
	while (v0n != v1p && CircularDist(P, S.v01, S.v00) < max_turns) {
		const Vertex v0nn = Circular(P, v0n + 1);
		const Vertex v1pp = Circular(P, v1p - 1);

		TRAVERSE_LOG("v0nn %lld v1pp %lld", v0nn, v1pp);

		// if the inputs are symmetric up to a flat edge which has an even number of vertices
		if (v1pp == v0n || v0nn == v1p) {
			TRAVERSE_LOG("overlapping at v1pp %lld v0nn %lld", v1pp, v0nn);
			break;
		}

		// Meeting, if the path has been symmetric until now the last edge will also be symmetric
		// updating this side of the bounds and moving to the other direction
		if (v0n == v1p) {
			TRAVERSE_LOG("meeting at v0n == v1p == %lld", v0n);
			break;
		}

		const int d0n = direction(P.col(v0n), P.col(v0nn));
		const int d1p = direction(P.col(v1p), P.col(v1pp));
		TRAVERSE_LOG("direction d0n %lld", d0n);
		TRAVERSE_LOG("direction d1p %lld", d1p);

		// This next edge is not symmetric
		if (sym_axis_flip[sym_axis][d0n] != d1p) {
			TRAVERSE_LOG("terminating not symmetric flip(d0n) = %d", sym_axis_flip[sym_axis][d0n]);
			break;
		}

		// moving to the next points
		v0n = v0nn;
		v1p = v1pp;

		// extending the subpath bounds
		S.v00 = v0n;
		S.v10 = v1p;
	}

	// walking along the negative direction from v0 and positive from v1
	// same loop as above (it could be moved to the same function which takes a direction parameter)
	Vertex v1n = v1, v0p = v0;
	while (v1n != v0p && CircularDist(P, S.v01, S.v00) < max_turns) {
		const Vertex v1nn = Circular(P, v1n + 1);
		const Vertex v0pp = Circular(P, v0p - 1);

		TRAVERSE_LOG("v1nn %lld v0pp %lld", v1nn, v0pp);

		if (v0pp == v1n || v1nn == v0p) {
			// PF_LOGF("overlapping v0p %lld v0pp %lld v1n %lld v11 %lld", v0p, v0pp, v1n, v1nn);
			break;
		}

		// Meeting at point belonging to the symmetric axis
		if (v1n == v0p) {
			TRAVERSE_LOG("meeting at v1p == v0n == %lld", v1p);
			break;
		}

		const int d0p = direction(P.col(v0p), P.col(v0pp));
		const int d1n = direction(P.col(v1n), P.col(v1nn));
		TRAVERSE_LOG("symmetries d0p %d d1n %d", d0p, d1n);

		// Edge is not symmetric
		if (sym_axis_flip[sym_axis][d1n] != d0p) {
			TRAVERSE_LOG("terminating not symmetric flip(d1n) = %d", sym_axis_flip[sym_axis][d1n]);
			break;
		}

		// moving to the next points
		v0p = v0pp;
		v1n = v1nn;

		// extending the subpath bounds
		S.v01 = v0p;
		S.v11 = v1n;
	}

	TRAVERSE_LOG("symmetric subpath v00 %lld v01 %lld v10 %lld v11 %lld", S.v00, S.v01, S.v10, S.v11);
}

// we only want to detect meaningful symmetries, thus consider subpaths as symmetric if there
// are at least two consecutive identical turns.
// If the supath is sequential (v00 == v10 || v01 == v11) the criteria has to be satisfied
// at least once along the *joined* boundary.
// 
// if v00 == v10 && v10 == v11:
//	    true
// if v00 overlaps v10: (joint0)
//      check v01 -> v11
// if v11 overlaps v01: (joint1)
//      check v10 -> v00
bool check_criteria_segment(const mat2x&P, const std::vector<int>& convexity, Vertex v0, Vertex v1, int direction);
bool check_criteria(const mat2x& P, const std::vector<int>& C, SubPath& s, bool circular) {
	//CRITERIA_LOG("checking subpath v00 %lld v01 %lld v10 %lld v11 %lld", s.v00, s.v01, s.v10, s.v11);

	// the shrinking could have resulted in subpaths of different length
	Vertex d0, d1;
	if (circular) {
		d0 = CircularDist(P, s.v01, s.v00);
		d1 = CircularDist(P, s.v10, s.v11);
	} else {
		d0 = (s.v00 - s.v01);
		d1 = (s.v11 - s.v10);
	}

	if (d0 != d1) {
		return false;
	}

	// This should really be handled during the traversal, but it's 3:29AM and it does the trick
	if (d0 > P.cols() / 2 || d1 > P.cols() / 2) {
		return false;
	}

	int last_turn = 0;
	const bool joint0 = is_subpath_joint0(P, s, circular);
	const bool joint1 = is_subpath_joint1(P, s, circular);

	if (joint0 && joint1) {
		return true;
	}

	if (joint0) {
		return check_criteria_segment(P, C, s.v01, s.v11, +1);
	}

	if (joint1) {
		return check_criteria_segment(P, C, s.v10, s.v00, +1);
	}

	return check_criteria_segment(P, C, s.v01, s.v00, +1) &&
			check_criteria_segment(P, C, s.v10, s.v11, +1);
}

// if either of the endpoints in the symmetric subpaths lies on a flat vertex, it is moved
// to the nearest corner in the appropriate direction
void shrink_to_corners (const mat2x& P, const std::vector<int>& C, SubPath& s, bool circular) {
	SubPath s0 = s;

	// shrinking the subpaths to the nearest corner being careful that v* can be valid
	// while still corresponding to a flat vertex. 
	// The shrinking should also be "symmetric", that is applied to both subpaths simoultaneously
	if (!is_subpath_joint0(P, s)) {
#if SHRINK_SYMMETRY_TO_CORNER_SAME_CONVEXITY
		while (C[s.v00] != C[s.v10]) {
			s.v00 = PathUtils::next_transition_vertex(C, s.v00, -1);
			s.v10 = PathUtils::next_transition_vertex(C, s.v10, +1);
		}
#else
		if (C[s.v00] == 0 || C[s.v10] == 0) {
			s.v00 = PathUtils::next_transition_vertex(C, s.v00, -1);
			s.v10 = PathUtils::next_transition_vertex(C, s.v10, +1);
		}
#endif
	}

	if (!is_subpath_joint1(P, s)) {
#if SHRINK_SYMMETRY_TO_CORNER_SAME_CONVEXITY
		while (C[s.v01] != C[s.v11]) {
			s.v01 = PathUtils::next_transition_vertex(C, s.v01, +1);
			s.v11 = PathUtils::next_transition_vertex(C, s.v11, -1);
		}
#else
		if (C[s.v01] == 0 || C[s.v11] == 0) {
			s.v01 = PathUtils::next_transition_vertex(C, s.v01, +1);
			s.v11 = PathUtils::next_transition_vertex(C, s.v11, -1);
		}
#endif
	}
			
	//PF_LOGF("shrink v00 %lld(%lld) v01 %lld(%lld) v10 %lld(%lld) v11 %lld(%lld)", s.v00, s0.v00, s.v01, s0.v01, s.v10, s0.v10, s.v11, s0.v11);
}

// if v0 and v1 are bounds of a symmetric region, the points inbetween are always valid 
// and circularity is not needed
bool check_criteria_segment(const mat2x& P, const std::vector<int>& convexity, Vertex v0, Vertex v1, int dir) {
	// I don't see why this could happen, but in case the code below would run an extra cycle
	if (v0 == v1) {
		return false;
	}

	int last_turn = 0;
	Vertex v_last_turn = -1;
	int turns = 1;

	// starting to count inflections from the point after, which avoid various noisy matches in organic
	// models. This doesn't remove subpaths which are symmetric w.r.t. to a single midpoint
	Vertex v = Circular(convexity, v0 + dir);

	bool includes_flat_edge = false;
	while (v != v1) {
		if (convexity[v] != 0) {
			if (convexity[v] == last_turn) {
				// considering this a flat edge only if it's longer than one pixel
				const double d = (P.col(v) - P.col(v_last_turn)).norm();
						
				if (d > 1. + PF_EPS) {
					includes_flat_edge = true;
				}
			}

			last_turn = convexity[v];
			v_last_turn = v;
			++turns;
		}

		v = Circular(convexity, v + dir);
	}

	if (!includes_flat_edge || turns <= 6) {
		return false;
	}

	return true;
}

bool are_candidate(vec2 p0, vec2 p1, int axis) {
	switch (axis) {
	case AXIS_VERTICAL:
		return abs(p0.y() - p1.y()) < PF_EPS;
	case AXIS_HORIZONTAL:
		return abs(p0.x() - p1.x()) < PF_EPS;
	case AXIS_SEC_QUAD13: {
		vec2 d01 = (p1 - p0).normalized();
		double dev0 = atan2(d01.y(), d01.x());
		return abs(dev0 + M_PI_4) < PF_EPS ||
			abs(dev0 - (3. / 4) * M_PI) < PF_EPS;
	}

	case AXIS_SEC_QUAD02: {
		vec2 d01 = (p1 - p0).normalized();
		double dev0 = atan2(d01.y(), d01.x());
		return abs(dev0 - M_PI_4) < PF_EPS ||
			abs(dev0 + (3. / 4) * M_PI) < PF_EPS;
	}
	default:
		break;
	}
	return false;
}

void order_symmetry_vertices_to_canonical(
	const mat2x& B,
	std::vector<Symmetry::SubPath>& S,
	bool circular
) {
	for (size_t i = 0; i < S.size(); ++i) {
		auto& s = S[i];

		// reordering the symmetry points consistently with the orientation of the boundary
		bool joint0 = Symmetry::is_subpath_joint0(B, s, circular);
		bool joint1 = Symmetry::is_subpath_joint1(B, s, circular);

		// making v00 and v10 the consecutive vertices
		if (joint1 && !joint0) {
			swap(s.v00, s.v01);
			swap(s.v10, s.v11);
		}

		// making the direction +1 from v11 -> v01
		if (CircularDist(B, s.v01, s.v00) > CircularDist(B, s.v00, s.v01)) {
			swap(s.v01, s.v11);
			swap(s.v00, s.v10);
		}
	}
}

std::vector<SubPath> find_longest(const mat2x& P, bool circular) {
	std::vector<int> convexity;
	PathUtils::compute_convexities(P, convexity);

	std::vector<SubPath> S;

	// Finding any two points which are potentially symmetric and walking along boundary
	// in both directions respectively checking that the direction to the next match
	// the symmetric version in the opposite point.
	for (Vertex i = 0; i < P.cols(); ++i) {
		for (Vertex j = i + 2; j < P.cols(); ++j) {
			// finding the next point after two points with two convexity
			if (!check_criteria_segment(P, convexity, i, j, +1) || 
				convexity[i] != convexity[j]) {
				continue;
			}

			for (int k = 0; k < AXIS_COUNT; ++k) {
				if (are_candidate(P.col(i), P.col(j), k)) {
					//PF_LOGF("candidates %lld %lld", i, j);

					SubPath s;
					find_symmetric_subpath(P, i, j, s, k, circular);

					if (s) {
						shrink_to_corners(P, convexity, s, circular);

						// exists already 
						if (find(S.begin(), S.end(), s) != S.end()) {
							continue;
						}

						// filters sequences which would be captured by find_shortest
						if (!check_criteria(P, convexity, s, circular)) {
							continue;
						}

						if (s) {
							// PF_LOGF("found symmetric subpath v00 %lld v01 %lld v10 %lld v11 %lld", s.v00, s.v01, s.v10, s.v11);
							S.emplace_back(std::move(s));
						}
					}
				}
			}
		}
	}

	// reorder the vertices of each symmetry so that they are ordered consistently (v01 -> v00 -> v10 -> v11)
	order_symmetry_vertices_to_canonical(P, S, circular);

	return S;
}

Vertex find_next_transition(const std::vector<int>& C, Vertex v, bool circular) {
	Vertex vn = v;
	do {
		if (!circular && v == C.size() - 1) {
			return -1;
		}
				
		vn = Circular(C, vn + 1);
				
		if (C[vn] != 0) {
			return vn;
		}
	} while (vn != v);

	return -1;
}

// For every 4 consecutive points which have the following concavities
// -1 +1 +1 -1 or +1 -1 -1 +1
// and the distance between the first and second and third and fourth point matches than
// it is considered one.
// To detect 45 degrees ones 
std::vector<SubPath> find_shortest(const mat2x& P, bool circular) {
	std::vector<int> C;
	PathUtils::compute_convexities(P, C);
			
	std::vector<SubPath> S;

	for (Vertex v = 0; v < P.cols(); ++v) {
		if (!circular && (v == 0 || v < P.cols() - 3)) {
			continue;
		}

		if (C[v] != 0) {
			Vertex vn = find_next_transition(C, v, circular);
			Vertex vnn = find_next_transition(C, vn, circular);
			Vertex vnnn = find_next_transition(C, vnn, circular);

			if (vn == -1 || vnn == -1 || vnnn == -1) {
				continue;
			}

			// axis aligned symmetry
			if (C[v] == C[vnnn] && C[vn] == C[vnn]) {
				double d0 = (P.col(v) - P.col(vn)).norm();
				double d1 = (P.col(vnnn) - P.col(vnn)).norm();

				if (abs(d0 - d1) < PF_EPS) {
					// PF_LOGF("found short symmetry %lld %lld %lld %lld", v, vn, vnn, vnnn);

					SubPath s;
					s.v01 = v;
					s.v00 = vn;
					s.v10 = vnn;
					s.v11 = vnnn;
										
					if (abs(P.col(vn).x() - P.col(vnn).x()) < PF_EPS) {
						s.axis = Symmetry::AXIS_VERTICAL;
					}
					else if (abs(P.col(vn).y() - P.col(vnn).y()) < PF_EPS) {
						s.axis = Symmetry::AXIS_HORIZONTAL;
					}
					else {
						PF_LOGF("unknown relationship in symmetric transition v00 %lld v01 %lld v10 %lld v11 %lld",
							s.v00, s.v01, s.v10, s.v11);
						PF_ABORT;
					}
							
					S.emplace_back(std::move(s));
				}
			}

			if (!circular && v >= P.cols() - 2) {
				continue;
			}

			Vertex v0 = v;
			Vertex v1 = Circular(P, v + 1);
			Vertex v2 = Circular(P, v + 2);

			// 45 degrees symmetries
			if (C[v0] == C[v2] && C[v0] == -C[v1]) {
				SubPath s;
				s.v00 = v0;
				s.v10 = v2;
				s.v01 = PathUtils::next_transition_vertex(C, v0, -1);
				s.v11 = PathUtils::next_transition_vertex(C, v2, +1);

				// filtering out only those for which the distance to the next transitions
				// is greater than two, which means that there won't be an edge crossing over the
				// symmetric region and we *reasonably* expect consistent behavior
				if ((P.col(s.v00) - P.col(s.v01)).norm() < 2. ||
					(P.col(s.v10) - P.col(s.v11)).norm() < 2.) {
					continue;
				}

				if (are_candidate(P.col(s.v00), P.col(s.v10), AXIS_SEC_QUAD02)) {
					s.axis = AXIS_SEC_QUAD02;
				} else if (are_candidate(P.col(s.v00), P.col(s.v10), AXIS_SEC_QUAD13)) {
					s.axis = AXIS_SEC_QUAD13;
				} else {
					// something is not axis aligned or the convexities are wrong
					PF_ABORT;
				}

				S.emplace_back(std::move(s));
			}
		}
	}

	return S;
}

void order_vertices(mat2x& R, std::vector<SubPath>& S, bool circular = false) {
	for (int i = 0; i < S.size(); ++i) {
		auto& s = S[i];

		// reordering the symmetry points consistently with the orientation of the boundary
		bool joint0 = Symmetry::is_subpath_joint0(R, s, circular);
		bool joint1 = Symmetry::is_subpath_joint1(R, s, circular);

		// making v00 and v10 the consecutive vertices
		if (joint1 && !joint0) {
			swap(s.v00, s.v01);
			swap(s.v10, s.v11);
		}

		// making the direction +1 from v11 -> v01
		if (CircularDist(R, s.v01, s.v00) > CircularDist(R, s.v00, s.v01)) {
			swap(s.v01, s.v11);
			swap(s.v00, s.v10);
		}
	}
}

// s.v00 matches/overlaps s.v10
bool is_subpath_joint0(const mat2x& P, const SubPath& s, bool circular) {
	const Vertex d = circular ? CircularDist(P, s.v00, s.v10) : (s.v10 - s.v00);
	return d <= 1;
}

// s.v01 matches/overlaps s.v11
bool is_subpath_joint1(const mat2x& P, const SubPath& s, bool circular) {
	const Vertex d = circular ? CircularDist(P, s.v11, s.v01) : (s.v01 - s.v11);
	return d <= 1;
}

// all symmetry endpoints differ
bool is_subpath_disjoint(const mat2x& P, const SubPath& s, bool circular) {
	return !is_subpath_joint0(P, s, circular) && !is_subpath_joint1(P, s, circular);
}

bool is_point_symmetric(const mat2x& P, const SubPath& s, const Vertex v, bool circular) {
	if (!is_subpath_joint0(P, s, circular) && !is_subpath_joint1(P, s, circular)) {
		return false;
	}
	return true;
}

vec2i subpath_bounds0(const mat2x& P, const SubPath& s) {
	if (CircularDist(P, s.v00, s.v01) > CircularDist(P, s.v01, s.v00)) {
		return vec2i(s.v01, s.v00);
	} else {
		return vec2i(s.v00, s.v01);
	}
}

vec2i subpath_bounds1(const mat2x& P, const SubPath& s) {
	if (CircularDist(P, s.v10, s.v11) > CircularDist(P, s.v11, s.v10)) {
		return vec2i(s.v11, s.v10);
	}
	else {
		return vec2i(s.v10, s.v11);
	}
}

Vertex calculate_length(const mat2x& P, const std::vector<int>& turns, const SubPath& s, bool circular) {
	PF_ASSERT(P.cols() == turns.size());

	// PF_LOG_DIVIDER

	Vertex vp = s.v01;
	Vertex vn = s.v11;
	while (CircularDist(P, vp, s.v01) > 
		CircularDist(P, vn, s.v01)) {
		if (!circular) {
			if (vp - 1 < 0 || vn + 1 >= P.cols()) {
				break;
			}
		}

		Vertex vpp = Circular(P, vp - 1);
		Vertex vnn = Circular(P, vn + 1);

		// PF_LOGF("vpp %d turn %d vnn %d turn %d", vpp, turns[vpp], vnn, turns[vnn]);

		if (turns[vpp] != turns[vnn]) {
			break;
		}

		vp = vpp;
		vn = vnn;
	}

	// PF_LOGF("symmetric region v01 %d v00 %d v10 %d v11 %d length %d vp %d vn %d", s.v01, s.v00, s.v10, s.v11, CircularDist(P, vp, vn), vp, vn);
	return CircularDist(P, vp, vn) + 1;
}

NAMESPACE_END(Symmetry)
NAMESPACE_END(polyfit)