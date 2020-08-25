// polyvec
#include <polyvec/polygon-tracer/detect_parallel_edges_relative_proximity.hpp>
#include <polyvec/regularity/check_opposite_points_winding_number.hpp>
#include <polyvec/regularity/are_points_facing_inside.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/path.hpp>

// libc++
#include <algorithm>

using namespace Eigen;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

// assumes the two intervals are ordered
bool test_overlap_interval_1d(
	const double a0,
	const double a1,
	const double b0,
	const double b1,
	vec2& i
) {
	if (b0 > a1 || a0 > b1) {
		return false;
	}
	
	i(0) = max(a0, b0);
	i(1) = min(a1, b1);
	return (i(1) - i(0)) > PF_EPS_LARGE;
}

int detect_parallel_edges_relative_proximity_axis_aligned(
	const mat2x& P,      // A list of points
	mat4xi& PL,          // pair of parallel straight edges, each column stores (src0, dst0, src1, dst1)
	const double min_feature_ratio,
	const double proximity_scale,
	const double min_length_ratio,
	const bool circular
) {
	int count = 0;

	vector<int> C;
	PathUtils::compute_convexities(P, C);

	vec2 bb_min(INFINITY, INFINITY);
	vec2 bb_max(-INFINITY, -INFINITY);

	for (Index i = 0; i < P.cols(); ++i) {
		bb_min = bb_min.cwiseMin(P.col(i));
		bb_max = bb_max.cwiseMax(P.col(i));
	}

	vec2 figure_size = bb_max - bb_min;

	PL.resize(4, 0);

	for (int i = 0; i < (int)C.size(); ++i) {
		if (!circular && i == C.size() - 1) {
			continue;
		}

		// corner to corner
		if (C[i] == 0) {
			continue;
		}

		const Vertex inext = PathUtils::next_transition_vertex(C, i, +1, circular);
		if (inext == -1) {
			continue;
		}

		const vec2 p0 = P.col(i);
		const vec2 p1 = CircularAt(P, inext);

		// skipping length 1 segments
		if (max(abs(p0(0) - p1(0)), abs(p0(1) - p1(1))) < 1. + PF_EPS) {
			continue;
		}

		for (int j = inext; j < (int)C.size(); ++j) {
			if (!circular && j == C.size() - 1) {
				continue;
			}

			// I don't understand how this can be true? (todo)
			if (i == j) {
				continue;
			}

			if (C[j] == 0) {
				continue;
			}

			const Vertex jnext = PathUtils::next_transition_vertex(C, j, +1, circular);
			if (jnext == -1) {
				continue;
			}

			const vec2 p2 = P.col(j);
			const vec2 p3 = CircularAt(P, jnext);
			
			// Test relative length of the two segments
			const double i_length = (p1 - p0).norm();
			const double j_length = (p3 - p2).norm();
			const double length_ratio = min(i_length, j_length) / max(i_length, j_length);
			if (length_ratio < min_length_ratio) {
				continue;
			}

			// skipping length 1 segments
			if (max(abs(p2(0) - p3(0)), abs(p2(1) - p3(1))) < 1. + PF_EPS) {
				continue;
			}

			// Are the two points facing each other?
            if (!polyvec::are_points_facing_inside(P, i, j) &&
                !polyvec::are_points_facing_inside(P, i, jnext) &&
                !polyvec::are_points_facing_inside(P, inext, j) &&
                !polyvec::are_points_facing_inside(P, inext, jnext)) {
                continue;
            }

			vec2 interval;

			// x-axis
			if (test_overlap_interval_1d(
				min(p0(0), p1(0)), max(p0(0), p1(0)),
				min(p2(0), p3(0)), max(p2(0), p3(0)),
				interval
			)) {
				PF_ASSERT(abs(p0(1) - p1(1)) < PF_EPS);
				PF_ASSERT(abs(p2(1) - p3(1)) < PF_EPS);

				const double overlap = abs(interval(0) - interval(1)); // abs unnecessary?
				const double proximity = abs(p0(1) - p2(1));

				// Test feature length with respect to overlap axis
				const double figure_ratio = min(i_length, j_length) / figure_size.x();
				if (figure_ratio < min_feature_ratio) {
					continue;
				}

				Index i00 = p0(0) < p1(0) ? i : inext;
				Index i01 = p0(0) < p1(0) ? inext : i;
				Index i10 = p2(0) < p3(0) ? j : jnext;
				Index i11 = p2(0) < p3(0) ? jnext : j;

				bool is_configuration_invalid = false;

				// Checking convexity of the lower bound
				if (abs(P(0, i00) - P(0, i10)) < PF_EPS_LARGE) {
					is_configuration_invalid = C[i00] != -1 || C[i10] != -1;
				} else if (P(0, i00) > P(0, i10)) {
					is_configuration_invalid = C[i00] != -1;
				} else {
					is_configuration_invalid = C[i10] != -1;
				}

				// Checking convexity of the upper bound
				if (!is_configuration_invalid) {
					if (abs(P(0, i01) - P(0, i11)) < PF_EPS_LARGE) {
						is_configuration_invalid = C[i01] != -1 || C[i11] != -1;
					}
					else if (P(0, i01) < P(0, i11)) {
						is_configuration_invalid = C[i01] != -1;
					}
					else if (P(0, i01) > P(0, i11)) {
						is_configuration_invalid = C[i11] != -1;
					}
				}

				if (!is_configuration_invalid && overlap > proximity_scale * proximity - PF_EPS) {
					bool duplicate_found = false;
					for (Index k = 0; k < PL.cols(); ++k) {
						if (PL(0, k) == j) {
							duplicate_found = true;
							break;
						}
					}

					if (!duplicate_found) {
						MatrixUtils::append(PL, vec4i(i, inext, j, jnext));
					}
				}
			}

			// y-axis
			if (test_overlap_interval_1d(
				min(p0(1), p1(1)), max(p0(1), p1(1)),
				min(p2(1), p3(1)), max(p2(1), p3(1)),
				interval
			)) {
				PF_ASSERT(abs(p0(0) - p1(0)) < PF_EPS);
				PF_ASSERT(abs(p2(0) - p3(0)) < PF_EPS);

				const double overlap = abs(interval(0) - interval(1)); // abs unnecessary?
				const double proximity = abs(p0(0) - p2(0));

				// Test feature length with respect to overlap axis
				const double figure_ratio = min(i_length, j_length) / figure_size.y();
				if (figure_ratio < min_feature_ratio) {
					continue;
				}

				Index i00 = p0(1) < p1(1) ? i : inext;
				Index i01 = p0(1) < p1(1) ? inext : i;
				Index i10 = p2(1) < p3(1) ? j : jnext;
				Index i11 = p2(1) < p3(1) ? jnext : j;

				bool is_configuration_invalid = false;

				// Checking convexity of the lower bound
				if (abs(P(1, i00) - P(1, i10)) < PF_EPS_LARGE) {
					is_configuration_invalid = C[i00] != -1 || C[i10] != -1;
				}
				else if (P(1, i00) > P(1, i10)) {
					is_configuration_invalid = C[i00] != -1;
				}
				else {
					is_configuration_invalid = C[i10] != -1;
				}

				// Checking convexity of the upper bound
				if (!is_configuration_invalid) {
					if (abs(P(1, i01) - P(1, i11)) < PF_EPS_LARGE) {
						is_configuration_invalid = C[i01] != -1 || C[i11] != -1;
					}
					else if (P(1, i01) < P(1, i11)) {
						is_configuration_invalid = C[i01] != -1;
					}
					else if (P(1, i01) > P(1, i11)) {
						is_configuration_invalid = C[i11] != -1;
					}
				}

				if (!is_configuration_invalid && overlap > proximity_scale * proximity - PF_EPS) {
					// Checking for duplicate, better to do it at the end?
					bool duplicate_found = false;
					for (Index k = 0; k < PL.cols(); ++k) {
						if (PL(0, k) == j) {
							duplicate_found = true;
							break;
						}
					}

					if (!duplicate_found) {
						MatrixUtils::append(PL, vec4i(i, inext, j, jnext));
					}
				}
			}
		}
	}

	return PL.cols();
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)