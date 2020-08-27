// polyvec
#include <polyvec/polygon-tracer/align_opposite_corners.hpp>
#include <polyvec/regularity/check_opposite_points_winding_number.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace polyvec;
using namespace Eigen;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int test_symmetry_1 (
	const mat2x& P,
	const int v
) {
	const vec2 p_prev = CircularAt(P, v - 1);
	const vec2 p = P.col(v);
	const vec2 p_next = CircularAt(P, v + 1);

	const double angle0 = AngleUtils::deviation_from_horizontal((p - p_prev).normalized());
	const double angle1 = AngleUtils::deviation_from_horizontal((p_next - p).normalized());

	const bool same_length =
		((p - p_prev).squaredNorm() - (p - p_next).squaredNorm()) < PF_EPS;

	const bool same_angle = abs(angle0 - angle1) < PF_EPS;
	return same_length && same_angle;
}

int test_symmetry_2 (
	const mat2x& P,
	const int v,
	const int d
) {
	const vec2 p0 = CircularAt(P, v - d);
	const vec2 p1 = CircularAt(P, v);
	const vec2 p2 = CircularAt(P, v + d * 1);
	const vec2 p3 = CircularAt(P, v + d * 2);

	if (abs(p1(0) - p2(0)) > PF_EPS && abs(p1(1) - p2(1)) > PF_EPS) {
		return false;
	}

	const double angle0 = AngleUtils::deviation_from_horizontal((p3 - p2).normalized());
	const double angle1 = AngleUtils::deviation_from_horizontal((p1 - p0).normalized());

	const bool same_angle = abs(angle0 - angle1) < PF_EPS;
	return same_angle;
}

int enumerate_moves (
	const mat2x& B,
	const mat2x& P,
	const std::vector<int>& C,
	const int v,
	const int c,
	int* V
) {
	// If any of the two edges is axis aligned, then we absolutely do not wish to move it!
	const vec2 p_prev = CircularAt(P, c - 1);
	const vec2 p = P.col(c);
	const vec2 p_next = CircularAt(P, c + 1);
	
	PF_VERBOSE_F("p_prev %f %f", p_prev(0), p_prev(1));
	PF_VERBOSE_F("p      %f %f", p(0), p(1));
	PF_VERBOSE_F("p_next %f %f", p_next(0), p_next(1));

	if (abs(p_prev(0) - p(0)) < PF_EPS || abs(p_prev(1) - p(1)) < PF_EPS ||
		abs(p_next(0) - p(0)) < PF_EPS || abs(p_next(1) - p(1)) < PF_EPS) {
		V[0] = v;
		return 1;
	}

	if (test_symmetry_1(P, Circular(P, c - 1)) || test_symmetry_1(P, Circular(P, c + 1)) ||
		test_symmetry_2(P, Circular(P, c - 1), -1) || test_symmetry_2(P, Circular(P, c + 1), +1)) {
		V[0] = v;
		return 1;
	}

	if (C[v] == 0) { // midpoint
		if (CircularAt(C, v - 1) != -1 || CircularAt(C, v + 1) != -1) {
			V[0] = v;
			return 1;
		}

		V[0] = Circular(C, v - 1);
		V[1] = v;
		V[2] = Circular(C, v + 1);
		return 3;
	} else if (C[v] == 1) {
		V[0] = v;
		return 1;
	} else {
        if (CircularAt(C, v - 1) == 0 && CircularAt(C, v - 2) == -1) {
            V[0] = v;
            V[1] = Circular(C, v - 1);
            V[2] = Circular(C, v - 2);
            return 3;
        }

        if (CircularAt(C, v + 1) == 0 && CircularAt(C, v + 2) == -1) {
            V[0] = v;
            V[1] = Circular(C, v + 1);
            V[2] = Circular(C, v + 2);
            return 3;
        }

		// only two possible moves
		if (CircularAt(C, v - 1) == -1) {
			V[0] = Circular(C, v - 1);
			V[1] = v;
			return 2;
		}

		if (CircularAt(C, v + 1) == -1) {
			V[0] = v;
			V[1] = Circular(C, v + 1);
			return 2;
		}

        V[0] = v;
        V[1] = Circular(C, v - 1);
        V[2] = Circular(C, v + 1);
        return 3;
	}
}

// Are we introducing an edge which could not have existed ?
bool test_edge_accuracy(
	const mat2x& B,
	const int vsrc,
	const int vdst
) {
	mat2 e;
	e.col(0) = B.col(vsrc);
	e.col(1) = B.col(vdst);

	const vec2 d_error = PathUtils::distance_bounds_from_points(B, e, vec2i(vsrc, vdst));

	PF_VERBOSE_F("test edge %d - %d accuracy %f %f", vsrc, vdst, d_error.minCoeff(), d_error.maxCoeff());

	return ErrorMetrics::accuracy_within_bounds(d_error.minCoeff(), d_error.maxCoeff(), e.col(0) - e.col(1));
}

// Are we breaking existing parallel edges ?
bool test_edge_alignment(
	const vec2 e0_old,
	const vec2 e0_new,
	const vec2 e1_old,
	const vec2 e1_new
) {
	const double dev0_old = AngleUtils::deviation_from_horizontal(e0_old);
	const double dev0_new = AngleUtils::deviation_from_horizontal(e0_new);
	const double dev1_old = AngleUtils::deviation_from_horizontal(e1_old);
	const double dev1_new = AngleUtils::deviation_from_horizontal(e1_new);

    PF_VERBOSE_F("dev0_old %f dev0_new %f dev1_old %f dev1_new %f", dev0_old, dev0_new, dev1_old, dev1_new);

	const double diff_old = abs(dev0_old - dev1_old);
	const double diff_new = abs(dev0_new - dev1_new);
	
    PF_VERBOSE_F("diff_old %f diff_new %f", diff_old, diff_new);

	if (diff_old < PF_EPS_MEDIUM && diff_new > PF_EPS_MEDIUM) {
		return false;
	} 

	if (diff_old )
	
	return true;
}

// Does it even make sense to attempt a regularization if the two edges
// are completely misaligned ?
bool test_opposite_corners_false_positive (
	const vec2 p0_prev,
	const vec2 p0,
	const vec2 p0_next,
	const vec2 p1_prev,
	const vec2 p1,
	const vec2 p1_next
) {
	const bool is_inflection = AngleUtils::have_opposite_convexity(p0_prev, p0, p1, p1_next);
	const double angle0 = AngleUtils::spanned_shortest(p0, p1, p1_next);
	const double angle1 = AngleUtils::spanned_shortest(p0_prev, p0, p1);
	const double angle_diff = abs(angle0 - angle1);

    PF_VERBOSE_F("is inflection %d angle0 %f angle1 %f angle difference %f", is_inflection, PF_DEG(angle0), PF_DEG(angle1), PF_DEG(angle_diff));

	// We expect to match closures, not spurious concave corners which are nearly flat
	const double angle0p = AngleUtils::spanned_shortest(p0_prev, p0, p0_next);
	const double angle1p = AngleUtils::spanned_shortest(p1_prev, p1, p1_next);
	if (angle0p > PF_RAD(135) || angle1p > PF_RAD(135)) {
		return false;
	}

	// We should match things which are already roughly aligned
	if (angle_diff > PF_RAD(50)) {
		return false;
	}

	// Let's also not break symmetries
	const double dev00 = AngleUtils::deviation_from_horizontal((p0_prev - p0).normalized());
	const double dev01 = AngleUtils::deviation_from_horizontal((p0_next - p0).normalized());
	const double dev10 = AngleUtils::deviation_from_horizontal((p1_prev - p1).normalized());
	const double dev11 = AngleUtils::deviation_from_horizontal((p1_next - p1).normalized());
	PF_VERBOSE_F("test symmetry dev00 %f dev01 %f dev10 %f dev11 %f", PF_DEG(dev00), PF_DEG(dev01), PF_DEG(dev10), PF_DEG(dev11));
	if (abs(dev00 - dev01) < PF_EPS || abs(dev10 - dev11) < PF_EPS) {
		return false;
	}

	// Also, if they are already on the same axis, let's not break it
	if (abs(p0(0) - p1(0)) < PF_EPS || abs(p0(1) - p1(1)) < PF_EPS) {
		return false;
	}

	return true;
}

// Returns true if the move is between two vertices which are only symmetric w.r.t. to 
// themselves. If any of the two vertices is symmetric w.r.t. any  other vertex than the move is
// not allowed.
// EXCEPT if the two vertices which v0 and v1 are symmetric to are also symmetric with respect to one another.
// Technically v0 and v1 could be symmetric to more than one set of vertices, we do not
// currently consider this scenario.
// ----------------------------------------------------------------------------
bool would_move_break_symmetry(
    const mat2x& B,
    const vecXi& V,
    const int v0,
    const int v1,
    const polyfit::Regularity::RegularityInformation& RE,
    const bool circular
) {
    int v0_other_symmetry = -1;
    int v1_other_symmetry = -1;

    for(auto& sym : RE.vertex_symmetries()) {        

        // Symmetry match does not contain this two corners
        if (!sym.matches(v0) && !sym.matches(v1)) {
            continue; 
        }

        // Points are symmetric w.r.t. each other
        if (sym.matches(v0, v1)) {
            continue;
        }

        if (v0_other_symmetry < 0 && sym.matches(v0)) {
            v0_other_symmetry = sym.v0 == v0 ? v1 : v0;
            continue;
        }

        if (v1_other_symmetry < 0 && sym.matches(v1)) {
            v1_other_symmetry = sym.v0 == v1 ? v1 : v0;
            continue;
        }
          
        // This is the scenario not handled mentioned in the comments
        // PF_ABORT; 
    }

    // Good!
    if (v0_other_symmetry == -1 && v1_other_symmetry == -1) {
        return false;
    }

    // Are the other matches symmetric w.r.t. each other. We are assuming that the move on v0 and v1 will
    // preserve this symmetry.
    for (auto& sym : RE.vertex_symmetries()) {        
        if (sym.matches(v0_other_symmetry, v1_other_symmetry)) {
            break;
        }

        if (sym.matches(v0_other_symmetry) || sym.matches(v1_other_symmetry)) {
            return true;
        }
    }

    return false;
}

void align_opposite_corners(
	const mat2x& B,
	std::vector<BoundaryGraph::Edge>& E,
	polyfit::Regularity::RegularityInformation&    RE,
    vecXi& V,
	const bool circular
) {
	PF_LOG_DIVIDER;
	PF_VERBOSE_S("align opposite corners");

	mat2x P;
	BoundaryGraph::trace_to_points(B, V, P);

	// please.. mat2xi
	std::vector<vec2i> candidates;

	// ehi you... the compute convexities is bugged
	std::vector<int> C;
	PathUtils::compute_convexities(P, C);
	
	for (Index i = 0; i < V.size(); ++i) {
		double d_min_sq = INFINITY;
		Index j_min = -1;

        // Distance of 2
		for (Index dist = 2; dist < V.size() - 1; ++dist) {
            const Index j = Circular(V, i + dist);
			// considering only inflections
			if (C[i] != -1 || C[j] != -1) {
				continue;
			}

			if (Regularity::check_opposite_points_winding_number(P, i, j)) {
				const double d_sq = (P.col(i) - P.col(j)).squaredNorm();

                // This is not the ideal way to handle this scenario, but in practice it works well.
                // We don't want to necessarily match the closest corners, but rather the most continuous ones
                // this would require reworking 5this main loop.
                // If two corners have roughly the same distance, let's pick the sharpest one
                // THIS IS AWFUL.
                const double delta_d_sq = abs(d_sq - d_min_sq);

                if (d_sq < d_min_sq && delta_d_sq >= 2.) {
					d_min_sq = d_sq;
					j_min = j;
				} 
                else if (delta_d_sq < 2.) {
                    const double angle0 = AngleUtils::spanned_shortest(
                        CircularAt(P, j_min - 1), CircularAt(P, j_min), CircularAt(P, j_min + 1)
                    );

                    const double angle1 = AngleUtils::spanned_shortest(
                        CircularAt(P, j - 1), CircularAt(P, j), CircularAt(P, j + 1)
                    );

                    if (angle1 < angle0) {
                        d_min_sq = d_sq;
                        j_min = j;
                    }
                }

				PF_VERBOSE_F("opposite corners %d - %d distance %f", V(i), V(j), d_sq);
			}
		}

		if (j_min != -1) {
			// Making sure that the candidate pairs are unique in both directions..
			// this could be done above
			bool append = true;
			for (size_t j = 0; j < candidates.size(); ++j) {
				PF_VERBOSE_F("compare %d %d - %d %d",
					V(candidates[j](0)), V(candidates[j](1)), V(i), V(j_min));
				if (candidates[j](1) == j_min || candidates[j](0) == j_min || candidates[j](1) == i) {
					const double d_other_sq = (P.col(candidates[j](0)) - P.col(candidates[j](1))).squaredNorm();

					if (d_min_sq < d_other_sq) {
						candidates[j] = vec2i(i, j_min);
					}

					append = false;
					break;
				}
			}

			if (append) {
				candidates.emplace_back(vec2i(i, j_min));
			}
		}
	}

	std::vector<int> C_B; 
    PathUtils::compute_convexities(B, C_B);

	// If you think this is a bunch of hacks, you are clearly not familiar with perceptual flow control
	for (Index i = 0; i < candidates.size(); ++i) {
		const int c0 = candidates[i](0);
		const int c1 = candidates[i](1);
		const int v0 = V(c0);
		const int v1 = V(c1);
		PF_LOG_DIVIDER;
		PF_VERBOSE_F("candidate %d(%d) %d(%d)", c0, v0, c1, v1);

        if (would_move_break_symmetry(B, V, c0, c1, RE, circular)) {
            continue;
        }

		// We only want to align corners which are ambiguous, that is they are in a 
		// 1 pixel concavity |_| where there could be a midpoint.
        int V0[3] = { -1, -1, -1 };
		int V0_count = enumerate_moves(B, P, C_B, v0, c0, V0);

        int V1[3] = { -1, -1, -1 };
		int V1_count = enumerate_moves(B, P, C_B, v1, c1, V1);

		// no candidates moves
		if (V0_count == 1 && V1_count == 1) {
			PF_VERBOSE_S("no candidate moves");
			continue;
		}

		PF_VERBOSE_F("V0 candidates %d", V0_count);
		PF_VERBOSE_F("V1 candidates %d", V1_count);

		// We need to decide which consecutive edges should be considered.
		// If we have a nearly flat edge, let's look up the point before
		const double angle_prev_opt0 = AngleUtils::spanned_shortest(B.col(CircularAt(V, c0 - 2)), B.col(CircularAt(V, c0 - 1)), B.col(v0));
		const double angle_prev_opt1 = AngleUtils::spanned_shortest(B.col(v0), B.col(CircularAt(V, c0 + 1)), B.col(CircularAt(V, c0 + 2)));
		const double angle_next_opt0 = AngleUtils::spanned_shortest(B.col(v1), B.col(CircularAt(V, c1 + 1)), B.col(CircularAt(V, c1 + 2)));
		const double angle_next_opt1 = AngleUtils::spanned_shortest(B.col(CircularAt(V, c1 - 2)), B.col(CircularAt(V, c1 - 1)), B.col(v1));

		const double angle_flat_thr = M_PI - PF_RAD(15);
		const bool angle_prev_opt0_is_flat = angle_prev_opt0 > angle_flat_thr;
		const bool angle_prev_opt1_is_flat = angle_prev_opt1 > angle_flat_thr;
		const bool angle_next_opt0_is_flat = angle_next_opt0 > angle_flat_thr;
		const bool angle_next_opt1_is_flat = angle_next_opt1 > angle_flat_thr;

		const bool d_prev_opt0_is_close = ShortestDist(B, CircularAt(V, c0 - 1), V(c0)) <= 2;
		const bool d_prev_opt1_is_close = ShortestDist(B, V(c0), CircularAt(V, c0 + 1)) <= 2;
		const bool d_next_opt0_is_close = ShortestDist(B, CircularAt(V, c1 + 1), V(c1)) <= 2;
		const bool d_next_opt1_is_close = ShortestDist(B, V(c1), CircularAt(V, c1 - 1)) <= 2;

		int d_prev_opt0 = (angle_prev_opt0_is_flat && d_prev_opt0_is_close) ? 2 : 1;
		int d_prev_opt1 = (angle_prev_opt1_is_flat && d_prev_opt1_is_close) ? 2 : 1;
		int d_next_opt0 = (angle_next_opt0_is_flat && d_next_opt0_is_close) ? 2 : 1;
		int d_next_opt1 = (angle_next_opt1_is_flat && d_next_opt1_is_close) ? 2 : 1;
        PF_VERBOSE_F("d_prev_opt0 %d d_next_opt0 %d", d_prev_opt0, d_next_opt0);
        PF_VERBOSE_F("d_prev_opt1 %d d_next_opt1 %d", d_prev_opt1, d_next_opt1);
		
        PF_VERBOSE_F("combination 0: %d - %d - %d - %d", CircularAt(V, c0 - d_prev_opt0), v0, v1, CircularAt(V, c1 + d_next_opt0));
        PF_VERBOSE_F("combination 1: %d - %d - %d - %d", CircularAt(V, c0 + d_prev_opt1), v0, v1, CircularAt(V, c1 - d_next_opt1));

		const vec2 p0_prev_opt0 = B.col(CircularAt(V, c0 - d_prev_opt0));
		const vec2 p0_prev_opt1 = B.col(CircularAt(V, c0 + d_prev_opt1));
		const vec2 p1_next_opt0 = B.col(CircularAt(V, c1 + d_next_opt0));
		const vec2 p1_next_opt1 = B.col(CircularAt(V, c1 - d_next_opt1));

		vec2 p0 = B.col(v0);
		vec2 p1 = B.col(v1);

		const double dev_reference = AngleUtils::deviation_from_horizontal((p1 - p0).normalized());
        PF_VERBOSE_F("dev_reference %f", PF_DEG(dev_reference));

		const double dev0_opt0 = AngleUtils::spanned_shortest(p0_prev_opt0, p0, p1);
		const double dev1_opt0 = AngleUtils::spanned_shortest(p0, p1, p1_next_opt0);
		const double dev0_opt1 = AngleUtils::spanned_shortest(p0_prev_opt1, p0, p1);
		const double dev1_opt1 = AngleUtils::spanned_shortest(p0, p1, p1_next_opt1);

        PF_VERBOSE_F("dev0_opt0 %f dev1_opt0 %f", PF_DEG(dev0_opt0), PF_DEG(dev1_opt0));
        PF_VERBOSE_F("dev0_opt1 %f dev1_opt1 %f", PF_DEG(dev0_opt1), PF_DEG(dev1_opt1));

		//const double dev_opt0_max = max(abs(dev_reference - dev0_opt0), abs(dev_reference - dev1_opt0));
		//const double dev_opt1_max = max(abs(dev_reference - dev0_opt1), abs(dev_reference - dev1_opt1));
		const double dev_opt0 = dev0_opt0 + dev1_opt0;
		const double dev_opt1 = dev0_opt1 + dev1_opt1;
        PF_VERBOSE_F("dev_opt0 %f dev_opt1 %f", dev_opt0, dev_opt1);

		bool use_opt0 = dev_opt0 > dev_opt1;
		const vec2 p0_prev = use_opt0 ? p0_prev_opt0 : p0_prev_opt1;
		const vec2 p0_next = use_opt0 ? p0_prev_opt1 : p0_prev_opt0;
		const vec2 p1_next = use_opt0 ? p1_next_opt0 : p1_next_opt1;
		const vec2 p1_prev = use_opt0 ? p1_next_opt1 : p1_next_opt0;
        PF_VERBOSE_F("use_opt0 %d", use_opt0);

		// Does it even make sense to attempt this regularization ?
		if (!test_opposite_corners_false_positive(p0_prev, p0, p0_next, p1_prev, p1, p1_next)) {
            PF_VERBOSE_S("discarding false positive");
			continue;
		}

		// picking the best combination
		int v0_min = -1;
		int v1_min = -1;
		double error_min = INFINITY;

		for (int j = 0; j < V0_count; ++j) {
			for (int k = 0; k < V1_count; ++k) {
				const int v0j = V0[j];
				const int v1k = V1[k];
				vec2 p0 = B.col(v0j);
				vec2 p1 = B.col(v1k);

				const double angle0 = AngleUtils::spanned_shortest(p0_prev, p0, p1);
				const double angle1 = AngleUtils::spanned_shortest(p0, p1, p1_next);

				double error = INFINITY;

				// If the result is axis-aligned, it's perfect!!!
				const bool same_axis = abs(p1(1) - p0(1)) < PF_EPS || abs(p1(0) - p0(0)) < PF_EPS;
				const bool axis_aligned0 = abs(p0_prev(0) - p0(0)) < PF_EPS || abs(p0_prev(1) - p0(1)) < PF_EPS;
				const bool axis_aligned1 = abs(p1(0) - p1_next(0)) < PF_EPS || abs(p1(1) - p1_next(1)) < PF_EPS;

				if (same_axis && (axis_aligned0 || axis_aligned1)) {
					error = 0.;
				}
				// otherwise let's go with optimal continuity
				// which means that in the case of inflections we prefer having the
				// highest minimum angle possible. 
				// In case the two opposite points have the same convexity, we try	
				// to minimize the difference between the two.
				// Generally, any solution with no inflections is better than one with inflections
				else {
					const bool has_inflection = AngleUtils::have_opposite_convexity(p0_prev, p0, p1, p1_next);
                    const double curved_penalty = 1.;
                    const double inflection_penalty = 1.;

					if (has_inflection) {
						error = curved_penalty + inflection_penalty + (M_PI - min(angle0, angle1)) / M_PI;
					} else {
                        // We always prefer straight lines to curved ones
                        const double angle_difference = abs(angle0 - angle1);
                        if (angle_difference < PF_EPS && abs(angle0 - M_PI) < PF_EPS && abs(angle1 - M_PI) < PF_EPS) {
                            error = 0.;
                        } else {
                            error = curved_penalty + angle_difference / M_PI;
                        }
					}
				}

				if (error < error_min) {
					error_min = error;
					v0_min = v0j;
					v1_min = v1k;
				}
                PF_VERBOSE_F("pair %d - %d error %f", v0j, v1k, error);
			}
		}

		// checking if the new edges that would be inserted are valid
		const bool test_exists_e00 = test_edge_accuracy(B, CircularAt(V, c0 - 1), v0_min);
		const bool test_exists_e01 = test_edge_accuracy(B, v0_min, CircularAt(V, c0 + 1));
		const bool test_exists_e10 = test_edge_accuracy(B, CircularAt(V, c1 - 1), v1_min);
		const bool test_exists_e11 = test_edge_accuracy(B, v1_min, CircularAt(V, c1 + 1));

		if (!test_exists_e00 || !test_exists_e01 || !test_exists_e10 || !test_exists_e11) {
            PF_VERBOSE_F("fail accuracy %d - %d", v0_min, v1_min);
			continue;
		}

		// are we breaking a previous alignment ?
		// mamma mia... can we not normalize() what about detecting the regularities *before*?
		const vec2 dir00_old = (B.col(V(c0)) - B.col(CircularAt(V, c0 - 1))).normalized();
		const vec2 dir00_new = (B.col(v0_min) - B.col(CircularAt(V, c0 - 1))).normalized();
		const vec2 dir10_old = (B.col(CircularAt(V, c1 + 1)) - B.col(V(c1))).normalized();
		const vec2 dir10_new = (B.col(CircularAt(V, c1 + 1)) - B.col(v1_min)).normalized();

		const vec2 dir01_old = (B.col(CircularAt(V, c0 + 1)) - B.col(V(c0))).normalized();
		const vec2 dir01_new = (B.col(CircularAt(V, c0 + 1)) - B.col(v0_min)).normalized();
		const vec2 dir11_old = (B.col(V(c1)) - B.col(CircularAt(V, c1 - 1))).normalized();
		const vec2 dir11_new = (B.col(v1_min) - B.col(CircularAt(V, c1 - 1))).normalized();

		const bool test_alignment_e0 = test_edge_alignment(dir00_old, dir00_new, dir10_old, dir10_new);
		const bool test_alignment_e1 = test_edge_alignment(dir01_old, dir01_new, dir11_old, dir11_new);

		if (!test_alignment_e0 || !test_alignment_e1) {
            PF_VERBOSE_F("fail alignment %d - %d", v0_min, v1_min);
			continue;
		}

		// relocating the two corners 
		PF_ASSERT(v0_min != -1 && v1_min != -1);
		V(c0) = v0_min;
		V(c1) = v1_min;

		// Have we introduced a degenerate edge which should be removed ?
		if (use_opt0 && d_prev_opt0 == 2 && CircularDist(V, Circular(V, c0 - 1), c0) == 1) {
			EraseOrdered(V, { Circular(V, c0 - 1) });
		}

		if (use_opt0 && d_next_opt0 == 2 && CircularDist(V, c1, Circular(V, c1 + 1)) == 1) {
			EraseOrdered(V, { Circular(V, c1 + 1) });
		}
		
		if (!use_opt0 && d_prev_opt1 == 2 && CircularDist(V, c0, Circular(V, c0 + 1)) == 1) {
			EraseOrdered(V, { Circular(V, c0 + 1) });
		}
		
		if (!use_opt0 && d_next_opt1 == 2 && CircularDist(V, Circular(V, c1 - 1), c1) == 1) {
			EraseOrdered(V, { Circular(V, c1 - 1) });
		}
	}

	return;

	// Removing any degenerate edge which we introduced in the previous apss
	int removed;
	do {
		removed = 0;

		for (Index i = 0; i < V.size(); ++i) {
			if (!circular && i == V.size() - 1) {
				continue;
			}

			const int i0 = i;
			const int i1 = Circular(V, i + 1);
			const double len_sq = (B.col(V(i0)) - B.col(V(i1))).squaredNorm();

			if (len_sq < 1. + PF_EPS) {
				PF_VERBOSE_F("test vertex %d (length^2 %f) against candidates", i0, len_sq);

				for (size_t j = 0; j < candidates.size(); ++j) {
					// The candidate corner can only be one of the two
					PF_ASSERT(!((candidates[j](0) == i0 && candidates[j](1) == i1) ||
						        (candidates[j](0) == i1 && candidates[j](1) == i0)));

					if (candidates[j](0) == i0 || candidates[j](1) == i0) {
						EraseOrdered(V, Circular(V, i + 1));
						++removed;
						break;
					}
					
					if (candidates[j](0) == i1 || candidates[j](1) == i1) {
						EraseOrdered(V, i);
						++removed;
						break;
					}
				}
			}
		}
	} while (removed > 0);
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)