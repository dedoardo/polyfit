// polyvec
#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/is_symmetric_path_valid.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/utils/matrix.hpp> 

// libc++
#include <cstdlib> // min, max

#define SIMPLER_VALID_LOGIC 1

#define TEST_OUTER_INFLECTIONS 1
#define TEST_CONCAVE_BOUNDARY_INFLECTION 1

// defines the length of consecutive flat boundaries which are a clear inflection
// in the raster. This shouldn't be changed.
#define MIN_LENGTH_INFLECTION 2. 
#define MIN_LENGTH_INFLECTION_SQ (MIN_LENGTH_INFLECTION * MIN_LENGTH_INFLECTION) - PF_EPS

using namespace std;
using namespace Eigen;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

bool is_symmetric_path_valid(
	const mat2x& B,                   
	const vecXi& P_old,
	const vecXi& P_new,
	const std::vector<Vertex>& S_ord, 
	const bool circular
) {
	// well, nothing to do here
	if (!P_new.size()) {
		return false;
	}
	
	bool is_fit_valid = true;
	
	// todo: this are redundant but a couple of methods work directly with a list
	// of points
	mat2x PP_old, PP_new;

	BoundaryGraph::trace_to_points(B, P_old, PP_old);
	BoundaryGraph::trace_to_points(B, P_new, PP_new);

	// todo: please... get rid of this allocations
	vector<int> C_old, C_new;
	PathUtils::compute_convexities(PP_old, C_old);
	PathUtils::compute_convexities(PP_new, C_new);

	// Counting the number of inflections in the older path
	const int n_inflections_old = AngleUtils::count_inflections(PP_old, circular);
	const int n_inflections_new = AngleUtils::count_inflections(PP_new, circular);
	PF_VERBOSE_F("test inflections old %d new %d", n_inflections_old, n_inflections_new);
	
	if (n_inflections_new > n_inflections_old) {
		is_fit_valid = false;
	}

	// If the inflection happens inside the symmetric region, we should permit it.
	// todo: this handles only the first and last corner for disjoint symmetries
	if (!is_fit_valid) {
#if !SIMPLER_VALID_LOGIC
		const Vertex v_beg = FindIndex(P_new, (int)S_ord.front());
		const Vertex v_end = FindIndex(P_new, (int)S_ord.back());
		PF_ASSERT(v_beg >= 0 && v_end >= 0);

				// Checking the raster at the two opposite corners
		const double a_beg = PathUtils::angle_at_corner(B, P_new, v_beg, circular);
		const double a_end = PathUtils::angle_at_corner(B, P_new, v_end, circular);

		// attempting to pair two opposite 90 degree concave corners
		// if the symmetry is not disjoint this is tipically a clear closure
		const bool same_convexity = C_new[v_beg] == C_new[v_end];
		const bool is_angle_90 = fabsf(a_beg - M_PI_2) < PF_EPS && fabsf(a_end - M_PI_2) < PF_EPS;

		// The two segments should at least be 90 degrees
		bool is_length_noise = false;

		if (!circular && v_beg == 0) {
			is_length_noise = true;
		}
		else if (!circular && v_end == P_new.size() - 1) {
			is_length_noise = true;
		}
		else {
			const double len00_sq = (CircularAt(PP_new, v_beg - 1) - CircularAt(PP_new, v_beg)).squaredNorm();
			const double len01_sq = (CircularAt(PP_new, v_beg + 1) - CircularAt(PP_new, v_beg)).squaredNorm();
			const double len10_sq = (CircularAt(PP_new, v_end - 1) - CircularAt(PP_new, v_end)).squaredNorm();
			const double len11_sq = (CircularAt(PP_new, v_end + 1) - CircularAt(PP_new, v_end)).squaredNorm();
			const double len0_sq = min(len00_sq, len01_sq);
			const double len1_sq = min(len10_sq, len11_sq);

			if (len0_sq < MIN_LENGTH_INFLECTION_SQ || len1_sq < MIN_LENGTH_INFLECTION_SQ) {
				is_length_noise = true;
			}
		}

		PF_LOGF("test inflection allow same-convexity %d is-angle-90 %d is-length-noise %d", same_convexity, is_angle_90, is_length_noise);

		if (same_convexity && is_angle_90 && !is_length_noise) {
			is_fit_valid = true;
		}
#endif
	}

#if TEST_OUTER_INFLECTIONS
	// This test is fundamentally the same as the one above as it aims to prevent the regularization
	// to introduce inflections. 
	// It is separated from the one which counts the inflections because inflections are currently
	// defined on edges. In cases where two+ consecutive edges are already inflections adding another concave corner
	// between the endpoints won't increase the inflection count while it will be visually disruptive.
	// Since the only points which would have been affected by the symmetry are the ones before after the symmetric region,
	// we check that the convexity on those has been preserved after the regularization.
	{
	const Index v_prev_old = PathUtils::find_closest_corner_before_strict(B, P_old, S_ord.front());
	const Index v_prev_new = PathUtils::find_closest_corner_before_strict(B, P_new, S_ord.front());
	const Index v_next_old = PathUtils::find_closest_corner_after_strict(B, P_old, S_ord.back());
	const Index v_next_new = PathUtils::find_closest_corner_after_strict(B, P_new, S_ord.back());
	if (C_old[v_prev_old] != C_new[v_prev_new] || C_old[v_next_old] != C_new[v_next_new]) {
		is_fit_valid = false;
	}
	}

	// The test before is more result-driven to be honest. The number of concave corners is not necessarily
	// increased, but they are relocated in such a way that it visually hurts.
	// I think that some of this tests could be tested more effectively, but we need to handle a similar
	// scenario as above where there is no relocation but a concave corner is introduced between an original edge
	{
		Index test_boundary_corners[] = {
			S_ord.front(),
			S_ord.back()
		};

		for (int i = 0; is_fit_valid && i < PV_LEN(test_boundary_corners); ++i) {
			const Index test_boundary_corner = test_boundary_corners[i];

			const Index v_prev_old = PathUtils::find_closest_corner_before_strict(B, P_old, test_boundary_corner);
			const Index v_next_old = PathUtils::find_closest_corner_after_strict(B, P_old, test_boundary_corner);
			const Index v_prev_new = FindIndex(P_new, P_old(v_prev_old));
			const Index v_next_new = FindIndex(P_new, P_old(v_next_old));
			if (v_prev_new == -1 || v_next_new == -1) {
				continue;
			}

			const Index d_old = CircularDist(P_old, v_prev_old, v_next_old);
			const Index d_new = CircularDist(P_new, v_prev_new, v_next_new);

			if (d_old && d_new > d_old) {
				const int dir = DirectionBetween(P_new, v_prev_new, v_next_new);

				for (Index j = 1; j < d_new; ++j) {
					const Index v_new = Circular(P_new, v_prev_new + j * dir);

					// Skipping points inside the symmetric region. This could be tested more nicely 
					// with distances rather than a linear search (todo), but it's not too bad as we
					// are only looking for new points inserted before the end of the
					// symmetric region (tipically just 1).
					if (FindIndex(S_ord, (Vertex)P_new(v_new)) == -1 && C_new[v_new] == -1) {
						is_fit_valid = false;
						break;
					}
				}
			}

#if TEST_CONCAVE_BOUNDARY_INFLECTION
			const Index v_boundary_old = FindIndex(P_old, (int)test_boundary_corner);
			const Index v_boundary_new = FindIndex(P_new, (int)test_boundary_corner);
			//fprintf(stderr, "found path vertex %lld corner %d\n", v_boundary_new, P_new(v_boundary_new));

			if (v_boundary_old == -1 && v_boundary_new != -1 && C_new[v_boundary_new] == -1) {
				PF_ASSERT(P_new(v_boundary_new) == test_boundary_corner);
				is_fit_valid = false;
				break;
			}
#endif
		}
	}

#endif

	// Regardless of whether the fit has failed, we need to make sure that no degenerate edges have been introduced
	// outside the symmetric region
	const int n_degenerate_old = AngleUtils::count_visually_degenerate_edges(B, P_old, circular);
	const int n_degenerate_new = AngleUtils::count_visually_degenerate_edges(B, P_new, circular);

    PF_VERBOSE_F("comparing degenerate edges old %d new %d", n_degenerate_old, n_degenerate_new);

	if (n_degenerate_new > n_degenerate_old) {
		PF_VERBOSE_S ("failing symmetry");
		is_fit_valid = false;
	}

	return is_fit_valid;
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)