// polyvec
#include <polyvec/polygon-tracer/pick_optimal_symmetry_side.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/utils/matrix.hpp>

// libc++
#include <algorithm>

using namespace std;

#define LOG_SUBPATH_EDGES 0
#define LOG_SUBPATH_ENERGIES 1

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int pick_optimal_symmetry_side(
	const mat2x& B,          // raster boundary
	const vecXi& P,          // polygonal approximation
	const std::vector<BoundaryGraph::Edge>& E, // edges
	const std::vector<Vertex> PE, // polygon edges
	const vec2i& sym0,       // (first, last) points in the raster for symmetry 0
	const vec2i& sym1,       // (first, last) points in the raster for symmetry 1
	std::vector<Vertex>& S0, // points of P contained in sym0
	std::vector<Vertex>& S1, // points of P contained in sym1
	const bool circular
) {
	PF_VERBOSE_F("pick optimal symmetric side sym0(%d %d) sym1(%d %d)\n",
		sym0(0), sym0(1), sym1(0), sym1(1));

	S0.clear();
	S1.clear();

	// todo: to avoid the extra allocation, S0 and S1 should be used to store the indices 
	// with respect to the path and convertex to the raster indices just before returning.
	std::vector<Vertex> Psub0, Psub1;

	for (int j = 0; j < P.size(); ++j) {
		if (PathUtils::contains_closed(B, sym0(0), sym0(1), P(j), circular)) {
			S0.emplace_back(P(j));
		}

		if (PathUtils::contains_closed(B, sym1(0), sym1(1), P(j), circular)) {
			S1.emplace_back(P(j));
		}

		const BoundaryGraph::Edge& e = E[PE[j]];

		if (PathUtils::contains_closed(B, sym0(0), sym0(1), e.v0, circular) &&
			PathUtils::contains_closed(B, sym0(0), sym0(1), e.v1, circular)) {
			Psub0.emplace_back(j);
		}

		if (PathUtils::contains_closed(B, sym1(0), sym1(1), e.v0, circular) &&
			PathUtils::contains_closed(B, sym1(0), sym1(1), e.v1, circular)) {
			Psub1.emplace_back(j);
		}
	}

	// nothing to do here as both edges are starting outside the symmetric region
	if (S0.size() < 2 && S1.size() < 2) {
		return 0;
	}

	// To decide which one of the two sides is more symmetric, look at the number of vertices
	// in both paths. 
	// If they differ, it's a good indicator for smoothness and/or regularities
	if (S0.size() != S1.size()) {
		//return S0.size() > S1.size() ? -1 : +1;
	}

	// This is horrendous, but in order to do this correctly I need to think a little bit
	// more how to go about it and it __literally__ makes a difference only in one model
	double E0 = 0.;
	for (int i = 0; i < Psub0.size(); ++i) {
		// the piece of code below won't work with non-circular paths
		assert(circular);
		
		const Vertex Ei_prev = CircularAt(PE, Psub0[i] - 1);
		const Vertex Ei = PE[Psub0[i]];
		const Vertex Ei_next = CircularAt(PE, Psub0[i] + 1);

		const BoundaryGraph::Edge e_prev = E[Ei_prev];
		const BoundaryGraph::Edge e = E[Ei];
		const BoundaryGraph::Edge e_next = E[Ei_next];

#if LOG_SUBPATH_EDGES
		PF_LOG_DIVIDER;
		PF_LOGF("e_prev v0 %lld v1 %lld", e_prev.v0, e_prev.v1);
		PF_LOGF("e v0 %lld v1 %lld", e.v0, e.v1);
		PF_LOGF("e_next v0 %lld v1 %lld", e_next.v0, e_next.v1);
#endif

		E0 += ErrorMetrics::calculate_inner(
			e.d_min, e.d_max,
			B.col(e_prev.v0), B.col(e_prev.v1), B.col(e_next.v0), B.col(e_next.v1)
		);
	}

#if LOG_SUBPATH_EDGES
	PF_LOG_DIVIDER;
	PF_LOG_DIVIDER;
#endif

	double E1 = 0.;
	for (int i = 0; i < Psub1.size(); ++i) {
		assert(circular);

		const Vertex Ei_prev = CircularAt(PE, Psub1[i] - 1);
		const Vertex Ei = PE[Psub1[i]];
		const Vertex Ei_next = CircularAt(PE, Psub1[i] + 1);

		const BoundaryGraph::Edge e_prev = E[Ei_prev];
		const BoundaryGraph::Edge e = E[Ei];
		const BoundaryGraph::Edge e_next = E[Ei_next];

#if LOG_SUBPATH_EDGES
		PF_LOG_DIVIDER;
		PF_LOGF("e_prev v0 %lld v1 %lld", e_prev.v0, e_prev.v1);
		PF_LOGF("e v0 %lld v1 %lld", e.v0, e.v1);
		PF_LOGF("e_next v0 %lld v1 %lld", e_next.v0, e_next.v1);
#endif

		E1 += ErrorMetrics::calculate_inner(
			e.d_min, e.d_max,
			B.col(e_prev.v0), B.col(e_prev.v1), B.col(e_next.v0), B.col(e_next.v1)
		);
	}

#if LOG_SUBPATH_ENERGIES
	PF_LOG_DIVIDER;
	PF_VERBOSE_F("subpath 0 energy %f", E0);
    PF_VERBOSE_F("subpath 1 energy %f", E1);
	PF_LOG_DIVIDER;
#endif

    // Any of the two paths will do if they are exactly the same
	if (E0 == E1) {
		return 1;
	}

	return E0 < E1 ? -1 : +1;
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)