// polyvec
#include <polyvec/polygon-tracer/iterative-global.hpp>
#include <polyvec/polygon-tracer/iterative-local.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/polygon-tracer/is_symmetric_path_valid.hpp>
#include <polyvec/polygon-tracer/pick_optimal_symmetry_side.hpp>
#include <polyvec/polygon-tracer/align_inflections_to_half_space.hpp>
#include <polyvec/polygon-tracer/remove_edges_crossing_implicit_corners.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/core/log.hpp>

// libc++
#include <algorithm>
#include <numeric>

#define PREVENT_SYMMETRY_OVERLAP 0

// logging
#define LOG_SYMMETRIC_REGIONS 0
#define LOG_TRACE_ITERATIONS 0
#define LOG_EDGE_DELETE_LISTS 0

#define HANDLE_NESTED_SYMMETRIES 1
#define ORDER_BIG_TO_SMALL 1

using namespace std;
using namespace polyvec;
using namespace Eigen;

#define ON_ATTEMPT ++attempts; if (on_iter) { on_iter(B, P, E); }

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void iterative_global_on_attempt_default(const mat2x&, const std::vector<BoundaryGraph::Edge>&, const int) { }
void iterative_global_on_result_default(const mat2x&, const vecXi&, const int) { }

bool is_symmetry_overlapping_remaining(
	const mat2x& B,
	const std::vector<Symmetry::SubPath>& S,
	const size_t i,
	const bool circular
) {
	const auto& s_test = S[i];

	for (size_t i = 0; i < S.size(); ++i) {
		const auto& s = S[i];

		// skipping entire model symmetries
		const bool joint0 = Symmetry::is_subpath_joint0(B, s, circular);
		const bool joint1 = Symmetry::is_subpath_joint1(B, s, circular);

		if (joint0 && joint1) {
			continue;
		}

		if (PathUtils::contains_open(B, s.v01, s.v00, s_test.v00, circular) ||
			PathUtils::contains_open(B, s.v01, s.v00, s_test.v10, circular) ||
			PathUtils::contains_open(B, s.v10, s.v11, s_test.v00, circular) ||
			PathUtils::contains_open(B, s.v10, s.v11, s_test.v10, circular)) {
			return true;
		}
	}

	return false;
}

// TIntervalFunctor: std::pair<Vertex, Vertex>(const Symmetry::SubPath&)
template <typename TIntervalFunctor>
void add_nested_symmetries_for_interval(const mat2x& B, const std::vector<Symmetry::SubPath>& raster_symmetries, bool circular, TIntervalFunctor&& interval_functor, std::vector<std::vector<int>>& nested_symmetries)
{
	std::vector<int> symmetry_order_by_begin(raster_symmetries.size());
	std::iota(symmetry_order_by_begin.begin(), symmetry_order_by_begin.end(), 0);
	std::sort(symmetry_order_by_begin.begin(), symmetry_order_by_begin.end(), [&](int i1, int i2) 
	{
		auto intvl1 = std::forward<TIntervalFunctor>(interval_functor)(raster_symmetries[i1]);
		auto intvl2 = std::forward<TIntervalFunctor>(interval_functor)(raster_symmetries[i2]);
		return intvl1.first < intvl2.first;
	});

	auto add_symmetries_in_interval = [&](int interval_lower, int interval_upper, int self, std::vector<int>& output)
	{
		// find another symmetry that is nested in this symmetry
		auto search_it = std::lower_bound(symmetry_order_by_begin.begin(), symmetry_order_by_begin.end(), interval_lower, [&](int j, const Vertex p)
		{
			auto intvl = std::forward<TIntervalFunctor>(interval_functor)(raster_symmetries[j]);
			return intvl.first < p;
		});

		bool potential_candidate = true;
		std::pair<int, int> intvl;

		auto wrap_it = [&](std::vector<int>::iterator it)
		{
			if (it == symmetry_order_by_begin.end())
				it = symmetry_order_by_begin.begin();

			intvl = std::forward<TIntervalFunctor>(interval_functor)(raster_symmetries[*it]);
			potential_candidate = PathUtils::contains_closed(B, interval_lower, interval_upper, intvl.first, circular);
			return it;
		};

		search_it = wrap_it(search_it);
		auto it = search_it;
		while (potential_candidate)
		{
			if (*it != self) // don't compare a symmetry to itself
			{
				auto& check_sym = raster_symmetries[*it];
				//check if this symmetry is nested inside the interval
				//only need to check the upper bound, we already checked the lower bound
				if (PathUtils::contains_closed(B, interval_lower, interval_upper, intvl.second, circular))
					output.push_back(*it);
			}
			++it;
			it = wrap_it(it);
			if (it == search_it)
				break;
		}
	};

	for (int i = 0; i < raster_symmetries.size(); ++i)
	{
		auto& s = raster_symmetries[i];

		auto& nested = nested_symmetries[i];

		add_symmetries_in_interval(s.v01, s.v00, i, nested);
		add_symmetries_in_interval(s.v10, s.v11, i, nested);		
	}
}

//Returns for each symmetry a list of other symmetries that are nested inside it
std::vector<std::vector<int>> find_nested_symmetries(const mat2x& B, const std::vector<Symmetry::SubPath>& raster_symmetries, bool circular)
{	
	std::vector<std::vector<int>> nested_symmetries(raster_symmetries.size()); //for each symmetry, a list of other symmetries nested inside it
	
	add_nested_symmetries_for_interval(B, raster_symmetries, circular, [](const Symmetry::SubPath& s) { return std::make_pair(s.v01, s.v00); }, nested_symmetries);
	add_nested_symmetries_for_interval(B, raster_symmetries, circular, [](const Symmetry::SubPath& s) { return std::make_pair(s.v10, s.v11); }, nested_symmetries);

	for (int i = 0; i < raster_symmetries.size(); ++i)
	{				
		auto& nested = nested_symmetries[i];		
		// re-establish original order
		std::sort(nested.begin(), nested.end());

		// remove duplicates
		for (auto it = nested.begin(); it != nested.end(); )
		{
			auto next = std::next(it);
			if (next == nested.end())
				break;
			if (*it == *next)
				it = nested.erase(it);
			else
				++it;
		}

		if (nested.size() > 0)
		{			
			auto& s = raster_symmetries[i];
			PF_VERBOSE_F("Symmetries nested in (%d %d - %d %d):", s.v01, s.v00, s.v10, s.v11);
			for (auto j : nested)
			{
				auto& s2 = raster_symmetries[j];
				PF_VERBOSE_F("(%d %d - %d %d)", s2.v01, s2.v00, s2.v10, s2.v11);
			}
		}
	}
	return nested_symmetries;
}

// Assumes that the path goes through the vertex
double calculate_angle_at_vertex(
    const mat2x& B,
    const vecXi& P,
    const int v
) {
    for (Index i = 0; i < P.size(); ++i) {
        if (P(i) == v) {
            return AngleUtils::spanned_shortest(B.col(CircularAt(P, i - 1)), B.col(P(i)), B.col(CircularAt(P, i + 1)));
        }
    }

    return INFINITY;
}

enum SYMMETRIC_ATTEMPTS {
    SYMMETRIC_ATTEMPT_SIDE0_SAME_ANGLE = 0,
    SYMMETRIC_ATTEMPT_SIDE1_SAME_ANGLE,
    SYMMETRIC_ATTEMPT_SIDE0_ANY_ANGLE,
    SYMMETRIC_ATTEMPT_SIDE1_ANY_ANGLE,
    SYMMETRIC_ATTEMPT_COUNT
};

int iterative_global(
	ShortestPath::State& G_state,
	const mat2x& B,
	const std::vector<Vertex>& M,
	const std::vector<BoundaryGraph::Edge>& original_E,
	vecXi& P,
	std::vector<Symmetry::SubPath>& raster_symmetries,
	const bool circular,
    IterationCallback on_iter,
    const std::vector<bool>& corners_C0
) {
	int attempts = 0;

	auto E = original_E;

	//
	std::vector<int> C;
	PathUtils::compute_convexities(B, C);

	// Sorting symmetries by priority
	std::sort ( raster_symmetries.begin(), raster_symmetries.end(), [&](const Symmetry::SubPath & lhs, const Symmetry::SubPath & rhs ) -> bool {
        const int len_rhs = CircularDist(B, rhs.v01, rhs.v00);
        const int len_lhs = CircularDist(B, lhs.v01, lhs.v00);
        if (len_rhs == len_lhs) {
            return (B.col(rhs.v01) - B.col(rhs.v11)).squaredNorm() < (B.col(lhs.v01) - B.col(lhs.v11)).squaredNorm();
        }

#if ORDER_BIG_TO_SMALL
		return len_lhs > len_rhs;
#else
		return len_lhs < len_rhs;
#endif	
	});

    for (auto s : raster_symmetries) {
        PF_VERBOSE_F("Symmetry %d %d %d %d length %d dist %f", s.v01, s.v00, s.v10, s.v11, CircularDist(B, s.v01, s.v00), (B.col(s.v01) - B.col(s.v00)).squaredNorm());
    }

#if HANDLE_NESTED_SYMMETRIES
	// Find symmetry nesting
	auto nested_symmetries = find_nested_symmetries(B, raster_symmetries, circular);
	std::vector<bool> symmetry_handled(raster_symmetries.size(), false);
#endif

	// holds the indices of the boundary graph vertices of the current polygon, they are used to lookup
	// error information without having to recompute it.
	std::vector<Vertex> PE;
    
    // subP_<symmetry_side>_<first/second half>
    vector<Index> subP_0_0, subP_1_1, subP_0_1, subP_1_0, subP_0, subP_1;

    // We make 4 attempts
    // 1. Preferred symmetric side, same angle outgoing from the symmetric region
    // 2. Preferred symmetric side
    // 3. Other side, same angle outgoing from the symmetric region
    // 4. Other side.
    vector<size_t> E_delete_attempts[SYMMETRIC_ATTEMPT_COUNT];
    unordered_set<BoundaryGraph::EdgeID> E_keep_attempts[2];

	std::function<void(int)> apply_symmetry = [&](int symmetry_index)
	{		
		const Symmetry::SubPath& s = raster_symmetries[symmetry_index];
#if HANDLE_NESTED_SYMMETRIES
		if (symmetry_handled[symmetry_index])
			return;
		
        symmetry_handled[symmetry_index] = true;
        
        for (auto nested : nested_symmetries[symmetry_index])
			apply_symmetry(nested);
#endif

        PF_VERBOSE_F("Applying symmetry %d %d %d %d", s.v01, s.v00, s.v10, s.v11);

        // Finding the optimal side and getting the subpaths to symmetrize
        PF_ASSERT(BoundaryGraph::trace_to_edges(E, P, PE, circular));
        const int optimal_side = pick_optimal_symmetry_side(B, P, E, PE, vec2i(s.v01, s.v00), vec2i(s.v10, s.v11), subP_0_0, subP_1_1, circular);
        if (optimal_side == 0) {
            return;
        }

        PathUtils::order(B, subP_0_0, s.v01);
        PathUtils::order(B, subP_1_1, s.v10);

        // Finding the target symmetric subpath
        subP_0_1.clear();
        for (size_t i = 0; i < subP_0_0.size(); ++i) {
            subP_0_1.emplace_back(Circular(B, s.v11 - CircularDist(B, s.v01, subP_0_0[subP_0_0.size() - i - 1])));
        }

        subP_1_0.clear();
        for (size_t i = 0; i < subP_1_1.size(); ++i) {
            subP_1_0.emplace_back(Circular(B, s.v01 + CircularDist(B, subP_1_1[subP_1_1.size() - i - 1], s.v11)));
        }

        // Concatenating them in a single path
        subP_0.clear();
        subP_0.insert(subP_0.end(), subP_0_0.begin(), subP_0_0.end());
        subP_0.insert(subP_0.end(), subP_0_1.begin(), subP_0_1.end());
        subP_1.clear();
        subP_1.insert(subP_1.end(), subP_1_0.begin(), subP_1_0.end());
        subP_1.insert(subP_1.end(), subP_1_1.begin(), subP_1_1.end());

        PF_VERBOSE_F("symmetry (%lld %lld) (%lld %lld) optimal side %d", s.v01, s.v00, s.v10, s.v11, optimal_side);
        for (Index i = 0; i < subP_0.size(); ++i) {
            PF_VERBOSE_F("subpath_0[%d] = %d", i, subP_0[i]);
        }
        for (Index i = 0; i < subP_1.size(); ++i) {
            PF_VERBOSE_F("subpath_1[%d] = %d", i, subP_1[i]);
        }

        // Which edges should be preserved in the first attempt
        E_keep_attempts[0].clear();
        for (size_t i = 0; i < subP_0.size() - 1; ++i) {
            const auto e_id = BoundaryGraph::make_edge_id(subP_0[i], subP_0[i + 1]);
            PF_ASSERT(E_keep_attempts[0].find(e_id) == E_keep_attempts[0].end());
            E_keep_attempts[0].insert(e_id);
        }

        E_keep_attempts[1].clear();
        for (size_t i = 0; i < subP_1.size() - 1; ++i) {
            const auto e_id = BoundaryGraph::make_edge_id(subP_1[i], subP_1[i + 1]);
            PF_ASSERT(E_keep_attempts[1].find(e_id) == E_keep_attempts[1].end());
            E_keep_attempts[1].insert(e_id);
        }

        // Now finding the list of edges which should be deleted for each attempt
        for (int i = 0; i < SYMMETRIC_ATTEMPT_COUNT; ++i) { E_delete_attempts[i].clear(); }

        const double angle_optimal_0 = calculate_angle_at_vertex(B, P, subP_0.front());
        const double angle_optimal_1 = calculate_angle_at_vertex(B, P, subP_1.back());
        PF_VERBOSE_F("angle_optimal_0 %f angle_optimal_1 %f", angle_optimal_0, angle_optimal_1);

        for (size_t i = 0; i < E.size(); ++i) {
            // Completely outside the symmetric region
            if (!PathUtils::contains_closed(B, s.v01, s.v00, E[i].v0, circular) &&
                !PathUtils::contains_closed(B, s.v01, s.v00, E[i].v1, circular) &&
                !PathUtils::contains_closed(B, s.v10, s.v11, E[i].v0, circular) &&
                !PathUtils::contains_closed(B, s.v10, s.v11, E[i].v1, circular)) {
                continue;
            }

            // Preserving all the edges which are incoming/outgoing from the last vertices in the symmetric region
            const bool keep_0 = E[i].v0 == subP_0_0.back() || E[i].v0 == subP_0_1.back() || E[i].v1 == subP_0_0.front() || E[i].v1 == subP_0_1.front();
            const bool keep_1 = E[i].v0 == subP_1_0.back() || E[i].v0 == subP_1_1.back() || E[i].v1 == subP_1_0.front() || E[i].v1 == subP_1_1.front();

            bool same_angle_0 = false, same_angle_1 = false;
            if (keep_0 && subP_0.size() >= 2) {
                const double angle = E[i].v0 == subP_0.back() ?
                    AngleUtils::spanned_shortest(B.col(subP_0[subP_0.size() - 2]), B.col(E[i].v0), B.col(E[i].v1)) :
                    AngleUtils::spanned_shortest(B.col(E[i].v0), B.col(E[i].v1), B.col(subP_0[1]));
                same_angle_0 = abs(angle - angle_optimal_0) < PF_EPS;
            }

            if (keep_1 && subP_1.size() >= 2) {
                const double angle = E[i].v0 == subP_1.back() ?
                    AngleUtils::spanned_shortest(B.col(subP_1[subP_1.size() - 2]), B.col(E[i].v0), B.col(E[i].v1)) :
                    AngleUtils::spanned_shortest(B.col(E[i].v0), B.col(E[i].v1), B.col(subP_1[1]));
                same_angle_1 = abs(angle - angle_optimal_1) < PF_EPS;
            }

            const auto e_id = BoundaryGraph::make_edge_id(E[i].v0, E[i].v1);
            // Is it one of the edges which should be preserved?
            if (E_keep_attempts[0].count(e_id) == 0) {
                if (keep_0 && !same_angle_0) {
                    E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_SAME_ANGLE].emplace_back(i);
                } else if (!keep_0) {
                    if (PathUtils::contains_closed(B, subP_0.front(), subP_0.back(), E[i].v0, circular) ||
                        PathUtils::contains_closed(B, subP_0.front(), subP_0.back(), E[i].v1, circular)) {
                        E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_ANY_ANGLE].emplace_back(i);
                        E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_SAME_ANGLE].emplace_back(i);
                    }
                }
            }

            if (E_keep_attempts[1].count(e_id) == 0) {
                if (keep_1 && !same_angle_1) {
                    E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_SAME_ANGLE].emplace_back(i);
                }
                else if (!keep_1) {
                    if (PathUtils::contains_closed(B, subP_1.front(), subP_1.back(), E[i].v0, circular) ||
                        PathUtils::contains_closed(B, subP_1.front(), subP_1.back(), E[i].v1, circular)) {
                        E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_ANY_ANGLE].emplace_back(i);
                        E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_SAME_ANGLE].emplace_back(i);
                    }
                }
            }
        }

        for (int attempt = 0; attempt < SYMMETRIC_ATTEMPT_COUNT; ++attempt) {
            // Figuring out which data is used for this attempt
            vector<size_t>* E_delete = nullptr;
            vector<Index>* subP = nullptr;
            bool test_same_angle = false;

            switch (attempt) {
            case SYMMETRIC_ATTEMPT_SIDE0_ANY_ANGLE:
                E_delete = optimal_side > 0 ? &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_ANY_ANGLE] : &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_ANY_ANGLE];
                subP = optimal_side > 0 ? &subP_0 : &subP_1;
                break;
            case SYMMETRIC_ATTEMPT_SIDE0_SAME_ANGLE:
                E_delete = optimal_side > 0 ? &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_SAME_ANGLE] : &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_SAME_ANGLE];
                subP = optimal_side > 0 ? &subP_0 : &subP_1;
                test_same_angle = true;
                break;
            case SYMMETRIC_ATTEMPT_SIDE1_ANY_ANGLE:
                E_delete = optimal_side > 0 ? &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_ANY_ANGLE] : &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_ANY_ANGLE];
                subP = optimal_side > 0 ? &subP_1 : &subP_0;
                break;
            case SYMMETRIC_ATTEMPT_SIDE1_SAME_ANGLE:
                E_delete = optimal_side > 0 ? &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE1_SAME_ANGLE] : &E_delete_attempts[SYMMETRIC_ATTEMPT_SIDE0_SAME_ANGLE];
                subP = optimal_side > 0 ? &subP_1 : &subP_0;
                test_same_angle = true;
                break;
            default:
                PF_ABORT;
            }

            for (Index j = 0; j < subP->size(); ++j) {
                PF_VERBOSE_F("%d", subP->at(j));
            }

            for (Index j = 0; j < E_delete->size(); ++j) {
                PF_VERBOSE_F("delete %d %d", E[E_delete->at(j)].v0, E[E_delete->at(j)].v1);
            }

            vector<BoundaryGraph::Edge> E_new = E;
            EraseOrdered(E_new, *E_delete);

            vecXi P_new;
            minimum(G_state, B, M, E_new, P_new, circular);

            if (is_symmetric_path_valid(B, P, P_new, *subP, circular)) {
                if (test_same_angle) {
                    const double angle_0 = calculate_angle_at_vertex(B, P_new, subP->front());
                    const double angle_1 = calculate_angle_at_vertex(B, P_new, subP->back());

                    if (abs(angle_0 - angle_1) > PF_EPS_MEDIUM) {
                        continue;
                    }
                }

                PF_VERBOSE_F("Accept symmetry %d %d %d %d attempt %d", s.v01, s.v00, s.v10, s.v11, attempt);
                P = P_new;
                E = E_new;
                break;
            }
        }
	};

	// If two paths are already symmetric, they are not explicitly skipped in order to prune the graph
	for (int i = 0; i < raster_symmetries.size(); ++i) {	
		apply_symmetry(i);
	}

	return attempts;
}

int find_edge_indices(
	const std::vector<BoundaryGraph::Edge>& E, // List of edges where the index is computed from
	const std::vector<BoundaryGraph::Edge>& V, // Values for which the respective indices need to be computed
	std::vector<size_t>& I                     // Where the indices are written
) {
	I.clear();

	for (size_t i = 0; i < V.size(); ++i) {
		for (size_t j = 0; j < E.size(); ++j) {
			if (E[j].v0 == V[i].v0 && E[j].v1 == V[i].v1) {
				I.emplace_back(j);
				break;
			}
		}
	}

	return V.size() == I.size();
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)