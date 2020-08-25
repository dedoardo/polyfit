// polyvec
#include <polyvec/polygon-tracer/prune_edges_force_vertex.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/core/log.hpp>

// Reports a pair of vertices for each edge removed
#define LOG_PRUNED_EDGES 1

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int prune_edges_force_vertex(
	const mat2x& B,
	std::vector<BoundaryGraph::Edge>& E,
	const Eigen::Index v,
	const bool circular
) {
	auto E_it = E.begin();
	int removed = 0;

#if LOG_PRUNED_EDGES
	PF_LOG_DIVIDER;
	PF_LOGF("prune edges crossing %lld", v);
#endif

	while (E_it != E.end()) {
		if (PathUtils::contains_open(B, E_it->v0, E_it->v1, v, circular)) {
#if LOG_PRUNED_EDGES
			PF_LOGF("remove %lld -> %lld", E_it->v0, E_it->v1);
#endif

			E_it = E.erase(E_it);
			++removed;
		} else {
			++E_it;
		}
	}

	return removed;
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)
