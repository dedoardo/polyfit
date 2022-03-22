#include <polyvec/regularity/regularity_action.hpp>

#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/core/log.hpp>

using namespace polyfit;
using namespace PolygonTracer;

void PolygonTracer::apply_regularities(
	ShortestPath::State& G_state,
	const mat2x& B,
	const std::vector<Vertex>& M,
	std::vector<BoundaryGraph::Edge>& E,
	vecXi& P,
	const bool circular,
	IterationCallback on_iter,
	const std::vector<bool>& corners_C0,
	const std::vector<RegularityAction*>& regularity_actions)
{
	struct PolygonState
	{
		vecXi P;
		mat2x PP;
		int inconsitent_edges;
        int inflections;
	};

	PolygonState current, old;

	old.P = P;
	BoundaryGraph::trace_to_points(B, P, old.PP);
	old.inconsitent_edges = AngleUtils::count_visually_inconsistent_edges(old.PP, circular);
    old.inflections = AngleUtils::count_inflections(old.PP, circular);

	for (auto action : regularity_actions) {
		auto E_delete = action->get_edges_to_delete(B, old.P, E, circular);

		if (E_delete.empty())
			continue;

		// Holds the edges which have been temporarily removed from the graph
		std::vector<BoundaryGraph::Edge> E_rm;
		E_rm.reserve(E_delete.size());
        for (auto i : E_delete)
            E_rm.push_back(E[i]);

        PF_VERBOSE_F("erasing %lld edges temporarily", E_delete.size());
		std::sort(E_delete.begin(), E_delete.end());
		EraseOrdered(E, E_delete);

		// Attempting to find a more regular polygon
		bool fail = !PolygonTracer::minimum(G_state, B, M, E, current.P, circular, corners_C0);

		// In case we found a solution, checking if it reduced the quality of the model ?
		if (!fail) {
			BoundaryGraph::trace_to_points(B, current.P, current.PP);
			current.inconsitent_edges = AngleUtils::count_visually_inconsistent_edges(current.PP, circular);
            current.inflections = AngleUtils::count_inflections(current.PP, circular);

            PF_VERBOSE_F("comparing visual quality old %d - new %d", old.inconsitent_edges, current.inconsitent_edges);
			if (current.inconsitent_edges > old.inconsitent_edges && !action->can_introduce_degenerate_edges()) {
				fail = true;
			}

            if (current.inflections > old.inflections && !action->can_introduce_inflections()) {
                fail = true;
            }

            if (!action->is_polygon_acceptable(B, old.P, current.P, circular)) {
                fail = true;
            }
		}

		// calling the callback before reinserting the edges in case of failure
		if (on_iter) {
			on_iter(B, P, E);
		}

		// In case of failure, restoring the older connectivty
		if (fail) {
            PF_VERBOSE_S("failed to apply regularity");
			current = old;
			E.insert(E.end(), E_rm.begin(), E_rm.end());
            PF_VERBOSE_F("restore older path and insert %d edges", E_rm.size());
		}
		else
			old = current;
	}

	P = old.P;
}