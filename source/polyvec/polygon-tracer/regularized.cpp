// polyvec
#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/regularized.hpp>
#include <polyvec/polygon-tracer/iterative-global.hpp>
#include <polyvec/polygon-tracer/detect_parallel_edges_relative_proximity.hpp>
#include <polyvec/polygon-tracer/remove_edges_crossing_implicit_corners.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/polygon-tracer/align_opposite_corners.hpp>
#include <polyvec/polygon-tracer/preserve_axis_aligned.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/regularity/replace_flat_nodes_with_edges.hpp>
#include <polyvec/regularity/collapse_opposite_short_segments.hpp>
#include <polyvec/options.hpp>
#include <polyvec/regularity/regularity_action.hpp>
#include <polyvec/regularity/add_local_symmetries.hpp>
#include <polyvec/regularity/add_parallel_edges.hpp>
#include <polyvec/regularity/add_90_degree_edges.hpp>
#include <polyvec/regularity/continuations.hpp>
#include <polyvec/regularity/continuations_2018.hpp>
#include <polyvec/core/constants.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

#define ON_ATTEMPT ++attempts; if (on_iter) { on_iter(B, P, E); }

int regularized (
    ShortestPath::State& G_state,
    const mat2x& B,
    const std::vector<Vertex>& M,
    std::vector<BoundaryGraph::Edge>& E,
	polyfit::Regularity::RegularityInformation&    RE,
    vecXi& P,
	std::vector<Symmetry::SubPath>& raster_symmetries,
    std::vector<Symmetry::SubPath>& raster_symmetries_local,
    const bool circular,
    IterationCallback on_iter,
    const std::vector<bool>& corners_C0,
    const std::vector<bool>& B_flat,
    const std::unordered_set<int>& kopf_junctions
) {
    // Forcing the polygon to include 2(+) x 2(+) corners
#if POLYVEC_REGULARIZE
    remove_edges_crossing_implicit_corners(B, E, circular);
#endif
    
    int attempts = 0;

    // ---
    if (!PolygonTracer::minimum(G_state, B, M, E, P, circular, corners_C0, B_flat)) {
        return 0;
    }
    ON_ATTEMPT;

#if !POLYVEC_REGULARIZE
	return 1;
#endif

	preserve_axis_aligned_edges(G_state, B, M, E, P, circular);

	std::vector<RegularityAction*> regularity_actions;
	add_local_symmetries(B, raster_symmetries_local, B_flat, circular, regularity_actions);
	add_parallel_edges(B, circular, regularity_actions);
    apply_regularities(G_state, B, M, E, P, circular, on_iter, corners_C0, regularity_actions);

	for (auto a : regularity_actions)
		delete a;
	regularity_actions.clear();

    add_90_degree_edges(B, P, circular, regularity_actions);
    apply_regularities(G_state, B, M, E, P, circular, on_iter, corners_C0, regularity_actions);

    for (auto a : regularity_actions)
        delete a;
    regularity_actions.clear();


    ON_ATTEMPT;
    // Solve with global symmetries
    iterative_global(G_state, B, M, E, P, raster_symmetries, circular, on_iter, corners_C0);
    ON_ATTEMPT;	

    // We now have an almost final approximation. We need to guarantee that the successive
    // postprocessing operations do not break existing regularities. Such operations also take
    // the existing regularity graph.
	RE = polyfit::Regularity::RegularityInformation();
	RE.update(B, P, raster_symmetries, raster_symmetries_local, circular);
	
	Regularity::replace_flat_nodes_with_edges(B, P, B_flat, RE, raster_symmetries, raster_symmetries_local);

	// Find continuations	
	RE.add_continuations(B, P, E, circular, true);	

	RE.clear_parallels();
	RE.update(B, P, raster_symmetries, raster_symmetries_local, circular);

    //polyvec::find_continuation_candidates_2018(B, P, RE, circular, true);

    return 1;
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)