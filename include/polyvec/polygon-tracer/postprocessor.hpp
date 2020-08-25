#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/curve-tracer/spline_types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/regularity/symmetry.hpp>

// libc++
#include <vector>

NAMESPACE_BEGIN(polyvec)

class FitClassifier;

/*
    Performs visual improvement operations after a symmetric graph, polygon and curve have
    been fit to the given boundaries.
*/
class PostProcessor {
public:
    PostProcessor(
        const polyfit::mat2x&                      B,
        const std::vector<polyfit::Vertex>&        M,
        std::vector<polyfit::BoundaryGraph::Edge>& E,
        polyfit::vecXi&                            P,
        polyfit::mat2x&                            PP,
        std::vector<polyfit::Symmetry::SubPath>& raster_symmetries,
        std::vector<polyfit::Symmetry::SubPath>&   raster_symmetries_local,
		polyfit::Regularity::RegularityInformation& regularity,
        std::vector<TangentFitType>&               tangent_fits,
        CurvePrimitiveSequence&                    curves,
        const char*                                classifier_uri,
        const int polygon_id = -1
    );

    void retrace_with_C0_corners();

    void collapse_short_segments();
    int  collapse_short_segments_iteration();    

    const std::vector<bool>& get_immovable_corners() const;
    const std::vector<polyfit::BoundaryGraph::Edge>& get_original_edges() const;
    const std::vector<size_t>& get_permanently_deleted_edges() const;

	//Tries to merge adjacent curves to produce more continuous results
	void merge_curves();

private:
	void trace_curves_from_polygon(FitClassifier& classifier);

    // Compares _P and P for degenerate edges/inflections
    bool _should_accept_trace() const;
    bool _relocate_short_edge(const int v_src, const int v_dst, const int v_force);
    bool _attempt_remove_edge_force_corner ( const int v_src, const int v_dst, const int v_force);
    void _reset_immovable_corners();
    void _set_immovable_corners();

    std::unordered_map<polyfit::BoundaryGraph::EdgeID, Eigen::Vector2d> _accuracy_map;
    const bool                                 circular;
    const polyfit::mat2x&                      B;
    const std::vector<polyfit::Vertex>&        M; // Yo, remove this
    const int                                  polygon_id;
    std::vector<polyfit::BoundaryGraph::Edge>& E;
    polyfit::vecXi&                            P;
    polyfit::mat2x&                            PP;
    std::vector<polyfit::Symmetry::SubPath>&   raster_symmetries;
    std::vector<polyfit::Symmetry::SubPath>&   raster_symmetries_local;
	polyfit::Regularity::RegularityInformation& regularity;
    std::vector<TangentFitType>&               tangent_fits;
    CurvePrimitiveSequence&                    curves;
    const std::string                          classifier_uri;

    // Permanent
    // Since each retracing operations starts from the original list of edges, the indices remain valid
    polyfit::ShortestPath::State              _G_state;
    std::vector<polyfit::BoundaryGraph::Edge> _E_original;
    std::vector<size_t>                       _E_delete_permanent;

    // This is the list of edges and its delete list in the current iteration
    polyfit::vecXi                            _P;
    std::vector<polyfit::BoundaryGraph::Edge> _E;
	polyfit::Regularity::RegularityInformation _RE;
    std::vector<size_t>                       _E_delete;
    std::vector<size_t>                       _E_delete_all; // _E_delete + _E_delete_permanent
    std::vector<size_t> _E_delete_temp;

    std::vector<bool> _immovable;
};

NAMESPACE_END(polyvec)