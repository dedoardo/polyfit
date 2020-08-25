#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/regularity/symmetry.hpp>

// Eigen
#include <Eigen/Core>

// libc++
#include <unordered_map>
#include <vector>
#include <memory>

NAMESPACE_BEGIN(polyvec)

// Forward declarations
class FitClassifier;

//
class ImageBoundary {
public:
    ImageBoundary(
        const Eigen::Matrix2Xi& points, 
        const int polygon_id,
        const Eigen::Vector4d& color
    );

    // Returns true if the boundary vertex 'v' is flat, if that is the case, edges will not be
    // constructed starting/from or to it
    bool is_flat(const int v) const;

    // Constructs the boundary graph allowing for flat vertices
    void initialize_graph_with_flat_vertices(const Eigen::VectorXi& B_flat);

    // Fits a polygon with regularization
    void trace_polygon(
        const std::unordered_set<int>& kopf_junctions = std::unordered_set<int>()
    );
    
    // Fits a polygon without regularization
    void trace_polygon_simple();

    void compute_regularities();

    // Resets the graph to its original state preserving the fitting results
    void reset_graph();

    // Adds an edge to the boundary graph, validating the indices
    void add_edge(const polyfit::BoundaryGraph::Edge& e);

    void remove_edges_crossing_vertices(const std::vector<int> V);

    double bounding_box_area()const;

    // Prunes the boundary graph so that the polygon will go through the specified corners.
    // It does not trace a new polygon.
    // The subpath *will* exist regardless of the current state of the boundary graph, but
    // a polygon is not guaranteed to exist.
    // It validates the subpath indices against the boundary
    void force_graph_subpath(const Eigen::VectorXi& subpath);

    // Returns the list of vertices in the polygon contained in [v_src, v_dst]
    Eigen::VectorXi get_polygon_subpath(const int v_src, const int v_dst) const;

    // Same as get_polygon_subpath but returns the index of the first and last vertex in the closed boundary interval [v_src, v_dst]
    std::pair<int, int> get_polygon_subpath_bounds(const int v_src, const int v_dst) const;

    // Returns the sum of the edge cost for all polygon edges contained within [v_first, v_last]
    double calculate_subpath_energy(const int v_first, const int v_last) const;

    // Returns the index of the polygon vertex matching the specified boundary vertex or -1 if it doesn't exist.
    int find_polygon_vertex_for_boundary_vertex(const int v)const;

    // assumes that v is a C^0 corner
    int find_line_before_polygon_vertex(const int v)const;
    int find_line_after_polygon_vertex(const int v)const;

    // Obtains a tangent classification for each corner. A polygon must have been traced
    // No curves are fit
    void classify_corners(FitClassifier& classifier);

    // Initializes the curves with their initial guesses
    void fit_curves_initial_guess();

    void flip_coordinate_system(const polyfit::BoundaryGraph::EdgeID edge);

    std::pair<int, int> get_primitives_interval_open   (const int v_first, const int v_last, const int direction) const;
    std::pair<int, int> get_primitives_interval_closed (const int v_first, const int v_last, const int direction) const;

    std::vector<CurvePrimitive> clone_curves_at_corner_reversed (const int v_P) const;

    // angle at *polygon* vertex
    double angle_at_vertex(const int v);

    // -----------------------------------------------
    // Getters
    int id() const { return _id; }
    const Eigen::Vector4d& color() const { return _color; }

    // Input raster
    const Eigen::Matrix2Xi& raster_points() const { return _points; }
    const Eigen::Matrix2Xd& raster_points_as_double()const { return _points_as_double; }

    // Boundary
    const Eigen::Matrix2Xd& boundary_points() const { return _polygon.B; }
    const polyfit::BoundaryGraph::AccuracyMap& accuracy_map()const { return _accuracy_map; }
    const std::vector<polyfit::BoundaryGraph::Edge>& edges() const { return _polygon.E; }
    std::vector<polyfit::BoundaryGraph::Edge>& edges() { return _polygon.E; }
    const std::vector<polyfit::BoundaryGraph::Edge>& edges_original() const { return _E_original; }
    std::vector<size_t>& E_delete() { return _E_delete; }
    const std::vector<int>& convexities() const { return _boundary_convexities; } // inside 
    std::vector<polyfit::Symmetry::SubPath>& raster_symmetries() { return _polygon.raster_symmetries; }
    std::vector<polyfit::Symmetry::SubPath>& raster_symmetries_local() { return _polygon.raster_symmetries_local; }
    const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries() const { return _polygon.raster_symmetries; }

    // Polygon
    Eigen::VectorXi& polygon_vertices() { return _polygon.P; }
    const Eigen::VectorXi& polygon_vertices() const { return _polygon.P; }
    Eigen::Matrix2Xd& polygon_points() { return _polygon.PP; }
    const Eigen::Matrix2Xd& polygon_points() const { return _polygon.PP; }
    std::vector<Eigen::Index>& midpoints() { return _polygon.M; }
    const std::vector<Eigen::Index>& midpoints() const { return _polygon.M; }
    std::vector<bool>& midpoints_bits() { return _polygon.midpoints; }
    const std::vector<bool>& midpoints_bits() const { return _polygon.midpoints; }
    const std::vector<int>& polygon_convexities_inside() const { return _convexities_inside; }
    const std::vector<int>& polygon_convexities_outside() const { return _convexities_outside; }
    const polyfit::Regularity::RegularityInformation& regularity_graph() const { return _polygon.RE; }
	polyfit::Regularity::RegularityInformation& regularity_graph() { return _polygon.RE; }

    // Smooth fits
    std::vector<TangentFitType>& tangents_fits() { return _tangents; }
    const std::vector<TangentFitType>& tangents_fits() const { return _tangents; }
    CurvePrimitiveSequence& spline() { return _spline; }
    const CurvePrimitiveSequence& spline()const { return _spline; }
    const std::vector<polyvec::AttemptInfo>& fitting_attempts() const { return _fitting_attempts; }

private:
    void _initialize_curve_fitter();
    void _update_polygon();
    void _update_convexities();
    bool _is_polygon_data_consistent_and_exists();
    void _initialize_graph_if_necessary();

    // Initialization 
    const int _id;
    const Eigen::Vector4d _color;

    // Raster
    Eigen::Matrix2Xi _points;
    Eigen::Matrix2Xd _points_as_double;
    std::vector<int> _boundary_convexities; // w.r.t. _polygon.B
    std::vector<bool> _B_flat;

    // Polygon
    PolygonData                               _polygon;
    std::vector<polyfit::BoundaryGraph::Edge> _E_original;
    polyfit::BoundaryGraph::AccuracyMap                _accuracy_map;
    polyfit::ShortestPath::State              _shortest_path_state;
    std::vector<size_t>                       _E_delete;
    
    //
    std::vector<int>                          _convexities_inside;
    std::vector<int>                          _convexities_outside;

    // Curves
    bool _curves_dirty;
    std::unique_ptr<CurveSequenceFitter> _spline_fitter;
    CurvePrimitiveSequence _spline;
    std::vector<TangentFitType> _tangents;
    std::vector<polyvec::AttemptInfo> _fitting_attempts;

    std::vector<bool> _edge_coordinate_systems;
};

NAMESPACE_END(polyvec)
