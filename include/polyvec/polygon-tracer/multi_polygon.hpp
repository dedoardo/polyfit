#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/mc/segment_connectivity.hpp>
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/polygon-tracer/image_boundary.hpp>

// libc++
#include <string>
#include <utility>
#include <vector>

// Forward declarations
namespace polyfit {
namespace mc {
class RasterImageConnectivity;
class SegmentConnectivity;
}
}

NAMESPACE_BEGIN(polyvec)

struct JunctionInfo
{
  // unique global-id
  int vertex;

  // assuming that a junction
  // cannot correspond to two points on a boundary. This is a valid assumption
  // if the image segmentation follows [Kopf11]
  std::vector<polyfit::BoundaryGraph::EdgeID> boundary_edges;
  
  // index of the neighboring regions (== valence of the junction)
  std::vector<int> boundary_ids;

  std::vector<int> polygon_vertices;
};

// Note that these are EdgeID computed from the polygon vertices, not
// the boundary vertices
struct CoupledEdge {
    polyfit::BoundaryGraph::EdgeID id_0;
    polyfit::BoundaryGraph::EdgeID id_1;
};

// A shared boundary between two regions (boundary
// This is only defined 
struct OpenSegmentInfo {
    int id;

    int boundary_src;
    int boundary_dst;
    int direction_src = 0;
    int direction_dst = 0;

    // Boundary point bounds for the two regions
    std::pair<int, int> B_bounds_src;
    std::pair<int, int> B_bounds_dst;

    Eigen::Matrix2Xi points;

    // Asserts that the information is not contracticting 
    void validate() const;
};

// If I could figure out how to use mc::MeshConnectivity...
struct HalfEdge {
    unsigned int id = -1;
    unsigned int id_prev = -1;
    unsigned int id_next = -1;
    unsigned int id_twin = -1;

    polyfit::BoundaryGraph::EdgeID edge;
    int boundary;

    const std::string str()const;
};

struct HalfEdgeConnectivity {
    int boundary(const int he);
    int vertex_from(const int he);
    int vertex_to(const int he);
    int next(const int he);
    int prev(const int he);
    int twin(const int he);

    std::vector<HalfEdge> half_edges;
};

// suffixes are important
struct ImageToBoundaryEdgeMap {
    using ID = polyfit::BoundaryGraph::EdgeID;

    ImageToBoundaryEdgeMap(const int cols, const std::vector<ImageBoundary>& boundaries);

    int image_vertex_id_from_point (const Eigen::Vector2d& p) const;

    // Preserves the order, skipping midpoints
    ID make_boundary_edge_id(
        const int boundary,
        int v0_B, 
        int v1_B,
        // This should be set to true if you are creating a boundary edge *not* from image points (i.e. you are 
        // aware of the existance of midpoints)
        const bool include_midpoints = false
    ) const;

    // Creates a unique edge id for a boundary point, setting invert_edge if the edge stored in the image has to be inverted
    ID make_image_edge_id(
        const int boundary, 
        const int v0_B, 
        const int v1_B,
        bool& invert_edge
    )const;

    // Just creates an image edge id from the points, who cares where they are coming from
    ID make_image_edge_id(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1
    )const;

    // Adds an half-edge to the map
    void add_boundary_half_edge (
        const int boundary, // boundary the half-edge is incident to
        const ID& e_B       // Boundary half-edge
    );

    // queries the twin of an half-edge
    ID twin(
        const ID& e_B,
        const int boundary_src,
        const int boundary_dst
    ) const;

    int to_boundary_vertex_from (
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const int boundary // for validation
    ) const;

    int to_boundary_vertex_to(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        const int boundary // for validation
    ) const;

    void get_incident_boundaries(
        const Eigen::Vector2d& p0,
        const Eigen::Vector2d& p1,
        int& boundary_0,
        int& boundary_1
    ) const;

private:
    // The boundary stuff is not necessarily but helpful for validation
    struct EdgePair {
        ID he;
        int boundary = -1;
        ID he_twin;
        int boundary_twin = -1;
    };

    int cols;
    std::unordered_map<ID, EdgePair> _pairs;
    const std::vector<ImageBoundary>& _boundaries;
};

class MultiPolygonTracer {
public:
    MultiPolygonTracer(
        polyfit::mc::RasterImageConnectivity& conn_raster,
        polyfit::mc::SegmentConnectivity& conn_segments,
        const std::vector<Eigen::Vector4d>& boundary_colors,
        const char* classifier_uri);

    // Operations
    void trace();

    // Getters
    std::vector<int> get_boundaries_ordered_by_area() const;
    const std::vector<ImageBoundary>& boundaries() const;
    const std::vector<JunctionInfo>& junctions() const;
    const std::vector<Eigen::Matrix2Xi> segments()const;
    const std::vector<HalfEdge>& half_edges()const;

    // Accounting for midpoints by scaling all the points by 2
    int boundary_vertex_id_from_pos(const Eigen::Vector2d pos) const;
    Eigen::Vector2d boundary_vertex_pos_from_id(const int id) const;

private:
    //
    // Algorithm steps
    //
    void _initialize_image_boundaries ();
    void _initialize_segments();
    void _build_image_to_boundary_map();
    void _remove_non_shared_midpoints();
    void _extract_junctions();
    void _extract_polygon_forcing_junctions();
    void _preserve_thickness_at_kopf_junctions();
    void _merge_polygons_and_compute_open_segments_info();
    void _construct_regularity_graph();
    void _classify_polygon_corners_and_junctions();
    void _merge_tangent_classifications_along_open_segments();
    void _downgrade_tangent_classifications_at_junctions();
    void _calculate_coordinate_system_directions_along_open_segments();
    void _construct_curves_initial_guess();
    void _initialize_curve_fitter();
    void _couple_curves_along_open_segments();
    void _merge_consecutive_parallel_lines_except_at_junctions();
    void _solve_curves();
    void _post_process();
    void _build_polygon_half_edges_and_fill_junctions_map();
    void _merge_curves_at_junctions();

private:
    void _report_info();

    polyfit::BoundaryGraph::EdgeID make_image_edge_id(
        const Eigen::Vector2d& v0,
        const Eigen::Vector2d& v1
    );

    std::pair<int, int> _get_boundary_subpath_for_segment(
        const int segment_id,
        const int boundary_id
    ) const;

    std::pair<int, int> _get_polygon_subpath_for_segment(
        const int b_src,
        const int b_dst,
        const int region,
        Eigen::VectorXi* P = nullptr
    ) const;

    //
    // Members
    //
private:
    // List of unique image boundaries
    std::vector<ImageBoundary>    _boundaries;
    std::vector<bool>             _boundaries_closed;
    std::vector<std::vector<int>> _boundary_junctions;
    std::vector<int>              _polygon_to_boundary;
    int                           _outer_boundary;

    // Image topology
    polyfit::mc::RasterImageConnectivity& _conn_raster;
    polyfit::mc::SegmentConnectivity&     _conn_segments;

    // File pointing to the serialized classifier
    std::string _classifier_uri;

    // Tangent fit classifier
    FitClassifierRandomForest _classifier;

    // E_{image} -> (HE_{boundary0}, HE_{boundary1})
    std::unique_ptr<ImageToBoundaryEdgeMap> _I2B_map;

    // List of unique junctions
    std::vector<JunctionInfo> _junctions;
    std::unordered_map<int, polyfit::BoundaryGraph::EdgeID> _junctions_to_half_edges;
    std::unordered_set<int> _kopf_junctions;
    std::vector<std::vector<int>> _kopf_junctions_per_boundary;

    // Colors for each polygon
    std::vector<Eigen::Vector4d> _colors;

    std::vector<OpenSegmentInfo> _open_segments;

    // Is the segment closed?
    std::vector<bool> _segment_closed;
    std::vector<Eigen::Matrix2Xi>  _segment_points;

    // Bitmap of [polygon][edge] which is set to 1 if the coordinate system for this region is to be
    // inverted in order to couple the parameters
    std::vector<std::vector<bool>> _inverted_edge_coord_systems;

    HalfEdgeConnectivity _he_polygon;

    // Primitives for which the degrees of freedom need to be reduced. Reducing the degrees of
    // freedom has to be done after initializing the curve fitter but the information is tipically computed before hand.
    std::vector<std::pair<CurvePrimitive*, DofOptions::Type>> _primitives_to_constraint;

    // Curve solver interface
    std::unique_ptr<CurveFitter> _curve_fitter;
};

NAMESPACE_END(polyvec)