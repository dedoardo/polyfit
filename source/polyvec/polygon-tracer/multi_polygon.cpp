// Polyvec
#include <polyvec/polygon-tracer/multi_polygon.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/polygon-tracer/regularized.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/regularity/replace_flat_nodes_with_edges.hpp>
#include <polyvec/curve-tracer/spline_types.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
#include <polyvec/mc/segment_connectivity.hpp>
#include <polyvec/mc/junction.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/directions.hpp>
#include <polyvec/polygon-tracer/postprocessor.hpp>
#include <polyvec/core/constants.hpp>

// debug remove
#include <polyvec/io/pdf.hpp>
#include "../../../apps/dev/drawing.hpp"

// C++ STL
#include <unordered_set>
#include <memory>
#include <numeric>

using namespace polyfit;
using namespace polyvec;
using namespace std;
using namespace Eigen;

NAMESPACE_BEGIN(polyvec)

void OpenSegmentInfo::validate() const {
    PF_ASSERT(direction_src != 0 && direction_dst != 0);
    PF_ASSERT(direction_src == -direction_dst);
    PF_ASSERT(boundary_src != boundary_dst);
}

const std::string HalfEdge::str()const {
    const auto eid = BoundaryGraph::unpack_edge_id(edge);
    return StringUtils::fmt("half-edge id %d boundary %d edge %d-%d", id, boundary, eid(0), eid(1));
}

int HalfEdgeConnectivity::boundary(const int he) {
    return half_edges[he].boundary;
}

int HalfEdgeConnectivity::vertex_from(const int he) {
    return BoundaryGraph::unpack_edge_id(half_edges[he].edge)(0);
}

int HalfEdgeConnectivity::vertex_to(const int he) {
    return BoundaryGraph::unpack_edge_id(half_edges[he].edge)(1);
}

int HalfEdgeConnectivity::next(const int he) {
    return half_edges[he].id_next;
}

int HalfEdgeConnectivity::prev(const int he) {
    return half_edges[he].id_prev;
}

int HalfEdgeConnectivity::twin(const int he) {
    return half_edges[he].id_twin;
}

ImageToBoundaryEdgeMap::ImageToBoundaryEdgeMap(const int cols, const std::vector<ImageBoundary>& boundaries) 
    : cols(cols), _boundaries(boundaries) { }

int ImageToBoundaryEdgeMap::image_vertex_id_from_point(const Eigen::Vector2d& p) const {
    const Eigen::Vector2i pos_i = (p * 2).cast<int>();
    return pos_i.x() + pos_i.y() * (cols * 2);
}

ImageToBoundaryEdgeMap::ID ImageToBoundaryEdgeMap::make_boundary_edge_id(
    const int boundary,
    int v0_B,
    int v1_B,
    const bool include_midpoints
) const {
    const Eigen::Matrix2Xd& pts = _boundaries[boundary].boundary_points();
    const vector<bool>& midpoints = _boundaries[boundary].midpoints_bits();
    PF_ASSERT(!midpoints[v0_B] || !midpoints[v1_B]);

    const bool are_aa = (pts.col(v0_B) - pts.col(v1_B)).cwiseAbs().minCoeff() < PF_EPS;
    if (!include_midpoints && midpoints[v0_B] && are_aa) {
        v0_B = Circular(_boundaries[boundary].boundary_points(), v0_B - 1);
    }

    if (!include_midpoints && midpoints[v1_B] && are_aa) {
        v1_B = Circular(_boundaries[boundary].boundary_points(), v1_B + 1);
    }

    return BoundaryGraph::make_edge_id(v0_B, v1_B);
}

ImageToBoundaryEdgeMap::ID ImageToBoundaryEdgeMap::make_image_edge_id(
    const int boundary,
    const int v0_B,
    const int v1_B,
    bool& invert_edge
) const {
    PF_ASSERT(boundary >= 0);
    const int v0_I = image_vertex_id_from_point(_boundaries[boundary].boundary_points().col(v0_B));
    const int v1_I = image_vertex_id_from_point(_boundaries[boundary].boundary_points().col(v1_B));
    invert_edge = v0_I > v1_I;
    return BoundaryGraph::make_edge_id(min(v0_I, v1_I), max(v0_I, v1_I));
}

ImageToBoundaryEdgeMap::ID ImageToBoundaryEdgeMap::make_image_edge_id(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1
)const {
    const int v0_I = image_vertex_id_from_point(p0);
    const int v1_I = image_vertex_id_from_point(p1);
    return BoundaryGraph::make_edge_id(min(v0_I, v1_I), max(v0_I, v1_I));
}

void ImageToBoundaryEdgeMap::add_boundary_half_edge(
    const int boundary,
    const ID& e_B
) {
    PF_ASSERT(boundary >= 0);
    const vec2i v_B = BoundaryGraph::unpack_edge_id(e_B);
    bool invert_edge;
    const ID e_I = make_image_edge_id(boundary, v_B(0), v_B(1), invert_edge);
    //PF_DEV_F("Add %d %d -> %lld boundary %d", v_B(0), v_B(1), e_I, boundary);
    EdgePair& edge_pair = _pairs[e_I];
    if (invert_edge) {
        PF_ASSERT(edge_pair.boundary_twin == -1);
        edge_pair.boundary_twin = boundary;
        edge_pair.he_twin = e_B;
    } else {
        PF_ASSERT(edge_pair.boundary == -1);
        edge_pair.boundary = boundary;
        edge_pair.he = e_B;
    }
}

ImageToBoundaryEdgeMap::ID ImageToBoundaryEdgeMap::twin(
    const ID& e_B,
    const int boundary_src,
    const int boundary_dst
) const {
    PF_ASSERT(boundary_src >= 0 && boundary_dst >= 0);
    const vec2i v_B = BoundaryGraph::unpack_edge_id(e_B);
    bool invert_edge;
    const ID e_I = make_image_edge_id(boundary_src, v_B(0), v_B(1), invert_edge);
    const auto find_pair = _pairs.find(e_I);
    if (find_pair == _pairs.end()) {
        return 0;
    }

    const EdgePair& edge_pair = find_pair->second;

    if (boundary_dst == edge_pair.boundary) {
        return edge_pair.he;
    } else if (boundary_dst == edge_pair.boundary_twin) {
        return edge_pair.he_twin;
    } else {
        return 0;
    }
}

#define FLIP(v, b) ((b) ? ((v) == 0 ? 1 : 0) : (v))
int ImageToBoundaryEdgeMap::to_boundary_vertex_from(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const int boundary
) const {
    PF_ASSERT(boundary >= 0);
    const int v0_I = image_vertex_id_from_point(p0);
    const int v1_I = image_vertex_id_from_point(p1);
    const bool invert_edge = v0_I > v1_I;
    const ID e_I = BoundaryGraph::make_edge_id(min(v0_I, v1_I), max(v0_I, v1_I));
    const EdgePair& edge_pair = _pairs.at(e_I);
    return boundary == edge_pair.boundary ?
        BoundaryGraph::unpack_edge_id(edge_pair.he)(FLIP(0, invert_edge)) :
        BoundaryGraph::unpack_edge_id(edge_pair.he_twin)(FLIP(1, invert_edge));
}

int ImageToBoundaryEdgeMap::to_boundary_vertex_to(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    const int boundary
) const {
    PF_ASSERT(boundary >= 0);
    const int v0_I = image_vertex_id_from_point(p0);
    const int v1_I = image_vertex_id_from_point(p1);
    const bool invert_edge = v0_I > v1_I;
    const ID e_I = BoundaryGraph::make_edge_id(min(v0_I, v1_I), max(v0_I, v1_I));
    const EdgePair& edge_pair = _pairs.at(e_I);
    return boundary == edge_pair.boundary ?
        BoundaryGraph::unpack_edge_id(edge_pair.he)(FLIP(1, invert_edge)) :
        BoundaryGraph::unpack_edge_id(edge_pair.he_twin)(FLIP(0, invert_edge));
}

void ImageToBoundaryEdgeMap::get_incident_boundaries(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1,
    int& boundary_0,
    int& boundary_1
) const {
    const int v0_I = image_vertex_id_from_point(p0);
    const int v1_I = image_vertex_id_from_point(p1);
    const ID e_I = BoundaryGraph::make_edge_id(min(v0_I, v1_I), max(v0_I, v1_I));
    const EdgePair& edge_pair = _pairs.at(e_I);
    boundary_0 = edge_pair.boundary;
    boundary_1 = edge_pair.boundary_twin;
}

MultiPolygonTracer::MultiPolygonTracer(
    mc::RasterImageConnectivity& conn_raster,
    mc::SegmentConnectivity& conn_segments,
    const std::vector<Eigen::Vector4d>& boundary_colors,
    const char*                  classifier_uri
) : 
    _conn_raster(conn_raster), 
    _conn_segments(conn_segments),
    _colors(boundary_colors),
    _classifier_uri(classifier_uri) {

    // Load the trained model
    _classifier.load_from_file(_classifier_uri);
}

void MultiPolygonTracer::trace() {
    init_pipeline();
    _initialize_image_boundaries();
    _initialize_segments(); 
    _build_image_to_boundary_map();
    _remove_non_shared_midpoints();
    _extract_junctions();
    _extract_polygon_forcing_junctions();
    _preserve_thickness_at_kopf_junctions();
    _merge_polygons_and_compute_open_segments_info();
    _calculate_coordinate_system_directions_along_open_segments();
    _build_polygon_half_edges_and_fill_junctions_map();
    _classify_polygon_corners_and_junctions();
    _merge_tangent_classifications_along_open_segments();
    _downgrade_tangent_classifications_at_junctions();
    _construct_curves_initial_guess();
    _merge_consecutive_parallel_lines_except_at_junctions();
    _initialize_curve_fitter();
    _couple_curves_along_open_segments();
    _solve_curves();
    _post_process();
    _merge_curves_at_junctions();
}

/*
    Creates a list of unique image boundaries, so that closed image regions do not 
    reference separate boundaries
*/
void MultiPolygonTracer::_initialize_image_boundaries() {
    // 
    // Reset containers
    //
    const int n_polygons = _conn_raster.n_polygons();
    _boundaries.clear();
    _polygon_to_boundary.clear();
    _polygon_to_boundary.resize(n_polygons, -1);

    // Compute a list of sets for each boundary which have duplicated entries 
    vector<unordered_set<int>> boundary_ids(n_polygons);
    mat2xi points;
    for (int i = 0; i < n_polygons; ++i) {
        points = _conn_raster.get_polygon_points(i);

        unordered_set<int> polygon_id;
        polygon_id.reserve(points.cols());
        for (Index j = 0; j < points.cols(); ++j) {
            polygon_id.insert(boundary_vertex_id_from_pos(points.col(j).cast<double>()));
        }
        boundary_ids[i] = move(polygon_id);
    }

    // The color is that of the outer boundary
    vector<Vector4d> colors;
    colors.resize(n_polygons);

    vector<int> hole_polygon_ids;
    int outer_polygon_id;
    for (int i = 0; i < _conn_raster.n_regions(); ++i) {
        _conn_raster.region_to_polygon(i, outer_polygon_id, hole_polygon_ids);

        for (int j = 0; j < (int)boundary_ids.size(); ++j) {
            if (outer_polygon_id == j ) continue;

            if (boundary_ids[outer_polygon_id] == boundary_ids[j]) {
                colors[j] = _colors[outer_polygon_id];
            }
        }

        colors[outer_polygon_id] = _colors[outer_polygon_id];
    }

    // Finding the duplicates and creating the polygon -> ImageRegion map
    for (int i = 0; i < n_polygons; ++i) {
        if (boundary_ids[i].empty()) {
            continue;
        }

        const int boundary_id = _boundaries.size();
        for (int j = i + 1; j < n_polygons; ++j) {
            if (boundary_ids[i] == boundary_ids[j]) {
                boundary_ids[j].clear();
                _polygon_to_boundary[j] = boundary_id;
            }
        }

        _polygon_to_boundary[i] = boundary_id;
        _boundaries.emplace_back(_conn_raster.get_polygon_points(i), boundary_id, colors[i]);
    }

    _outer_boundary = -1;
    for (int i = 0; i < n_polygons; ++i) {
        if (_conn_raster.polygon_to_region(i) == mc::RasterImageConnectivity::infinite_region) {
            _outer_boundary = _polygon_to_boundary[i];
        }
    }


    PF_VERBOSE_F( "Total polygons %d", _boundaries.size() );
    for (int i = 0; i < (int)_polygon_to_boundary.size(); ++i) {
        PF_VERBOSE_F("polygon %d -> %d", i, _polygon_to_boundary[i]);
    }

    // todo: this is a bit hidden ....
    _primitives_to_constraint.clear();
}

void MultiPolygonTracer::_initialize_segments() {
    const int n_segments = _conn_segments.n_segments();
    _segment_closed.resize(n_segments);
    _segment_points.resize(n_segments);
    for (int segment = 0; segment < n_segments; ++segment) {
        // Find all open segments
        bool is_closed;
        _conn_segments.get_segment_points(segment, _segment_points[segment], is_closed);
        _segment_closed[segment] = is_closed;

        if (_segment_closed[segment]) {
            continue;
        }
    }
}

void MultiPolygonTracer::_remove_non_shared_midpoints() {
    _boundaries_closed.clear();
    _boundaries_closed.resize(_boundaries.size(), true);
    
    Eigen::Matrix2Xi segment_pts;
    for (size_t i = 0; i < _conn_segments.n_segments(); ++i) {
        bool is_segment_closed;
        _conn_segments.get_segment_points(i, segment_pts, is_segment_closed);

        if (is_segment_closed) {
            continue;
        }

        int region_0, region_1;
        _conn_segments.get_incident_polygons_for_segment(i, region_0, region_1);
        _boundaries_closed[_polygon_to_boundary[region_0]] = false;
        _boundaries_closed[_polygon_to_boundary[region_1]] = false;
    }

    unordered_map<int, int> midpoint_counter;
    for (const ImageBoundary& boundary : _boundaries) {
        if (_boundaries_closed[boundary.id()]) {
            continue;
        }

        for (Index v : boundary.midpoints()) {
            const Vector2d pm = boundary.boundary_points().col(v);
            const int v_id = _I2B_map->image_vertex_id_from_point(pm);
            ++midpoint_counter[v_id];
        }
    }

    for (ImageBoundary& boundary : _boundaries) {
        if (_boundaries_closed[boundary.id()]) {
            continue;
        }

        vector<Index> midpoints_new;
        for (Index v : boundary.midpoints()) {
            const Vector2d pm = boundary.boundary_points().col(v);
            const int v_id = _I2B_map->image_vertex_id_from_point(pm);
            const int occurrences = midpoint_counter[v_id];
            PF_ASSERT(occurrences <= 2);

            if (occurrences == 1) {
                PF_VERBOSE_F("Removing unique midpoint boundary %d vertex %d", boundary.id(), v);
                boundary.midpoints_bits()[v] = false;
            } else {
                midpoints_new.emplace_back(v);
            }
        }

        boundary.midpoints() = move(midpoints_new);
    }
}

/*
    Add each boundary edge and half-edge to the Image -> Boundary map
*/
void MultiPolygonTracer::_build_image_to_boundary_map() {
    _I2B_map.reset(new ImageToBoundaryEdgeMap(_conn_raster.pixel_vertex_dims().x(), _boundaries));

    const size_t n_boundaries = _boundaries.size();
    for (size_t boundary_id = 0; boundary_id < n_boundaries; ++boundary_id) {
        const auto& pts = _boundaries[boundary_id].boundary_points();
        const auto& midpoints = _boundaries[boundary_id].midpoints_bits();

        for (Index vp_B = 0; vp_B < pts.cols(); ++vp_B) {
            if (midpoints[vp_B]) {
                // Careful of not looping around and considering the same vertex twice.
                if (vp_B == pts.cols() - 1) {
                    break;
                }

                Index vn_B = Circular(pts, vp_B + 1);
                const ImageToBoundaryEdgeMap::ID he_B = _I2B_map->make_boundary_edge_id((int)boundary_id, vp_B, vn_B, true);
                _I2B_map->add_boundary_half_edge((int)boundary_id, he_B);
                vp_B = Circular(pts, vp_B + 1);
            }

            Index vn_B = Circular(pts, vp_B + 1);
            if (midpoints[vn_B]) {
                const ImageToBoundaryEdgeMap::ID he_B = _I2B_map->make_boundary_edge_id((int)boundary_id, vp_B, vn_B, true);
                _I2B_map->add_boundary_half_edge((int)boundary_id, he_B);
                vn_B = Circular(pts, vn_B + 1);
            }

            PF_ASSERT(abs(1. - (pts.col(vn_B) - pts.col(vp_B)).squaredNorm()) < PF_EPS);
            const ImageToBoundaryEdgeMap::ID he_B = _I2B_map->make_boundary_edge_id((int)boundary_id, vp_B, vn_B);
            _I2B_map->add_boundary_half_edge((int)boundary_id, he_B);
        }
    }
}

// Obtain a list of joints to region/continuity mappings from the per-vertex junction information
// If a junction has only two neighboring regions it means that the other region is the outmost region
// (i.e. the background)
void MultiPolygonTracer::_extract_junctions() {
    // Getting the junction type for each vertex in the image
    std::vector<mc::JunctionType> junctions_per_vertex;
    mc::find_junctions(_conn_raster, junctions_per_vertex);
    std::vector<int> junction_vertices;
    
    // Create a set of kopf junctions to resolve boundary <-> image
    for (size_t i = 0; i < junctions_per_vertex.size(); ++i) {
        if (junctions_per_vertex[i] != mc::JUNCTION_INVALID &&
            junctions_per_vertex[i] != mc::JUNCTION_KOPF) {
            junction_vertices.emplace_back((int)i);
        }

        if (junctions_per_vertex[i] == mc::JUNCTION_KOPF) {
            const int v_id = boundary_vertex_id_from_pos(_conn_raster.pixel_vertex_index_to_pos(i).cast<double>());
            PF_ASSERT(_kopf_junctions.find(v_id) == _kopf_junctions.end());
            _kopf_junctions.insert(v_id);
        }
    }

    _kopf_junctions_per_boundary.clear();
    _kopf_junctions_per_boundary.resize(_boundaries.size());
    for (const ImageBoundary& boundary : _boundaries) {
        const auto& pts = boundary.boundary_points();
        for (int i = 0; i < pts.cols(); ++i) {
            const int vertex_id = boundary_vertex_id_from_pos(pts.col(i));
            if (_kopf_junctions.count(vertex_id)) {
                _kopf_junctions_per_boundary[boundary.id()].emplace_back(i);
            }
        }
    }

    // This loop should be merged with the one above...
    for (int i = 0; i < (int)junction_vertices.size(); ++i) {
        JunctionInfo junction;

        // changing the id from pixel-unique to boundary-unique 
        junction.vertex = boundary_vertex_id_from_pos(_conn_raster.pixel_vertex_index_to_pos(junction_vertices[i]).cast<double>());
        junction.boundary_edges.resize(_boundaries.size(), 0x0);

        for (int boundary_id = 0; boundary_id < _boundaries.size(); ++boundary_id) {
            const auto& pts = _boundaries[boundary_id].boundary_points();
            for (int j = 0; j < pts.cols(); ++j) {
                const int vertex = boundary_vertex_id_from_pos(pts.col(j));
                if (vertex == junction.vertex) {
                    const auto e_id = BoundaryGraph::make_edge_id(j, Circular(pts, j + 1));
                    junction.boundary_edges[boundary_id] = e_id;
                    junction.boundary_ids.emplace_back(boundary_id);
                    break;
                }
            }
        }

        // Otherwise it shouldn't be a junction
        PF_ASSERT(junction.boundary_ids.size() > 1);
        _junctions.emplace_back(junction);
    }
}

// Construct the graph for each region, forcing the paths to go through the junction corners
void MultiPolygonTracer::_extract_polygon_forcing_junctions() {
    vecXi B_flat;

    _boundary_junctions.reserve(_boundaries.size());

    for (int boundary_id = 0; boundary_id < _boundaries.size(); ++boundary_id) {
        ImageBoundary& boundary = _boundaries[boundary_id];

        // Boundary data
        auto& E = boundary.edges();
        const auto& B = boundary.boundary_points();
        const auto& C = boundary.convexities();

        // For every junction, find if the current region has a vertex there and remove
        // all the edges which are crossing the junction .
        std::vector<int> junction_vertices;
        B_flat.resize(0);
        for (size_t i = 0; i < _junctions.size(); ++i) {
            if (_junctions[i].boundary_edges[boundary_id] == 0) {
                continue;
            }

            const int boundary_vertex = BoundaryGraph::unpack_edge_id(_junctions[i].boundary_edges[boundary_id])(0);
            if (boundary_vertex != -1) {
                junction_vertices.emplace_back(boundary_vertex);
                if (boundary.is_flat(boundary_vertex)) {
                    MatrixUtils::append(B_flat, boundary_vertex);
                    PF_VERBOSE_F("Boundary %d has flat vertex %d", boundary_id, boundary_vertex);
                }
            }
        }

        // Initialize boundary graph permitting edges from/to junction vertices even if flat
        PF_VERBOSE_F("Boundary %d has %d junctions", boundary_id, B_flat.size());
        boundary.initialize_graph_with_flat_vertices(B_flat);
        boundary.E_delete().clear();

        for (size_t i = 0; i < E.size(); ++i) {
            for (size_t j = 0; j < junction_vertices.size(); ++j) {
                const int junction_vertex = junction_vertices[j];
                if (PathUtils::contains_open(boundary.boundary_points(), boundary.edges()[i].v0, boundary.edges()[i].v1, junction_vertex, true)) {
                    PF_VERBOSE_F("Removing edge %d - %d because it contains junction vertex %d", boundary.edges()[i].v0, boundary.edges()[i].v1, junction_vertex);
                    boundary.E_delete().emplace_back(i);
                    break;
                }
            }
        }
        PF_VERBOSE_F("Tracing boundary %d, removing %d edges", boundary_id, boundary.E_delete().size());
        EraseOrdered(boundary.edges(), boundary.E_delete());

        boundary.trace_polygon();
        _boundary_junctions.emplace_back(move(junction_vertices));
    }
}

// Since a kopf junction cannot be valence 3 or 4, each junction can refere at most to
// 2 different unique segments
struct KopfJunctionSegments {
    // Image edges, we need to store the points and not the id to preserve the order
    // pp_0 and pp_1 are the junction vertices.
    Vector2d pp_0, pn_0;
    Vector2d pp_1, pn_1;

    // Debugging
    int first_segment = -1;
    int second_segment = -1;
};

// Third attempt, I am really tired
void MultiPolygonTracer::_preserve_thickness_at_kopf_junctions() {
#if 0
    {
        DevicePDF* pdf = new DevicePDF("D:/data/polyvec/out/kopf-junctions.svg", 1, 1);
        for (ImageBoundary& boundary : _boundaries) {
            draw_raster(boundary.boundary_points());
        }

        for (const int v_kopf_I : _kopf_junctions) {
            const auto p = boundary_vertex_pos_from_id(v_kopf_I);
            draw::point(p, .5, Style::fill(colors::red));
            draw::text(p, to_string(v_kopf_I), draw::font_pdf / 3, Style::text());
        }

        pdf->draw(0, 0);
        delete pdf;
    }
#endif
    
    // This stuff should be stored somewhere else, it's being recomputed each time
    Eigen::Matrix2Xi segment_points;
    bool is_segment_closed;

    // junction -> image edges
    unordered_map<int, KopfJunctionSegments> kopf_junctions_to_segments;
    
    // Get the segments edges for each kopf junction
    for (int segment_id = 0; segment_id < _conn_segments.n_segments(); ++segment_id) {
        _conn_segments.get_segment_points(segment_id, segment_points, is_segment_closed);

        for (Index i = 0; i < segment_points.cols(); ++i) {
            const int vertex_id = boundary_vertex_id_from_pos(segment_points.col(i).cast<double>());
            if (_kopf_junctions.count(vertex_id)) {
                const Vector2d pp = segment_points.col(i).cast<double>();
                const Vector2d pn = segment_points.col((i + 1) % segment_points.cols()).cast<double>();

                KopfJunctionSegments& segments = kopf_junctions_to_segments[vertex_id];
                if (segments.first_segment == -1) {
                    segments.first_segment = segment_id;
                    segments.pp_0 = pp;
                    segments.pn_0 = pn;
                } else if (segments.second_segment == -1) {
                    segments.second_segment = segment_id;
                    segments.pp_1 = pp;
                    segments.pn_1 = pn;
                } else {
                    PF_ERROR_F("Kopf junction %d has valence > 2", vertex_id);
                    PF_ABORT;
                }
            }
        }
    }

    // Guess what.. the regularization is dependent on whether the polygon is passing through it in the first place
    // 
    // bit set of points that each boundary is traversing.
    vector<unordered_map<int, int>> points_crossed_B(_boundaries.size());
    for (ImageBoundary& boundary : _boundaries) {
        unordered_map<int, int> points;
        for (Index i = 0; i < boundary.polygon_vertices().size(); ++i) {
            const int v_I = _I2B_map->image_vertex_id_from_point(boundary.polygon_points().col(i));
            points[boundary.polygon_vertices()(i)] = i;
        }
        points_crossed_B[boundary.id()] = move(points);
    }

    // for every kopf junction we should have two valid edges, which could potentially correspond to the same
    // segment, but that's not an issue.
    vector<vector<int>> boundary_vertices_to_regularize(_boundaries.size());
    for (const pair<int, KopfJunctionSegments>& kv : kopf_junctions_to_segments) {
        const int v_kopf_I = kv.first;
        const KopfJunctionSegments& segments = kv.second;

        int boundary_00, boundary_01;
        _I2B_map->get_incident_boundaries(segments.pp_0, segments.pn_0, boundary_00, boundary_01);
        int v_00 = -1, v_01 = -1;
        if (boundary_00 != -1) {
            v_00 = _I2B_map->to_boundary_vertex_from(segments.pp_0, segments.pn_0, boundary_00);
        }

        if (boundary_01 != -1) {
            v_01 = _I2B_map->to_boundary_vertex_from(segments.pp_0, segments.pn_0, boundary_01);
        }

        int boundary_10, boundary_11;
        _I2B_map->get_incident_boundaries(segments.pp_1, segments.pn_1, boundary_10, boundary_11);
        int v_10 = -1, v_11 = -1;
        if (boundary_10 != -1) {
            v_10 = _I2B_map->to_boundary_vertex_from(segments.pp_1, segments.pn_1, boundary_10);
        }
        if (boundary_11 != -1) {
            v_11 = _I2B_map->to_boundary_vertex_from(segments.pp_1, segments.pn_1, boundary_11);
        }

        // We want to regularize boundaries only if their polygons are going through them
        const bool segment_0_crossing =
            (boundary_00 == -1 ? false : (points_crossed_B[boundary_00].count(v_00) > 0)) ||
            (boundary_01 == -1 ? false : (points_crossed_B[boundary_01].count(v_01) > 0));
        const bool segment_1_crossing =
            (boundary_10 == -1 ? false : (points_crossed_B[boundary_10].count(v_10) > 0)) ||
            (boundary_11 == -1 ? false : (points_crossed_B[boundary_11].count(v_11) > 0));

        if (!segment_0_crossing && !segment_1_crossing) {
            continue;
        }

        // Figuring out which of the two sides we need to regularize
        PF_VERBOSE_F("Kopf junction %d", v_kopf_I);
        PF_VERBOSE_F("boundary_00 %d vertex %d", boundary_00, v_00);
        PF_VERBOSE_F("boundary_01 %d vertex %d", boundary_01, v_01);
        PF_VERBOSE_F("boundary_10 %d vertex %d", boundary_10, v_10);
        PF_VERBOSE_F("boundary_11 %d vertex %d", boundary_11, v_11);
        PF_VERBOSE_F("segment_0_crossing %d segment_1_crossing %d", segment_0_crossing, segment_1_crossing);

        bool regularize_0 = false, regularize_1 = false;
        // If only one of two segment boundaries is going through the junction then we are regularizing the other one
        if (segment_0_crossing && !segment_1_crossing) {
            regularize_1 = true;
        }
        else if (segment_1_crossing && !segment_0_crossing) {
            regularize_0 = true;
        }
        // Which one of the two segments is the outside and which one is the inside, if the segments correspond to the different 
        // boundaries than we can use the area as a measure of which shape "prevails" over the other
        else if (kv.second.first_segment != kv.second.second_segment) {
            double segment_0_area = 0., segment_1_area = 0.;
            segment_0_area = max(segment_0_area, boundary_00 != -1 ? _boundaries[boundary_00].bounding_box_area() : 0.);
            segment_0_area = max(segment_0_area, boundary_01 != -1 ? _boundaries[boundary_01].bounding_box_area() : 0.);
            segment_1_area = max(segment_1_area, boundary_10 != -1 ? _boundaries[boundary_10].bounding_box_area() : 0.);
            segment_1_area = max(segment_1_area, boundary_11 != -1 ? _boundaries[boundary_11].bounding_box_area() : 0.);

            if (abs(segment_0_area - segment_1_area) < PF_EPS_MEDIUM) {
                regularize_0 = regularize_1 = true;
            }else {
                regularize_0 = segment_0_area > segment_1_area;
                regularize_1 = segment_1_area > segment_0_area;
            }
        }
        // If the segments are part of the same boundary, we use some vodoo ifs to figure out which side should be regularized
        // We prefer to regularize the side which has the same number of options, if both sides have the same number of corners
        // then we regularize both
        else {
            PF_ASSERT((boundary_00 != -1) || (boundary_01 != -1));
            PF_ASSERT((boundary_10 != -1) || (boundary_11 != -1));
            const vector<int>& C_0 = _boundaries[boundary_00 != -1 ? boundary_00 : boundary_01].convexities();
            const vector<int>& C_1 = _boundaries[boundary_10 != -1 ? boundary_10 : boundary_11].convexities();

            const int v_0 = v_00 != -1 ? v_00 : v_01;
            const int v_1 = v_10 != -1 ? v_10 : v_11;

            const int corners_near_0 = (CircularAt(C_0, v_0 - 1) != 0 ? 1 : 0) + (CircularAt(C_0, v_0 + 1) != 0 ? 1 : 0);
            const int corners_near_1 = (CircularAt(C_1, v_1 - 1) != 0 ? 1 : 0) + (CircularAt(C_1, v_1 + 1) != 0 ? 1 : 0);
            if (corners_near_0 == corners_near_1) {
                regularize_0 = regularize_1 = true;
            } else {
                regularize_0 = corners_near_0 > corners_near_1;
                regularize_1 = corners_near_1 > corners_near_0;
            }
        }

        PF_VERBOSE_F("regularize_0 %d regularize_1 %d", regularize_0, regularize_1);
        if (regularize_0) {
            if (boundary_00 != -1) {
                boundary_vertices_to_regularize[boundary_00].emplace_back(v_00);
            }
            if (boundary_01 != -1) {
                boundary_vertices_to_regularize[boundary_01].emplace_back(v_01);
            }
        }
        if (regularize_1) {
            if (boundary_10 != -1) {
                boundary_vertices_to_regularize[boundary_10].emplace_back(v_10);
            }
            if (boundary_11 != -1) {
                boundary_vertices_to_regularize[boundary_11].emplace_back(v_11);
            }
        }
    }

    // Time to regularize the boundary graphs!
    unordered_set<int> is_kopf_junction_for_this_boundary;
    for (ImageBoundary& boundary : _boundaries) {
        const auto& B = boundary.boundary_points();
        const auto& C = boundary.convexities();
        auto& E = boundary.edges();

        is_kopf_junction_for_this_boundary.clear();

        for (const int v_kopf : boundary_vertices_to_regularize[boundary.id()]) {
            is_kopf_junction_for_this_boundary.insert(v_kopf);
        }

        boundary.E_delete().clear();
        for (size_t i = 0; i < E.size(); ++i) {
            for (int v_kopf : boundary_vertices_to_regularize[boundary.id()]) {
                const int v_kopf_prev = Circular(B, v_kopf - 1);
                const int v_kopf_next = Circular(B, v_kopf + 1);

                if (!PathUtils::contains_closed(B, E[i].v0, E[i].v1, v_kopf, true)) {
                    continue;
                }

                if (E[i].v0 == v_kopf_prev && E[i].v1 != v_kopf_next) {
                    const vec2 dir_junction = B.col(v_kopf) - B.col(v_kopf_prev);
                    const vec2 dir_poly = B.col(E[i].v1) - B.col(v_kopf_prev);
                    const double junction_angle = AngleUtils::spanned_shortest_between(dir_junction, dir_poly);
                    if (junction_angle > PF_RAD(44)) {
                        continue;
                    }
                }

                if (E[i].v1 == v_kopf_next && E[i].v0 != v_kopf_prev) {
                    const vec2 dir_junction = B.col(v_kopf) - B.col(v_kopf_next);
                    const vec2 dir_poly = B.col(E[i].v0) - B.col(v_kopf_next);
                    const double junction_angle = AngleUtils::spanned_shortest_between(dir_junction, dir_poly);
                    if (junction_angle > PF_RAD(44)) {
                        continue;
                    }
                }

                boundary.E_delete().emplace_back(i);
                break;
            }
        }

        EraseOrdered(E, boundary.E_delete());

        auto add_edge_safe = [&](const BoundaryGraph::Edge& e) { if (e.v0 != e.v1) { E.emplace_back(e); } };

        for (int v_kopf : boundary_vertices_to_regularize[boundary.id()]) {
            const int v_kopf_prev = Circular(B, v_kopf - 1);
            const int v_kopf_next = Circular(B, v_kopf + 1);
            add_edge_safe(BoundaryGraph::Edge::make_pixel_diagonal(v_kopf_prev, v_kopf_next));

            if (C[v_kopf_prev] == 0 || C[v_kopf_next] == 0) {
                if (C[v_kopf_prev] == 0) {
                    const int v_kopf_prev_prev = PathUtils::next_transition_vertex(boundary.convexities(), v_kopf_prev, -1);
                    const int v_kopf_prev_prev_I = boundary_vertex_id_from_pos(B.col(v_kopf_prev_prev));
                    if (!is_kopf_junction_for_this_boundary.count(v_kopf_prev_prev_I)) {
                        add_edge_safe(BoundaryGraph::Edge::make_flat(v_kopf_prev_prev, v_kopf_prev));
                        PF_VERBOSE_F("Adding edge %d %d case 0", E.back().v0, E.back().v1);
                    } else {
                        add_edge_safe(BoundaryGraph::Edge::make_flat(Circular(B, v_kopf_prev_prev + 1), v_kopf_prev));
                        PF_VERBOSE_F("Adding edge %d %d case 1", E.back().v0, E.back().v1);
                    }
                }

                if (C[v_kopf_next] == 0) {
                    const int v_kopf_next_next = PathUtils::next_transition_vertex(boundary.convexities(), v_kopf_next, +1);
                    add_edge_safe(BoundaryGraph::Edge::make_flat(v_kopf_next, v_kopf_next_next));
                    PF_VERBOSE_F("Adding edge %d %d case 2", E.back().v0, E.back().v1);
                    add_edge_safe(BoundaryGraph::Edge::make_flat(v_kopf_next, Circular(B, v_kopf_next_next - 1)));
                    PF_VERBOSE_F("Adding edge %d %d case 3", E.back().v0, E.back().v1);
                }
            }
        }
        
#if 0
        if (boundary.id() == 1) {
            DevicePDF* pdf = new DevicePDF("D:/data/polyvec/out/graph.svg", 1, 1);
            draw_raster(B);
            for (size_t i = 0; i < boundary.edges().size(); ++i) {
                const auto& e = boundary.edges()[i];
                draw::line(B.col(e.v0), B.col(e.v1), Style::outline(colors::forest_green, .5));
                PF_VERBOSE_F("%d -> %d", e.v0, e.v1);
            }
            draw_raster_indices(B);

            pdf->draw(0, 0);
            delete pdf;
        }
#endif

        boundary.trace_polygon();
    }
}

void MultiPolygonTracer::_merge_polygons_and_compute_open_segments_info() {
    // Subpaths
    vecXi subP0, subP1;
    
    // Holds the minimum subpath vertices remapped on the other boundary
    vecXi subP_mapped;
    
    _open_segments.clear();
    _open_segments.reserve(_conn_segments.n_segments());

    // For the polygons that we want to retrace we need to reset the graph.
    std::vector<bool> has_graph_been_reset(_boundaries.size(), false);

    for (int segment_id = 0; segment_id < _conn_segments.n_segments(); ++segment_id) {
        if (_segment_closed[segment_id]) {
            continue;
        }

        OpenSegmentInfo segment;
        segment.id = segment_id;
        
        bool is_segment_closed;
        _conn_segments.get_segment_points(segment_id, segment.points, is_segment_closed);

        // Get the two regions which are sharing the boundary
        int region_0, region_1;
        _conn_segments.get_incident_polygons_for_segment(segment.id, region_0, region_1);
        const int boundary_0 = _polygon_to_boundary[region_0];
        const int boundary_1 = _polygon_to_boundary[region_1];
        PF_VERBOSE_F("merge segment %d between %d %d", segment.id, boundary_0, boundary_1);

        // Find the two respective subpaths and partial energies
        const auto B_bounds_0 = _get_boundary_subpath_for_segment(segment.id, boundary_0);
        const auto B_bounds_1 = _get_boundary_subpath_for_segment(segment.id, boundary_1);
        const pair<int, int> v_bounds0 = _get_polygon_subpath_for_segment(B_bounds_0.first, B_bounds_0.second, boundary_0, &subP0);
        const pair<int, int> v_bounds1 = _get_polygon_subpath_for_segment(B_bounds_1.first, B_bounds_1.second, boundary_1, &subP1);
        PF_ASSERT(v_bounds0.first != -1 && v_bounds0.second != -1);
        PF_ASSERT(v_bounds1.first != -1 && v_bounds1.second != -1);

        const double E0 = _boundaries[boundary_0].calculate_subpath_energy(v_bounds0.first, v_bounds0.second);
        const double E1 = _boundaries[boundary_1].calculate_subpath_energy(v_bounds1.first, v_bounds1.second);
        PF_ASSERT(E0 >= 0.);
        PF_ASSERT(E1 >= 0.);

        // Preserving the volume if possible, but preferring lower energy by default
        bool prefer_0;
        if (subP0.size() == 2 && subP1.size() > 2) {
            prefer_0 = false;
        } else if (subP0.size() > 2 && subP1.size() == 2) {
            prefer_0 = true;
        } else {
            prefer_0 = E0 < E1;
        }

        if (prefer_0) {
            segment.boundary_src = boundary_0;
            segment.boundary_dst = boundary_1;
            segment.B_bounds_src = B_bounds_0;
            segment.B_bounds_dst = B_bounds_1;
            segment.direction_src = 1;
            segment.direction_dst = -1;
        } else {
            segment.boundary_src = boundary_1;
            segment.boundary_dst = boundary_0;
            segment.B_bounds_src = B_bounds_1;
            segment.B_bounds_dst = B_bounds_0;
            segment.direction_src = -1;
            segment.direction_dst = +1;
        }
        PF_VERBOSE_F("Preferring boundary %d", segment.boundary_src);

        // Replace the other region's subpath with the one with minimum energy
        // If they have the same energy, either of them will do (they could still be different!)
        vecXi& subP_src = prefer_0 ? subP0 : subP1;
        vecXi& subP_dst = prefer_0 ? subP1 : subP0;
        const int boundary_idx_src = prefer_0 ? boundary_0 : boundary_1;
        const int boundary_idx_dst = prefer_0 ? boundary_1 : boundary_0;
        ImageBoundary& boundary_src = _boundaries[boundary_idx_src];
        ImageBoundary& boundary_dst = _boundaries[boundary_idx_dst];

        PF_VERBOSE_F("Polygon %d energy %f", boundary_0, E0);
        PF_VERBOSE_F("Polygon %d energy %f", boundary_1, E1);

        const mat2x& pts_src = boundary_src.boundary_points();
        const mat2x& pts_dst = boundary_dst.boundary_points();

        subP_mapped.resize(subP_src.size());
        for (Index i = 0; i < subP_src.size(); ++i) {
            int vertex_src_mapped;
            if (i == 0) {
                auto e_id_B = _I2B_map->make_boundary_edge_id(boundary_idx_src, subP_src(i), Circular(pts_src, subP_src(i) + 1), true);
                auto e_id_B_twin = _I2B_map->twin(e_id_B, boundary_idx_src, boundary_idx_dst);
                if (e_id_B_twin == 0) {
                    e_id_B = _I2B_map->make_boundary_edge_id(boundary_idx_src, subP_src(i), Circular(pts_src, subP_src(i) + 2), true);
                    e_id_B_twin = _I2B_map->twin(e_id_B, boundary_idx_src, boundary_idx_dst);
                    PF_ASSERT(e_id_B_twin != 0);
                }

                subP_mapped(i) = BoundaryGraph::unpack_edge_id(e_id_B_twin)(1);
            } else if (i == subP_src.size() - 1) {
                auto e_id_B = _I2B_map->make_boundary_edge_id(boundary_idx_src, Circular(pts_src, subP_src(i) - 1), subP_src(i), true);
                auto e_id_B_twin = _I2B_map->twin(e_id_B, boundary_idx_src, boundary_idx_dst);
                if (e_id_B_twin == 0) {
                    e_id_B = _I2B_map->make_boundary_edge_id(boundary_idx_src, Circular(pts_src, subP_src(i) - 2), subP_src(i), true);
                    e_id_B_twin = _I2B_map->twin(e_id_B, boundary_idx_src, boundary_idx_dst);
                    PF_ASSERT(e_id_B_twin != 0);
                }

                subP_mapped(i) = BoundaryGraph::unpack_edge_id(e_id_B_twin)(0);
            } else {
                auto e_id_B = _I2B_map->make_boundary_edge_id(boundary_idx_src, subP_src(i), Circular(pts_src, subP_src(i) + 1), true);
                auto e_id_B_twin = _I2B_map->twin(e_id_B, boundary_idx_src, boundary_idx_dst);
                if (e_id_B_twin == 0) {
                    e_id_B = _I2B_map->make_boundary_edge_id(boundary_idx_src, subP_src(i), Circular(pts_src, subP_src(i) + 2), true);
                    e_id_B_twin = _I2B_map->twin(e_id_B, boundary_idx_src, boundary_idx_dst);
                    PF_ASSERT(e_id_B_twin != 0);
                }

                subP_mapped(i) = BoundaryGraph::unpack_edge_id(e_id_B_twin)(1);
            }
            PF_VERBOSE_F("mapped corner %d(%d) -> %d(%d)", subP_src(i), boundary_idx_src, subP_mapped(i), boundary_idx_dst);
        }

        subP_mapped.reverseInPlace();
        
        // the path is already shared as it's going through all the junctions
        if (subP_mapped.size() == 0 || subP_dst.size() == 0) {
            continue;
        }

        if (!has_graph_been_reset[boundary_0]) {
            _boundaries[boundary_0].reset_graph();
            has_graph_been_reset[boundary_0] = true;
        }

        if (!has_graph_been_reset[boundary_1]) {
            _boundaries[boundary_1].reset_graph();
            has_graph_been_reset[boundary_1] = true;
        }

        // Ordered interval limits
        boundary_dst.force_graph_subpath(subP_mapped);

        // The same way we have forced the subpath in the destination graph we need to do it in the source graph
        // In the case the source region is the destination region for another segment the path might change.
        boundary_src.force_graph_subpath(subP_src);

        _open_segments.emplace_back(move(segment));
    }

    // Obtaining new polygons
    for (ImageBoundary& boundary : _boundaries) {
        if (has_graph_been_reset[boundary.id()]) {
            // forcing flat edges to go through the flat vertex
            PF_VERBOSE_F("Retracing boundary %d", boundary.id());
            boundary.remove_edges_crossing_vertices(_boundary_junctions[boundary.id()]);
            boundary.trace_polygon_simple();
        }
    };
}

/*
    In a polygon with multiple colored regions, continuations are only detected between 
    polygons facing inside the same region or between shapes themselves.

    Detecting regularities between multiple regions is very unstable.
*/
void MultiPolygonTracer::_construct_regularity_graph() {
    for (ImageBoundary& boundary : _boundaries) {
        boundary.regularity_graph().reset();
		boundary.regularity_graph().update(
            boundary.boundary_points(),
            boundary.polygon_vertices(),
            boundary.raster_symmetries(),
            boundary.raster_symmetries_local(),
            true
        );
    }
}

/*
    For each junction and incident region determines the continuity of the fit in order to 
    compute additional fitting objectives.
    It also classifies the non-shared polygon corners.
*/
void MultiPolygonTracer::_classify_polygon_corners_and_junctions() {
    for (ImageBoundary& boundary : _boundaries) {
        boundary.classify_corners(_classifier);
    }
}

/*
    Making sure that corners which are shared between two segments have the same tangent classification.
    We use the classifications from the side which has been preferred during the polygon 
*/
void MultiPolygonTracer::_merge_tangent_classifications_along_open_segments() {
    vecXi subpath_src, subpath_dst;

    for (const OpenSegmentInfo& segment : _open_segments) {
        segment.validate();
        subpath_src = _boundaries[segment.boundary_src].get_polygon_subpath(segment.B_bounds_src.first, segment.B_bounds_src.second);
        subpath_dst = _boundaries[segment.boundary_dst].get_polygon_subpath(segment.B_bounds_dst.first, segment.B_bounds_dst.second);

        // If the polygon merging and the segment info extraction worked correctly, the two paths will be identical
        PF_ASSERT(subpath_src.size() == subpath_dst.size());
        const int n_points = subpath_src.size();

        // We also expect the tangent fit classification to have already happened
        ImageBoundary& boundary_src = _boundaries[segment.boundary_src];
        ImageBoundary& boundary_dst = _boundaries[segment.boundary_dst];

        subpath_dst.reverseInPlace();

        PF_VERBOSE_F("SEGMENT %d (%d -> %d)", segment.id, segment.boundary_src, segment.boundary_dst);
        for (Index i = 0; i < subpath_src.size(); ++i) {
            PF_VERBOSE_F("%d %d", subpath_src(i), subpath_dst(i));
        }

        // We want to exclude junctions
        if (n_points < 3) {
            continue;
        }

        for (int i = 1; i < n_points - 1; ++i) {
            const int v_src = subpath_src(i);
            const int v_dst = subpath_dst(i);
            //PF_DEV_F("Polygon %d vertex %d from %s to %s (copied from polygon %d vertex %d)", segment.boundary_dst, v_dst, 
            //    tangent_fit_type_to_string(boundary_dst.tangents_fits()[v_dst]), tangent_fit_type_to_string(boundary_src.tangents_fits()[v_src]),
            //    segment.boundary_src, v_src);
            boundary_dst.tangents_fits()[v_dst] = boundary_src.tangents_fits()[v_src];
        }
    }
}

/*
    Makes sures that at every junction at least one of the regions has been classified as G^0.
*/
void MultiPolygonTracer::_downgrade_tangent_classifications_at_junctions() {

	struct AngleInfo
	{
		double angle;
		int boundary_id;
		int v_polygon;

		AngleInfo(double angle, int boundary_id, int v_polygon)
			: angle(angle), boundary_id(boundary_id), v_polygon(v_polygon)
		{ }

		bool operator<(const AngleInfo& other) const
		{
			if (angle != other.angle)
				return angle < other.angle;
			return boundary_id < other.boundary_id;
		}
	};

    std::vector<AngleInfo> angles; // only angles of smooth fit types

    for (JunctionInfo& junction : _junctions) {
        PF_VERBOSE_F("Junction vertex %d", junction.vertex);

		// Gather triplets and pairs of regions where we need to downgrade one
		std::vector<std::tuple<int, int, int>> need_to_downgrade_one_of_three;
		std::vector<std::pair<int, int>> need_to_downgrade_one_of_two;

		struct RegionInfo
		{			
			const ImageBoundary& boundary;
			const int v_boundary;
			const int v_polygon;
			TangentFitType fit_type;
		};

		auto get_region_info = [&](int boundary_id)
		{
			const ImageBoundary& boundary = _boundaries[boundary_id];

			PF_ASSERT(junction.boundary_edges[boundary_id] != 0);
			const int v_boundary = BoundaryGraph::unpack_edge_id(junction.boundary_edges[boundary_id])(0);
			PF_ASSERT(v_boundary != -1);
			const int v_polygon = boundary.find_polygon_vertex_for_boundary_vertex(v_boundary);
			PF_ASSERT(v_polygon != -1);

			return RegionInfo
			{
				boundary, v_boundary, v_polygon, boundary.tangents_fits()[v_polygon]
			};
		};

        junction.polygon_vertices.clear();

		// incoming to the junction vertex
		auto junction_halfedge = _junctions_to_half_edges.at(junction.vertex);
		
        angles.clear();
		for (int i = 0; i < junction.boundary_ids.size(); ++i) {
			auto boundary_id = junction.boundary_ids[i];
			while (_he_polygon.boundary(junction_halfedge) != boundary_id)
				junction_halfedge = _he_polygon.twin(_he_polygon.next(junction_halfedge));

			//        \  b  /
			//         \   /
			//           j
			//     left  |  right
			//           |
			//

			auto boundary = get_region_info(boundary_id);
			auto left_id = _he_polygon.boundary(_he_polygon.twin(_he_polygon.next(junction_halfedge)));
			auto right_id = _he_polygon.boundary(_he_polygon.twin(junction_halfedge));

			auto left_neighbor = get_region_info(left_id);
			auto right_neighbor = get_region_info(right_id);			
			
			std::cout << "Boundary " << boundary_id << " (vertex " << junction.vertex << ", left: " << left_id << ", right: " << right_id << std::endl;

			const vec2 pp = CircularAt(boundary.boundary.polygon_points(), boundary.v_polygon - 1);
			const vec2 p = boundary.boundary.polygon_points().col(boundary.v_polygon);
			const vec2 pn = CircularAt(boundary.boundary.polygon_points(), boundary.v_polygon + 1);

			// A region with smooth transition cannot have two smooth neighbors
			if (boundary.fit_type != TANGENT_FIT_CONSTANT && left_neighbor.fit_type != TANGENT_FIT_CONSTANT && right_neighbor.fit_type != TANGENT_FIT_CONSTANT)
			{
				std::cout << "One of the following need to be downgraded around vertex " << junction.vertex << ": " << left_id << ", " << boundary_id << ", " << right_id << std::endl;
				need_to_downgrade_one_of_three.emplace_back(left_id, boundary_id, right_id);
			}
			
			// Two smooth neighbors are only allowed across axis-aligned edges			
			bool is_axis_aligned = (pn - p).cwiseAbs().minCoeff() < PF_EPS && (pn - p).squaredNorm() > 4 - PF_EPS;
			if (!is_axis_aligned && boundary.fit_type != TANGENT_FIT_CONSTANT && left_neighbor.fit_type != TANGENT_FIT_CONSTANT)
			{
				std::cout << "One of the following need to be downgraded around vertex " << junction.vertex << ": " << left_id << ", " << boundary_id << std::endl;
				need_to_downgrade_one_of_two.emplace_back(boundary_id, left_id);
			}
        
            //PF_DEV_F("Boundary id %d vertex %d pp %f %f p %f %f pn %f %f", boundary_id, v_polygon, pp(0, 0), pp(1, 0), p(0, 0), p(1, 0), pn(0, 0), pn(1, 0));

            junction.polygon_vertices.emplace_back(boundary.v_polygon);
			if(boundary.fit_type != TANGENT_FIT_CONSTANT)
				angles.emplace_back(AngleUtils::spanned_shortest(pp, p, pn), boundary_id, boundary.v_polygon);
        }

		// sort from acute to obtuse
		std::sort(angles.begin(), angles.end());
		for (auto it = angles.begin(); it != angles.end(); ++it)
		{
			if (need_to_downgrade_one_of_three.empty() && need_to_downgrade_one_of_two.empty())
				break; // no more action necessary

			//downgrade this angle
			ImageBoundary& boundary = _boundaries[it->boundary_id];
			PF_VERBOSE_F("Downgrading polygon %d vertex %d", boundary.id(), it->v_polygon);
			boundary.tangents_fits()[it->v_polygon] = TANGENT_FIT_CONSTANT;

			need_to_downgrade_one_of_two.erase(std::remove_if(
				need_to_downgrade_one_of_two.begin(), 
				need_to_downgrade_one_of_two.end(), 
				[&](const std::pair<int, int>& p) 
			{
				return p.first == it->boundary_id || p.second == it->boundary_id;
			}), need_to_downgrade_one_of_two.end());

			need_to_downgrade_one_of_three.erase(std::remove_if(
				need_to_downgrade_one_of_three.begin(),
				need_to_downgrade_one_of_three.end(),
				[&](const std::tuple<int, int, int>& p)
			{
				return std::get<0>(p) == it->boundary_id || std::get<1>(p) == it->boundary_id || std::get<2>(p) == it->boundary_id;
			}), need_to_downgrade_one_of_three.end());
		}		        
    }
}

/*
    For each polygon shared by two boundaries, calculates which curve parameterization's need to 
    invert their front/back coordinate systems in order to couple the opposite curve parameters.

    Since the open segments are not overlapping, flipping either of the coordinate systems is identical.
*/
void MultiPolygonTracer::_calculate_coordinate_system_directions_along_open_segments() {
    for (OpenSegmentInfo& segment : _open_segments) {
        segment.validate();
        ImageBoundary& boundary_src = _boundaries[segment.boundary_src];
        ImageBoundary& boundary_dst = _boundaries[segment.boundary_dst];

        const auto v_bounds_src = boundary_src.get_polygon_subpath_bounds(segment.B_bounds_src.first, segment.B_bounds_src.second);
        const auto v_bounds_dst = boundary_dst.get_polygon_subpath_bounds(segment.B_bounds_dst.first, segment.B_bounds_dst.second);

        PF_VERBOSE_F("Calculate coordinate system sign for boundaries %d %d", segment.boundary_src, segment.boundary_dst);
        PF_VERBOSE_F("boundary_src has bounds %d %d", v_bounds_src.first, v_bounds_src.second);
        PF_VERBOSE_F("boundary_dst has bounds %d %d", v_bounds_dst.first, v_bounds_dst.second);

        const int n_vertices_src = CircularDist(boundary_src.polygon_vertices(), v_bounds_src.first, v_bounds_src.second);
        const int n_vertices_dst = CircularDist(boundary_dst.polygon_vertices(), v_bounds_dst.first, v_bounds_dst.second);
        PF_ASSERT(n_vertices_src == n_vertices_dst);

        for (int i = 0; i < n_vertices_src; ++i) {
            const int v_dst = Circular(boundary_dst.polygon_points(), (Index)(v_bounds_dst.first + i));
            const int v_dst_next = Circular(boundary_dst.polygon_points(), (Index)(v_dst + 1));
            PF_VERBOSE_F("Flipping edge %d %d boundary %d", v_dst, v_dst_next, boundary_dst.id());
            boundary_dst.flip_coordinate_system(BoundaryGraph::make_edge_id(v_dst, v_dst_next));
        }
    }
}

/*
    Initializes the curves types from the tangent classifications
*/
void MultiPolygonTracer::_construct_curves_initial_guess() {
    for (ImageBoundary& boundary : _boundaries) {
        boundary.fit_curves_initial_guess();
    }
}

void MultiPolygonTracer::_initialize_curve_fitter() {
    const vec2i dims = _conn_raster.pixel_vertex_dims();
    const double aabb_diagonal = dims.cast<double>().norm();
    _curve_fitter.reset(new CurveFitter(aabb_diagonal, true));

    // Insert all the curves from all the region in the optimizer
    for (ImageBoundary& boundary : _boundaries) {
        auto& prims = boundary.spline().primitives;
		auto is_edge_important = get_important_or_axis_aligned_edges(boundary.polygon_points(), boundary.regularity_graph());
		auto edge_angles = find_edge_angles_from_parallel(boundary.polygon_points(), boundary.regularity_graph());
		addCurveSequenceToFitter(*_curve_fitter, boundary.spline(), true, 0, prims.size(), boundary.polygon_points().cols(), is_edge_important, edge_angles, true);
    }    

    // Adding the constraints recorded during the previous steps
    for (size_t i = 0; i < _primitives_to_constraint.size(); ++i) {
        _primitives_to_constraint[i].first->curve->reduce_degrees_of_freedom(_primitives_to_constraint[i].second);
    }
}

/*
    Compiles a list of true/false values for each edge and polygon which indicates
    which coordinate systems need to be inverted in order for the coupling to be successful.

    Note that at this point the subpaths lying in shared segments are assumed to be the same.
*/
void MultiPolygonTracer::_couple_curves_along_open_segments() {
    for (const OpenSegmentInfo& segment : _open_segments) {
        segment.validate();

        PF_VERBOSE_F("Segment %d is shared between boundary %d and %d", segment.id, segment.boundary_src, segment.boundary_dst);

        // src/dst is not relevant here, we just want to couple the segments
        ImageBoundary& boundary_0 = _boundaries[segment.boundary_src];
        ImageBoundary& boundary_1 = _boundaries[segment.boundary_dst];

        VectorXi subP_0 = boundary_0.get_polygon_subpath(segment.B_bounds_src.first, segment.B_bounds_src.second);
        VectorXi subP_1 = boundary_1.get_polygon_subpath(segment.B_bounds_dst.first, segment.B_bounds_dst.second);

        // Making sure that the two segments are identical and have the same corner classifications
        if (segment.direction_src < 0) {
            subP_0.reverseInPlace();
        }
        if (segment.direction_dst < 0) {
            subP_1.reverseInPlace();
        }
        
        PF_ASSERT(subP_0.size() == subP_1.size());
        PF_VERBOSE_F("Coupling boundaries %d %d on segment %d", segment.boundary_src, segment.boundary_dst, segment.id);

        for (Index i = 1; i < subP_0.size() - 1; ++i) {
            const int v_0 = subP_0(i);
            const int v_1 = subP_1(i);
            const TangentFitType tangent_fit_0 = boundary_0.tangents_fits()[v_0];
            const TangentFitType tangent_fit_1 = boundary_1.tangents_fits()[v_1];
            PF_VERBOSE_F("Expected same tangent fits boundary %d vertex %d tangents %s", boundary_0.id(), v_0, tangent_fit_type_to_string(tangent_fit_0));
            PF_VERBOSE_F("Expected same tangent fits boundary %d vertex %d tangents %s", boundary_1.id(), v_1, tangent_fit_type_to_string(tangent_fit_1));

            PF_ASSERT(boundary_0.tangents_fits()[v_0] == boundary_1.tangents_fits()[v_1]);
        }

        PF_VERBOSE_F("Open Segment between %d and %d", segment.boundary_src, segment.boundary_dst);
        for (Index i = 0; i < subP_0.size(); ++i) {
            PF_VERBOSE_F("Couple corners %d - %d", subP_0(i), subP_1(i));
        }

        vector<CurvePrimitive>& curves_0 = boundary_0.spline().primitives;
        vector<CurvePrimitive>& curves_1 = boundary_1.spline().primitives;
        const auto curve_interval_0 = boundary_0.get_primitives_interval_closed(subP_0(0), subP_0(subP_0.size() - 1), segment.direction_src);
        const auto curve_interval_1 = boundary_1.get_primitives_interval_closed(subP_1(0), subP_1(subP_1.size() - 1), segment.direction_dst);

        PF_VERBOSE_F("curve_interval_0 %d %d", curve_interval_0.first, curve_interval_0.second);
        PF_VERBOSE_F("curve_interval_1 %d %d", curve_interval_1.first, curve_interval_1.second);
        
        const int n_curves_0 = segment.direction_src > 0 ?
            CircularDist(curves_0, curve_interval_0.first, curve_interval_0.second) :
            CircularDist(curves_0, curve_interval_0.second, curve_interval_0.first);

        const int n_curves_1 = segment.direction_dst > 0 ?
            CircularDist(curves_1, curve_interval_1.first, curve_interval_1.second) :
            CircularDist(curves_1, curve_interval_1.second, curve_interval_1.first);

        if (n_curves_0 != n_curves_1) {
            PF_VERBOSE_F("Mismatching curves curve_interval_0: %d %d (%d) curve_interval_1: %d %d (%d)",
                curve_interval_0.first, curve_interval_0.second, segment.direction_src,
                curve_interval_1.first, curve_interval_1.second, segment.direction_dst);
            PF_ABORT;
        }

        PF_ASSERT(n_curves_0 == n_curves_1);
        PF_VERBOSE_F("Number of primitives to couple %d", n_curves_0);

        for (int i = 0; i < n_curves_0; ++i) {
            const int curve_idx_0 = Circular(curves_0, curve_interval_0.first + i * segment.direction_src);
            const int curve_idx_1 = Circular(curves_1, curve_interval_1.first + i * segment.direction_dst);
            auto& curve_0 = curves_0[curve_idx_0].curve;
            auto& curve_1 = curves_1[curve_idx_1].curve;

            if (i > 0 && i < n_curves_0 - 1) {
                PF_ASSERT(curve_0->get_curve()->get_type() == curve_1->get_curve()->get_type());
            }

            PF_VERBOSE_F("Coupling curves %d(%d) - %d(%d)", curve_idx_0, boundary_0.id(), curve_idx_1, boundary_1.id());
            
            _curve_fitter->make_overlap_start_to_end(curve_0.get(), curve_1.get());
            _curve_fitter->make_overlap_end_to_start(curve_0.get(), curve_1.get());
        }
    }
}

void MultiPolygonTracer::_merge_consecutive_parallel_lines_except_at_junctions() {
    vector<vector<bool>> junctions_per_boundary(_boundaries.size());
    for (const ImageBoundary& boundary : _boundaries) {
        junctions_per_boundary[boundary.id()].resize(boundary.polygon_vertices().size());
    }

    for (const JunctionInfo& junction : _junctions) {
        PF_ASSERT(junction.boundary_ids.size() == junction.polygon_vertices.size());
        for (size_t i = 0; i < junction.boundary_ids.size(); ++i) {
            junctions_per_boundary[junction.boundary_ids[i]][junction.polygon_vertices[i]] = true;
        }
    }
    
    for (Index i = 0; i < _boundaries.size(); ++i) {
        PF_VERBOSE_F("Boundary %d has %d junctions", i, junctions_per_boundary[i].size());
        for (size_t j = 0; j < junctions_per_boundary[i].size(); ++j) {
            if (junctions_per_boundary[i][j]) {
                PF_VERBOSE_F("Boundary %d vertex %d is unmergable", i, j);
            }
        }
    }

    for (size_t i = 0; i < _boundaries.size(); ++i) {
        int merge_from = 0;
        int merge_to = (int)_boundaries[i].spline().primitives.size();
        merge_consecutive_parallel_lines(_boundaries[i].spline(), true, merge_from, merge_to, junctions_per_boundary[i]);
    }

    // We cannot merge lines across junctions, but to enforce G^1 continuity we need to fix the tangents before
    for (ImageBoundary& boundary : _boundaries) {
        auto& primitives = boundary.spline().primitives;
        for (int prim_id_prev = 0; prim_id_prev < primitives.size(); ++prim_id_prev) {
            const int prim_id_next = Circular(primitives, prim_id_prev + 1);
            auto line_prev = dynamic_cast<GlobFitLineParametrization*>(primitives[prim_id_prev].curve.get());
            auto line_next = dynamic_cast<GlobFitLineParametrization*>(primitives[prim_id_next].curve.get());
           
            if (!line_prev || !line_next) {
                continue;
            }

            if (!can_merge_consecutive_lines(line_prev->get_curve()->dposdt(0.), line_next->get_curve()->dposdt(0.))) {
                continue;
            }

            // If this is triggered it means that the merging above failed in a case where it should have succeded
            const bool is_prev_junction = junctions_per_boundary[boundary.id()][primitives[prim_id_prev].corner];
            const bool is_next_junction = junctions_per_boundary[boundary.id()][primitives[prim_id_next].corner];
            PF_ASSERT(is_prev_junction || is_next_junction);

            // If the merge is at the entrance of the junction we fix the source tangent
            if (is_prev_junction) {
                _primitives_to_constraint.emplace_back(std::make_pair(&primitives[prim_id_prev], DofOptions::FIX_FRONT_TANGENT));
            }

            // If at the outside we fix the end tangent
            if (is_next_junction) {
                _primitives_to_constraint.emplace_back(std::make_pair(&primitives[prim_id_next], DofOptions::FIX_BACK_TANGENT));
            }
        }
    }
}

void print_primitives_cps(const vector<CurvePrimitive>& prims) {
    for (size_t i = 0; i < prims.size(); ++i) {
        cout << "Primitive " << i << endl;
        const auto& curve = prims[i].curve->get_curve().get();
        switch (curve->get_type()) {
        case GLOBFIT_CURVE_BEZIER: {
            std::cout << ((BezierCurve*)curve)->get_control_points() << endl;
            break;
        }
        case GLOBFIT_CURVE_LINE: {
            cout << ((GlobFitCurve_Line*)curve)->get_points() << endl;
            cout << "Params " << curve->get_params().transpose() << endl;
            cout << "Length " << ((GlobFitCurve_Line*)curve)->dposdt(0.).norm() << endl;
            cout << dynamic_pointer_cast<GlobFitLineParametrization>(prims[i].curve)->get_front_coordinate_system().matrix().col(0).transpose() << endl;
            cout << dynamic_pointer_cast<GlobFitLineParametrization>(prims[i].curve)->get_back_coordinate_system().matrix().col(0).transpose() << endl;
            break;
        }
        }
    }
}

/*
    Fits smooth splines in each polygonal region, coupling shared segments and introducing
    additional objectives at junctions to prevent gaps.
*/
void MultiPolygonTracer::_solve_curves() {
    //PF_DEV_S("Curves before optimization 12");
    //print_primitives_cps(_boundaries[12].spline().primitives);
    _curve_fitter->solve_with_stiffening();
    //PF_DEV_S("Curves after optimization 12");
    //print_primitives_cps(_boundaries[3].spline().primitives);
}

void MultiPolygonTracer::_post_process() {
    for (ImageBoundary& boundary : _boundaries) {
        PF_VERBOSE_F("Postprocessing boundary %d", boundary.id());

        PostProcessor pp(
            boundary.boundary_points(),
            boundary.midpoints(),
            boundary.edges(),
            boundary.polygon_vertices(),
            boundary.polygon_points(),
            boundary.raster_symmetries(),
            boundary.raster_symmetries_local(),
            boundary.regularity_graph(),
            boundary.tangents_fits(),
            boundary.spline(),
            _classifier_uri.c_str(),
            boundary.id()
        );

        //pp.collapse_short_segments();
        //pp.merge_curves();
    }
}

polyfit::BoundaryGraph::EdgeID MultiPolygonTracer::make_image_edge_id(
    const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p1
) {
    const int v0 = boundary_vertex_id_from_pos(p0);
    const int v1 = boundary_vertex_id_from_pos(p1);
    return BoundaryGraph::make_edge_id(min(v0, v1), max(v0, v1));
}

std::pair<int, int> MultiPolygonTracer::_get_boundary_subpath_for_segment(
    const int segment_id,
    const int boundary_id
) const {
    const auto& P = _boundaries[boundary_id].boundary_points();
    const auto& segment_pts = _segment_points[segment_id];
    
    // Obtain image space edge IDs from the segment points
    const int B_first = _I2B_map->to_boundary_vertex_from(segment_pts.col(0).cast<double>(), segment_pts.col(1).cast<double>(), boundary_id);
    const int B_last = _I2B_map->to_boundary_vertex_to(segment_pts.col(segment_pts.cols() - 2).cast<double>(), segment_pts.col(segment_pts.cols() - 1).cast<double>(), boundary_id);

    // Getting the vertices of the queried boundary
    PF_ASSERT(B_first != -1 && B_last != -1);
    int B_first_next = Circular(P, B_first + 1);
    if ((P.col(B_first) - P.col(B_first_next)).norm() < .5 + PF_EPS) {
        B_first_next = Circular(P, B_first_next + 1);
    }

    const int B_first_next_I = _I2B_map->image_vertex_id_from_point(P.col(B_first_next));
    const int S_first_next_I = _I2B_map->image_vertex_id_from_point(segment_pts.col(1).cast<double>());
    const int S_limit_direction = B_first_next_I == S_first_next_I ? +1 : -1;
    const int S_limit_src = S_limit_direction > 0 ? B_first : B_last;
    const int S_limit_dst = S_limit_direction > 0 ? B_last : B_first;

    PF_STATUS_F("The ordered limits for polygon %d subpath for segment %d are %d %d", boundary_id, segment_id, S_limit_src, S_limit_dst);
    return std::pair<int, int>(S_limit_src, S_limit_dst);
}

std::pair<int, int> MultiPolygonTracer::_get_polygon_subpath_for_segment(
    const int S_limit_src,
    const int S_limit_dst,
    const int boundary_id,
    Eigen::VectorXi* P // this should be computed from the returned value...
) const {
    // Reset output
    if (P) {
        P->resize(0);
    }

    const ImageBoundary& boundary = _boundaries[boundary_id];

    // Finding the first point the polygon subpath. They are not guaranteed to overlap
    int subP_src = -1;
    int subP_src_dist = INT_MAX;

    for (Index i = 0; i < boundary.polygon_vertices().size(); ++i) {
        PF_VERBOSE_F("Polygon vertex %d = %d", i, boundary.polygon_vertices()(i));
        if (!PathUtils::contains_closed(boundary.boundary_points(), S_limit_src, S_limit_dst, boundary.polygon_vertices()(i), true)) {
            continue;
        }

        const int dist = abs(CircularDist(boundary.boundary_points(), S_limit_src, boundary.polygon_vertices()(i)));
        if (dist < subP_src_dist) {
            subP_src = i;
            subP_src_dist = dist;
        }
    }

    // The subpath doesn't contain any vertex
    if (subP_src == -1) {
        return make_pair(-1, -1);
    }

    // Gathering all the vertices inside the subpath 
    int subP_dst = subP_src;
    size_t counter = 0;
    do {
        if (P) {
            MatrixUtils::append(*P, boundary.polygon_vertices()(subP_dst));
        }
        subP_dst = (int)Circular(boundary.polygon_vertices(), (Index)subP_dst + 1);
    } while (PathUtils::contains_closed(boundary.boundary_points(), S_limit_src, S_limit_dst, boundary.polygon_vertices()(subP_dst), true) &&
        ++counter < boundary.polygon_vertices().size());
    return make_pair(subP_src, (int)Circular(boundary.polygon_vertices(), (Index)subP_dst - 1));
}

void MultiPolygonTracer::_report_info() {
    for (size_t i = 0; i < _boundaries.size(); ++i) {
        cout << "-----------------------------------------------" << endl;
        cout << "Polygon " << i << " has " << _boundaries[i].polygon_vertices().size() << " points" << endl;
        cout << _boundaries[i].polygon_vertices().transpose() << endl;
    }

    for (size_t i = 0; i < _open_segments.size(); ++i) {
        const auto& segment = _open_segments[i];
        cout << "-----------------------------------------------" << endl;
        cout << "Segment " << segment.id << endl;
        cout << "boundary_src " << segment.boundary_src << endl;
        cout << "boundary_dst " << segment.boundary_dst << endl;
        cout << "B_bounds_src " << segment.B_bounds_src.first << " " << segment.B_bounds_src.second << endl;
        cout << "B_bounds_dst " << segment.B_bounds_dst.first << " " << segment.B_bounds_dst.second << endl;
    }

    for (size_t i = 0; i < _boundaries.size(); ++i) {
        const ImageBoundary& boundary = _boundaries[i];

        cout << "-----------------------------------------------" << endl;
        cout << "Spline " << i << endl;
        for (Index i = 0; i < boundary.spline().primitives.size(); ++i) {
            const CurvePrimitive& prim = boundary.spline().primitives[i];
            const auto edge_src = BoundaryGraph::unpack_edge_id(prim.fitting_info.edge_src);
            const auto edge_dst = BoundaryGraph::unpack_edge_id(prim.fitting_info.edge_dst);
            cout << "> Corner " << prim.corner << endl;
            cout << "Type " << curve_type_to_string(prim.curve->get_curve()->get_type()) << endl;
            cout << "Edge src " << edge_src(0) << "-" << edge_src(1) << endl;
            cout << "Edge dst " << edge_dst(0) << "-" << edge_dst(1) << endl;
        }
    }
}

void MultiPolygonTracer::_build_polygon_half_edges_and_fill_junctions_map() {
    unordered_map<BoundaryGraph::EdgeID, vec2i> he_map;

    _he_polygon.half_edges.clear();
    _junctions_to_half_edges.clear();
    for (JunctionInfo& junction : _junctions) {
        _junctions_to_half_edges[junction.vertex] = -1;
    }

    unsigned int he_id_offset = 0;
    for (size_t boundary_id = 0; boundary_id < _boundaries.size(); ++boundary_id) {
        const ImageBoundary& boundary = _boundaries[boundary_id];
        for (size_t i = 0; i < boundary.polygon_vertices().size(); ++i) {
            HalfEdge he;
            he.boundary = boundary_id;
            he.edge = BoundaryGraph::make_edge_id(i, Circular(boundary.polygon_vertices(), i + 1));
            he.id_prev = he_id_offset + Circular(boundary.polygon_vertices(), (Index)i - 1);
            he.id = he_id_offset + i;
            he.id_next = he_id_offset + Circular(boundary.polygon_vertices(), (Index)i + 1);

            const int v = boundary_vertex_id_from_pos(boundary.polygon_points().col(i));
            const int vn = boundary_vertex_id_from_pos(CircularAt(boundary.polygon_points(), i + 1));
            if (_junctions_to_half_edges.find(vn) != _junctions_to_half_edges.end()) {
                _junctions_to_half_edges[vn] = he.id; // overwriting, it doesnt' matter
            }

            const BoundaryGraph::EdgeID edge_id_unique = BoundaryGraph::make_edge_id(min(v, vn), max(v, vn));
            if (he_map.find(edge_id_unique) == he_map.end()) {
                he_map.insert(pair<BoundaryGraph::EdgeID, vec2i>(edge_id_unique, vec2i(he.id, -1)));
            } else {
                he_map[edge_id_unique](1) = he.id;
            }

            PF_VERBOSE_F("Boundary %d edge %d %d -> %d (%d %d)", boundary_id, i, Circular(boundary.polygon_vertices(), i + 1), edge_id_unique, v, vn);
            _he_polygon.half_edges.emplace_back(move(he));
        }

        he_id_offset += boundary.polygon_vertices().size();
    }

    for (const auto& kv : he_map) {
        PF_ASSERT(kv.second(0) != -1);
        if (kv.second(0) != -1 && kv.second(1) != -1) {
            _he_polygon.half_edges[kv.second(0)].id_twin = kv.second(1);
            _he_polygon.half_edges[kv.second(1)].id_twin = kv.second(0);
        }
    }
}

// Exists only between G0/G1, G0/G0 do not need to be processed.
struct JunctionSegment {
    int id;
    bool outgoing;
    int boundary_G0;
    BoundaryGraph::EdgeID edge_G0;
    int boundary_G1;
    BoundaryGraph::EdgeID edge_G1;
};

/*
    Do not run the optimization after calling this method.
*/
void MultiPolygonTracer::_merge_curves_at_junctions() {
    vector<int> he_at_junction; // Keeps track of all the half-edges which are incident on the junction vertex
    vector<CurvePrimitive> curves_to_insert; // Temporary container for the curves which needs to be replace in the G^0 junctions
    vector<CurvePrimitive> curves_split_prev, curves_split_next; // Temporary containers used in case the curves need to be split
    unordered_map<int, vector<CurvePrimitive>> he_to_cloned_curves; // Maps half-edge ids to the curves which will be substituted to the existing lines

    // Starting from any polygon edge incoming to the junction
    for (const pair<int, BoundaryGraph::EdgeID>& junction_to_he : _junctions_to_half_edges) {
        PF_ASSERT(junction_to_he.second >= 0 && junction_to_he.second < _he_polygon.half_edges.size());
        
        // Getting all the half_edges around the junction with twin ordering
        he_at_junction.clear();

        int he_next = junction_to_he.second;

        do {
            he_at_junction.emplace_back(he_next);
            he_next = _he_polygon.half_edges[_he_polygon.half_edges[he_next].id_next].id_twin;
        } while (he_next != junction_to_he.second);

        int G0_boundaries = 0;
        int he_G0 = -1;
        for (const int he : he_at_junction) {
            if (_boundaries[_he_polygon.boundary(he)].tangents_fits()[_he_polygon.vertex_to(he)] == TANGENT_FIT_CONSTANT) {
                ++G0_boundaries;
                he_G0 = he;
            }
        }

        const int valence = he_at_junction.size();
        PF_VERBOSE_F("Junction %d has valence %d with %d G0 boundaries", junction_to_he.first, valence, G0_boundaries);

        if (junction_to_he.first == 2306) {
            printf("moo");
        }

        // For all G^1 boundaries, get the list of curves which needs to be copied on the respective G0 boundary and
        // use the half-edge id to reference them
        he_to_cloned_curves.clear();
        for (const int he : he_at_junction) {
            const int boundary_G1 = _he_polygon.boundary(he);
            if (_boundaries[boundary_G1].tangents_fits()[_he_polygon.vertex_to(he)] == TANGENT_FIT_CONSTANT) {
                continue;
            }

            const int boundary_0 = _he_polygon.boundary(_he_polygon.twin(he));
            const int boundary_1 = _he_polygon.boundary(_he_polygon.twin(_he_polygon.next(he)));
            PF_ASSERT(boundary_G1 != boundary_0);
            PF_ASSERT(boundary_G1 != boundary_1);

            const int v_junction_0 = _he_polygon.vertex_from(_he_polygon.twin(he));
            const int v_junction_1 = _he_polygon.vertex_to(_he_polygon.twin(_he_polygon.next(he)));
            const TangentFitType tangent_fit_0 = _boundaries[boundary_0].tangents_fits()[v_junction_0];
            const TangentFitType tangent_fit_1 = _boundaries[boundary_1].tangents_fits()[v_junction_1];
            PF_VERBOSE_F("boundary_G1 %d (%s) has neighbors %d (%s) and %d (%s)", boundary_G1, _he_polygon.half_edges[he].str().c_str(),
                boundary_0, tangent_fit_type_to_string(tangent_fit_0), boundary_1, tangent_fit_type_to_string(tangent_fit_1));

            auto cloned_curves = _boundaries[boundary_G1].clone_curves_at_corner_reversed(_he_polygon.vertex_to(he));
            PF_ASSERT(cloned_curves.size() <= 2);

            const int he_G0_prev = tangent_fit_0 == TANGENT_FIT_CONSTANT ? _he_polygon.twin(he) : -1;
            const int he_G0_next = tangent_fit_1 == TANGENT_FIT_CONSTANT ? _he_polygon.twin(_he_polygon.next(he)) : -1;

            if (he_G0_prev != -1) {
                PF_ASSERT(he_to_cloned_curves.find(he_G0_prev) == he_to_cloned_curves.end());
            }

            if (he_G0_next != -1) {
                PF_ASSERT(he_to_cloned_curves.find(he_G0_next) == he_to_cloned_curves.end());
            }

            // If both boundaries are G^0 we need to split the curve
            vec2 corner_pos_new;
            if (he_G0_prev != -1 && he_G0_next != -1) {
                curves_split_prev.clear();
                curves_split_next.clear();

                const vec2 corner_pos = boundary_vertex_pos_from_id(junction_to_he.first);

                // todo: we should split at the intersection of the line and the curve (or any other generic sequence of primitives)
                // todo: we should be more careful in the presence of valence 4 junctions
                if (cloned_curves.size() == 2 && cloned_curves[0].curve->get_curve()->get_type() == GLOBFIT_CURVE_BEZIER) {
                    const double corner_pos_t = cloned_curves[0].curve->get_curve()->project(corner_pos);
                    corner_pos_new = cloned_curves[0].curve->get_curve()->pos(corner_pos_t);
                    pair<CurvePrimitive, CurvePrimitive> curves_split = cloned_curves[0].split(corner_pos_t);
                    curves_split_next.emplace_back(move(curves_split.first));
                    curves_split_prev.emplace_back(move(curves_split.second));
                    curves_split_prev.emplace_back(move(cloned_curves[1]));
                    PF_VERBOSE_F("Junction %d case 0", junction_to_he.first);
                }
                else if (cloned_curves.size() == 2 && cloned_curves[1].curve->get_curve()->get_type() == GLOBFIT_CURVE_BEZIER) {
                    const double corner_pos_t = cloned_curves[1].curve->get_curve()->project(corner_pos);
                    corner_pos_new = cloned_curves[1].curve->get_curve()->pos(corner_pos_t);
                    pair<CurvePrimitive, CurvePrimitive>  curves_split = cloned_curves[1].split(corner_pos_t);
                    curves_split_next.emplace_back(move(cloned_curves[0]));
                    curves_split_next.emplace_back(move(curves_split.first));
                    curves_split_prev.emplace_back(move(curves_split.second));
                    PF_VERBOSE_F("Junction %d case 1", junction_to_he.first);
                }
                else {
                    PF_ASSERT(cloned_curves.size() == 1);
                    PF_ASSERT(cloned_curves[0].curve->get_curve()->get_type() == GLOBFIT_CURVE_BEZIER);
                    const double corner_pos_t = cloned_curves[0].curve->get_curve()->project(corner_pos);
                    corner_pos_new = cloned_curves[0].curve->get_curve()->pos(corner_pos_t);
                    pair<CurvePrimitive, CurvePrimitive> curves_split = cloned_curves[0].split(corner_pos_t);
                    curves_split_next.emplace_back(move(curves_split.first));
                    curves_split_prev.emplace_back(move(curves_split.second));
                    PF_VERBOSE_F("Junction %d case 2", junction_to_he.first);
                }

                PF_ASSERT(!curves_split_prev.empty() && !curves_split_next.empty());
                he_to_cloned_curves[he_G0_prev] = move(curves_split_prev);
                he_to_cloned_curves[he_G0_next] = move(curves_split_next);
            }
            // Otherwise no need to copy
            else if (he_G0_prev != -1) {
                PF_VERBOSE_F("Storing %d curves in half-edge %s", cloned_curves.size(), _he_polygon.half_edges[he_G0_prev].str().c_str());
                he_to_cloned_curves[he_G0_prev] = move(cloned_curves);
            }
            else if (he_G0_next != -1) {
                PF_VERBOSE_F("Storing %d curves in half-edge %s", cloned_curves.size(), _he_polygon.half_edges[he_G0_next].str().c_str());
                he_to_cloned_curves[he_G0_next] = move(cloned_curves);
            }
            // All smooth, nothing to do 
            else {
                continue;
            }
        }

        // Reconstructing the junction curves from the G^0 boundaries
        for (const int he_prev : he_at_junction) {
            const int boundary_G0 = _he_polygon.boundary(he_prev);
            if (_boundaries[boundary_G0].tangents_fits()[_he_polygon.vertex_to(he_prev)] != TANGENT_FIT_CONSTANT) {
                continue;
            }

            const int he_next = _he_polygon.next(he_prev);

            const bool replace_prev = he_to_cloned_curves.find(he_prev) != he_to_cloned_curves.end();
            const bool replace_next = he_to_cloned_curves.find(he_next) != he_to_cloned_curves.end();

            // Sorrounded by G0 boundaries
            if (!replace_prev && !replace_next) {
                continue;
            }

            // Removing the curves from the G^0 boundary making sure that 
            vector<size_t> prims_delete;
            prims_delete.reserve(2);
            int insert_index = -1;
            if (replace_prev) {
                insert_index = _boundaries[boundary_G0].find_line_before_polygon_vertex(_he_polygon.vertex_to(he_prev));
                prims_delete.emplace_back((size_t)insert_index);
            }
            if (replace_next) {
                prims_delete.emplace_back((size_t)_boundaries[boundary_G0].find_line_after_polygon_vertex(_he_polygon.vertex_to(he_prev)));

                if (insert_index == -1) {
                    insert_index = (int)prims_delete.back();
                }
            }
            EraseOrdered(_boundaries[boundary_G0].spline().primitives, prims_delete);

            // Putting all the curves to be inserted into a single array
            curves_to_insert.clear();
            if (replace_prev) {
                curves_to_insert.insert(curves_to_insert.end(), make_move_iterator(he_to_cloned_curves[he_prev].begin()), make_move_iterator(he_to_cloned_curves[he_prev].end()));
            }

            if (replace_next) {
                curves_to_insert.insert(curves_to_insert.end(), make_move_iterator(he_to_cloned_curves[he_next].begin()), make_move_iterator(he_to_cloned_curves[he_next].end()));
            }

            for (auto& c : curves_to_insert) {
                c.corner = _he_polygon.vertex_to(he_prev);
                PF_VERBOSE_F("Inserting curve %d", c.corner);
            }

            _boundaries[boundary_G0].spline().primitives.insert(_boundaries[boundary_G0].spline().primitives.begin() + insert_index,
                make_move_iterator(curves_to_insert.begin()), make_move_iterator(curves_to_insert.end()));

            PF_DEV_F("Junction %d relocating C^0 transition boundary %d vertex %d", junction_to_he.first, boundary_G0, _he_polygon.vertex_to(he_prev));

            // If we replaced only a single part of the curve, we need to relocate the next endpoint (line) to the previous curve, 
            // essentially making the transition really C0
            const int v_C0 = _he_polygon.vertex_to(he_prev);
            if (replace_prev && !replace_next) {
                const int curve_id_after = _boundaries[boundary_G0].find_line_after_polygon_vertex(v_C0);
                auto line_param_after_junction = dynamic_pointer_cast<GlobFitLineParametrization>(_boundaries[boundary_G0].spline().primitives[curve_id_after].curve);
                const auto bezier_param_before_junction = dynamic_pointer_cast<GlobFitBezierAngleBasedParametrization>(CircularAt(_boundaries[boundary_G0].spline().primitives, curve_id_after - 1).curve);
                PF_ASSERT(line_param_after_junction);
                PF_ASSERT(bezier_param_before_junction);
                auto points = dynamic_pointer_cast<GlobFitCurve_Line>(line_param_after_junction->get_curve())->get_points();
                points.col(0) = dynamic_pointer_cast<BezierCurve>(bezier_param_before_junction->get_curve())->get_control_points().col(3);
                dynamic_pointer_cast<GlobFitCurve_Line>(line_param_after_junction->get_curve())->set_points(points);
            }

            if (replace_next && !replace_prev) {
                const int curve_id_before = _boundaries[boundary_G0].find_line_before_polygon_vertex(v_C0);
                auto line_param_before_junction = dynamic_pointer_cast<GlobFitLineParametrization>(_boundaries[boundary_G0].spline().primitives[curve_id_before].curve);
                const auto bezier_param_after_junction = dynamic_pointer_cast<GlobFitBezierAngleBasedParametrization>(CircularAt(_boundaries[boundary_G0].spline().primitives, curve_id_before + 1).curve);
                PF_ASSERT(line_param_before_junction);
                PF_ASSERT(bezier_param_after_junction);
                auto points = dynamic_pointer_cast<GlobFitCurve_Line>(line_param_before_junction->get_curve())->get_points();
                points.col(1) = dynamic_pointer_cast<BezierCurve>(bezier_param_after_junction->get_curve())->get_control_points().col(0);
                dynamic_pointer_cast<GlobFitCurve_Line>(line_param_before_junction->get_curve())->set_points(points);
            }
        }
    }
}

int MultiPolygonTracer::boundary_vertex_id_from_pos(const Eigen::Vector2d pos) const {
    const int cols = _conn_raster.pixel_vertex_dims().x() * 2;
    const Eigen::Vector2i pos_i = (pos * 2).cast<int>();
    return pos_i.x() + pos_i.y() * cols;
}

Eigen::Vector2d MultiPolygonTracer::boundary_vertex_pos_from_id(const int id) const {
    const int cols = _conn_raster.pixel_vertex_dims().x() * 2;
    const Eigen::Vector2i pos_i(id % cols, id / cols);
    return pos_i.cast<double>() * .5;
}

/*
    Sorts the image boundaries by area and removes the outmost one
*/
std::vector<int> MultiPolygonTracer::get_boundaries_ordered_by_area() const {
    vector<int> regions_sorted_by_area = _conn_raster.regions_sorted_by_area();
    vector<int> boundaries_sorted_by_area;
    boundaries_sorted_by_area.reserve(_boundaries.size());

    vector<int> holes;
    int outer_boundary;
    for (const int region : regions_sorted_by_area) {
        _conn_raster.region_to_polygon(region, outer_boundary, holes);
        const int outer_boundary_id = _polygon_to_boundary[outer_boundary];

        //if (outer_boundary_id != _outer_boundary) {
            boundaries_sorted_by_area.emplace_back(outer_boundary_id);
        //}
    }

    return boundaries_sorted_by_area;   
}

const std::vector<ImageBoundary>& MultiPolygonTracer::boundaries() const {
    return _boundaries;
}

const std::vector<JunctionInfo>& MultiPolygonTracer::junctions() const {
    return _junctions;
}

const std::vector<Eigen::Matrix2Xi> MultiPolygonTracer::segments()const {
    return _segment_points;
}

const std::vector<HalfEdge>& MultiPolygonTracer::half_edges()const {
    return _he_polygon.half_edges;
}

NAMESPACE_END(polyvec)