 // Polyvec
#include <polyvec/polygon-tracer/image_boundary.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/polygon-tracer/regularized.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/regularity/continuations.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/core/constants.hpp>

// libc++
#include <utility>
#include <algorithm>

using namespace polyfit;
using namespace Eigen;
using namespace std;

NAMESPACE_BEGIN(polyvec)

ImageBoundary::ImageBoundary(
    const Eigen::Matrix2Xi& points, 
    const int polygon_id,
    const Vector4d& color) :
    _points(points), 
    _points_as_double(points.cast<double>()), 
    _curves_dirty(true), 
    _id(polygon_id),
    _color(color) {
    PF_ASSERT(_id >= 0);

    BoundaryGraph::create_fitting_path(_points_as_double, _polygon.B, _polygon.M, true);
    PathUtils::compute_convexities(_polygon.B, _boundary_convexities);

    _polygon.midpoints.clear();
    _polygon.midpoints.resize(_polygon.B.cols(), false);
    for (size_t i = 0; i < _polygon.M.size(); ++i) {
        const Index v_p = Circular(_polygon.B, _polygon.M[i] - 1);
        if ((_polygon.B.col(v_p) - _polygon.B.col(_polygon.M[i])).squaredNorm() < .5 + PF_EPS) {
            _polygon.midpoints[_polygon.M[i]] = true;
        }
    }
}

bool ImageBoundary::is_flat(const int v) const {
    return _boundary_convexities[v] == 0;
}

void ImageBoundary::initialize_graph_with_flat_vertices(const Eigen::VectorXi& B_flat) {
    for (Index i = 0; i < B_flat.size(); ++i) {
        PF_ASSERT(is_flat(B_flat(i)));
    }

    _B_flat.clear();
    _B_flat.resize(_polygon.B.cols(), false);
    for (Index i = 0; i < B_flat.size(); ++i) {
        _B_flat[B_flat(i)] = true;
    }

    BoundaryGraph::connect_valid_edges(_polygon.B, _polygon.M, _polygon.E, true, &_accuracy_map, _B_flat, true);
    BoundaryGraph::connect_valid_edges(_polygon.B, _polygon.M, _E_original, true, &_accuracy_map, _B_flat, false);
}

void ImageBoundary::trace_polygon(
    const std::unordered_set<int>& kopf_junctions
) {
    _initialize_graph_if_necessary();
    _polygon.raster_symmetries = Symmetry::find_longest(_polygon.B, true);
    _polygon.raster_symmetries_local = Symmetry::find_shortest(_polygon.B, true);
    PolygonTracer::regularized(
        _shortest_path_state, _polygon.B, _polygon.M, _polygon.E, _polygon.RE, _polygon.P, 
        _polygon.raster_symmetries,  _polygon.raster_symmetries_local, true,
        nullptr, std::vector<bool>(), _B_flat, kopf_junctions
    );
    _update_polygon();
    PF_VERBOSE_F("Continuations %d", _polygon.RE.continuations().size());
}

void ImageBoundary::trace_polygon_simple() {
    _initialize_graph_if_necessary();

    if (!PolygonTracer::minimum(_shortest_path_state, _polygon.B, _polygon.M, _polygon.E, _polygon.P, true)) {
        PF_DEV_S("Failed to find polygon");
    }

    _update_polygon();
}

void ImageBoundary::reset_graph() {
    _polygon.E = _E_original;
}

void ImageBoundary::add_edge(const BoundaryGraph::Edge& e) {
    PF_ASSERT(e.v0 >= 0 && e.v0 < _polygon.B.cols());
    PF_ASSERT(e.v1 >= 0 && e.v1 < _polygon.B.cols());
    _polygon.E.emplace_back(e);
}

void ImageBoundary::remove_edges_crossing_vertices(const std::vector<int> V) {
    _E_delete.clear();
    for (size_t i = 0; i < _polygon.E.size(); ++i) {
        for (const int v : V) {
            if (PathUtils::contains_open(_polygon.B, _polygon.E[i].v0, _polygon.E[i].v1, v, true)) {
                _E_delete.emplace_back(i);
                break;
            }
        }
    }

    EraseOrdered(_polygon.E, _E_delete);
}

double ImageBoundary::bounding_box_area()const {
    return (_polygon.B.row(0).maxCoeff() - _polygon.B.row(0).minCoeff()) *
        (_polygon.B.row(1).maxCoeff() - _polygon.B.row(1).minCoeff());
}

void ImageBoundary::force_graph_subpath ( const Eigen::VectorXi& subpath ) {
    BoundaryGraph::force_graph_subpath(_polygon.B, _polygon.E, subpath, _accuracy_map);
}

Eigen::VectorXi ImageBoundary::get_polygon_subpath(const int v_src, const int v_dst) const {
    std::pair<int, int> v_bounds = get_polygon_subpath_bounds(v_src, v_dst);
    
    if (v_bounds.first == -1 || v_bounds.second == -1) {
        return VectorXi();
    }

    VectorXi subP;
    const int n_vertices = CircularDist(_polygon.P, v_bounds.first, v_bounds.second);
    for (int i = 0; i <= n_vertices; ++i) {
        const int v = Circular(_polygon.P, v_bounds.first + i);
        MatrixUtils::append(subP, v);
    }

    return subP;
}

std::pair<int, int> ImageBoundary::get_polygon_subpath_bounds(const int v_src, const int v_dst) const {
    // Finding the first point the polygon subpath. They are not guaranteed to overlap
    int subP_src = -1;
    int subP_src_dist = INT_MAX;

    for (Index i = 0; i < _polygon.P.size(); ++i) {
        if (!PathUtils::contains_closed(_polygon.B, v_src, v_dst, _polygon.P(i), true)) {
            continue;
        }

        const int dist = abs(CircularDist(_polygon.B, v_src, _polygon.P(i)));
        if (dist < subP_src_dist) {
            subP_src = i;
            subP_src_dist = dist;
        }
    }

    // The subpath doesn't contain any vertex
    if (subP_src == -1) {
        return std::pair<int, int>(-1, -1);
    }

    // Gathering all the vertices inside the subpath 
    int subP_dst = subP_src;
    size_t counter = 0;
    do {
        subP_dst = (int)Circular(_polygon.P, (Index)subP_dst + 1);
    } while (PathUtils::contains_closed(_polygon.B, v_src, v_dst, _polygon.P(subP_dst), true) &&
        ++counter < _polygon.P.size());

    return std::pair<int, int>(subP_src, (int)Circular(_polygon.P, subP_dst - 1));
}

double ImageBoundary::calculate_subpath_energy(const int v_first, const int v_last) const {
    PF_ASSERT(v_first >= 0 && v_first < _polygon.P.size());
    PF_ASSERT(v_last >= 0 && v_last < _polygon.P.size());

    if (v_first == -1 || v_last == -1) {
        return 0.;
    }

    double E = 0.;

    // If a single point is in the subpath, just return the angle
    if (v_first == v_last) {
        const double angle = AngleUtils::spanned_shortest(CircularAt(_polygon.PP, v_first - 1), _polygon.PP.col(v_first), CircularAt(_polygon.PP, v_first + 1));
        E += ErrorMetrics::smoothness(angle);
        return E;
    }

    int v = v_first;
    while (v != v_last) {
        const vec2 p0 = CircularAt(_polygon.PP, (Index)v - 1);
        const vec2 p1 = CircularAt(_polygon.PP, (Index)v);
        const vec2 p2 = CircularAt(_polygon.PP, (Index)v + 1);
        const vec2 p3 = CircularAt(_polygon.PP, (Index)v + 2);

        // This should be an assert, but the accuracy map is not tied to the boundary graph and might be late in being updated...
        const BoundaryGraph::EdgeID edge_id = BoundaryGraph::make_edge_id(_polygon.P(v), CircularAt(_polygon.P, (Index)v + 1));
        vec2 d_error;
        if (_accuracy_map.find(edge_id) == _accuracy_map.end()) {
            d_error = GeomRaster::distance_bounds_from_points_with_slack(_polygon.B, _polygon.PP.col(v), CircularAt(_polygon.PP, v + 1), _polygon.P(v), CircularAt(_polygon.P, v + 1));
        } else {
            d_error = _accuracy_map.at(edge_id);
        }

        E += BoundaryGraph::approximate_edge_cost(p0, p1, p2, p3, d_error.minCoeff(), d_error.maxCoeff());
        v = (int)Circular(_polygon.PP, (Index)v + 1);
    }

    return E;
}

int ImageBoundary::find_polygon_vertex_for_boundary_vertex(const int v)const {
    PF_ASSERT(v >= 0 && v < _polygon.B.cols());
    for (Index i = 0; i < _polygon.P.size(); ++i) {
        if (_polygon.P(i) == v) {
            return i;
        }
    }

    return -1;
}

int ImageBoundary::find_line_before_polygon_vertex(const int polygon_v) const {
    PF_ASSERT(_tangents[polygon_v] == TANGENT_FIT_CONSTANT);
    const BoundaryGraph::EdgeID edge_before_corner = BoundaryGraph::make_edge_id(Circular(_polygon.P, polygon_v - 1), polygon_v);
    for (size_t i = 0; i < _spline.primitives.size(); ++i) {
        if (_spline.primitives[i].corner == polygon_v &&
            _spline.primitives[i].fitting_info.edge_src == _spline.primitives[i].fitting_info.edge_dst &&
            _spline.primitives[i].fitting_info.edge_src == edge_before_corner) {
            return (int)i;
        }
    }

    PF_ABORT;
    return -1;
}

int ImageBoundary::find_line_after_polygon_vertex(const int polygon_v) const {
    PF_ASSERT(_tangents[polygon_v] == TANGENT_FIT_CONSTANT);
    const BoundaryGraph::EdgeID edge_after_corner = BoundaryGraph::make_edge_id(polygon_v, Circular(_polygon.P, polygon_v + 1));
    for (size_t i = 0; i < _spline.primitives.size(); ++i) {
        if (_spline.primitives[i].corner == polygon_v &&
            _spline.primitives[i].fitting_info.edge_src == _spline.primitives[i].fitting_info.edge_dst &&
            _spline.primitives[i].fitting_info.edge_dst == edge_after_corner) {
            return (int)i;
        }
    }

    PF_ABORT;
    return -1;
}

void ImageBoundary::classify_corners(FitClassifier& classifier) {
    _initialize_curve_fitter();
    _spline_fitter->classify_evolutionary_simple(classifier, _tangents, 
        [&](const CurveSequenceFitter::EvolutionaryFittingState & state, const CurveSequenceFitter::FittingAttempt * fits) {
            _fitting_attempts.emplace_back();
            _fitting_attempts.back().attempt[0] = fits[0];
            _fitting_attempts.back().attempt[1] = fits[1];
            _fitting_attempts.back().attempt[2] = fits[2];
            _fitting_attempts.back().state = state;
    });
}

// Assumes a tangent fit classification for each polygon corner has been obtained
void ImageBoundary::fit_curves_initial_guess() {
    PF_ASSERT(_spline_fitter.get());
    _spline = _spline_fitter->fit_initial_guess(_tangents);
}

// The function will fail if you attempt to flip the same edge two times.
void ImageBoundary::flip_coordinate_system(const BoundaryGraph::EdgeID edge) {
    PF_ASSERT(!_edge_coordinate_systems.empty()); // Assuming we have fit a polygon
    const vec2i v_edge = BoundaryGraph::unpack_edge_id(edge);
    PF_ASSERT(v_edge(0) >= 0 && v_edge(0) < _polygon.P.size());
    PF_ASSERT(v_edge(1) >= 0 && v_edge(1) < _polygon.P.size());
    PF_ASSERT(_edge_coordinate_systems[v_edge(0)] == false);
    _edge_coordinate_systems[v_edge(0)] = true;
}

std::pair<int, int> ImageBoundary::get_primitives_interval_open(const int v_first, const int v_last, const int direction) const {
    const auto interval_closed = get_primitives_interval_closed(v_first, v_last, direction);
    auto interval_open = interval_closed;
    while (_spline.primitives[interval_open.first].corner == v_first) {
        interval_open.first = Circular(_spline.primitives, interval_open.first + direction);
    }

    while (_spline.primitives[interval_open.second].corner == v_last) {
        interval_open.second = Circular(_spline.primitives, interval_open.second - direction);
    }

    if (interval_closed.first == interval_open.second && interval_closed.second == interval_open.first) {
        return pair<int, int>(interval_closed.first, interval_closed.first);
    }

    return interval_open;
}

std::pair<int, int> ImageBoundary::get_primitives_interval_closed(const int v_first, const int v_last, const int direction) const {
    PF_ASSERT(v_first >= 0 && v_first < _polygon.P.size());
    PF_ASSERT(v_last >= 0 && v_last < _polygon.P.size());

    const int v_first_ordered = direction > 0 ? v_first : v_last;
    const int v_last_ordered = direction > 0 ? v_last : v_first;

    // Find the last primitive fitting the corner
    int curve_begin = distance(_spline.primitives.begin(), find_if(_spline.primitives.begin(), _spline.primitives.end(),
        [v_first_ordered](const CurvePrimitive& prim) -> bool {
            return prim.corner == v_first_ordered;
    }));
    
    int curve_end = distance(_spline.primitives.begin(), find_if(_spline.primitives.begin(), _spline.primitives.end(),
        [v_last_ordered](const CurvePrimitive& prim) -> bool {
            return prim.corner == v_last_ordered;
        }));

    while (CircularAt(_spline.primitives, curve_begin + 1).corner == v_first_ordered) {
        ++curve_begin;
    }

    while (CircularAt(_spline.primitives, curve_end - 1).corner == v_last_ordered) {
        --curve_end;
    }

    if (direction > 0) {
        return pair<int, int>(curve_begin, curve_end);
    } else {
        return pair<int, int>(curve_end, curve_begin);
    }
}

std::vector<CurvePrimitive> ImageBoundary::clone_curves_at_corner_reversed(const int v_P) const {
    int v_first= distance(_spline.primitives.begin(), find_if(_spline.primitives.begin(), _spline.primitives.end(),
        [v_P](const CurvePrimitive& prim) { return prim.corner == v_P; }));
    while (CircularAt(_spline.primitives, v_first - 1).corner == v_P) {
        v_first = Circular(_spline.primitives, v_first - 1);
    }

    std::vector<CurvePrimitive> ret;
    while (_spline.primitives[v_first].corner == v_P) {
        CurvePrimitive prim;
        prim.curve = shared_ptr<GlobFitCurveParametrization>(_spline.primitives[v_first].curve->clone());
        prim.curve->reverse();
        ret.emplace_back(move(prim));
        PF_VERBOSE_F("Boundary %d clone curve %d at corner %d", _id, v_first, v_P);
        v_first = Circular(_spline.primitives, v_first + 1);
    }
    reverse(ret.begin(), ret.end());
    return ret;
}

double ImageBoundary::angle_at_vertex(const int v) {
    PF_ASSERT(_polygon.P.size() > 2);
    return AngleUtils::spanned_shortest(CircularAt(_polygon.PP, v - 1), _polygon.PP.col(v), CircularAt(_polygon.PP, v + 1));
}

// Assumes a polygonal fit has been obtained
void ImageBoundary::_initialize_curve_fitter() {
    PF_ASSERT(_is_polygon_data_consistent_and_exists());

    if (_curves_dirty) {
        PF_VERBOSE_F("Regularities %d continuations %d", _polygon.RE.size(), _polygon.RE.continuations().size());
        _spline_fitter.reset(new CurveSequenceFitter(
            _polygon.B, 
            _polygon.PP, 
            _polygon.PV, 
            _polygon.RE, 
            _id, 
            _edge_coordinate_systems));
        _curves_dirty = false;
    }
}

void ImageBoundary::_update_polygon() {
    BoundaryGraph::trace_to_points(_polygon.B, _polygon.P, _polygon.PP);

    _polygon.PV.resize(_polygon.P.size());
    for (int v = 0; v < _polygon.P.size(); ++v) {
        _polygon.PV[v] = _polygon.P(v);
    }

    _update_convexities();
    _curves_dirty = true;
    _edge_coordinate_systems.clear();
    _edge_coordinate_systems.resize(_polygon.P.size(), false);

#if POLYVEC_REGULARIZE
    _polygon.RE.reset();	
	_polygon.RE.update(_polygon.B, _polygon.P, _polygon.raster_symmetries, _polygon.raster_symmetries_local, true);
	_polygon.RE.add_continuations(_polygon.B, _polygon.P, true);    
#endif
}

// Since the convexities are computed on the polygon corners, it assumes that a polygon exists
void ImageBoundary::_update_convexities() {
    PF_ASSERT(_is_polygon_data_consistent_and_exists());

    PathUtils::compute_convexities(_polygon.PP, _convexities_inside);
    _convexities_outside = _convexities_inside;
    for (size_t i = 0; i < _convexities_outside.size(); ++i) {
        _convexities_outside[i] = -_convexities_outside[i];
    }
}

bool ImageBoundary::_is_polygon_data_consistent_and_exists() {
    return _polygon.P.size() > 0 &&
        _polygon.PP.cols() == _polygon.P.size() &&
        _polygon.PV.size    () == _polygon.P.size();
}

void ImageBoundary::_initialize_graph_if_necessary() {
    if (_polygon.E.empty()) {
        BoundaryGraph::connect_valid_edges(_polygon.B, _polygon.M, _polygon.E, true, &_accuracy_map);
        _E_original = _polygon.E;
        _B_flat.clear();
        _B_flat.resize(_polygon.B.cols(), false);
    }
}

NAMESPACE_END(polyvec)