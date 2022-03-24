// Polyvec
#include <polyvec/polygon-tracer/postprocessor.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/polygon-tracer/regularized.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/options.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/curve-tracer/bezier_merging.hpp>
#include <polyvec/regularity/find_symmetric_primitives.hpp>
#include <polyvec/regularity/continuations.hpp>

// libc++
#include <unordered_set>

#define FANCY_SYMMETRIC_LOGIC 1

using namespace Eigen;
using namespace polyvec;
using namespace polyfit;
using namespace std;

NAMESPACE_BEGIN(polyvec)

PostProcessor::PostProcessor(
    const mat2x&                      B,
    const std::vector<Vertex>&        M,
    std::vector<BoundaryGraph::Edge>& E,
    vecXi&                            P,
    mat2x&                            PP,
	std::vector<polyfit::Symmetry::SubPath>&   raster_symmetries,
    std::vector<polyfit::Symmetry::SubPath>&   raster_symmetries_local,
    polyfit::Regularity::RegularityInformation&    RE,
    std::vector<TangentFitType>&      tangent_fits,
    CurvePrimitiveSequence&           curves,
    const char*                       classifier_uri,
    const int polygon_id
) :
    circular(true),
    B(B), M(M), E(E), P(P), PP(PP), raster_symmetries(raster_symmetries), raster_symmetries_local(raster_symmetries_local),
    regularity(RE), tangent_fits(tangent_fits), curves(curves),
    classifier_uri(classifier_uri), polygon_id(polygon_id) {
    BoundaryGraph::connect_valid_edges(B, M, _E_original, true, &_accuracy_map);
}

/*
    Finds all the polygon edges which are within 2 pixels and for which one of the two corners
    has been identified as C^0 and attempts to remove them.
    The edge connecting the two pixels is removed and the shortest path is rerun, making sure that 
    none of the regularities is broken.
*/
struct Collapsable {
    Collapsable(const int v_src, const int v_dst, const int v_force) :
        v_src(v_src), v_dst(v_dst), v_force(v_force) { }

    int v_src;
    int v_dst;
    int v_force;
};

void PostProcessor::retrace_with_C0_corners() {
    std::vector<bool> corners_C0(B.cols(), false);
    std::vector<size_t> corners;

    for (Index i = 0; i < P.size(); ++i) {
        if (tangent_fits[i] == TANGENT_FIT_CONSTANT) {
            corners_C0[P(i)] = true;
            corners.emplace_back(P(i));
        }
    }

    // Actually forcing the path to go through the corners
    _E = _E_original;
    _E_delete.clear();

    for (size_t i = 0; i < _E.size(); ++i) {
        const vec2 p = B.col(_E[i].v0);
        const vec2 pn = B.col(_E[i].v1);
        
        for (size_t j = 0; j < corners.size(); ++j) {
            if (PathUtils::contains_open(B, _E[i].v0, _E[i].v1, corners[j], circular)) {
                _E_delete.emplace_back(i);
                break;
            }
        }
    }

    EraseOrdered(_E, _E_delete);
    PolygonTracer::regularized(_G_state, B, M, _E, regularity, P, raster_symmetries, raster_symmetries_local, circular, nullptr, corners_C0);
    BoundaryGraph::trace_to_points(B, P, PP);
}

void PostProcessor::collapse_short_segments() {
	FitClassifierRandomForest classifier;
	classifier.load_from_file(classifier_uri.c_str());

    collapse_short_segments_iteration();
    trace_curves_from_polygon(classifier);
}

int PostProcessor::collapse_short_segments_iteration() {
    double max_edge_length_sq = Options::get()->post_processor_max_short_edge_length;
    max_edge_length_sq *= max_edge_length_sq;

    // Let's gather all the edges which need to be collapsed
    std::vector<Collapsable> to_collapse;
    
    std::vector<int> C;
    PathUtils::compute_convexities(B, C);

    for (Index i = 0; i < P.size(); ++i) {
        const int v = P(i);
        const int vn = CircularAt(P, i + 1);

        const vec2 p = B.col(v);
        const vec2 pn = B.col(vn);
        const double edge_length_sq = (p - pn).squaredNorm();

        const bool is_inflection = AngleUtils::have_opposite_convexity(
            B.col(CircularAt(P, i - 1)), p, pn, B.col(CircularAt(P, i + 2))
        );

        if (is_inflection) {
            continue;
        }
        
        if (edge_length_sq > max_edge_length_sq + PF_EPS) {
            continue;
        }

        PF_ASSERT(tangent_fits.size() == P.size());
        const TangentFitType fit_src = tangent_fits[i];
        const TangentFitType fit_dst = CircularAt(tangent_fits, i + 1);

        if (fit_src != TANGENT_FIT_CONSTANT && fit_dst != TANGENT_FIT_CONSTANT) {
            continue;
        }

        //
        const int is_diagonal = CircularDist(B, v, vn) == 3;
        if (is_diagonal && fit_src == TANGENT_FIT_CONSTANT && fit_dst == TANGENT_FIT_CONSTANT) {
            const int cn = CircularAt(C, vn - 1);
            const int cp = CircularAt(C, v + 1);
            if (cn == -1) {
                to_collapse.emplace_back(v, vn, Circular(C, vn - 1));
            }
            else if (cp == -1) {
                to_collapse.emplace_back(v, vn, Circular(C, v + 1));
            }
        }
        // Figuring out where to place the corner, if only one of the corners is C^0 then that is where it will 
        // be placed. If both are C^0 and it's symmetric then into a midpoint if it exists, otherwise we need 
        // attempt both corners.
        else if (fit_src == TANGENT_FIT_CONSTANT && fit_dst == TANGENT_FIT_CONSTANT) {
            const bool midpoint_exists = CircularDist(P, v, vn) == 2 && GeomRaster::are_points_axis_aligned(p, pn);

            if (midpoint_exists) {
                    to_collapse.emplace_back(v, vn, (int)Circular(B, v + 1));
            }
            else {
                const int dist_prev = CircularDist(B, PathUtils::next_transition_vertex(C, v, -1, circular), v);
                const int dist_next = CircularDist(B, vn, PathUtils::next_transition_vertex(C, vn, +1, circular));

                if (dist_prev > dist_next) {
                    to_collapse.emplace_back(v, vn, v);
                }
                else if (dist_prev < dist_next) {
                    to_collapse.emplace_back(v, vn, vn);
                }
            }
        }
        else if (fit_src == TANGENT_FIT_CONSTANT) {
            to_collapse.emplace_back(v, vn, v);
        }
        else if (fit_dst == TANGENT_FIT_CONSTANT) {
            to_collapse.emplace_back(v, vn, vn);
        }
        else {
            continue;
        }
    }

    _set_immovable_corners();
    //return 0;

    int collapsed = 0;
    for (size_t i = 0; i < to_collapse.size(); ++i) {
        // Are we collapsing onto an immovable corners
        bool collapsing_regularity_safely =
            (!_immovable[to_collapse[i].v_dst] && _immovable[to_collapse[i].v_src] && to_collapse[i].v_src == to_collapse[i].v_force) ||
            (!_immovable[to_collapse[i].v_src] && _immovable[to_collapse[i].v_dst] && to_collapse[i].v_dst == to_collapse[i].v_force) || 
            (!_immovable[to_collapse[i].v_src] && !_immovable[to_collapse[i].v_dst]);

        // but not both should be immovable
        if (_immovable[to_collapse[i].v_src] && _immovable[to_collapse[i].v_dst]) {
            collapsing_regularity_safely = false;
        }

        // Are we collapsing a regularity or collapsing it into the very same vertex
        PF_DEV_F("Attempting to collapse (%d %d) -> %d is_safe %d", to_collapse[i].v_src, to_collapse[i].v_dst, to_collapse[i].v_force, collapsing_regularity_safely);
        if (collapsing_regularity_safely) {
            if (_relocate_short_edge(to_collapse[i].v_src, to_collapse[i].v_dst, to_collapse[i].v_force)) {
                ++collapsed;
            }
        }

        if (collapsed == 2) {
            break;
        }
    }

    return collapsed;
}

/*
    Run steps which have not been run during post-processing (i.e. curve fitting)
*/
void PostProcessor::trace_curves_from_polygon(FitClassifier& classifier) {
    // come on...
    std::vector<Vertex> PV(P.size());
    for (size_t i = 0; i < PV.size(); ++i) {
        PV[i] = P(i);
    }

    CurveSequenceFitter fitter(B, PP, PV, regularity, polygon_id);
    curves = fitter.fit_evolutionary_simple(classifier, tangent_fits, nullptr);
}

namespace {
    bool is_degenerate(const mat2x& B, const vecXi& P, const int v0, const int v1, const bool circular) {
        if (!circular && (v0 == 0 || v1 == P.size() - 1)) {
            return false;
        }

        const vec2 p0p = B.col(CircularAt(P, v0 - 1));
        const vec2 p0 = B.col(CircularAt(P, v0));
        const vec2 p1 = B.col(CircularAt(P, v1));
        const vec2 p1n = B.col(CircularAt(P, v1 + 1));

        if (!AngleUtils::have_opposite_convexity(p0p, p0, p1, p1n)) {
            return false;
        }

        if ((p0 - p1).squaredNorm() > 1. + PF_EPS) {
            return false;
        }

        return true;
    }
}

bool PostProcessor::_should_accept_trace() const {
    unsigned int degenerate_edges_before = 0;
    for (Index i = 0; i < P.size(); ++i) {
        degenerate_edges_before += is_degenerate(B, P, i, Circular(P, i + 1), circular) ? 1 : 0;
    }

    unsigned int degenerate_edges_after = 0;
    for (Index i = 0; i < _P.size(); ++i) {
        degenerate_edges_after += is_degenerate(B, _P, i, Circular(_P, i + 1), circular) ? 1 : 0;
    }

    if (degenerate_edges_after > degenerate_edges_before) {
        return false;
    }

    // Count symmetries
    unsigned int symmetries_before = regularity.vertex_symmetries().size();    

    return true;
}

// let's not retrace just move the point to a more smooth location.
bool PostProcessor::_relocate_short_edge(const int v_src, const int v_dst, const int v_force) {   
    // Not enough degrees of freedom
    if (P.size() < 6) {
        return false;
    }
    
    PF_VERBOSE_F("Postprocessor graph has %d edges", E.size());
    PF_VERBOSE_F("Relocating %d %d -> %d", v_src, v_dst, v_force);

    // v_force can either be v_src or v_dst or the midpoint
    const bool collapsing_to_midpoint = CircularDist(B, v_src, v_dst) == 2 && Circular(B, v_src + 1) == v_force;
    PF_ASSERT(v_force == v_src || v_force == v_dst || collapsing_to_midpoint);

    if (collapsing_to_midpoint && find(M.begin(), M.end(), v_force) == M.end()) {
        return false;
    }

    // Replace the first half of the subpath
    const bool replace_to_src = collapsing_to_midpoint || v_force == v_dst;
    const bool replace_from_dst = collapsing_to_midpoint || v_force == v_src;

    // I am sorry
    std::vector<int> C;
    BoundaryGraph::trace_to_points(B, P, PP);
    PathUtils::compute_convexities(PP, C);

    // We need to find the polygon indices for this vertices
    int v_src_poly = -1;
    for (Index i = 0; i < P.size(); ++i) {
        if (P(i) == v_src) {
            v_src_poly = i;
            break;
        }
    }
    int v_dst_poly = Circular(P, v_src_poly + 1);

    // This has been collapsed in previous edits
    if (P(v_dst_poly) != v_dst) {
        return false;
    }

    // Making sure that we are not breaking existing symmetries
#if FANCY_SYMMETRIC_LOGIC
    bool would_break_symmetry = false;
    for (const polyfit::Regularity::Symmetry& r : regularity.vertex_symmetries()) {
        const int dist_to_src = CircularDist(P, v_dst_poly, r.v0);
        const int dist_from_dst = CircularDist(P, r.v1, v_src_poly);
        if (min(dist_to_src, dist_from_dst) >= 2) {
            continue;
        }
        
        if (Circular(P, r.v0 + 1) == r.v1 || Circular(P, r.v0 - 1) == r.v1) {
            // We want to preserve the same angle
            const int v_prev_0 = PathUtils::opposite(P.size(), r.v1, r.v0);
            const int v_prev_1 = PathUtils::opposite(P.size(), r.v0, r.v1);

            bool are_symmetric = false;
            for (const polyfit::Regularity::Symmetry& s : regularity.vertex_symmetries(v_prev_0)) {
                if ((s.v0 == v_prev_0 && s.v1 == v_prev_1) ||
                    (s.v1 == v_prev_0 && s.v0 == v_prev_1)) {
                    are_symmetric = true;
                    break;
                }
            }

            PF_DEV_F("angle_0 %d %d %d angle_1 %d %d %d", v_prev_0, r.v0, r.v1, v_prev_1, r.v1, r.v0);
            const float angle_0 = AngleUtils::spanned_shortest(B.col(P(v_prev_0)), B.col(P(r.v0)), B.col(P(r.v1)));
            const float angle_1 = AngleUtils::spanned_shortest(B.col(P(v_prev_1)), B.col(P(r.v1)), B.col(P(r.v0)));
            PF_DEV_F("Vertex symmetry %d %d angle_0 %f s %f", r.v0, r.v1, angle_0, angle_1);

            int c_src = C[v_src_poly];
            int c_dst = C[v_dst_poly];
            int c_next = replace_to_src ? C[Circular(P, v_src_poly - 1)] : C[Circular(P, v_dst_poly + 1)];

            // The purpose of the collapse is to remove visually disruptive inflections, making sure we are removing an inflection
            if (abs(angle_0 - angle_1) < PF_EPS && !are_symmetric && c_src == c_dst && c_next == c_dst) {
                would_break_symmetry = true;
            }
        }
    }

    if (would_break_symmetry) {
        return false;
    }
#endif

    PF_ASSERT(P(v_dst_poly) == v_dst);

    vec3i fix_src, fix_dst;
    if (v_src == v_force) {
        fix_src = vec3i(Circular(P, v_src_poly - 2), Circular(P, v_src_poly - 1), Circular(P, v_src_poly));
        fix_dst = vec3i(Circular(P, v_dst_poly + 2), Circular(P, v_dst_poly +3), Circular(P, v_dst_poly + 4));
    } else if (v_dst == v_force) {
        fix_src = vec3i(Circular(P, v_src_poly - 4), Circular(P, v_src_poly - 3), Circular(P, v_src_poly - 2));
        fix_dst = vec3i(Circular(P, v_dst_poly), Circular(P, v_dst_poly + 1), Circular(P, v_dst_poly + 2));
    } else {
        fix_src = vec3i(Circular(P, v_src_poly - 3), Circular(P, v_src_poly - 2), Circular(P, v_src_poly - 1));
        fix_dst = vec3i(Circular(P, v_dst_poly + 1), Circular(P, v_dst_poly + 2), Circular(P, v_dst_poly + 3));
    }

    PF_VERBOSE_F("fix_src %d %d %d", P(fix_src(0)), P(fix_src(1)), P(fix_src(2)));
    PF_VERBOSE_F("fix_dst %d %d %d", P(fix_dst(0)), P(fix_dst(1)), P(fix_dst(2)));

    vector<bool> corners_C0(B.cols(), false);
    if (collapsing_to_midpoint) {
        corners_C0[P(fix_src(2))] = true;
        corners_C0[P(fix_dst(0))] = true;
    } else {
        corners_C0[v_force] = true;
    }

    vector<BoundaryGraph::Node> N;
    BoundaryGraph::construct_nodes(B, E, M, N);
    vector<ShortestPath::Node> G;
    BoundaryGraph::construct_trace_graph(B, E, M, N, G, true, corners_C0);

    Index hn_src = -1, hn_dst = -1;
    for (Index i = 0; i < G.size(); ++i) {
        const auto& n = N[G[i].v];

        if (hn_src != -1 && hn_dst != -1) {
            break;
        }

        if (n.v0 == P(fix_src(0)) && n.v1 == P(fix_src(1)) && n.v2 == P(fix_src(2))) {
            PF_ASSERT(hn_src == -1);
            hn_src = i;
        }

        if (n.v0 == P(fix_dst(0)) && n.v1 == P(fix_dst(1)) && n.v2 == P(fix_dst(2))) {
            PF_ASSERT(hn_dst == -1);
            hn_dst = i;
        }
    }

    // The path has changed since the last time
    if (hn_src == -1 || hn_dst == -1) {
        return false;
    }

    vector<Index> subP;
    ShortestPath::find(_G_state, G, hn_src, hn_dst, subP);

    for (Index v : subP) {
        PF_VERBOSE_F("%d", N[G[v].v].v1);
    }

    int v_remove = fix_src(2);
    const int v_remove_last = P(fix_dst(1));
    while (P(v_remove) != v_remove_last) {
        PF_VERBOSE_F("%d != %d", P(v_remove), v_remove_last);
        EraseOrdered(P, { v_remove });
        v_remove = Circular(P, v_remove);
        PF_ASSERT(P.size() > 0);
    }
    
    PF_VERBOSE_F("Remove index %d (%d)", v_remove, P(v_remove));
    for (size_t i = 1; i < subP.size() - 1; ++i) {
        MatrixUtils::insert_at(P, v_remove + i - 1, N[G[subP[i]].v].v1);
    }

    PF_VERBOSE_F("Total length %d", P.size());
    for (Index i = 0; i < P.size(); ++i) {
        PF_VERBOSE_F("P %d", P(i));
    }

    BoundaryGraph::trace_to_points(B, P, PP);
    regularity.reset();
    regularity.update(B, P, raster_symmetries, raster_symmetries_local, true);

    auto candidates = Regularity::find_continuation_candidates(B, P, regularity, circular, false);
    Regularity::pick_continuations(B, P, candidates, regularity, circular, E);

    return true;
}

/*
    Removes the edge connecting v_src and v_dst and forces the path to go through v_force.
    Returns true if the polygon has been updated without introducing new degenerate edges or inflections.
    If the operation has succeded the list of edges which have been delete is added to the permanent list of deleted 
    edges.
*/
bool PostProcessor::_attempt_remove_edge_force_corner ( const int v_src, const int v_dst, const int v_force ) {
    _E = _E_original;
    _E_delete.clear();
    _E_delete_temp.clear();

    // This whole thing needs to be rewritten anyway
    PF_VERBOSE_F("Remove edge %d %d", v_src, v_dst);
    int v_poly = -1;
    for (Index i = 0; i < P.size(); ++i) {
        if (P(i) == v_src) {
            v_poly = i;
            break;
        }
    }
    PF_ASSERT(v_poly != -1);
    const vec3i v_fix_before(CircularAt(P, v_poly - 5), CircularAt(P, v_poly - 4), CircularAt(P, v_poly - 3));
    const vec3i v_fix_after(CircularAt(P, v_poly + 2), CircularAt(P, v_poly + 3), CircularAt(P, v_poly + 4));
    PF_VERBOSE_F("Fix before (%d %d %d)", v_fix_before(0), v_fix_before(1), v_fix_before(2));
    PF_VERBOSE_F("Fix after (%d %d %d)", v_fix_after(0), v_fix_after(1), v_fix_after(2));

    for (size_t i = 0; i < _E.size(); ++i) {
        if (_E[i].v0 == v_src && _E[i].v1 == v_dst) {
            _E_delete.emplace_back(i);
            continue;
        }

        if (PathUtils::contains_open(B, _E[i].v0, _E[i].v1, v_force, true)) {
            _E_delete.emplace_back(i);
            continue;
        }

        if ((_E[i].v0 == v_fix_before(0) && _E[i].v1 == v_fix_before(1)) ||
            (_E[i].v0 == v_fix_before(1) && _E[i].v1 == v_fix_before(2)) ||
            (_E[i].v0 == v_fix_after(0) && _E[i].v1 == v_fix_after(1)) ||
            (_E[i].v0 == v_fix_after(1) && _E[i].v1 == v_fix_after(2))) {
            continue;
        }

        if (PathUtils::contains_open(B, v_fix_after(0), v_fix_after(2), _E[i].v0, true) ||
            PathUtils::contains_open(B, v_fix_after(0), v_fix_after(2), _E[i].v1, true) ||
            PathUtils::contains_open(B, v_fix_before(0), v_fix_before(2), _E[i].v0, true) ||
            PathUtils::contains_open(B, v_fix_before(0), v_fix_before(2), _E[i].v1, true)) {
            _E_delete_temp.emplace_back(i);
        }
    }

    _E_delete_all.clear();
    _E_delete_all.insert(_E_delete_all.end(), _E_delete.begin(), _E_delete.end());
    _E_delete_all.insert(_E_delete_all.end(), _E_delete_permanent.begin(), _E_delete_permanent.end());
    //_E_delete_all.insert(_E_delete_all.end(), _E_delete_temp.begin(), _E_delete_temp.end());
    sort(_E_delete_all.begin(), _E_delete_all.end());
    EraseDuplicatesOrdered(_E_delete_all);
    EraseOrdered(_E, _E_delete_all);

    ShortestPath::State G_state;
    PolygonTracer::regularized(G_state, B, M, _E, _RE, _P, raster_symmetries, raster_symmetries_local, true);

    if (_P.size() > 0 && _should_accept_trace()) {
        P = _P;
        regularity = _RE;
        BoundaryGraph::trace_to_points(B, P, PP);

        // Making the current delete list permanent
        _E_delete_permanent.insert(_E_delete_permanent.end(), _E_delete.begin(), _E_delete.end());
        return true;
    } 
    return false;
}

void PostProcessor::_reset_immovable_corners() {
    _immovable.resize(B.cols(), false);
}

/*
    Setting every corner which is a critical part of a regularity to 0
*/
void PostProcessor::_set_immovable_corners() {
    double max_edge_length_sq = Options::get()->post_processor_max_short_edge_length;
    max_edge_length_sq *= max_edge_length_sq;

    PF_VERBOSE_F("Setting immovable corners with %d regularities", (int)regularity.size());

    _reset_immovable_corners();
	for (auto& r : regularity.parallels()) {
		PF_VERBOSE_S("Regularity parallel");
		if (r.connected_00_11) {
			_immovable[P(r.v00)] = true;
			_immovable[P(r.v11)] = true;
			PF_VERBOSE_F("Making %d immovable", r.v00);
			PF_VERBOSE_F("Making %d immovable", r.v11);
		}

		if (r.connected_01_10) {
			_immovable[P(r.v01)] = true;
			_immovable[P(r.v10)] = true;
			PF_VERBOSE_F("Making %d immovable", r.v01);
			PF_VERBOSE_F("Making %d immovable", r.v10);
		}
	}

	for (auto& r : regularity.continuations()) {
		PF_VERBOSE_S("Regularity continuation");
		//_immovable[P(r.data.continuation.v0_prev)] = true;
		//_immovable[P(r.data.continuation.v0)] = true;
		//_immovable[P(r.data.continuation.v1)] = true;
		//_immovable[P(r.data.continuation.v1_next)] = true;
	}	

    // We also do not want to break axis-aligned lines
    for (Index i = 0; i < P.size(); ++i) {
        const int v = P(i);
        const int vn = CircularAt(P, i + 1);
        const vec2 p = B.col(v);
        const vec2 pn = B.col(vn);
        
        if (!GeomRaster::are_points_axis_aligned(p, pn)) {
            continue;
        }

        const double len_sq = (p - pn).squaredNorm();
        if (len_sq < max_edge_length_sq + PF_EPS) {
            continue;
        }

        _immovable[v] = true;
        _immovable[vn] = true;
    }

    // Aswell as axis-aligned angles
    for (Index i = 0; i < P.size(); ++i) {
        const int v = P(i);
        const int vn = CircularAt(P, i + 1);
        const int vnn = CircularAt(P, i + 2);

        const double angle_vn = AngleUtils::spanned_shortest(B.col(v), B.col(vn), B.col(vnn));
        const double len_sq = min((B.col(v) - B.col(vn)).squaredNorm(), (B.col(vn) - B.col(vnn)).squaredNorm());
        if (abs(angle_vn - M_PI_2) < PF_EPS && len_sq > 1. + PF_EPS) {
            _immovable[v] = true;
            _immovable[vn] = true;
            _immovable[vnn] = true;
        }
    }

    // We shouldn't be moving points which are along the axis of symmetry
    // That is exactly a midpoint!
    for (Index i = 0; i < M.size(); ++i) {
        _immovable[M[i]] = true;
    }
}

const std::vector<bool>& PostProcessor::get_immovable_corners() const {
    return _immovable;
}

const std::vector<BoundaryGraph::Edge>& PostProcessor::get_original_edges() const {
    return _E_original;
}

const std::vector<size_t>& PostProcessor::get_permanently_deleted_edges() const {
    return _E_delete_permanent;
}

void PostProcessor::merge_curves()
{
	polyfit::BezierMerging::MergeRecursivelyOptions options;	
	options.allow_merging_lines = true;
	options.distance_tolerance = 0.25;
	options.alpha_max = 0.8;
	options.per_bezier_const = 1.5;

	auto symmetric_curves = Regularity::find_symmetric_primitives(regularity.edge_symmetries(), curves, P.size());

	std::vector<Eigen::Matrix2Xd> in_curves, out_curves;
	std::vector<std::vector<int>> out2in;

	for (auto& p : curves.primitives)
	{
		auto c = p.curve->get_curve();
		auto bezier = dynamic_pointer_cast<BezierCurve>(c);
		auto line = dynamic_pointer_cast<GlobFitCurve_Line>(c);

		if (bezier)
		{
			in_curves.push_back(bezier->get_control_points());
			continue;
		}

		if (line)
		{
			in_curves.push_back(line->get_points());
			continue;
		}

		throw std::runtime_error("There is a primitive whose type was not anticipated during curve merging.");
	}

	polyfit::BezierMerging::merge_curves(
		in_curves, true, options, symmetric_curves, out_curves, out2in);

	CurvePrimitiveSequence new_sequence;
	std::vector<TangentFitType> new_tangent_fits;
	for (int i = 0; i < out_curves.size(); ++i)
	{
		new_sequence.primitives.emplace_back();
		if (out_curves[i].cols() == 2)
			new_sequence.primitives.back().curve = std::make_shared<GlobFitCurveParametrization>(std::make_shared<GlobFitCurve_Line>(out_curves[i]));
		else if (out_curves[i].cols() == 4)
			new_sequence.primitives.back().curve = std::make_shared<GlobFitCurveParametrization>(std::make_shared<BezierCurve>(out_curves[i]));
		else
			throw std::runtime_error("Curve merging produced a curve with " + std::to_string(out_curves[i].cols()) + " control points");

		// dummy type
		new_tangent_fits.push_back(TANGENT_FIT_CONSTANT);
	}

	curves = std::move(new_sequence);
	tangent_fits = std::move(new_tangent_fits);
}

NAMESPACE_END(polyvec)