// polyvec
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/core/options.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/regularity/symmetry.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/core/constants.hpp>

// libc++
#include <unordered_map>
#include <algorithm>

#define ALLOW_VERTICES_NEXT_TO_CORNERS 0

using namespace Eigen;
using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(BoundaryGraph)

Edge Edge::make_pixel_diagonal(const int v0, const int v1) {
    BoundaryGraph::Edge e;
    e.v0 = v0;
    e.v1 = v1;
    e.d_min = 0.;
    e.d_max = .5;
    e.flags = 0x0;
    return e;
}

Edge Edge::make_flat(const int v0, const int v1) {
    BoundaryGraph::Edge e;
    e.v0 = v0;
    e.v1 = v1;
    e.d_min = e.d_max = 0.;
    e.flags = 0x0;
    return e;
}

EdgeID make_edge_id(const int c0, const int c1) {
    PF_ASSERT(c0 <= INT32_MAX);
    PF_ASSERT(c1 <= INT32_MAX);
    const uint64_t c0_as_u64 = (c0 & 0xffffffff);
    const uint64_t c1_as_u64 = (c1 & 0xffffffff);
    return (c0_as_u64 << 32) | (c1_as_u64);
}

vec2i unpack_edge_id(const EdgeID id) {
    vec2i c;
    c(0) = (int)((id >> 32) & 0xffffffff);
    c(1) = (int)(id & 0xffffffff);
    return c;
}

// Inserts midpoints on small flat segments 
void create_fitting_path(const mat2x& R, mat2x& P, std::vector<Vertex>& M, bool circular) {
	P.resize(2, 0);
	M.clear();
			
	const Vertex voff = circular ? 0 : 1;
	Vertex vstart = 0;

	// flat points are set to 0
	std::vector<int> convexity;
	PathUtils::compute_convexities(R, convexity);

	if (R.cols() > 3) {
		// Starting from the first non-flat vertex
		for (Vertex i = voff; i < R.cols() - voff; ++i) {
			if (convexity[i] != 0) {
				vstart = i;
				break;
			}
		}
	}

	bool R_ccw;
	WindingNumber::compute_orientation(R, R_ccw);
	const int orientation = R_ccw ? -1 : +1;

	for (Vertex i = 0; i < R.cols(); ++i) {
		const Vertex v = Circular(R, vstart + i);

		// finding the next convex corner
		int v_cvx = 0;
		if (circular || i < R.cols() - 1) {
			v_cvx = convexity[v];
		}

		if (i == 0 || abs((R.col(v) - P.col(P.cols() - 1)).norm()) > PF_EPS) {
			//PF_LOGF("appending vertex %lld pos %f %f", v, R.col(v).x(), R.col(v).y());
			MatrixUtils::append(P, R.col(v));
		}

		if (v_cvx != 0) {
			for (Vertex j = 1; j <= PV_MIDPOINT_EDGE_LENGTH; ++j) {
				if (!circular && v + j >= R.cols()) {
					break;
				}

				const Vertex vn = Circular(R, v + j);
				if (convexity[vn] == convexity[v]) {
					// The points should be symmetric, same length and convexity expected
					int vp_cvx = 0, vnn_cvx = 0;
					double vp_dist, vnn_dist;

					Vertex k = 1;
					while (vp_cvx == 0 && k < R.cols()) {
						if (!circular && v - k < 0) {
							break;
						}

						vp_cvx = convexity[Circular(R, v - k)];
						vp_dist = (CircularAt(R, v - k) - R.col(v)).norm();
						++k;
					}

					k = 1;
					while (vnn_cvx == 0 && k < R.cols()) {
						if (!circular && vn + k >= R.cols()) {
							break;
						}

						vnn_cvx = convexity[Circular(R, vn + k)];
						vnn_dist = (CircularAt(R, vn + k) - R.col(vn)).norm();
						++k;
					}

					bool same_cvx = vp_cvx == vnn_cvx && vp_cvx != 0;
					bool same_length = abs(vnn_dist - vp_dist) < 1e-4;

					if (same_cvx && same_length) {
						M.emplace_back(P.cols());
						MatrixUtils::append(P, .5 * (R.col(v) + R.col(vn)));
					}

					break;
				}

				if (convexity[vn] != 0) {
					break;
				}
			}
		}
	}
}

// returns true if the any of the edge endpoints is within the symmetric region, 
// but not symmetric. 
bool is_edge_violating_axis_aligned_edge(
    const mat2x & P, 
    const std::vector<int>& C,
    const Edge & e, 
    const Vertex a0, 
    const Vertex a1, 
    const bool circular,
    const double P_bb_diagonal
) {
	// overlaps
	if (e.v0 == a0 && e.v1 == a1) {
		return false;
	}

	// midpoints
    //const double length_visible_corner = max(1, std::floor(.05 * P_bb_diagonal));

    const Vertex a0_prev = PathUtils::next_transition_vertex(C, a0, -1, circular);
    const Vertex a1_next = PathUtils::next_transition_vertex(C, a1, +1, circular);
	double dist_a0a1 = (P.col(a0) - P.col(a1)).norm();
    const double side_length_0 = (P.col(a0_prev) - P.col(a0)).norm();
    const double side_length_1 = (P.col(a1_next) - P.col(a1)).norm();
    const bool old_test = side_length_0 < dist_a0a1 + PF_EPS || side_length_1 < dist_a0a1 + PF_EPS;
    const bool new_test = min(side_length_0, side_length_1) < dist_a0a1 * P_bb_diagonal + PF_EPS;

    if (dist_a0a1 < 2. + 1e-4 &&
        old_test) {
        //new_test) {
        const vec2 pmid = .5 * (P.col(a0) + P.col(a1));

        if ((P.col(e.v1) - pmid).norm() < 1e-4 ||
            (P.col(e.v0) - pmid).norm() < 1e-4) {
            return false;
        }
    }

	// crossing
	if (PathUtils::contains_closed(P, e.v0, e.v1, a0, circular)) {
		return e.v1 != a0;
	}

	if (PathUtils::contains_closed(P, e.v0, e.v1, a1, circular)) {
		return e.v0 != a1;
	}

	return false;
}
bool is_edge_crossing_subpath_asymmetric(const mat2x& P, const Edge& e, const Symmetry::SubPath& s, bool circular) {
	// crossing
	switch (s.axis) {
	case Symmetry::AXIS_HORIZONTAL:
	case Symmetry::AXIS_VERTICAL:
		return (e.v0 == s.v01 && e.v1 == s.v10) ||
			(e.v0 == s.v00 && e.v1 == s.v11) ||
			(e.v0 == s.v00 && e.v1 != s.v10) ||
			(e.v0 != s.v00 && e.v1 == s.v10);

	case Symmetry::AXIS_SEC_QUAD02:
	case Symmetry::AXIS_SEC_QUAD13: {
		break;
	}
	default:
		PF_ABORT;
	}

	return false;
}

vec2i calc_subpath_region_bounds(const std::vector<int>& C, const Symmetry::SubPath& s, bool circular) {
	vec2i bounds;
	bounds(0) = (int)PathUtils::end_of_convex_region(C, s.v01, -1, circular);
	bounds(1) = (int)PathUtils::end_of_convex_region(C, s.v11, +1, circular);
	return bounds;
}

bool is_edge_in_subpath_region(const mat2x& P, const Edge& e, const vec2i bounds, bool circular) {
	return PathUtils::contains_closed(P, bounds(0), bounds(1), e.v0, circular) ||
		PathUtils::contains_closed(P, bounds(0), bounds(1), e.v1, circular);
}

// Adding the edge if it's not a single pixel long to the subgraph
void append_edge_to_subpath_region_graph_if_not_degenerate(const mat2x& P, vector<ShortestPath::Node>& G, const Edge& e) {
	if ((P.col(e.v0) - P.col(e.v1)).norm() < 1. + PF_EPS) {
		//SYM_LOG("discarding degenerate edge %lld %lld", e.v0, e.v1);
		return;
	}

	Vertex nsrc = ShortestPath::find_or_add_node(G, e.v0);
	Vertex ndst = ShortestPath::find_or_add_node(G, e.v1);

	//SYM_LOG("adding edge %lld %lld nodes %lld(%lld) %lld(%lld)", e.v0, e.v1, nsrc, G[nsrc].v, ndst, G[ndst].v);

	// finally.. appending the neighbor
	G[nsrc].add_neighbor(ndst, 1.);

#ifdef PF_DEBUG
	G[nsrc].id = StringUtils::fmt("%lld -> %lld", e.v0, e.v1);
#endif
}

void remove_asymmetric_edges(const mat2x& B, const std::vector<Vertex>& M, std::vector<Edge>& E, bool circular) {
	std::vector<int> C;
	PathUtils::compute_convexities(B, C);

	std::vector<Symmetry::SubPath> S = Symmetry::find_shortest(B);

	// todo: this should be moved in the symmetry detection code or put back in the loop below
	for (int i = 0; i < S.size(); ++i) {
		PF_ASSERT(C[S[i].v00] == C[S[i].v10]);
	}

	// preference to convex then length
	std::sort(S.begin(), S.end(), [&](const Symmetry::SubPath& lhs, const Symmetry::SubPath& rhs) -> bool {
		if (C[lhs.v00] == +1 && C[rhs.v00] == -1) {
			return true;
		}
		if (C[lhs.v00] == -1 && C[rhs.v00] == +1) {
			return false;
		}
				
		return CircularDist(B, lhs.v00, lhs.v10) > CircularDist(B, rhs.v00, rhs.v10);
	});

	// temporary graph used to check if applying a certain symmetry (removing the edges)
	// leaves only paths with bad edges (In this stage we haven't constructed the graph nodes yet)
	std::vector<ShortestPath::Node> Gsub;

	// holds the edges to be removed in the current iteration
	std::vector<size_t> E_dlist;

	for (int i = 0; i < S.size(); ++i) {
		Gsub.clear();
		E_dlist.clear();

		PF_VERBOSE_F("candidate symmetry v00 %lld v01 %lld v10 %lld v11 %lld", S[i].v00, S[i].v01, S[i].v10, S[i].v11);

		// finds the firsts and last corner (before/after) the symmetric region
		// which defines the region where we expect to not have a degenerate path
		// these are the first two changes in convexity
		vec2i Gsub_bounds = calc_subpath_region_bounds(C, S[i], circular);
        PF_VERBOSE_F("convex region bounds %lld %lld", Gsub_bounds(0), Gsub_bounds(1));

		// creating the graph nodes for the bounds
		Gsub.emplace_back(Gsub_bounds(0));
		Gsub.emplace_back(Gsub_bounds(1));

		// Find all the edges to be removed while constructing the neighboring information
		// a map could be used here to only choose v00 v01 v10 v11
		for (int j = 0; j < E.size(); ++j) {
			// is this node within the region we are interested in? 
			if (!is_edge_in_subpath_region(B, E[j], Gsub_bounds, circular)) {
				continue;
			}

			// is it asymmetric ? 
			if (is_edge_crossing_subpath_asymmetric(B, E[j], S[i])) {
				E_dlist.emplace_back(j);
			}
			// let's keep it for now
			else {
				append_edge_to_subpath_region_graph_if_not_degenerate(B, Gsub, E[j]);
			}
		}

		// todo: check if this is how we want to approach this models
		bool skip = S[i].axis == Symmetry::AXIS_SEC_QUAD02 || S[i].axis == Symmetry::AXIS_SEC_QUAD13;

		// the subgraph containing what will be the valid edges in the region around the symmetry
		// should now contain only valid edges. Since append_edge_to_subpath_region_* only added valid
		// ones, any graph traversal should reach the destination
		if (!skip && !ShortestPath::find_breadth_first(Gsub, 0, 1)) {
            PF_VERBOSE_F("failed to find valid connection %d", 0);
			continue;
		}

        PF_VERBOSE_F("path found, removing %lld edges", E_dlist.size());
		for (int j = 0; j < E_dlist.size(); ++j) {
            PF_VERBOSE_F("removed edge %lld: %lld -> %lld", E_dlist[j], E[E_dlist[j]].v0, E[E_dlist[j]].v1);
		}

		Vertex d = CircularDist(B, S[i].v00, S[i].v10);
		if (d <= 2) {
			for (int j = 0; j < E_dlist.size(); ++j) {
				PF_LOGF("deleting edge %lld -> %lld", E[E_dlist[j]].v0, E[E_dlist[j]].v1);
			}
		}

		// If the path exists, we can safely delete the edges
		// this has to be done at every iteration to properly check for existing paths
		EraseOrdered(E, E_dlist);
	}
}

// Creates all edges which pass the accuracy threshold
void connect_valid_edges(
    const mat2x& P, 
    const std::vector<Vertex>& M, 
    std::vector<Edge>& edges, 
    bool circular,
    std::unordered_map<EdgeID, Eigen::Vector2d>* accuracy_map,
    const std::vector<bool>& _B_flat,
    const bool force_B_flat
) {
    std::vector<bool> B_flat = _B_flat;
    if (B_flat.empty()) {
        B_flat.resize(P.cols(), false);
    }

	std::vector<int> convexity;
	PathUtils::compute_convexities(P, convexity);

    unordered_map<int, vec2i> B_flat_force;
    assert(B_flat.size() == P.cols());
    for (Index i = 0; i < P.cols(); ++i) {
        if (B_flat[i]) {
            const int vertex_prev = PathUtils::next_transition_vertex(convexity, i, -1);
            const int vertex_next = PathUtils::next_transition_vertex(convexity, i, +1);

            if (CircularDist(P, vertex_prev, vertex_next) > 2) {
                B_flat_force[i] = vec2i(vertex_prev, vertex_next);
            }
            else {
                B_flat_force[i] = vec2i(-1, -1);
            }
        }
    }

    if (accuracy_map) {
        accuracy_map->clear();
    }

	// returns if a boundary point can act as a graph vertex
	auto is_vertex = [&](Vertex i, bool& out_is_midpoint, bool& out_is_flat) {
		out_is_midpoint = std::find(M.begin(), M.end(), i) != M.end();		

		if (out_is_midpoint)
			return true;

		if (B_flat[i])
			return true; //we allow outgoing edges from these vertices

		out_is_flat = AngleUtils::is_flat(CircularAt(P, i - 1), P.col(i), CircularAt(P, i + 1));

		if (out_is_flat)
		{
#if ALLOW_VERTICES_NEXT_TO_CORNERS
			// if this point is flat, we only allow it as a vertex if it is next to a corner
			return !AngleUtils::is_flat(CircularAt(P, i - 2), CircularAt(P, i - 1), CircularAt(P, i))
				|| !AngleUtils::is_flat(CircularAt(P, i), CircularAt(P, i + 1), CircularAt(P, i + 2));
#else
			// never allow flat vertices to be corners
			return false;
#endif			
		}

		return true;
	};

	for (Vertex i = 0; i < P.cols(); ++i) {
		const vec2 vi = P.col(i);

		bool is_i_midpoint, is_i_flat;
		if (!is_vertex(i, is_i_midpoint, is_i_flat))
			continue;		

		for (Vertex ioff = 1; ioff < P.cols(); ++ioff) {
			if (!circular && i + ioff >= P.cols()) {
				break;
			}
					
			const Vertex j = Circular(P, i + ioff);
			const vec2 vj = P.col(j); 

			// skipping connections to flat vertices	
			bool is_j_midpoint, is_j_flat;
			if (!is_vertex(j, is_j_midpoint, is_j_flat))
				continue;
			
			// prohibit connections between midpoints and their adjacent corners
			if ((is_i_midpoint || is_j_midpoint) && PathUtils::are_axis_aligned(vi, vj)) {
				continue;
			}

			// Avoid cutting corners between flat vertices
			if (is_i_flat && is_j_flat && ((vi - vj).cwiseAbs() - vec2(1.0, 1.0)).squaredNorm() < PF_EPS)
				continue;

            // Forcing flat vertices to be flat if possible
            if (force_B_flat && B_flat[i]) {
                const int vertex_next = B_flat_force[i](1);
                if (vertex_next != -1 && j != vertex_next && !B_flat[j] && CircularDist(P, i, j) > 2) {
                    continue;
                }
            }
            if (force_B_flat && B_flat[j]) {
                const int vertex_prev = B_flat_force[j](0);
                if (vertex_prev != -1 && i != vertex_prev && !B_flat[i] && CircularDist(P, i, j) > 2) {
                    continue;
                }
            }

            vec2 err = GeomRaster::distance_bounds_from_points_with_slack(P, P.col(i), P.col(j), i, j);
            const double phi = atan2((vj - vi).cwiseAbs().minCoeff(), (vj - vi).cwiseAbs().maxCoeff());
            const double d_eps = .1;
            const double d_max = 0.5 + (.5 + d_eps) * cos(phi) - .5 * sin(phi);

            if (err.maxCoeff() > d_max) {
                continue;
            }

            if (!ErrorMetrics::accuracy_within_bounds(err.minCoeff(), err.maxCoeff(), vj - vi)) {
                if (err.maxCoeff() >= 2. - PF_EPS) {
                    break;
                }
				continue;
			}

			PF_VERBOSE_F("edge %lld -> %lld error %f %f", i, j, err.x(), err.y());

			Edge edge;
			edge.d_min = err.minCoeff();
			edge.d_max = err.maxCoeff();
			edge.v0 = i;
			edge.v1 = j;
			edge.flags = 0x0;

			// precomputing the flags used to calculate the error
			if (!circular && i == 0) {
				edge.flags |= ErrorMetrics::EDGE_FIRST_IN_PATH;
			}

			if (!circular && i == P.cols() - 1) {
				edge.flags |= ErrorMetrics::EDGE_LAST_IN_PATH;
			}

			if (circular || i > 0 && i < P.cols() - 2) {
				const vec2 p0 = CircularAt(P, i - 1);
				const vec2 p1 = CircularAt(P, i);
				const vec2 p2 = CircularAt(P, i + 1);
				const vec2 p3 = CircularAt(P, i + 2);

				if (AngleUtils::have_opposite_convexity(p0, p1, p2, p3)) {
					edge.flags |= ErrorMetrics::EDGE_HAS_INFLECTION;
				}
			}

            if (accuracy_map) {
                const uint64_t edge_id = make_edge_id(i, j);
                PF_ASSERT(accuracy_map->find(edge_id) == accuracy_map->end());
                accuracy_map->insert(make_pair(edge_id, err));
            }

			edges.emplace_back(edge);
		}
	}
}

int construct_nodes (
    const mat2x& B,
    const std::vector<Edge>& E,
    const std::vector<Vertex>& M,
    std::vector<Node>& N,
    bool circular,
    const std::vector<bool>& _B_flat
) {
    std::vector<bool> B_flat = _B_flat;
    if (_B_flat.empty()) {
        B_flat.resize(B.cols(), false);
    }

	N.clear();

	// tracking neighbors for each point in both directions
	vector<vector<Vertex>> efrom;
	vector<vector<Vertex>> eto;
	efrom.resize(B.size());
	eto.resize(B.size());

	for (size_t i = 0; i < E.size(); ++i) {
		const Edge& e = E[i];

		efrom[e.v0].emplace_back(i);
		eto[e.v1].emplace_back(i);
	}

	bool is_ccw;
	WindingNumber::compute_orientation(B, is_ccw);

	const double flat_angle_thr = M_PI - PF_RAD(5);

	// transforming the list of midpoints into a bitmap indicating for every node whether it is a midpoint or not
	std::vector<bool> midpoint(B.size(), false);
	for (size_t i = 0; i < M.size(); ++i) {
		midpoint[M[i]] = true;
	}

	// convexities for each point
	std::vector<int> C;
	PathUtils::compute_convexities(B, C);

	// Keeping track of which small symmetric region each point belongs to 
	std::vector<Symmetry::SubPath> S_short = Symmetry::find_shortest(B, circular);
	std::vector<int> SYM(B.size(), -1);
	for (int i = 0; i < S_short.size(); ++i) {
		SYM[S_short[i].v00] = i;
		SYM[S_short[i].v10] = i;
	}

    // Assuming that 
    std::unordered_set<uint64_t> active_edges;
    for (size_t i = 0; i < E.size(); ++i) {
        active_edges.insert(make_edge_id(E[i].v0, E[i].v1));
    }

	// node = (i, j, k)
	for (Vertex j = 0; j < (Vertex)B.size(); ++j) {
		for (Vertex ii = 0; ii < (Vertex)eto[j].size(); ++ii) {
			for (Vertex kk = 0; kk < (Vertex)efrom[j].size(); ++kk) {
				const Vertex i = E[eto[j][ii]].v0;
				const Vertex k = E[efrom[j][kk]].v1;

				const vec2 vi = B.col(i);
				const vec2 vj = B.col(j);
				const vec2 vk = B.col(k);

				// is the angle flat ? 
				const double angle = AngleUtils::spanned_shortest(vi, vj, vk);
				bool is_flat = !B_flat[j] && angle > flat_angle_thr;
                bool is_null = angle < PF_EPS;
               
                if (is_null) {
                    continue;
                }

                const uint64_t edge_id = make_edge_id ( i, k );
                const bool is_straight_alternative_valid = active_edges.find(edge_id) != active_edges.end();
				if (!is_flat || !is_straight_alternative_valid) {
					Node n;
					n.v0 = i;
					n.v1 = j;
					n.v2 = k;
					n.e0 = eto[j][ii];
					n.e1 = efrom[j][kk];
					n.angle = isnan(angle) ? M_PI : angle;
                    //PF_DEV_F("Node %d %d %d", i, j, k);
					N.emplace_back(n);
				}
			}
		}
	}

	return (int)N.size();
}

inline bool point_on_symmetry(const mat2x& B, const std::vector<Symmetry::SubPath>& S, const Vertex v, bool circular) {
	return find_if(S.begin(), S.end(), [&B, v, circular](const Symmetry::SubPath & s) {
		return PathUtils::contains_open(B, s.v00, s.v10, v, circular);
		}) != S.end();
}

using EdgeKey = std::pair<Vertex, Vertex>;

int construct_trace_graph(const mat2x& B, 
	const std::vector<Edge>& E, 
	const std::vector<Vertex>& M, 
	const std::vector<Node>& N, 
	std::vector<ShortestPath::Node>& GN,
	bool circular,
    const std::vector<bool>& corners_c0) {

    PF_ASSERT(corners_c0.empty() || corners_c0.size() == B.cols());

	std::unordered_map <EdgeKey, std::vector<Vertex>, _StdHashVertexPair> nodes_in;
	std::unordered_map <EdgeKey, std::vector<Vertex>, _StdHashVertexPair> nodes_out;

	GN.clear();

	for (size_t i = 0; i < N.size(); ++i) {
		EdgeKey key0 = EdgeKey(N[i].v0, N[i].v1);
		EdgeKey key1 = EdgeKey(N[i].v1, N[i].v2);

		nodes_out[key0].emplace_back(i);
		nodes_in[key1].emplace_back(i);
		//PF_LOGF("node %lld %lld %lld key0 (%lld %lld) key1 (%lld %lld)", N[i].v0, N[i].v1, N[i].v2, key0.first, key0.second, key1.first, key1.second);
	}

	// transforming the list of midpoints into a bitmap indicating for every node whether it is a midpoint or not
	std::vector<bool> midpoint(B.cols(), false);
	for (size_t i = 0; i < M.size(); ++i) {
		midpoint[M[i]] = true;
	}

	// todo: this is done twice, it should be either passed as arguments or remove_asymmetric nodes integrated here
	// we avoid inserting asymmetric edges when connecting the neighbors.
	// first, consistently ordering the symmetries
	std::vector<Symmetry::SubPath> S_short = Symmetry::find_shortest(B, circular);
	std::vector<int> SYM(B.size(), -1);
	PF_VERBOSE_F("Symmetries %llu\n---------", S_short.size());
	for (int i = 0; i < S_short.size(); ++i) {
		PF_VERBOSE_F("v00 %lld v01 %lld v10 %lld v11 %lld", S_short[i].v00, S_short[i].v01, S_short[i].v10, S_short[i].v11);
		SYM[S_short[i].v00] = i;
		SYM[S_short[i].v10] = i;
	}

	GN.resize(N.size());

	for (size_t i = 0; i < N.size(); ++i) {
		const auto& n0 = N[i];
		ShortestPath::Node& node = GN[i];
		node.v = i;

#ifdef PF_DEBUG
		node.id = StringUtils::fmt("(%lld %lld %lld)", n0.v0, n0.v1, n0.v2);
#endif

		EdgeKey key1 = EdgeKey(N[i].v1, N[i].v2);

		if (nodes_out.find(key1) == nodes_out.end()) {
			continue;
		}

		// connecting to all the nodes which have the incoming edge shared with the outgoing
		// edge of n0.
		auto& nhbs = nodes_out[key1];
		for (Vertex nhb : nhbs) {
			const auto& n1 = N[nhb];			
                    
            const bool is_n0_C0 = corners_c0.empty() ? false : corners_c0[n0.v1];
            const bool is_n1_C0 = corners_c0.empty() ? false : corners_c0[n1.v1];

			node.add_neighbor(nhb, approx_edge_cost(
                B, midpoint, E, N[i], N[nhb], circular, is_n0_C0, is_n1_C0
            ));
		}
	}

	return (int)GN.size();
}

int force_graph_subpath(
    const mat2x& B,                 // Boundary fitting points
    std::vector<Edge>& E,           // Edges that will be pruned
    const vecXi& subpath,           // Subpath to be forced
    const AccuracyMap& accuracy_map // Accuracy for every edge in the subpath
) {
    if (subpath.size() == 0) {
        return 0;
    }

    std::vector<size_t> E_delete;

    // Can these two logics be merged?
    if (subpath.size() == 1) {
        for (size_t i = 0; i < E.size(); ++i) {
            if (PathUtils::contains_open(B, E[i].v0, E[i].v1, subpath(0), true)) {
                E_delete.emplace_back(i);
            }
        }

        EraseOrdered(E, E_delete);
    }
    else {
        const int B_src = subpath(0);
        const int B_dst = subpath(subpath.size() - 1);

        E_delete.clear();
        for (size_t i = 0; i < E.size(); ++i) {
            if (E[i].v0 == subpath(subpath.size() - 1) && E[i].v1 == subpath(0)) {
                continue;
            }

            // Removing all the edges within the region
            if (PathUtils::contains_closed(B, B_src, B_dst, E[i].v0, true) &&
                PathUtils::contains_closed(B, B_src, B_dst, E[i].v1, true)) {
                E_delete.emplace_back(i);
                PF_VERBOSE_F("Removing %d %d", E[i].v0, E[i].v1);
                continue;
            }

            // And those crossing the region endpoints
            if (PathUtils::contains_open(B, E[i].v0, E[i].v1, B_src, true) ||
                PathUtils::contains_open(B, E[i].v0, E[i].v1, B_dst, true)) {
                E_delete.emplace_back(i);
                PF_VERBOSE_F("Removing %d %d", E[i].v0, E[i].v1);
                continue;
            }
        }

        EraseOrdered(E, E_delete);

        if (subpath.size() > 1) {
            for (Index i = 0; i < subpath.size() - 1; ++i) {
                BoundaryGraph::Edge e;
                e.v0 = subpath(i);
                e.v1 = subpath(i + 1);

                const BoundaryGraph::EdgeID e_id = BoundaryGraph::make_edge_id(e.v0, e.v1);
                vec2 d_error;
                if (accuracy_map.find(e_id) == accuracy_map.end()) {
                    d_error = GeomRaster::distance_bounds_from_points_with_slack(B, B.col(e.v0), B.col(e.v1), e.v0, e.v1);
                } else {
                    d_error = accuracy_map.at(e_id);
                }

                e.d_min = d_error.minCoeff();
                e.d_max = d_error.maxCoeff();
                PF_VERBOSE_F("Adding %d %d", e.v0, e.v1);
                E.emplace_back(e);
            }
        }
    }

    return (int)E_delete.size();
}

// Some of the symmetry criteria are "duplicated" from the construct_nodes, but it is necessary to do so
// as a node could be locally symmetric, but when connected to the next node end up invalidating a local symmetry
bool trace(const mat2x& B,  // fitting path
	const std::vector<Node>& N, // fitting path nodes (unnecessary, but find_min_valence uses them)
	const std::vector<ShortestPath::Node>& GN, // construct_**_trace_graph
	std::vector<Vertex>& P, bool circular) {

	const Vertex vmin = find_min_valence(B, N);
	PF_LOGF("minimum valence vertex %lld", vmin);

	std::vector<Vertex> Nmin;
	for (size_t i = 0; i < N.size(); ++i) {
		if (N[i].v1 == vmin) {
			Nmin.emplace_back(i);
		}
	}

	std::vector<Vertex> PG;
	ShortestPath::State s;

	if (circular) {
		ShortestPath::find_cycle(s, GN, Nmin, PG);
	} else {
		ShortestPath::find(s, GN, Nmin[0], Nmin[0], PG);
	}

	if (PG.empty()) {
		PF_LOGS("no path found");
		return false;
	}

	PF_LOGF("found cycle cost %f", ShortestPath::cost(GN, PG, circular));

	for (int i = 0; i < PG.size(); ++i) {
		P.emplace_back(N[GN[PG[i]].v].v1);
	}

	return true;
}

void trace_to_points(const mat2x& B, const std::vector<Vertex>& V, mat2x& P) {
	P.resize(2, V.size());

	for (int i = 0; i < V.size(); ++i) {
		P.col(i) = B.col(V[i]);
	}
}

void trace_to_points(const mat2x& B, const vecXi& PV, mat2x& P) {
	P.resize(2, PV.size());

	for (int i = 0; i < PV.size(); ++i) {
		P.col(i) = B.col(PV(i));
	}
}
		
bool trace_symmetric(const mat2x& R, 
	std::vector<BoundaryGraph::Edge>& E, 
	std::vector<BoundaryGraph::Node>& N,
	const std::vector<Vertex>& M, 
	const std::vector<ShortestPath::Node>& GN, 
	std::vector<Vertex>& V, 
	bool circular, 
	OnTrace callback) {
	// Obtaining a first asymmetric initial solutio
	trace(R, N, GN, V, circular);

	// The boundary points are disconnected in the graph
	if (V.empty()) {
		return false;
	}

	int iter = 0;
	if (callback) {
		callback(R, V, E, Symmetry::SubPath(), iter++, true);
	}

	// Finding all the global symmetries
	std::vector<Symmetry::SubPath> S = Symmetry::find_longest(R, circular);
	PF_VERBOSE_F("found %lld symmetries", S.size());

	std::vector<int> C;
	PathUtils::compute_convexities(R, C);

	// Sorting them convex to concave and then shortest to longest
	std::sort(S.begin(), S.end(), [&](const Symmetry::SubPath & lhs, const Symmetry::SubPath & rhs) -> bool {
		if (C[lhs.v00] == +1 && C[rhs.v00] == -1) {
			return true;
		}
		if (C[lhs.v00] == -1 && C[rhs.v00] == +1) {
			return false;
		}

		return CircularDist(R, lhs.v00, lhs.v10) < CircularDist(R, rhs.v00, rhs.v10);
	});

	// Holds the edges that are about to be deleted in the current iteration
	std::vector<size_t> E_dlist;
			
	// These hold the vertices of V in the symmetric region being processed
	std::vector<Vertex> valid0, valid1;

	// Holds all the symmetric points ordered with respect to the boundary and the symmetry
	std::vector<Vertex> valid_ord;

	// If a symmetry cannot be applied, we want to preserve the edges which
	// are incoming or outgoing from the symmetric region.
	std::vector<Edge> E_boundary;

	// holds the polygon path points for the current iteration
	mat2x P, Psym;
	trace_to_points(R, V, P);

	int V_inflections = AngleUtils::count_inflections(P, circular);
	PF_VERBOSE_F("inflections %d", V_inflections);

	// Temporary data required for graph construction before deciding whether the new symmetric
	// path is accepted as an alternative.
	// The copies are not necessary since the individual elements will be removed, but this is less
	// prone to errors and the memory required is bounded by the input since the intent is to
	// prune the graph of all the invalid edges.
	vector<Edge> Esym;
	vector<Node> Nsym;
	vector<ShortestPath::Node> GNsym;
	vector<Vertex> Vsym;
	vector<int> Csym; // holds convexities

	// If two paths are already symmetric, they are not explicitly skipped in order to prune the graph
	for (int i = 0; i < S.size(); ++i) {
		Symmetry::SubPath s = S[i];
		PF_VERBOSE_F("trace iteration %d --------", iter);

		// reordering the symmetry points consistently with the orientation of the boundary
		bool joint0 = Symmetry::is_subpath_joint0(R, s, circular);
		bool joint1 = Symmetry::is_subpath_joint1(R, s, circular);

		// making v00 and v10 the consecutive vertices
		if (joint1 && !joint0) {
			swap(s.v00, s.v01);
			swap(s.v10, s.v11);
		}

		// making the direction +1 from v11 -> v01
		if (CircularDist(R, s.v01, s.v00) > CircularDist(R, s.v00, s.v01)) {
			swap(s.v01, s.v11);
			swap(s.v00, s.v10);
		}

		PF_VERBOSE_F("symmetry v01 %lld v00 %lld v10 %lld v11 %lld axis %d", s.v01, s.v00, s.v10, s.v11, s.axis);

		// Finding the vertices of the polygon which belong to the symmetric region
		valid0.clear();
		valid1.clear();
		for (int j = 0; j < V.size(); ++j) {
			if (PathUtils::contains_closed(R, s.v01, s.v00, V[j], circular)) {
				valid0.emplace_back(V[j]);
			}

			if (PathUtils::contains_closed(R, s.v10, s.v11, V[j], circular)) {
				valid1.emplace_back(V[j]);
			}
		}

#if 0
		for (int j = 0; j < valid0.size(); ++j) {
			draw::point(R.col(valid0[j]), 5., Style::fill(colors::red));
		}

		for (int j = 0; j < valid1.size(); ++j) {
			draw::point(R.col(valid1[j]), 5., Style::fill(colors::forest_green));
		}
#endif
				
		PF_VERBOSE_F("valid0 %llu valid1 %llu", valid0.size(), valid1.size());

		// To decide which one of the two sides is more symmetric, look at the number of vertices
		// in both paths. It's a good indicator for smoothness or regularities (likely exclusive to one
		// another). 
		// todo: if it fail, then pick the subpath which has the least cost
		if (valid0.size() == valid1.size()) {
			continue;
		}

		bool pick0 = valid0.size() > valid1.size();
		if (pick0) {
			valid1.clear();
			valid1.resize(valid0.size());
		}
		else {
			valid0.clear();
			valid0.resize(valid1.size());
		}

		PF_VERBOSE_F("pick0 %d", pick0);

		// Finding the set of corners which are symmetry w.r.t. the good side of the symmetry
		for (int j = 0; j < valid0.size(); ++j) {
			if (pick0) {
				Vertex v0 = valid0[j];
				Vertex v1 = Circular(R, s.v10 + CircularDist(R, v0, s.v00));
				PF_LOGF("symmetric point v0 %lld -> v1 %lld", v0, v1);
				valid1[j] = v1;
			}
			else {
				Vertex v1 = valid1[j];
				Vertex v0 = Circular(R, s.v01 + CircularDist(R, v1, s.v11));
				PF_LOGF("symmetric point v1 %lld -> v0 %lld", v1, v0);
				valid0[j] = v0;
			}
		}

		PF_ASSERT(valid0.size() == valid1.size());
		for (int j = 0; j < valid0.size(); ++j) {
			const vec2 v0 = R.col(valid0[j]);
			const vec2 v1 = R.col(valid1[j]);
			//draw::line(v0, v1, Style::outline(colors::red, 1.));
		}

		// having the points in valid0 and valid1 unordered, makes checking for consecutive
		// symmetric edges more complex. Putting the points into a single array and ordering 
		// them, this way a consecutive edge is expected to connect two consecutive points
		// in the ordered array. This is transparent to disjoint symmetries.
		valid_ord.clear();
		for (int j = 0; j < valid0.size(); ++j) {
			valid_ord.emplace_back(valid0[j]);
		}

		// adding the other points being without introducing duplicates in the case of joint symmetries
		for (int j = 0; j < valid1.size(); ++j) {
			// todo: can be replaced with a bunch of ifs
			if (find(valid_ord.begin(), valid_ord.end(), valid1[j]) == valid_ord.end()) {
				valid_ord.emplace_back(valid1[j]);
			}
		}

		for (int j = 0; j < valid_ord.size(); ++j) {
			PF_VERBOSE_F("valid unord %d", valid_ord[j]);
		}

		PathUtils::order(R, valid_ord, s.v01);
		PF_VERBOSE_F("ordered %lld symmetric points from %lld", valid_ord.size(), s.v01);
		for (int j = 0; j < valid_ord.size(); ++j) {
			PF_VERBOSE_F("valid ord %d", valid_ord[j]);
		}

		// Operating on the graph representation used for the shortest path is messy and
		// perfectly symmetric nodes might not exist as the raster boundary outside
		// the symmetric region is different on the two sides (otherwise it would be symmetric)
		// We first attempt to find a perfectly symmetric solution by removing all the 
		// asymmetric edges from the graph. 
		// If the shortest path fails or introduces a bad inflection, the requirement is relaxed
		// and all the edges which cross the boundaries of the symmetric region are reinserted in 
		// the graph. The edges are required to connect to the corner in the symmetric region.
		// If this also fails, then the symmetry is discarded and the graph restored to the original
		// form.
		// List of edges which are to be deleted in the current iteration
		E_dlist.clear();

		// List of edges which cross the boundary of the symmetric region
		E_boundary.clear();

		for (int j = 0; j < E.size(); ++j) {
			const auto& e = E[j];

			// Testing each corner to see if it is inside any of the two symmetric regions
			const bool contain00 = PathUtils::contains_closed(R, s.v01, s.v00, e.v0, circular);
			const bool contain01 = PathUtils::contains_closed(R, s.v01, s.v00, e.v1, circular);
			const bool contain10 = PathUtils::contains_closed(R, s.v10, s.v11, e.v0, circular);
			const bool contain11 = PathUtils::contains_closed(R, s.v10, s.v11, e.v1, circular);

			// If the edge is outside the symmetric region, it is left untouched
			if (!contain00 && !contain01 && !contain10 && !contain11) {
				continue;
			}

			// If only one of the edge endpoints is contained in one of the symmetric region
			// it is considered a boundary edge.
			const bool is_boundary = (contain01 != contain00) || (contain10 != contain11);

			// If the point is connecting two valid points in opposite symmetric regions, the criteria
			// for removal should follow that of regular edges
			const bool is_shared = (contain00 || contain01) && (contain10 || contain11);

			if (is_boundary && !is_shared) {
				// we are interested in preserving only the boundary vertices which connect to the first/last point 
				// in each symmetric region. This means that there should be no other point belonging to the symmetric
				// region in the open interval defined by the edge.
				// This has to be tested separately for the two symmetric regions, otherwise it would fail
				// in the presence of disjoint symmetries.
				bool discard = false;

				if (contain01 != contain00) {
					for (size_t k = 0; k < valid0.size(); ++k) {
						if (PathUtils::contains_open(R, e.v0, e.v1, valid0[k], circular)) {
							discard = true;
							break;
						}
					}
				}

				if (contain10 != contain11) {
					for (size_t k = 0; k < valid1.size(); ++k) {
						if (PathUtils::contains_open(R, e.v0, e.v1, valid1[k], circular)) {
							discard = true;
							break;
						}
					}
				}

				if (discard) {
					E_dlist.emplace_back(j);
					E_boundary.emplace_back(e);

					continue;
				}

#if 0
				draw::line(R.col(e.v0), R.col(e.v1), Style::outline(colors::forest_green, 2.));
#endif
			}
			else {
				// the symmetric edge is kept if it connects two consecutive vertices, otherwise it is marked for deletion
				// is the source endpoint part of the symmetric path?
				const auto& src_valid_it = find(valid_ord.begin(), valid_ord.end(), e.v0);

				if (src_valid_it == valid_ord.end() ) {
					E_dlist.emplace_back(j);
					continue;
				}

				// todo: circularity should be asserted more strictly ?
				// is the destination endpoint the next expected symmetric point
				const Vertex isrc = distance(valid_ord.begin(), src_valid_it);
				const bool is_dst_valid = CircularAt(valid_ord, isrc + 1) == e.v1;

				if (!is_dst_valid) {
					E_dlist.emplace_back(j);
					continue;
				}

#if 0
				draw::line(R.col(e.v0), R.col(e.v1), Style::outline(colors::talking_orange, 5.));
#endif
			}
		}

#if 0
		for (int j = 0; j < E_dlist.size(); ++j) {
			const vec2 p0 = R.col(E[E_dlist[j]].v0);
			const vec2 p1 = R.col(E[E_dlist[j]].v1);
			draw::line(p0, p1, Style::outline(colors::red, 1.));
		}

#endif

		//
		// 1st attempt: remove all the invalid edges and find a new polygon
		//
		PF_VERBOSE_F("refitting polygon with %lld edges", E.size() - E_dlist.size());
		Esym = E;
		EraseOrdered(Esym, E_dlist);
		BoundaryGraph::construct_nodes(R, Esym, M, Nsym, circular);
		BoundaryGraph::construct_trace_graph(R, Esym, M, Nsym, GNsym, circular);
		BoundaryGraph::trace(R, Nsym, GNsym, Vsym, circular);
				
		// How many times does this really happen?
		bool fit_failed = false;

		if (Vsym.empty()) {
			PF_VERBOSE_S("no path found");
			fit_failed = true;
		}
		else {
			// Has this path introduced any inflections?
			trace_to_points(R, Vsym, Psym);
			const int Vsym_inflections = AngleUtils::count_inflections(Psym, circular);

			if (Vsym_inflections > V_inflections) {
				PF_VERBOSE_S("bad symmetric fit, it would introduce an inflection");
				fit_failed = true;
			}

			fit_failed = false;
		}

		// good, applying the changes
		if (!fit_failed) {
			PF_VERBOSE_S("symmetric fit found, applying changes");
			E = std::move(Esym);
			N = std::move(Nsym);
			V = std::move(Vsym);

			if (callback) {
				callback(R, V, E, s, iter, true);
				++iter;
			}

			continue;
		}

		//
		// 2nd attempt: inserting boundary edges and fitting a new polygon
		//
		PF_VERBOSE_F("inserting %lld boundary edges", E_boundary.size());
		PF_VERBOSE_F("refitting polygon with %lld edges", E.size() - E_dlist.size() + E_boundary.size());
	}

	return true;
}

Vertex find_min_valence(const mat2x& B, const std::vector<Node>& N) {
	std::vector<int> valences(B.cols(), 0);

	for (size_t i = 0; i < N.size(); ++i) {
		++valences[N[i].v2];
	}

	Vertex vmin = -1;
	int vcount = std::numeric_limits<int>::max();
	for (Vertex i = 0; i < B.cols(); ++i) {
		//PF_LOGF("node %lld valence %d max %d", i, valences[i], vcount);
		if (valences[i] < vcount && valences[i] > 0) {
			vcount = valences[i];
			vmin = i;
		}
	}

	// No nodes
	assert_break(vmin != -1);
	return vmin;
}

double approx_edge_cost(
    const mat2x& B,
    const std::vector<bool>& M,
    const std::vector<Edge>& E,
    const Node& n0, 
    const Node& n1, 
    bool circular,
    const bool is_n0_C0,
    const bool is_n1_C0) {
	assert_break(n0.e1 == n1.e0);
	const vec2 p0 = B.col(n0.v0);
	const vec2 p1 = B.col(n0.v1);
	const vec2 p2 = B.col(n1.v1);
	const vec2 p3 = B.col(n1.v2);

	const int flags = AngleUtils::have_opposite_convexity(p0, p1, p2, p3) ? ErrorMetrics::EDGE_HAS_INFLECTION : 0x0;
	double e_accuracy, e_continuity, e_smoothness;
			
	if (VectorOptions::get()->accuracy_metric == AccuracyMetric::OneSided) {
		e_accuracy = ErrorMetrics::accuracy(E[n0.e1].d_min, E[n0.e1].d_max, flags);
	}
	else if (VectorOptions::get()->accuracy_metric == AccuracyMetric::Implicit) { 
		e_accuracy = ErrorMetrics::accuracy_implicit(E[n0.e1].d_min, E[n0.e1].d_max, flags);
	}else {
		PF_ABORT;
	}

	e_continuity = ErrorMetrics::continuity(n0.angle, n1.angle, flags);
	e_smoothness = ErrorMetrics::smoothness(n0.angle, n1.angle, flags);
	
    const bool is_short_edge = CircularDist(B, n0.v1, n1.v1) <= 2 && M[Circular(B, n0.v1 + 1)];
	double e_penalty = is_short_edge ? Options::get()->short_edge_penalty : 0.;
    e_penalty += Options::get()->edge_penalty;

	if (!is_n0_C0 && !is_n1_C0) {
		return e_penalty + e_accuracy + e_continuity + e_smoothness;
	}
	else if (is_n0_C0 && is_n1_C0) {
		return e_penalty + e_accuracy;
	}
	else if (is_n0_C0) {
		return e_penalty + e_accuracy + .5 * ErrorMetrics::smoothness(n1.angle);
	}
	else if (is_n1_C0) {
		return e_penalty + e_accuracy + .5 * ErrorMetrics::smoothness(n0.angle);
	}
	PF_ABORT;
}

double approximate_edge_cost (
    const vec2 p0,
    const vec2 p1,
    const vec2 p2,
    const vec2 p3,
    const double d_min,
    const double d_max
) {
    const int flags = AngleUtils::have_opposite_convexity(p0, p1, p2, p3) ? ErrorMetrics::EDGE_HAS_INFLECTION : 0x0;
    const double angle0 = AngleUtils::spanned_shortest(p0, p1, p2);
    const double angle1 = AngleUtils::spanned_shortest(p1, p2, p3);
    const double e_accuracy = ErrorMetrics::accuracy_implicit(d_min, d_max, flags);
    const double e_continuity = ErrorMetrics::continuity(angle0, angle1, flags);
    const double e_smoothness = ErrorMetrics::smoothness(angle0, angle1, flags);
    return e_accuracy + e_continuity + e_smoothness;
}

bool is_valid_midpoint_node(const mat2x& P, const std::vector<int>& C, const Vertex v0, const Vertex v1, const Vertex v2) {
	// if v1 is a midpoint, it is assumed to be flat
	PF_ASSERT(C[v1] == 0);

	const Vertex v1p = PathUtils::next_transition_vertex(C, v1, -1);
	const Vertex v1n = PathUtils::next_transition_vertex(C, v1, +1);
	if (PathUtils::next_transition_vertex(C, v1p, -1) == v0 && PathUtils::next_transition_vertex(C, v1n, +1) == v2) {
		return false;
	}

	const double l01 = (P.col(v0) - P.col(v1)).norm();
	const double l12 = (P.col(v1) - P.col(v2)).norm();
	if (abs(l01 - l12) > PF_EPS) {
		return false;
	}
			
	return true;
}

bool trace_to_edges(const std::vector<Edge>& E, const std::vector<Vertex>& V, std::vector<Vertex>& VE, bool circular) {
	const size_t ioff = circular ? 0 : 1;

	VE.clear();

	for (size_t i = 0; i < V.size() - ioff; ++i) {
		const Vertex v0 = V[i];
		const Vertex v1 = CircularAt(V, i + 1);
		const auto e_it = find_if(E.begin(), E.end(), [v0, v1](const Edge& e) {
			return e.v0 == v0 && e.v1 == v1;
		});

		if (e_it == E.end()) {
			PF_LOGF("failed to find edge %lld -> %lld", v0, v1);
			return false;
		}

		VE.emplace_back((Vertex)distance(E.begin(), e_it));
	}

	PF_ASSERT(VE.size() == V.size() - ioff);

	return true;
}

bool trace_to_edges(const std::vector<Edge>& E, const vecXi& V, std::vector<Vertex>& VE, bool circular) {
	const size_t ioff = circular ? 0 : 1;

	VE.clear();

	for (size_t i = 0; i < V.rows() - ioff; ++i) {
		const Vertex v0 = V(i);
		const Vertex v1 = V((i + 1) % V.rows());
		const auto e_it = find_if(E.begin(), E.end(), [v0, v1](const Edge & e) {
			return e.v0 == v0 && e.v1 == v1;
			});

		if (e_it == E.end()) {
			PF_LOGF("failed to find edge %lld -> %lld", v0, v1);
			return false;
		}

		VE.emplace_back((Vertex)distance(E.begin(), e_it));
	}

	PF_ASSERT(VE.size() == V.rows() - ioff);

	return true;
}

double angle_between_edges(
	const mat2x& B,
	const BoundaryGraph::Edge& e_prev,
	const BoundaryGraph::Edge& e_next
) {
	PF_ASSERT(e_prev.v1 == e_next.v0);

	const vec2 p0 = B.col(e_prev.v0);
	const vec2 p1 = B.col(e_prev.v1);
	const vec2 p2 = B.col(e_next.v1);

	return AngleUtils::spanned_shortest(p0, p1, p2);
}

NAMESPACE_END(BoundaryGraph)
NAMESPACE_END(polyfit)