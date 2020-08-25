// polyvec
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/utils/matrix.hpp>

// If activated, performs a test to see if a cycle exists before finding
// the shortest cycle in a graph
#define PRE_CYCLE_CHECK 1

using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int minimum(const mat2x& R, vecXi& P, mat2x& B, std::vector<BoundaryGraph::Edge>& E, const bool circular) {
	// Constructs a fitting path given the inputs boundary, placing midpoints on symmetric
	// segments on length <= 2
	std::vector<Vertex> M; // midpoints
	BoundaryGraph::create_fitting_path(R, B, M, circular);

	// For each pair of corners, if the line connecting the two samples is a valid approximation
	// of the underlying boundary, then it will be a candidate edge for the polygonal approximation
	BoundaryGraph::connect_valid_edges(B, M, E, circular);

	ShortestPath::State s;
	return minimum(s, B, M, E, P, circular);
}

// Implementation of cycle detection from https://www.geeksforgeeks.org/detect-cycle-in-a-graph/
bool has_cycle_from(int v, std::vector<bool>& visited, std::vector<bool>& is_on_stack, const std::vector<ShortestPath::Node>& G)
{
	if (!visited[v])
	{
		visited[v] = true;
		is_on_stack[v] = true;

		for (auto& e : G.at(v).neighbors)
		{
			if (!visited[e.to] && has_cycle_from(e.to, visited, is_on_stack, G))
				return true;
			else if (is_on_stack[e.to])
				return true;
		}
	}
	is_on_stack[v] = false;
	return false;
}

bool has_cycle(const std::vector<ShortestPath::Node>& G)
{
	std::vector<bool> visited(G.size(), false);
	std::vector<bool> is_on_stack(G.size(), false);
	for (int i = 0; i < G.size(); ++i)
		if (has_cycle_from(i, visited, is_on_stack, G))
			return true;
	return false;
}

int minimum(ShortestPath::State& s, 
    const mat2x& B, 
    const std::vector<Vertex>& M, 
    const std::vector<BoundaryGraph::Edge>& E, 
    vecXi& P, 
    const bool circular,
    const std::vector<bool>& corners_C0,
    const std::vector<bool>& B_flat
) {
    PF_ASSERT(circular);

    // We cannot construct an hypergraph with only two nodes, but if they connect the same two points,
    // we are safe to return the subpath connecting the two points;
    if (E.size() == 2 && E[0].v0 == E[1].v1 && E[0].v1 == E[1].v0) {
        P = vecXi(2);
        P(0) = E[0].v0;
        P(1) = E[0].v1;
        return 1;
    }

	// Forms nodes by connecting all triples of corners whose edges exist in the graph
	std::vector<BoundaryGraph::Node> N;
	if (!BoundaryGraph::construct_nodes(B, E, M, N, circular, B_flat)) {
		return 0;
	}

	// Each node can at most share one of the two edges with another node as
	// the corners are strictly ordered. If two nodes have an overlapping edge
	// they are considered neighbors, forming a secondary graph.
	std::vector<ShortestPath::Node> G;

	if (!BoundaryGraph::construct_trace_graph(B, E, M, N, G, circular, corners_C0)) {
		return 0;
	}

#if PRE_CYCLE_CHECK
	if (!has_cycle(G))
		return 0;
#endif

	// Increasing the likelihood of finding the shortest cycle by setting as
	// source/destination the vertex with minimum valence. If there are no nodes or 
    // no path is found, just keep attempting nodes randomly.
    bool path_found = false;
    Vertex vmin = BoundaryGraph::find_min_valence(B, N);
    
    P.resize(0);
    int attempt = 0;
    int max_attempts = 10;

    // Pick any vertex which is crossed by one of the nodes
    // todo: This needs to be absolutely reworked
    std::vector<int> vertex_tries;
    std::vector<bool> vertex_visited(B.cols(), false);
    for (const auto& n : N) {
        if (vertex_visited[n.v1]) {
            continue;
        }

        vertex_tries.emplace_back(n.v1);
        vertex_visited[n.v1] = true;
    }
    std::sort(vertex_tries.begin(), vertex_tries.end());

    do {
        // Finding the secondary (shortest path) graph components which are passing
        // through this point. 
        std::vector<Vertex> Vsrc, Vdst;
        for (size_t i = 0; i < N.size(); ++i) {
            if (circular && N[i].v1 == vmin) {
                Vsrc.emplace_back(i);
            }

            if (!circular && N[i].v1 == 0) {
                Vsrc.emplace_back(i);
            }

                if (!circular && N[i].v1 == B.cols() - 1) {
                Vdst.emplace_back(i);
            }

            if (Vsrc.size() && Vdst.size()) {
                PF_ASSERT(Vsrc.back() != Vdst.back());
            }
        }

        if (Vsrc.empty()) {
            vmin = vertex_tries[attempt];
            continue;
        }
        
        PF_ASSERT(!Vsrc.empty());

        // shortest path in G
        std::vector<Vertex> GP;
        ShortestPath::find_cycle(s, G, Vsrc, GP);

        // converting the graph node indices into boundary corners
        P.resize(GP.size());
        for (int i = 0; i < GP.size(); ++i) {
            P(i) = N[G[GP[i]].v].v1;
        }
        
        path_found = P.size() > 0;
    } while (P.size() == 0 && ++attempt < max_attempts);

    return P.size();
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)