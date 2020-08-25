// polyvec
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int symmetric(const mat2x& R, vecXi& P, mat2x& B, std::vector<BoundaryGraph::Edge>& E, const bool circular) {
	// Constructs a fitting path given the inputs boundary, placing midpoints on symmetric
	// segments on length <= 2
	std::vector<Vertex> M; // midpoints
	BoundaryGraph::create_fitting_path(R, B, M, circular);

	// For each pair of corners, if the line connecting the two samples is a valid approximation
	// of the underlying boundary, then it will be a candidate edge for the polygonal approximation
	BoundaryGraph::connect_valid_edges(B, M, E, circular);

	// Pruning locally asymmetric edges
	BoundaryGraph::remove_asymmetric_edges(B, M, E, circular);

	// Forms nodes by connecting all triples of corners whose edges exist in the graph
	std::vector<BoundaryGraph::Node> N;
	if (!BoundaryGraph::construct_nodes(B, E, M, N, circular)) {
		return 0;
	}

	// Each node can at most share one of the two edges with another node as
	// the corners are strictly ordered. If two nodes have an overlapping edge
	// they are considered neighbors, forming a secondary graph.
	std::vector<ShortestPath::Node> G;

	if (!BoundaryGraph::construct_trace_graph(B, E, M, N, G, circular)) {
		return 0;
	}

	// Increasing the likelihood of finding the shortest cycle by setting as
	// source/destination the vertex with minimum valence
	const Vertex vmin = BoundaryGraph::find_min_valence(B, N);

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

	// First, we find a baseline solution with no regularities which highlights
	// fundamental axis-aligned edges.
	ShortestPath::State s;

	// shortest path in G
	std::vector<Vertex> GP;
	if (circular) {
		ShortestPath::find_cycle(s, G, Vsrc, GP);
	} else {
		//ShortestPath::find(s, G, Vsrc, Vdst, GP);
	}

	// converting the graph node indices into boundary corners
	P.resize(GP.size());
	for (int i = 0; i < (int)GP.size(); ++i) {
		P(i) = (int)(N[G[GP[i]].v].v1);
	}

	return (int)P.cols();
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)