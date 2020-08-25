#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

#include <queue>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(ShortestPath)

struct Node {

	struct Edge
	{
		Vertex to;
		double cost;

		Edge() { }
		Edge(Vertex to, double cost)
			: to(to), cost(cost)
		{ }
	};

	Node(Vertex v = -1) : v(v) { }

	Vertex v;
	std::vector<Edge> neighbors;

	void add_neighbor(Vertex v, double cost)
	{
		neighbors.emplace_back(v, cost);
	}

#ifdef PF_DEBUG
	std::string id;
#endif
};

using PriorityQueueEl = std::pair<double, Vertex>;
using PriorityQueue = std::priority_queue<PriorityQueueEl, std::vector<PriorityQueueEl >, std::greater<PriorityQueueEl>>;
		
struct State {
	struct Visit {
		Vertex prev = -1;
		double cost = 1e16;
	};

	PriorityQueue queue;
	std::vector<Visit> visits;
	std::vector<bool> visited;

	void reset();
};

// sums up the cost of all the edges between P[0] and P[-1]
double cost(const std::vector<Node>& G, const std::vector<Vertex>& P, bool circular = true);
		
// The minimum cost path using Dijkstra's algorithm is stored in P
// Returns true if it exists, false otherwise
bool find(State& s, const std::vector<Node>& G, Vertex src, Vertex dst, std::vector<Vertex>& P);
		
// Multiple sources and multiple destination
// ok, let's make it a different function as long as we are sure that it works
bool find(State& s, const std::vector<Node>& G, const std::vector<Vertex>& src, const std::vector<Vertex>& dst, std::vector<Vertex>& P);

// Attempts (reliably) to find the shortest cycle 
// In theory, to find the shortest cycle we would have to run the shortest path from every node to
// itself and then pick the one with the least cost. 
// To avoid this expensive operation a shortest path is first found starting from the node with the
// minimum valence. Then, a minimum number of attempts are made to find the shortest path starting from the
// most distant node in the current path and updates it if a least cost one. If such path does not exists
// the remaining attemps are made starting from a random node in the path.
// In practice, 3 seems to be a good number of attempts to overcome the starting point bias
bool find_cycle(State& s, const std::vector<Node>& G, const std::vector<Vertex>& srcdst, std::vector<Vertex>& P);

// Finds a path between src and dst traversing G in a depth-first fashion. The first path to reach dst is returned
bool find_depth_first(const std::vector<Node>& G, Vertex src, Vertex dst, std::vector<Vertex>* P);

bool find_breadth_first(const std::vector<Node>& G, Vertex src, Vertex dst);

// utilities to operate on vector<Node>
// ---
Vertex find_or_add_node(std::vector<Node>& G, const Vertex v);

void print_graph(const std::vector<Node>& G, FILE *file);

NAMESPACE_END(ShortestPath)
NAMESPACE_END(polyfit)