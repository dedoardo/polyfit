#include <polyvec/shortest-path/dijkstra.hpp>

#include <polyvec/api.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/utils/matrix.hpp>

#include <cstdlib>
#include <algorithm>

#define DEBUG_CYCLE 0
#define DEBUG_TRAVERSE 0
#define DEBUG_DFS 0
#define DEBUG_BFS 0

#if DEBUG_TRAVERSE
#define TRAVERSE_LOG(fmt, ...) PF_LOGF(fmt, __VA_ARGS__)
#else
#define TRAVERSE_LOG(fmt, ...)
#endif

#if DEBUG_CYCLE
#define CYCLE_LOG(fmt, ...) PF_LOGF(fmt, __VA_ARGS__)
#else
#define CYCLE_LOG(fmt, ...)
#endif

#if DEBUG_DFS
#define DFS_LOG(fmt, ...) PF_LOGF(fmt, __VA_ARGS__)
#else
#define DFS_LOG(fmt, ...)
#endif

#if DEBUG_BFS
#define BFS_LOG(fmt, ...) PF_LOGF(fmt, __VA_ARGS__)
#else
#define BFS_LOG(fmt, ...)
#endif

using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(ShortestPath)

void State::reset() {
    queue = PriorityQueue();
	visits.clear();
	visited.clear();
}

double cost(const std::vector<Node>& G, const std::vector<Vertex>& P, bool circular) {
	double tot = 0.;
	for (size_t i = 0; i < P.size(); ++i) {
		const Node& n = G[P[i]];

		if (!circular && i == P.size() - 1) {
			break;
		}

		TRAVERSE_LOG("node %lld %s", P[i], n.id.c_str());

		int j = 0;
		for (; j < n.neighbors.size(); ++j) {
			TRAVERSE_LOG("nhb %lld expected %lld", n.nhbs[j], CircularAt(P, i + 1));

			auto& edge = n.neighbors[j];
			if (edge.to == CircularAt(P, i + 1)) {
				tot += edge.cost;
				break;
			}
		}

		assert_break(j < n.neighbors.size());
	}

	return tot;
}

bool find(State& s, const std::vector<Node>& G, Vertex src, Vertex dst, std::vector<Vertex>& P) {
	P.clear();

	if (G[src].neighbors.empty()) {
		return false;
	}

	// Clearing traverse structures
	s.reset();
	s.visits.resize(G.size(), State::Visit());
	s.visited.resize(G.size(), false);

	// circular source/destination
	State::Visit visit_dst;

	// initializing source node
	s.visits[src].cost = 0.;
	s.visits[src].prev = -1;
	s.queue.push(PriorityQueueEl(s.visits[src].cost, src));

	Vertex node_dst = -1;
	const bool circular = src == dst;
	unsigned int iters = 0;

	while (!s.queue.empty()) {
		PriorityQueueEl visit_el = s.queue.top();
		assert_break(visit_el.first >= 0.); // This should just not happen

		// next element
		const double   node_at_cost = visit_el.first;
		const Vertex inode_at = visit_el.second;
		const Node & node_at = G[inode_at];
		const State::Visit& node_at_visit_info = s.visits[inode_at];

		TRAVERSE_LOG("visiting %s", node_at.id.c_str());

		// removing it from the queue
		s.queue.pop();

		// Termination condition
		if (inode_at == dst) {
			// If it's circular and the dst_visit has never been touched then this is the first iteration
			if (!circular || (circular && visit_dst.prev != -1)) {
				TRAVERSE_LOG("Terminating circular %d visit-dst-cost %f visit-dst-prev %lld", circular, visit_dst.cost, visit_dst.prev);
				break;
			}
		}

		// Marking current as visited
		s.visited[inode_at] = true;

		TRAVERSE_LOG("node %s neighbors %d", node_at.id.c_str(), node_at.nhbs.size());

		// Neighboring nodes outgoing from the current
		for(auto& edge : node_at.neighbors) {
			const Vertex inode_to = edge.to;
			const Node& node_to = G[inode_to];

			double edge_cost = edge.cost;
			assert_break(edge_cost >= 0.);
			const double cost_so_far = node_at_visit_info.cost + edge_cost;

            //if (!circular) {
                //PF_DEV_F("%s -> %s cost %f", node_at.id.c_str(), node_to.id.c_str(), edge_cost);
            //}

			TRAVERSE_LOG("neighbor %s cost-to-here %f edge-cost %f cost-so-far %f",
				node_to.id.c_str(), node_at_visit_info.cost, edge_cost, cost_so_far);

			State::Visit visit_info_to;
			visit_info_to.cost = cost_so_far;
			visit_info_to.prev = inode_at;

			// updating cost if less than what there is right node
			State::Visit* visit = nullptr;

			if (inode_to == dst && circular) {
				visit = &visit_dst;
			} else {
				visit = &s.visits[inode_to];
			}

			if (cost_so_far < visit->cost) {
				*visit = visit_info_to;

				if (!s.visited[inode_to]) {
					s.queue.push(PriorityQueueEl(cost_so_far, inode_to));
				}
				else {
					TRAVERSE_LOG("node %s already visited", node_to.id.c_str());
				}
			}
		}
	}

	// did we reach the destination ?
	State::Visit* visit = circular ? &visit_dst : &s.visits[dst];
	if (visit->prev == -1) {
		TRAVERSE_LOG("shortest path failed, destination %s not reached", G[dst].id.c_str());
		return false;
	}
			
	TRAVERSE_LOG("destination prev %s", G[visit->prev].id.c_str());

	TRAVERSE_LOG("destination reached, reconstructing path");

	if (!circular) {
		P.emplace_back(dst);
	}

	P.emplace_back(visit->prev);

	Vertex inode_next = visit->prev;
	while (inode_next != -1) {
		P.emplace_back(s.visits[inode_next].prev);
		inode_next = s.visits[inode_next].prev;
	}
			
	if (P.back() == -1) {
		P.pop_back();
	}

	std::reverse(P.begin(), P.end());
	return true;
}

bool find(State& s, const std::vector<Node>& G, const std::vector<Vertex>& src, const std::vector<Vertex>& dst, std::vector<Vertex>& P) {
	if (src.empty() || dst.empty()) {
		TRAVERSE_LOG("no source or destination nodes");
		return false;
	}

	P.clear();

	// Clearing traverse structures
	s.reset();
	s.visits.resize(G.size(), State::Visit());
	s.visited.resize(G.size(), false);

	// circular source/destination
	std::vector<State::Visit> visit_dst(dst.size());

	// initializing source nodes
	for (int i = 0; i < src.size(); ++i) {
		s.visits[src[i]].cost = 0.;
		s.visits[src[i]].prev = -1;
		s.queue.push(PriorityQueueEl(s.visits[src[i]].cost, src[i]));
	}

	Vertex node_dst = -1;
	const bool circular = src == dst;
	unsigned int iters = 0;

	while (!s.queue.empty()) {
		PriorityQueueEl visit_el = s.queue.top();
		assert_break(visit_el.first >= 0.); // This should just not happen

		// next element
		const double   node_at_cost = visit_el.first;
		const Vertex inode_at = visit_el.second;
		const Node & node_at = G[inode_at];
		const State::Visit node_at_visit_info = s.visits[inode_at];

		TRAVERSE_LOG("visiting %s", node_at.id.c_str());

		// removing it from the queue
		s.queue.pop();

		// Termination condition
		bool dst_reached = false;
		for (int i = 0; i < dst.size(); ++i) {
			if (inode_at == dst[i]) {
				if (!circular || visit_dst[i].prev != -1) {
					TRAVERSE_LOG("terminating path circular %d visit-dst-cost %f visit-dst-prev %lld", circular, visit_dst[i].cost, visit_dst[i].prev);
					dst_reached = true;
					break;
				}
			}
		}

		if (dst_reached) {
			break;
		}

		// Marking current as visited
		s.visited[inode_at] = true;

		TRAVERSE_LOG("node %s neighbors %d", node_at.id.c_str(), node_at.nhbs.size());

		// Neighboring nodes outgoing from the current
		for(auto& edge : node_at.neighbors) {
			const Vertex inode_to = edge.to;
			const Node& node_to = G[inode_to];

			double edge_cost = edge.cost;
			assert_break(edge_cost >= 0.);
			const double cost_so_far = node_at_visit_info.cost + edge_cost;

			TRAVERSE_LOG("neighbor %s cost-to-here %f edge-cost %f cost-so-far %f",
				node_to.id.c_str(), node_at_visit_info.cost, edge_cost, cost_so_far);

			State::Visit visit_info_to;
			visit_info_to.cost = cost_so_far;
			visit_info_to.prev = inode_at;

			// updating cost if less than what there is right node
			State::Visit * visit = nullptr;
			for (int j = 0; j < dst.size(); ++j) {
				if (circular && inode_to == dst[j]) {
					visit = &visit_dst[j];
					break;
				}
				else {
					visit = &s.visits[inode_to];
					break;
				}
			}
					
			if (cost_so_far < visit->cost) {
				*visit = visit_info_to;

				if (!s.visited[inode_to]) {
					s.queue.push(PriorityQueueEl(cost_so_far, inode_to));
				}
				else {
					TRAVERSE_LOG("node %s already visited", node_to.id.c_str());
				}
			}
		}
	}

	// did we reach the destination ?
	// picking the destination node with the least cost
	State::Visit* visit = nullptr;
	for (int i = 0; i < dst.size(); ++i) {
		if (circular && visit_dst[i].prev != -1) {
			if (!visit || visit->cost > visit_dst[i].cost) {
				visit = &visit_dst[i];
			}
		}
		else if (s.visits[dst[i]].prev != -1) {
			if (!visit || visit->cost > s.visits[dst[i]].cost) {
				visit = &s.visits[dst[i]];
			}

			P.emplace_back(dst[i]);
		}
	}

	if (!visit) {
		TRAVERSE_LOG("no shortest path found");
		return false;
	}

	TRAVERSE_LOG("destination prev %s", G[visit->prev].id.c_str());
	TRAVERSE_LOG("destination reached, reconstructing path");

	P.emplace_back(visit->prev);

	Vertex inode_next = visit->prev;
	while (inode_next != -1) {
		P.emplace_back(s.visits[inode_next].prev);
		inode_next = s.visits[inode_next].prev;
	}

	if (P.back() == -1) {
		P.pop_back();
	}

	std::reverse(P.begin(), P.end());
	return true;
}
		
bool find_cycle(State& s, const std::vector<Node>& G, const std::vector<Vertex>& srcdst, std::vector<Vertex>& P) {
	CYCLE_LOG("find-cycle source nodes");
	for (int i = 0; i < srcdst.size(); ++i) {
		CYCLE_LOG("\tnode %lld id %s", srcdst[i], G[srcdst[i]].id.c_str());
	}

	// First, finding the cycle from srcdst which costs the least
	std::vector<Vertex> GP0, GP1;
	double GP0_cost = INFINITY;

	for (int i = 0; i < srcdst.size(); ++i) {
		const Vertex v = srcdst[i];

		if (::polyfit::ShortestPath::find(s, G, { v }, { v }, GP1)) {
			const double GP1_cost = cost(G, GP1);

			if (GP1_cost < GP0_cost) {
				GP0 = std::move(GP1);
				GP0_cost = GP1_cost;
				CYCLE_LOG("found a better path from node %lld id %s", v, G[v].id.c_str());
			}
		}
	}

	if (GP0.empty()) {
		return false;
	}
	
	CYCLE_LOG("found initial cycle length %lld", GP0.size());
	double e_tot = 0.;
	for (int i = 0; i < GP0.size(); ++i) {
		CYCLE_LOG("node %lld id %s", GP0[i], G[GP0[i]].id.c_str());

		auto ni = G[GP0[i]];

		for(auto& edge : ni.neighbors) {
			if (CircularAt(GP0, i + 1) == edge.to) {
				CYCLE_LOG("cost to next node %f", edge.cost);
				e_tot += edge.cost;
			}
		}
	}

	CYCLE_LOG("total path cost %f", e_tot);

	// not trying to find the cycle from the most distant node
	const int iters_min = 3;
	const int iters_max = 10;

	Vertex v = GP0[GP0.size() / 2];
	for (int iter = 0; iter < iters_max; ++iter) {
		CYCLE_LOG("find shortest cycle from node %lld id %s", v, G[v].id.c_str());
		assert_break(::polyfit::ShortestPath::find(s, G, v, v, GP1));
		const double GP1_cost = cost(G, GP1);
				
		CYCLE_LOG("found path cost %f length %llu", GP1_cost, GP1.size());
		for (int j = 0; j < GP1.size(); ++j) {
			CYCLE_LOG("node %lld id %s", GP1[j], G[GP1[j]].id.c_str());
		}

		if (GP1_cost < GP0_cost) {
			CYCLE_LOG("found a better path length %lld", GP1.size());
			CYCLE_LOG("old cost %f new cost %f ", GP0_cost, GP1_cost);

			for (int j = 0; j < GP1.size(); ++j) {
				CYCLE_LOG("node %lld id %s", GP1[j], G[GP1[j]].id.c_str());
			}

			v = GP1[GP1.size() / 2];
			GP0 = GP1;
			GP0_cost = GP1_cost;
		}
		// let's try some other random node in the path
		else {
			size_t v_rnd = (size_t)(rand() % GP0.size());
			CYCLE_LOG("attempting another random node %lld", v_rnd);
			v = GP0[v_rnd];
		}
	}

    // This is not necessary, but allows to diff the files with the tracing results without any knowledge of the content
    std::sort(GP0.begin(), GP0.end());

	P = std::move(GP0);
	return true;
}

bool dfs_traverse_rec(const std::vector<Node>& G, std::vector<bool>& visited, Vertex v, Vertex vdst, std::vector<Vertex>* P) {
	if (P) {
		P->emplace_back(v);
	}

	if (visited[v]) {
		return false;
	}

	visited[v] = true;

	if (v == vdst) {
		return true;
	}

	for(auto& edge : G[v].neighbors) {
		DFS_LOG("dfs %lld(%lld) nhb %lld(%lld)", v, G[v].v, G[v].nhbs[i], G[G[v].nhbs[i]].v);
		if (dfs_traverse_rec(G, visited, edge.to, vdst, P) && !visited[edge.to]) {
			return true;
		}
	}

	return false;
}

bool find_depth_first(const std::vector<Node>& G, Vertex src, Vertex dst, std::vector<Vertex>* P) {
	PF_ASSERT(src >= 0 && src < (Vertex)G.size());
	PF_ASSERT(dst >= 0 && dst < (Vertex)G.size());

	if (src == dst) {
		return true;
	}

	std::vector<bool> visited(G.size(), false);

	if (P) {
		P->clear();
	}

	DFS_LOG("dfs source %lld(%lld) destination %lld(%lld)", src, G[src].v, dst, G[dst].v);
	return dfs_traverse_rec(G, visited, src, dst, P);
}

bool find_breadth_first(const std::vector<Node>& G, Vertex src, Vertex dst) {
	PF_ASSERT(src >= 0 && src < (Vertex)G.size());
	PF_ASSERT(dst >= 0 && dst < (Vertex)G.size());

	BFS_LOG("bfs source %lld(%lld) destination %lld(%lld)", src, G[src].v, dst, G[dst].v);

	if (src == dst) {
		return true;
	}

	std::vector<bool> visited(G.size(), false);

	std::queue<Vertex> queue;
	queue.push(src);

	while (!queue.empty()) {
		Vertex v = queue.front();
		queue.pop();
				
		if (visited[v]) {
			continue;
		}

		BFS_LOG("bfs visit %lld(%lld)", v, G[v].v);

		visited[v] = true;

		if (v == dst) {
			return true;
		}

		for(auto& edge : G[v].neighbors) {
			if (!visited[edge.to]) {
				queue.push(edge.to);
			}
		}
	}

	return false;
}

Vertex find_or_add_node(std::vector<Node>& G, const Vertex v) {
	const auto& n_it = std::find_if(G.begin(), G.end(), [v](const Node & n) {
		return n.v == v;
	});
			
	if (n_it == G.end()) {
		G.emplace_back(v);
		return G.size() - 1;
	}
	else {
		return distance(G.begin(), n_it);
	}
}

void print_graph(const std::vector<Node>& G, FILE *file) {
int id = 0;
for (const Node & node :G ) {
	fprintf(file,"# node id \n%d\n", id);
	fprintf(file,"# node v \n%d\n", (int)node.v);
	fprintf(file,"# node n_neighbours \n%d\n", (int)node.neighbors.size());
	fprintf(file,"# node neighbours \n");
	for (int i = 0 ; i < node.neighbors.size() ; ++i ) {
	fprintf(file,"%d ", (int)node.neighbors[i].to);			
	}
fprintf(file,"\n");
	fprintf(file,"# node n costs \n");
	for (int i = 0 ; i < node.neighbors.size() ; ++i ) {
	fprintf(file,"%.5g ", node.neighbors[i].cost);							
	}			 
fprintf(file,"\n # ------------ \n");
++id;
} 		
}

NAMESPACE_END(ShortestPath)
NAMESPACE_END(polyfit)
