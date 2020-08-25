#pragma once

// Polyvec
#include <polyvec/core/types.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/regularity/symmetry.hpp>

// libc++
#include <unordered_map>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(BoundaryGraph)

struct Edge {
	Vertex v0 = -1;
	Vertex v1 = -1;
	double d_min = INFINITY;
	double d_max = INFINITY;
	int flags = 0x0; // ErrorMetrics::EDGE_*
    
    static Edge make_pixel_diagonal(const int v0, const int v1);
    static Edge make_flat(const int v0, const int v1);
};

using EdgeID = uint64_t;
EdgeID make_edge_id (const int c0, const int c1);
vec2i  unpack_edge_id(const EdgeID id);

using AccuracyMap = std::unordered_map<polyfit::BoundaryGraph::EdgeID, Eigen::Vector2d>;

struct Node {
	Vertex v0 = -1, v1 = -1, v2 = -1; // points -> R
	Vertex e0 = -1, e1 = -1;     // in/out edges -> E
	double angle = INFINITY;
};

// inserts midpoints in short flat symmetric edges
// R is a list of raster points
// B is a list of points which includes midpoints in short edges
// M is a list of indices to 'B' which indicates the midpoints
void create_fitting_path(
    const mat2x& R, 
    mat2x& B, 
    std::vector<Vertex>& M, 
    bool circular = true
);
		
// B and M are the fitting points obtained from `create_fitting_path`
// a list of edges which pass the accuracy threshold are passed to E
// inserts all connections (one direction) whose edges are withing the accuracy threshold.
// edges are inserted only between corners and midpoints
void connect_valid_edges(
    const mat2x& B, 
    const std::vector<Vertex>& M, 
    std::vector<Edge>& E, 
    bool circular = true,
    std::unordered_map<EdgeID, Eigen::Vector2d>* accuracy_map = nullptr, // (corner, corner) -> (d_min, d_max)
    const std::vector<bool>& _B_flat = std::vector<bool>(),
    const bool force_B_flat = false
);

// Loops through all symmetric regions in the image from convex to concave and then in order of magnitude
// Removes all the edges which don't traverse the region symmetrically as long as 
// no paths with bad flow are inserted.
// such paths contain axis-aligned (raster) edges of 1 pixel
void remove_asymmetric_edges(
    const mat2x& B, 
    const std::vector<Vertex>& M, 
    std::vector<Edge>& E, 
    bool circular = true
);

// Connects consecutive pairs of edges to form nodes, stores the angle in Node::angle and prunes
// flat nodes as long as there is a valid replacement.
// 
// if the node is flat, it is not inserted if the edge connecting the source and destination exists
// the node vertex will reference the vertex it is defined on (multiple nodes will reference the same)
// Returns the number of nodes 
int construct_nodes(
    const mat2x& B, 
    const std::vector<Edge>& E, 
    const std::vector<Vertex>& M, 
    std::vector<Node>& N, 
    bool circular = true,
    const std::vector<bool>& B_flat = std::vector<bool>()
);

// finds the shortest cycle (or path) in the graph defined by the nodes and edges N, E
// construct_trace_graph_* converts the edges and nodes obtained from the boundary a data structure
// that can be used to find the shortest-path or navigate
// Returns the number of nodes
int construct_trace_graph(
    const mat2x& B, 
    const std::vector<Edge>& E, 
    const std::vector<Vertex>& M, 
    const std::vector<Node>& N, 
    std::vector<ShortestPath::Node>& GN, 
    bool circular = true,
    const std::vector<bool>& corners_C0 = std::vector<bool>()
);

// Prunes the boundary graph so that the polygon will go through the specified corners.
// The subpath *will* exist regardless of the current state of the boundary graph, but
// a polygon is not guaranteed to exist.
// An accuracy_map needs to be provided which contains the accuracy values for the subpath edges only
// (it doesn't need to be the full accuracy map)
// If the subpath has only one vertex, no values needs to be provided in the accuracy map
// Returns the number of edges pruned
int force_graph_subpath(
    const mat2x& B,                 // Boundary fitting points
    std::vector<Edge>& E,           // Edges that will be pruned
    const vecXi& subpath,           // Subpath to be forced
    const AccuracyMap& accuracy_map // Accuracy for every edge in the subpath
);

// Finds the shortest cycle/path in the graph and returns the list of boundary points
bool trace(const mat2x& B, const std::vector<Node>& N, const std::vector<ShortestPath::Node>& GN, std::vector<Vertex>& P, bool circular = true);

// simply converts the list of indices into the list of respective points
void trace_to_points(const mat2x& B, const std::vector<Vertex>& V, mat2x& P);
void trace_to_points(const mat2x& B, const vecXi& PV, mat2x& P);

using OnTrace = std::function<void(const mat2x& R, const std::vector<Vertex>& V, const std::vector<Edge>& E, const Symmetry::SubPath& path,int iter, bool accepted)>;

// fits a polygonal approximation to the boundary B trying to preserve as many of the raster
// symmetries as possible without introducing additional inflections. The shortest path nodes GN
// have already been constructed in a previous call to construct_symmetric_trace_graph.
bool trace_symmetric(const mat2x& B,  // fitting path
	std::vector<Edge>& E,			  // valid edges (updated)
	std::vector<Node>& N,       // triples of corners obtained through construct_nodes (updated)
	const std::vector<Vertex>& M,     // midpoint vertices (flat)
	const std::vector<ShortestPath::Node>& GN, // nodes constructed through construct_symmetric_trace_graph
	std::vector<Vertex>& V,           // output polygonal path
	bool circular = true,             // true if the path has to be circular
	OnTrace callback = nullptr        // called every time a symmetry is successfully applied and the polygon updated
);

// Returns the *boundary* vertex (Node::v1) which has the minimum valence
Vertex find_min_valence(const mat2x& B, const std::vector<Node>& N);

// Returns true if the edge e connecting points in P doesn't exactly overlap the inner symmetry endpoints
bool is_edge_crossing_subpath_asymmetric(const mat2x& P, const Edge& e, const Symmetry::SubPath& s, bool circular = true);

bool is_edge_violating_axis_aligned_edge(
    const mat2x& P, 
    const std::vector<int>& C,
    const Edge& e, 
    const Vertex a0, 
    const Vertex a1, 
    const bool circular,
    const double P_bb_diagonal
);

// Returns the cost for an edge connecting two nodes with overlapping endpoints
double approx_edge_cost(
    const mat2x& B,
    const std::vector<bool>& M,
    const std::vector<Edge>& E,
	const Node& n0,
    const Node& n1,
    bool circular = true,
    bool is_n0_C0 = false,
    bool is_n1_C0 = false
);

// A version of approx_edge_cost which is usable without node instances
double approximate_edge_cost (
    const vec2 p0,
    const vec2 p1,
    const vec2 p2,
    const vec2 p3,
    const double e_min,
    const double e_max
);

// Given an array of convexities C for every point and assuming that v1 is a midpoint, it checks
// whether the edges connecting v0 -> v1 and v1 -> v2 are not crossing the |__| region connecting the midpoints
// It also returns false if the node is not symmetric
bool is_valid_midpoint_node(const mat2x& P, const std::vector<int>& C, const Vertex v0, const Vertex v1, const Vertex v2);

// FIlls VE with the list of edge indices given the list of polygon vertices
// Returns false if any of the edges in V is not contained in E
bool trace_to_edges(const std::vector<Edge>& E, const std::vector<Vertex>& V, std::vector<Vertex>& VE, bool circular = true);
bool trace_to_edges(const std::vector<Edge>& E, const vecXi& V, std::vector<Vertex>& VE, bool circular = true);

// Simply calculates the angle between the three points and aborts if the two edges are not neighbors
double angle_between_edges(
	const mat2x& B,                    // points used to construct the boundary graph
	const BoundaryGraph::Edge& e_prev, // previous edge
	const BoundaryGraph::Edge& e_next  // next edge
);

NAMESPACE_END(BoundaryGraph)
NAMESPACE_END(polyfit)