// polyvec
#include <polyvec/polygon-tracer/iterative-local.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/regularity/symmetry.hpp>
#include <polyvec/regularity/regularity_action.hpp>

// libc++
#include <cmath>
#include <algorithm>

#define LOG_PATH_ITERATIONS 0
#define SORT_BY_LENGTH_OF_SYMMETRIC_REGION 1

using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

//int iterative_local(const mat2x& R, 
//	vecXi& P, 
//	mat2x&raster, 
//	std::vector<BoundaryGraph::Edge>& E, 
//	const bool circular,
//	IterationCallback on_iter) {
//	// Constructs a fitting path given the inputs boundary, placing midpoints on symmetric
//	// segments on length <= 2
//	std::vector<Vertex> M; // midpoints
//	BoundaryGraph::create_fitting_path(R, raster, M, circular);
//
//	// For each pair of corners, if the line connecting the two samples is a valid approximation
//	// of the underlying boundary, then it will be a candidate edge for the polygonal approximation
//	BoundaryGraph::connect_valid_edges(raster, M, E, circular);
//
//	// Finding an initial solution
//	ShortestPath::State G_state;
//	if (!PolygonTracer::minimum(G_state, raster, M, E, P, circular)) {
//		return 0;
//	}
//
//	return 1 + iterative_local(G_state, raster, M, E, P, circular, on_iter);
//}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)