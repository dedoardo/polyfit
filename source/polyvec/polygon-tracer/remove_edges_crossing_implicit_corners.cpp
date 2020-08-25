// polyvec
#include <polyvec/polygon-tracer/remove_edges_crossing_implicit_corners.hpp>
#include <polyvec/polygon-tracer/find_2x2_corners.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace Eigen;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

int remove_edges_crossing_implicit_corners(
	const mat2x& B,                      // Raster points
	std::vector<BoundaryGraph::Edge>& E, // Edges to be pruned
	const bool circular
) {
	// Index of the edges which are to be deleted
	std::vector<size_t> E_dlist;

	// Finding the corners which should be removed
	vecXi C;
	find_2x2_corners(B, C, circular);

	// This could be simpler than O(|E| * |C|) but |C| is tipically low
	for (size_t i = 0; i < E.size(); ++i) {
		const BoundaryGraph::Edge& e = E[i];

		for (Index j = 0; j < C.size(); ++j) {
			if (PathUtils::contains_open(B, e.v0, e.v1, C(j), circular)) {
				E_dlist.emplace_back(i);
				break;
			}
		}
	}

	EraseOrdered(E, E_dlist);
	return (int)E_dlist.size();
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)