#ifndef POLYFIT_POLYGON_TRACER_SYMMETRIC_H_
#define POLYFIT_POLYGON_TRACER_SYMMETRIC_H_

#include <polyvec/api.hpp>
#include <polyvec/typedefs.hpp>
#include <polyvec/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Finds the cycle of minimum cost for the boundary B, pruning all edges which
	cross axis aligned edges of curvature extrema (3 edge symmetries).

	in:
		R is a clockwise oriented raster polyline (grid spacing)

	out:
		V are the vertices of the polygonal approximation (indexes B)
		B is the fitting polyline to which the polygon refers to
		E are the list of edges which connect other plausible polygon edges

	Returns 0 on failure
*/
int symmetric(const mat2x& R, vecXi& P, mat2x& B, std::vector<BoundaryGraph::Edge>& E, const bool circular);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)


#endif // POLYFIT_POLYGON_TRACER_SYMMETRIC_H_