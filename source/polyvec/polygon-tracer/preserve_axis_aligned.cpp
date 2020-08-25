#include <polyvec/polygon-tracer/preserve_axis_aligned.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace polyfit;
using namespace PolygonTracer;

void PolygonTracer::preserve_axis_aligned_edges(ShortestPath::State & G_state, const mat2x & B, const std::vector<Vertex>& M, std::vector<BoundaryGraph::Edge>& E, vecXi & P, const bool circular)
{
	// getting all the axis-aligned edges of length > 2
	std::vector<vec2i> B_aa;

	for (Vertex i = 0; i < P.rows(); ++i) {
		if (!circular && i == P.rows() - 1) {
			break;
		}

		const Vertex v = P(i);
		const Vertex vn = P((i + 1) % P.rows());
		const vec2 d = (B.col(vn) - B.col(v));

        if (d.squaredNorm() < 4. + PF_EPS) {
            continue;
        }

		if (std::abs(d.x()) < PF_EPS || std::abs(d.y()) < PF_EPS) {
			PF_VERBOSE_F("found axis-aligned edge %lld - %lld", v, vn);
			B_aa.emplace_back(vec2i(v, vn));
		}
	}

	// the graph construction guarantees that an axis-aligned edge will exactly
	// overlap the boundary.
	std::vector<size_t> E_dlist;

    std::vector<int> C;
    PathUtils::compute_convexities(B, C);

    const double B_bb_diagonal = GeomRaster::get_resolution_scaling_factor_32x32(B);

	for (size_t i = 0; i < E.size(); ++i) {
		for (size_t j = 0; j < B_aa.size(); ++j) {
			if (BoundaryGraph::is_edge_violating_axis_aligned_edge(B, C, E[i], B_aa[j](0), B_aa[j](1), circular, B_bb_diagonal)) {
				E_dlist.emplace_back(i);
				break;
			}
		}
	}

	EraseOrdered(E, E_dlist);

	// rerunning a shortest path, this should not fail but return the exact same result
	const vecXi Pold = P;
	if (!PolygonTracer::minimum(G_state, B, M, E, P, circular)) {
		PF_VERBOSE_F("shortest cycle with %lld edges should not fail", E.size());
		PF_ABORT;
	}

	// the shortest cycle heuristic is not deterministic
	// todo: check all the entries
	PF_ASSERT(Pold.cols() == P.cols());
}
