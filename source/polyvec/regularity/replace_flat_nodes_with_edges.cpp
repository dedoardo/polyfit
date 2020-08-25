// polyvec
#include <polyvec/regularity/replace_flat_nodes_with_edges.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/geometry/raster.hpp>

using namespace Eigen;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

void replace_flat_nodes_with_edges(
    const mat2x& raster,
    vecXi& polygon,
    const std::vector<bool>& raster_flat,
	polyfit::Regularity::RegularityInformation& regularities,
	std::vector<polyfit::Symmetry::SubPath>& raster_symmetries,
    std::vector<polyfit::Symmetry::SubPath>& raster_symmetries_local
) {
    PF_ASSERT(raster_flat.empty() || raster_flat.size() == raster.cols());

    bool changed;
    do {
        changed = false;
        for (Index i = 0; i < polygon.size(); ++i) {
            const vec2 pp = raster.col(CircularAt(polygon, i - 1));
            const vec2 p = raster.col(polygon(i));
            const vec2 pn = raster.col(CircularAt(polygon, i + 1));
            const bool is_perfectly_flat = abs(AngleUtils::spanned_shortest(pp, p, pn) - M_PI) < PF_EPS_MEDIUM;

            // Skipping vertices which are meant to be flat
            if (!raster_flat.empty() && raster_flat[polygon(i)]) {
                continue;
            }

            // please... use the accuracy map computed in BoundaryGraph::connect_valid_edges
            if (is_perfectly_flat) {
                const vec2 err = GeomRaster::distance_bounds_from_points_with_slack(raster, pp, pn, CircularAt(polygon, i - 1), CircularAt(polygon, i + 1));
                if (ErrorMetrics::accuracy_within_bounds(err.minCoeff(), err.maxCoeff(), pp - pn)) {
                    changed = true;
                    EraseOrdered(polygon, { i });
					regularities.reindex_after_vertex_deletion({ (int)i }, polygon.size() + 1);
                    break;
                }
            }
        }
    } while (changed);
	regularities.update(raster, polygon, raster_symmetries, raster_symmetries_local, true);
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)