#include <polyvec/regularity/add_local_symmetries.hpp>

#include <polyvec/core/log.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/utils/matrix.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

class LocalSymmetryRegularity
	: public RegularityAction
{
public:
	LocalSymmetryRegularity (
        const Symmetry::SubPath& symmetry,

        // Computed outside the symmetric region
        int v01_ext, int v11_ext)
		: symmetry(symmetry), v01_ext(v01_ext), v11_ext(v11_ext)
	{ }

	bool can_introduce_degenerate_edges() const { return false; }
    bool can_introduce_inflections() const { return true; }
    bool is_polygon_acceptable (const mat2x& raster, const vecXi& polygon_old, const vecXi& polygon_new, bool circular) const override { return true; }

	std::vector<size_t> get_edges_to_delete(const mat2x& raster, const vecXi& polygon, const std::vector<BoundaryGraph::Edge>& E, bool circular) const
	{
		std::vector<size_t> result;

		PF_VERBOSE_F("attempting to regularize local symmetry v01 %d v00 %d v10 %d v11 %d",
			symmetry.v01, symmetry.v00, symmetry.v10, symmetry.v11);

        const int l_sides = CircularDist(raster, symmetry.v01, symmetry.v00);
        assert(l_sides == CircularDist(raster, symmetry.v10, symmetry.v11));
        const int l_flat = CircularDist(raster, symmetry.v00, symmetry.v10);

        const bool is_length_2_or_above = (raster.col(symmetry.v00) - raster.col(symmetry.v10)).squaredNorm() > 4. - PF_EPS;
        const int smooth_kopf_0 = Circular(raster, symmetry.v00 + 1);
        const int smooth_kopf_1 = Circular(raster, symmetry.v10 - 1);

        std::vector<int> C;
        PathUtils::compute_convexities(raster, C);

        const double raster_bb_diagonal = GeomRaster::get_resolution_scaling_factor_32x32(raster);

		// Removing all edges violating  this symmetry
		for (int j = 0; j < E.size(); ++j) {
            const auto& e = E[j];

#if 0
            // Skipping edges which smooth our kopf junctions
            if (is_length_2_or_above &&
                (E[j].v0 == symmetry.v01 && E[j].v1 == smooth_kopf_0) ||
                (E[j].v0 == smooth_kopf_1 && E[j].v1 == symmetry.v11)) {
                continue;
            }
#endif

            if (BoundaryGraph::is_edge_violating_axis_aligned_edge(raster, C, e, symmetry.v00, symmetry.v10, circular, raster_bb_diagonal)) {
                result.emplace_back(j);
            }
		}

        PF_VERBOSE_F("Symmetry %d %d %d %d", symmetry.v01, symmetry.v00, symmetry.v10, symmetry.v11);
        for (size_t i = 0; i < result.size(); ++i) {
            PF_VERBOSE_F("Removing edge %d -> %d", E[result[i]].v0, E[result[i]].v1);
        }

        return result;
	}

private:
	Symmetry::SubPath symmetry;
    int v01_ext, v11_ext;
};

void add_local_symmetries(
	const mat2x& raster,
    std::vector<Symmetry::SubPath>& raster_symmetries_local,
    const std::vector<bool>& raster_flat_bits,
	const bool circular,
	std::vector<RegularityAction*>& regularity_actions) {

	std::vector<int> turns;
	PathUtils::compute_convexities(raster, turns);

	std::sort(raster_symmetries_local.begin(), raster_symmetries_local.end(), [&](const Symmetry::SubPath & lhs, const Symmetry::SubPath & rhs) -> bool {
		return Symmetry::calculate_length(raster, turns, lhs, circular) >
			Symmetry::calculate_length(raster, turns, rhs, circular);
	});

	for (auto& s : raster_symmetries_local)
	{
		if (s.axis != Symmetry::AXIS_HORIZONTAL &&
			s.axis != Symmetry::AXIS_VERTICAL) {
			continue;
		}

        int v01_ext = s.v01;
        int v11_ext = s.v11;

        PF_VERBOSE_F("Finding extended bounds for symmetry %d %d %d %d", s.v01, s.v00, s.v10, s.v11);
        while (turns[v01_ext] == turns[v11_ext] && CircularDist(raster, v01_ext, s.v00) < CircularDist(raster, v11_ext, s.v00)) {
            v01_ext = Circular(turns, v01_ext - 1);
            v11_ext = Circular(turns, v11_ext + 1);
        }
        v01_ext = Circular(turns, v01_ext + 1);
        v11_ext = Circular(turns, v11_ext - 1);
        PF_VERBOSE_F("Symmetry %d %d %d %d has extended bounds %d %d", s.v01, s.v00, s.v10, s.v11, v01_ext, v11_ext);

		regularity_actions.push_back(new LocalSymmetryRegularity(s, v01_ext, v11_ext));
	}
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)