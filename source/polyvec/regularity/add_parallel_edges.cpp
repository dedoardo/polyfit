#include <polyvec/regularity/add_parallel_edges.hpp>

#include <polyvec/core/log.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/options.hpp>
#include <polyvec/regularity/check_opposite_points_winding_number.hpp>

#include <iostream>

using namespace Eigen;

namespace std {
	template <> struct hash<polyfit::vec4i>
	{
		size_t operator()(const polyfit::vec4i & x) const
		{
			return (x(0) * 73856093) ^ (x(1) * 19349663) ^ (x(2) * 83492791) ^ (x(3) * 42673801);
		}
	};
}

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

void add_edges_bypassing_vertex(const mat2x& raster, bool circular, const std::vector<BoundaryGraph::Edge>& E, int vertex, std::vector<size_t>& delete_list)
{
    // Can we preemptively remove a degenerate edge if there are other alternatives?
    int n_alternative_options_before = 0;
    int n_alternative_options_after = 0;
    Index degenerate_before = -1;
    Index degenerate_after = -1;

    const int v_force_prev = Circular(raster, vertex - 1);
    const int v_force_next = Circular(raster, vertex + 1);

    for (size_t j = 0; j < E.size(); ++j) {
        const BoundaryGraph::Edge& e = E[j];

        if (PathUtils::contains_open(raster, e.v0, e.v1, vertex, circular))
            delete_list.emplace_back(j);

        if (e.v0 != vertex && e.v1 != vertex)
            continue;

        if (e.v0 == v_force_prev) {
            PF_ASSERT(degenerate_before == -1);
            degenerate_before = j;
        }

        if (e.v1 == v_force_next) {
            PF_ASSERT(degenerate_after == -1);
            degenerate_after = j;
        }

        n_alternative_options_before += (e.v1 == vertex) ? 1 : 0;
        n_alternative_options_after += (e.v0 == vertex) ? 1 : 0;
    }

    if (degenerate_before != -1 && n_alternative_options_before > 1)
        delete_list.emplace_back(degenerate_before);

    if (degenerate_after != -1 && n_alternative_options_after > 1)
        delete_list.emplace_back(degenerate_after);
}


class ParallelEdgeRegularity
	: public RegularityAction
{
public:

	ParallelEdgeRegularity(int i00, int i01, int i10, int i11, int dim, double length)
		: i00(i00), i01(i01), i10(i10), i11(i11), dim(dim), length(length)
	{
		PF_ASSERT(i00 != i01 && i10 != i11);
		PF_VERBOSE_F("New parallel edge pair along axis %d", dim);
	}
	
	bool can_introduce_degenerate_edges() const override { return true; }
    bool can_introduce_inflections()  const override { return true; }
    bool is_polygon_acceptable(const mat2x& raster, const vecXi& polygon_old, const vecXi& polygon_new, bool circular ) const override { return true; }

	std::vector<size_t> get_edges_to_delete(const mat2x& raster, const vecXi& polygon, const std::vector<BoundaryGraph::Edge>& E, bool circular) const
	{
		PF_VERBOSE_F("attempting to regularize parallel edge pair v01 %d v00 %d v10 %d v11 %d",
			i01, i00, i10, i11);

		std::vector<size_t> result;

		// Checking that the two edges are orthogonal to the overlap
		PF_ASSERT(abs(raster(1 - dim, i00) - raster(1 - dim, i01)) < PF_EPS);
		PF_ASSERT(abs(raster(1 - dim, i10) - raster(1 - dim, i11)) < PF_EPS);

		// Getting the values of the interval bounds properly ordered
		int i0_min, i0_max, i1_min, i1_max;
		i0_min = raster(dim, i00) < raster(dim, i01) ? i00 : i01;
		i0_max = raster(dim, i00) < raster(dim, i01) ? i01 : i00;
		i1_min = raster(dim, i10) < raster(dim, i11) ? i10 : i11;
		i1_max = raster(dim, i10) < raster(dim, i11) ? i11 : i10;

		// The same values, but ordered consistently w.r.t the polyline
		int i0_src, i0_dst, i1_src, i1_dst;
		i0_src = CircularDist(raster, i00, i01) < CircularDist(raster, i01, i00) ? i00 : i01;
		i0_dst = i0_src == i00 ? i01 : i00;
		i1_src = CircularDist(raster, i10, i11) < CircularDist(raster, i11, i10) ? i10 : i11;
		i1_dst = i1_src == i10 ? i11 : i10;

		const double b0_min = raster(dim, i0_min);
		const double b0_max = raster(dim, i0_max);
		const double b1_min = raster(dim, i1_min);
		const double b1_max = raster(dim, i1_max);
		PF_ASSERT(b0_max > b0_min && b1_max > b1_min); // this is just because I am tired

		// We are interested in configurations when the polygon is passing through 
		// a pair of endpoints (not necessarily aligned) and is *NOT* passing through
		// a pair of aligned interval endpoints
		bool contains0_min = false, contains0_max = false;
		bool contains1_min = false, contains1_max = false;
		bool contains0_mid = false, contains1_mid = false;

		for (int j = 0; j < polygon.size(); ++j) {
			contains0_min = polygon(j) == i0_min ? true : contains0_min;
			contains0_max = polygon(j) == i0_max ? true : contains0_max;
			contains1_min = polygon(j) == i1_min ? true : contains1_min;
			contains1_max = polygon(j) == i1_max ? true : contains1_max;
			contains0_mid = PathUtils::contains_open(raster, i0_src, i0_dst, polygon(j), circular) ? true : contains0_mid;
			contains1_mid = PathUtils::contains_open(raster, i1_src, i1_dst, polygon(j), circular) ? true : contains1_mid;
		}

		// Nothing to do
		if (contains0_min && contains0_max && contains1_min && contains1_max) {
			PF_VERBOSE_S("nothing to do, interval is already satisfied");
			return result;
		}

		// We want the corner to be inserted to be aligned, while the polygon should
		// go through the opposite
		const bool contains_min = contains0_min || contains1_min;
		const bool contains_max = contains0_max || contains1_max;
		const bool can_regularize_min = contains_max || contains_min;// && !contains_min;
		const bool can_regularize_max = contains_min || contains_max;// && !contains_max;

		// Not regularizing blinding, at least one of the two segments should already be
		// in place as an hint.
		//if (can_regularize_min && can_regularize_max) {
		//	PF_LOGS("skipping interval, too many things to regularize");
		//	continue;
		//}

		// More importantly, we want to increase parallelism and not break existing one.
		// If the polygon is already parallel, we skip the processing as regularizing any two
		// of the endpoints will break it.
		// We identify such segments by checking if they lay in the same quadrant
		vec2i ie0_min(-1, -1),
			ie0_max(-1, -1),
			ie1_min(-1, -1),
			ie1_max(-1, -1);

		for (int j = 0; j < polygon.size(); ++j) {
			if (!circular && j == polygon.size() - 1) {
				continue;
			}

			const vec2i ie(polygon(j), CircularAt(polygon, j + 1));

			if (PathUtils::contains_open(raster, ie(0), ie(1), i0_min, circular)) {
				PF_ASSERT(ie0_min.minCoeff() < 0);
				ie0_min = ie;
			}

			if (PathUtils::contains_open(raster, ie(0), ie(1), i0_max, circular)) {
				PF_ASSERT(ie0_max.minCoeff() < 0);
				ie0_max = ie;
			}

			if (PathUtils::contains_open(raster, ie(0), ie(1), i1_min, circular)) {
				PF_ASSERT(ie1_min.minCoeff() < 0);
				ie1_min = ie;
			}

			if (PathUtils::contains_open(raster, ie(0), ie(1), i1_max, circular)) {
				PF_ASSERT(ie1_max.minCoeff() < 0);
				ie1_max = ie;
			}
		}

		bool are_parallel_min = false, are_parallel_max = false;
		if (ie0_min.minCoeff() >= 0 && ie1_min.minCoeff() >= 0) {
			const vec2 e0_min = raster.col(ie0_min(1)) - raster.col(ie0_min(0));
			const vec2 e1_min = raster.col(ie1_min(1)) - raster.col(ie1_min(0));

			are_parallel_min = AngleUtils::lie_in_the_same_quadrant(e0_min, e1_min) ||
				AngleUtils::lie_in_the_same_quadrant(e0_min, -e1_min);
		}

		if (ie0_max.minCoeff() >= 0 && ie1_max.minCoeff() >= 0) {
			const vec2 e0_max = raster.col(ie0_max(1)) - raster.col(ie0_max(0));
			const vec2 e1_max = raster.col(ie1_max(1)) - raster.col(ie1_max(0));

			are_parallel_max = AngleUtils::lie_in_the_same_quadrant(e0_max, e1_max) ||
				AngleUtils::lie_in_the_same_quadrant(e0_max, -e1_max);
		}

		// Finding the vertex of the polygon where the two new corners should be introduced
		if (can_regularize_min && !are_parallel_min) {
			add_edges_bypassing_vertex(raster, circular, E, i0_min, result);
			add_edges_bypassing_vertex(raster, circular, E, i1_min, result);			
			PF_VERBOSE_F("forcing path through interval min %d %d", (int)i0_min, (int)i1_min);
		}

		if (can_regularize_max && !are_parallel_max) {
			add_edges_bypassing_vertex(raster, circular, E, i0_max, result);
			add_edges_bypassing_vertex(raster, circular, E, i1_max, result);
            PF_VERBOSE_F("forcing path through interval max %d %d", (int)i0_max, (int)i1_max);
		}

		// Another case we are interested in regularizing is when one of the edges is parallel and the
		// other contains a midpoint.
		const bool no_opposite_pair = !can_regularize_min && !can_regularize_max; // These avoid inserting duplicate indices
		if (no_opposite_pair && contains0_min && contains0_max && contains1_mid) {
			PF_ASSERT(!contains0_mid);
			add_edges_bypassing_vertex(raster, circular, E, i1_min, result);
			add_edges_bypassing_vertex(raster, circular, E, i1_max, result);			
            PF_VERBOSE_F("forcing path through interval 1 %d %d", (int)i1_min, (int)i1_max);
		}

		if (no_opposite_pair && contains1_min && contains1_max && contains0_mid) {
			PF_ASSERT(!contains1_mid);
			add_edges_bypassing_vertex(raster, circular, E, i0_min, result);
			add_edges_bypassing_vertex(raster, circular, E, i0_max, result);
            PF_VERBOSE_F("forcing path through interval 0 %d %d", (int)i0_min, (int)i0_max);
		}

        std::sort(result.begin(), result.end());
        EraseDuplicatesOrdered(result);
		return result;
	}

private:
	//Raster indices of parallel edge pairs
	int i00, i01, i10, i11;

	//The dimension in which the underlying edges extend
	int dim;

	double length;
};

// assumes the two intervals are ordered
bool test_overlap_interval_1d(
	const double a0,
	const double a1,
	const double b0,
	const double b1,
	vec2& i
) {
	if (b0 > a1 || a0 > b1) {
		return false;
	}

	i(0) = std::max(a0, b0);
	i(1) = std::min(a1, b1);
	return (i(1) - i(0)) > PF_EPS_LARGE;
}

void add_parallel_edges(
	const mat2x& raster,
	const bool circular,
	std::vector<RegularityAction*>& regularity_actions) {

	const double proximity_scale = Options::get()->regularity_raster_parallel_proximity_scale;
	const double min_length_ratio = Options::get()->regularity_raster_parallel_min_length_ratio;
	const double min_feature_ratio = Options::get()->regularity_raster_parallel_min_feature_ratio;

	std::vector<int> C;
	PathUtils::compute_convexities(raster, C);

	vec2 bb_min(INFINITY, INFINITY);
	vec2 bb_max(-INFINITY, -INFINITY);

	for (int i = 0; i < raster.cols(); ++i) {
		bb_min = bb_min.cwiseMin(raster.col(i));
		bb_max = bb_max.cwiseMax(raster.col(i));
	}

	vec2 figure_size = bb_max - bb_min;

	std::unordered_set<vec4i> seen_pairs;

	//TODO: To make it faster, first find axis-aligned edges, then test the overlap criterion
	for (int i = 0; i < (int)C.size(); ++i) {
		if (!circular && i == C.size() - 1) {
			continue;
		}

		// corner to corner
		if (C[i] == 0) {
			continue;
		}

		const Vertex inext = PathUtils::next_transition_vertex(C, i, +1, circular);
		if (inext == -1) {
			continue;
		}

		const vec2 p0 = raster.col(i);
		const vec2 p1 = CircularAt(raster, inext);

		// skipping length 1 segments
		if (std::max(std::abs(p0(0) - p1(0)), std::abs(p0(1) - p1(1))) < 1. + PF_EPS) {
			continue;
		}

		const double i_length = (p1 - p0).norm();		

		for (int j = inext; j < (int)C.size(); ++j) {
			if (!circular && j == C.size() - 1) {
				continue;
			}

			// I don't understand how this can be true? (todo)
			if (i == j) {
				continue;
			}

			if (C[j] == 0) {
				continue;
			}

			const Vertex jnext = PathUtils::next_transition_vertex(C, j, +1, circular);
			if (jnext == -1) {
				continue;
			}

			const vec2 p2 = raster.col(j);
			const vec2 p3 = CircularAt(raster, jnext);

			// Test relative length of the two segments			
			const double j_length = (p3 - p2).norm();
			const double length_ratio = std::min(i_length, j_length) / std::max(i_length, j_length);
			if (length_ratio < min_length_ratio) {
				continue;
			}

			// skipping length 1 segments
			if (std::max(abs(p2(0) - p3(0)), abs(p2(1) - p3(1))) < 1. + PF_EPS) {
				continue;
			}

			// Test feature length
			const double figure_ratio = std::min(i_length, j_length) / figure_size.x();
			if (figure_ratio < min_feature_ratio) {
				continue;
			}

			vec2 interval;

			for (int dim = 0; dim < 2; ++dim)
			{
				if (test_overlap_interval_1d(
					std::min(p0(dim), p1(dim)), std::max(p0(dim), p1(dim)),
					std::min(p2(dim), p3(dim)), std::max(p2(dim), p3(dim)),
					interval
				)) {
					PF_ASSERT(std::abs(p0(1 - dim) - p1(1 - dim)) < PF_EPS);
					PF_ASSERT(std::abs(p2(1 - dim) - p3(1 - dim)) < PF_EPS);

					const double overlap = std::abs(interval(0) - interval(1)); // abs unnecessary?
					const double proximity = std::abs(p0(1 - dim) - p2(1 - dim));					

					Eigen::Matrix2i ind;

					ind(0, 0) = p0(dim) < p1(dim) ? i : inext;
					ind(0, 1) = p0(dim) < p1(dim) ? inext : i;
					ind(1, 0) = p2(dim) < p3(dim) ? j : jnext;
					ind(1, 1) = p2(dim) < p3(dim) ? jnext : j;

					bool is_configuration_invalid = false;

					// Checking convexity of the lower bound
					for (int side = 0; side < 2; ++side)
					{
						if (is_configuration_invalid)
							continue;
						int relationFlip = side == 0 ? 1 : -1;
						if (std::abs(raster(dim, ind(0, side)) - raster(dim, ind(1, side))) < PF_EPS_LARGE) {
							is_configuration_invalid = C[ind(0, side)] != -1 || C[ind(1, side)] != -1;
						}
						else if (raster(dim, ind(0, side)) * relationFlip > raster(dim, ind(1, side)) * relationFlip) {
							is_configuration_invalid = C[ind(0, side)] != -1;
						}
						else {
							is_configuration_invalid = C[ind(1, side)] != -1;
						}
					}					

					vec4i pair(i, inext, j, jnext);
					if (!is_configuration_invalid && overlap > proximity_scale * proximity - PF_EPS && seen_pairs.find(pair) == seen_pairs.end()) {						
						seen_pairs.insert(pair);

                        // Are the two points facing each other?
                        if (!polyfit::Regularity::check_opposite_points_winding_number(raster, i, j) &&
                            !polyfit::Regularity::check_opposite_points_winding_number(raster, i, jnext) &&
                            !polyfit::Regularity::check_opposite_points_winding_number(raster, inext, j) &&
                            !polyfit::Regularity::check_opposite_points_winding_number(raster, inext, jnext)) {
                            continue;
                        }

						regularity_actions.push_back(new ParallelEdgeRegularity(i, inext, j, jnext, dim, i_length + j_length));
					}
				}
			}			
		}
	}
}

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)