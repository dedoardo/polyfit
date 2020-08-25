#include <polyvec/regularity/continuations.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/regularity/are_points_facing_inside.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>

#define NEW_PRIORITY 1
#define NEW_ANGLE_CRITERION 1

using namespace polyfit;
using namespace polyvec;
using namespace Regularity;

#define USE_CONTINUATION_LOGIC_2018 0

// Returns true if the regularity continuation is between degenerate edges
// These do not include the edges *after* the continuation, but only the continuation itself
// ----------------------------------------------------------------------------
bool has_degenerate_edges(
	const mat2x& P, const vecXi& V,
	const Continuation& r	
) {
	return GeomRaster::is_edge_degenerate(P.col(r.v0), P.col(r.v0_prev)) ||
		GeomRaster::is_edge_degenerate(P.col(r.v1), P.col(r.v1_next));
}

bool is_edge_accurate(const mat2x& B, const int vsrc, const int vdst)
{
    const vec2 psrc = B.col(vsrc);
    const vec2 pdst = B.col(vdst);
	const vec2 d_error = GeomRaster::distance_bounds_from_points_with_slack(B, psrc, pdst, vsrc, vdst);
	return ErrorMetrics::accuracy_within_bounds(d_error.minCoeff(), d_error.maxCoeff(), psrc - pdst);
}

std::vector<int> get_vertex_move_candidates(
	const mat2x& B,
	const vecXi& V,
	const int i,
	const bool circular,
	const RegularityInformation& reg,
	bool allow_move)
{
	std::vector<int> result = { 0 };

	if (!allow_move)
		return result;

	auto p = B.col(V(i));
	auto prev_id = CircularAt(V, i - 1);
	auto next_id = CircularAt(V, i + 1);
	auto prev = B.col(prev_id);
	auto next = B.col(next_id);

	bool is_x_axis_aligned = std::abs((p - prev).y()) < PF_EPS || std::abs((p - next).y()) < PF_EPS;
	bool is_y_axis_aligned = std::abs((p - prev).x()) < PF_EPS || std::abs((p - next).x()) < PF_EPS;

	bool is_self_symmetric = false;
	for (auto& s : reg.vertex_symmetries(i))
	{
		if (s.other(i) == i)
		{
			is_self_symmetric = true;
			break;
		}
	}

	int candidates[] = { -1, 1 };
	for (auto move : candidates)
	{
		if (is_self_symmetric)
			break;

		auto new_id = Circular(B, V(i) + move);
		auto new_p = B.col(new_id);
		auto diff = new_p - p;

		//preserve axis-alignment
		if (is_x_axis_aligned && std::abs(diff.y()) > PF_EPS)
			continue;

		if (is_y_axis_aligned && std::abs(diff.x()) > PF_EPS)
			continue;

		//make sure that the introduced edges fulfill accuracy
		if (!is_edge_accurate(B, prev_id, new_id) || !is_edge_accurate(B, new_id, next_id))
			continue;

		result.push_back(move);
	}

	return result;
}

// Returns if a polygon path between two vertices is longer than a given distance
bool is_subpath_longer_than(const mat2x& P, Eigen::Index polygon_index0, Eigen::Index polygon_index1, int direction, double min_distance)
{
	assert_break(polygon_index0 != polygon_index1);

	auto v = polygon_index0;
	double current_distance = 0;
	do
	{
		auto p0 = P.col(v);
		v = Circular(P, v + direction);
		auto p1 = P.col(v);

		current_distance += (p1 - p0).norm();
		// if the subpath is already longer than min_distance, we have succeeded
		if (current_distance > min_distance)
			return true;

	} while (v != polygon_index1);

	// we got to the end of the subpath with a path shorter than min_distance
	return false;
}

// Evaluates the distance from the line connecting the opposite points and the one connecting
// the first control points, returning true if the minimum of the two is within the
// expected threshold
// ----------------------------------------------------------------------------
void compute_distance_metrics(
	const mat2x& P,
	const vecXi& V,
	Continuation& r,
	const bool circular
) {
	const vec2 p0 = P.col(r.v0);
	const vec2 p1 = P.col(r.v1);
	const vec2 p0_prev = P.col(r.v0_prev);
	const vec2 p1_next = P.col(r.v1_next);
	const vec2 p_mid = .5 * (p0 + p1);

	vec2 ray_d = (p1 - p0).normalized();
	ray_d = vec2(ray_d(1), -ray_d(0));

	double ray_t_0 = INFINITY, ray_t_1 = INFINITY;
	double line_t_0 = INFINITY, line_t_1 = INFINITY;

	const bool intersect_0 = LineUtils::intersect(p_mid, p_mid + ray_d, p0_prev, p0, ray_t_0, line_t_0);
	const bool intersect_1 = LineUtils::intersect(p_mid, p_mid + ray_d, p1_next, p1, ray_t_1, line_t_1);
	if (intersect_0 && intersect_1) {
		const vec2 p_project0 = Num::lerp(p0_prev, p0, line_t_0);
		const vec2 p_project1 = Num::lerp(p1_next, p1, line_t_1);
		r.distance_midpoints_sq = (p_project0 - p_project1).squaredNorm();
	}
	else {
		r.distance_midpoints_sq = INFINITY;
	}
}

// ----------------------------------------------------------------------------
void compute_curvature_metrics(
	const mat2x& P,
	const vecXi& V,
	Continuation& r,
	const bool circular
) {
	// Continuation angle
	const vec2 p0 = P.col(r.v0);
	const vec2 p1 = P.col(r.v1);

	const int v0_next_expected = PathUtils::opposite(P.cols(), r.v0_prev, r.v0);
	const int v1_prev_expected = PathUtils::opposite(P.cols(), r.v1_next, r.v1);
	const vec2 p0_next = P.col(PathUtils::find_next_non_degenerate_corner(P, r.v0, v0_next_expected, circular));
	const vec2 p1_prev = P.col(PathUtils::find_next_non_degenerate_corner(P, r.v1, v1_prev_expected, circular));

	// Angle between the next polygon edge and the continuation edge
	r.angle_separation_0 = AngleUtils::spanned_shortest(p1, p0, p0_next);
	r.angle_separation_1 = AngleUtils::spanned_shortest(p0, p1, p1_prev);

	// Polygon angle
	r.angle_polygon_0 = M_PI - AngleUtils::spanned_shortest(P.col(r.v0_prev), p0, p0_next);
	r.angle_polygon_1 = M_PI - AngleUtils::spanned_shortest(P.col(r.v1_next), p1, p1_prev);

	// Angle between the curve tangent and the continuation edge
	r.angle_continuation_0 = AngleUtils::spanned_shortest(P.col(r.v0_prev), p0, p1);
	r.angle_continuation_1 = AngleUtils::spanned_shortest(P.col(r.v1_next), p1, p0);

	double angle_continuation_difference_00 = M_PI - AngleUtils::spanned_shortest(p0_next, p0, P.col(r.v0_prev));
	double angle_continuation_difference_01 = M_PI - AngleUtils::spanned_shortest(p1, p0, P.col(r.v0_prev));
	double angle_continuation_difference_10 = M_PI - AngleUtils::spanned_shortest(P.col(r.v1_next), p1, p1_prev);
	double angle_continuation_difference_11 = M_PI - AngleUtils::spanned_shortest(p0, p1, P.col(r.v1_next));
	r.angle_continuation_difference_0 = angle_continuation_difference_00 - angle_continuation_difference_01;
	r.angle_continuation_difference_1 = angle_continuation_difference_10 - angle_continuation_difference_11;

	const vec2 dir0 = (P.col(r.v1_next) - P.col(r.v1)).normalized();
	const vec2 dir1 = (P.col(r.v0_prev) - P.col(r.v0)).normalized();
	r.intersection_angle = acos(dir0.dot(dir1));
}

bool move_breaks_symmetry(
	int v0, int move_v0,
	int v1, int move_v1,
	const RegularityInformation& reg)
{
	//check if we are breaking a symmetry
	for (auto& sym : reg.vertex_symmetries(v0))
	{
		if (sym.other(v0) != v1)
			return true; // the vertex is symmetric to something else - do not move it
		//TODO: check if the move agrees with the symmetry axis
	}
	return false;
}

bool line_intersects_circle(
    const vec2 center,
    const double radius,
    const vec2 p0,
    const vec2 p1
) {
    double t = LineUtils::project_t(center, p0, p1);
    if (t < -PF_EPS) {
        t = 0.;
    }
    else if (t > 1. + PF_EPS) {
        t = 1.;
    }

    const vec2 p = Num::lerp(p0, p1, t);
    return (p - center).norm() < radius;
}

bool is_continuation_valid_2018 (
    const mat2x& P,
    const vecXi& V,
    Continuation& r,
    const bool circular,
    double resolution_scaling_factor,
    const RegularityInformation& reg
) {
    compute_distance_metrics(P, V, r, circular);
    compute_curvature_metrics(P, V, r, circular);

    if (has_degenerate_edges(P, V, r)) {
        return false;
    }

    if (r.move_v0 != 0)
    {
        if (move_breaks_symmetry(r.v0, r.move_v0, r.v1, r.move_v1, reg))
            return false;
    }

    if (r.move_v1 != 0)
    {
        if (move_breaks_symmetry(r.v1, r.move_v1, r.v0, r.move_v0, reg))
            return false;
    }

    if (std::min(r.angle_separation_0, r.angle_separation_1) < Options::get()->regularity_continuation_angle_separation_max - PF_EPS) {
        PF_VERBOSE_F("Fails angle separation");
        return false;
    }

    // Testing the continuation angle difference
    if (r.angle_continuation_difference_0 < -PF_RAD(5) || r.angle_continuation_difference_1 < -PF_RAD(5)) {
        PF_VERBOSE_F("Fails continuation difference");
        return false;
    }

    // Testing the intersection angle
    if (r.intersection_angle < PF_RAD(120)) {
        PF_VERBOSE_F("Fails intersection separation");
        return false;
    }

    const int i = r.v0;
    const int j = r.v1;
    const vec2 pi = P.col(i);
    const vec2 pj = P.col(j);
    const vec2 pi_prev = CircularAt(P, r.v0_prev);
    const vec2 pj_next = CircularAt(P, r.v1_next);
    if (AngleUtils::have_opposite_convexity(pi_prev, pi, pj, pj_next)) {
        const double angle_i = M_PI - AngleUtils::spanned_shortest(pi_prev, pi, pj);
        const double angle_j = M_PI - AngleUtils::spanned_shortest(pj_next, pj, pi);

        if (angle_i + angle_j > PF_RAD(50)) {
            return false;
        }
    }

    // I need to fill in 
    // angle_continuation_0/1
    // intersection_angle
    // distance_midpoints_sq
    bool is_valid_continuation = true;

    // TEST 1: Are a bunch of points around i j?
    const auto center = .5 * (pi + pj);
    const auto radius = .5 * (pi - pj).norm() * .9;
    for (int i = 0; i < V.size(); ++i) {
        const vec2 p0 = P.col(i);
        const vec2 p1 = CircularAt(P, i + 1);
        if (line_intersects_circle(center, radius, p0, p1)) {
            is_valid_continuation = false;
        }
    }

    if (is_valid_continuation) {
        return true;
    }

    // TEST2 : Same line continuation
    is_valid_continuation = true;
    const double angle_fw_i = AngleUtils::spanned_shortest(P.col(r.v0_prev), pi, pj);
    const double angle_fw_j = AngleUtils::spanned_shortest(P.col(r.v1_next), pj, pi);
    if (angle_fw_i < PF_RAD(159) || angle_fw_j < PF_RAD(159)) {
        is_valid_continuation = false;
    }
     
    if (!is_valid_continuation) {
        return false;
    }

    return true;
}

// For a possible continuation it computes the angles and distance measures and returns
// true if they satisfy the thresholds.
// ----------------------------------------------------------------------------
bool is_continuation_valid(
	const mat2x& P,
	const vecXi& V,
	Continuation& r,
	const bool circular,
	double resolution_scaling_factor,
	const RegularityInformation& reg
) {
	// Should we even match this?
	if (has_degenerate_edges(P, V, r)) {
		return false;
	}

	if (r.move_v0 != 0)
	{
		if (move_breaks_symmetry(r.v0, r.move_v0, r.v1, r.move_v1, reg))
			return false;
	}

	if (r.move_v1 != 0)
	{
		if (move_breaks_symmetry(r.v1, r.move_v1, r.v0, r.move_v0, reg))
			return false;
	}

	compute_distance_metrics(P, V, r, circular);
	compute_curvature_metrics(P, V, r, circular);

	PF_VERBOSE_F("Checking continuation %d %d %d %d\ncontinuation angles: %f %f, continuation difference: %f %f, polygon angle: %f %f, separation angle: %f %f, intersection angle: %f\nmidpoint distance: %f",
		r.v0_prev, r.v0, r.v1, r.v1_next,
		r.angle_continuation_0, r.angle_continuation_1, r.angle_continuation_difference_0, r.angle_continuation_difference_1, r.angle_polygon_0, r.angle_polygon_1, r.angle_separation_0, r.angle_separation_1,
		r.intersection_angle,
		r.distance_midpoints_sq)

		double max_distance_midpoints_sq = 2 * Options::get()->regularity_continuation_distance_max_32x32;
	max_distance_midpoints_sq *= resolution_scaling_factor;
	max_distance_midpoints_sq *= max_distance_midpoints_sq;

	// Test distance between midpoints projected on the tangents
	if (r.distance_midpoints_sq > max_distance_midpoints_sq) {
		return false;
	}

	// Test separation angle
	if (std::min(r.angle_separation_0, r.angle_separation_1) < Options::get()->regularity_continuation_angle_separation_max - PF_EPS) {
		PF_VERBOSE_F("Fails angle separation");
		return false;
	}

	// Testing the continuation angle
	if (r.angle_continuation_difference_0 < -PF_RAD(5) || r.angle_continuation_difference_1 < -PF_RAD(5)) {
		PF_VERBOSE_F("Fails continuation difference");
		return false;
	}

	if ((r.angle_continuation_0 > PF_EPS && r.angle_continuation_0 < Options::get()->regularity_continuation_angle_max) ||
		(r.angle_continuation_1 > PF_EPS && r.angle_continuation_1 < Options::get()->regularity_continuation_angle_max)) {
		PF_VERBOSE_F("Fails continuation angle");
		return false;
	}

#if NEW_ANGLE_CRITERION
	if (r.angle_polygon_0 < M_PI_4 - PF_EPS && r.angle_polygon_1 < M_PI_4 - PF_EPS) {
		PF_VERBOSE_F("Fails polygon angle");
		return false;
	}

	if (r.angle_polygon_0 < M_PI_4 + PF_EPS && r.angle_polygon_1 < M_PI_4 + PF_EPS) {
		//if it does not fulfill the stricter criterion, it has to fulfill circumference
		
		double min_subpath_length = 1.15 * (P.col(r.v0) - P.col(r.v1)).norm();
		if (!is_subpath_longer_than(P, r.v0, r.v1, +1, min_subpath_length) || !is_subpath_longer_than(P, r.v0, r.v1, -1, min_subpath_length))
		{
			PF_VERBOSE_S("Fails circumference");
			return false;
		}
	}

#else
	if (r.angle_polygon_0 < M_PI_4 + PF_EPS && r.angle_polygon_1 < M_PI_4 + PF_EPS) {
		PF_VERBOSE_F("Fails polygon angle");
		return false;
	}
#endif

	// Testing the flat polygon angle
	if (std::min(r.angle_polygon_0, r.angle_polygon_1) < Options::get()->regularity_continuation_angle_polygon_max) {
		PF_VERBOSE_F("Fails polygon angle 2");
		return false;
	}

	// Testing the intersection angle
	if (r.intersection_angle < Options::get()->regularity_continuation_angle_intersection_max) {
		PF_VERBOSE_F("Fails intersection separation");
		return false;
	}

	PF_VERBOSE_F("Passes");

	return true;
}

// Applying the continuations in order of preference
// ----------------------------------------------------------------------------
void Regularity::pick_continuations(
	const mat2x& B, vecXi& V,
	std::vector<Continuation>& candidates,
	polyfit::Regularity::RegularityInformation& regularity,
	const bool circular,
	std::vector<BoundaryGraph::Edge>& boundary_graph
) {

	auto distance_priority = [&](const Continuation& c)
	{
		return -c.distance_midpoints_sq;// -1 * (B.col(V(c.v0)) - B.col(V(c.v1))).norm();
	};

	auto angle_priority = [&](const Continuation& c)
	{
		auto ang0 = M_PI - c.angle_continuation_0;
		auto ang1 = M_PI - c.angle_continuation_1;

		auto is_inflection = std::abs(ang0 + ang1 - M_PI + c.intersection_angle) > PF_EPS;
		double inflection_penalty = (is_inflection ? 2 * (ang0 + ang1) : 0);

		return -(ang0 * ang0 + ang1 * ang1) - inflection_penalty;
		return -(ang0 + ang1) - inflection_penalty;
	};

	auto priority = [&](const Continuation& c)
	{		
		return distance_priority(c) + angle_priority(c) * 10;
	};

	sort(candidates.begin(), candidates.end(), [&](const Continuation& lhs, const Continuation& rhs) {
#if NEW_PRIORITY
		return priority(lhs) > priority(rhs);		
#else
		if (abs(lhs.distance_midpoint6s_sq - rhs.distance_midpoints_sq) < PF_EPS) {
			const double d_sq_opposite_lhs = (CircularAt(B, V(lhs.v0) + lhs.move_v0) - CircularAt(B, V(lhs.v1) + lhs.move_v1)).squaredNorm();
			const double d_sq_opposite_rhs = (CircularAt(B, V(rhs.v0) + rhs.move_v0) - CircularAt(B, V(rhs.v1) + rhs.move_v1)).squaredNorm();
			return d_sq_opposite_lhs < d_sq_opposite_rhs;
		}

		return lhs.distance_midpoints_sq < rhs.distance_midpoints_sq;	
#endif
	});

	const int INVALID = 0xffff;
	std::vector<int> vertex_movement(V.size(), INVALID);
	for (auto& r : candidates) {

		PF_VERBOSE_F("Continuation (%d -- %d + %d) - (%d + %d -- %d)", r.v0_prev, r.v0, r.move_v0, r.v1, r.move_v1, r.v1_next);
		PF_VERBOSE_F("Distance prio: %f, angle prio: %f", distance_priority(r), angle_priority(r));
		PF_VERBOSE_F("Angles: %f %f", r.angle_continuation_0, r.angle_continuation_1);

		if (vertex_movement[r.v0] != INVALID && vertex_movement[r.v0] != r.move_v0)
		{
			PF_VERBOSE_F("Already moved vertex");
			continue; //We already moved v0 differently
		}
		
		if (vertex_movement[r.v1] != INVALID && vertex_movement[r.v1] != r.move_v1)
		{
			PF_VERBOSE_F("Already moved vertex");
			continue; //We already moved v1 differently
		}

		auto target_v0 = Circular(B, V(r.v0) + r.move_v0);
		auto target_v1 = Circular(B, V(r.v1) + r.move_v1);

		if (target_v0 == CircularAt(V, r.v0 - 1) || target_v0 == CircularAt(V, r.v0 + 1))
		{
			PF_VERBOSE_F("Vertex V0 moves onto other vertex.");
			continue;
		}

		if (target_v1 == CircularAt(V, r.v1 - 1) || target_v1 == CircularAt(V, r.v1 + 1))
		{
			PF_VERBOSE_F("Vertex V1 moves onto other vertex.");
			continue;
		}

		auto cont0 = regularity.get_edge_continuation(r.v0_prev, r.v0);
		auto cont1 = regularity.get_edge_continuation(r.v1_next, r.v1);

		if (cont0 || cont1)
		{
			PF_VERBOSE_F("One of the edges already has a continuation.");
			continue;
		}			

		// Update the polygon
		if (vertex_movement[r.v0] == INVALID)
		{
			vertex_movement[r.v0] = r.move_v0;
			V(r.v0) = target_v0;
		}

		if (vertex_movement[r.v1] == INVALID)
		{
			vertex_movement[r.v1] = r.move_v1;
			V(r.v1) = target_v1;
		}

		regularity.add(std::move(r));
	}

	auto force_edge = [&](int v0, int v1)
	{
        if (CircularAt(V, v1 + 1) == V(v0))
            std::swap(v0, v1);

		v0 = V(v0);
		v1 = V(v1);

		const vec2 d_error = GeomRaster::distance_bounds_from_points_with_slack(B, B.col(v0), B.col(v1), v0, v1);

		vecXi subpath(2);
		BoundaryGraph::AccuracyMap amap;
		subpath << v0, v1;
		amap[BoundaryGraph::make_edge_id(v0, v1)] = d_error;
		BoundaryGraph::force_graph_subpath(B, boundary_graph, subpath, amap);
	};

	// fix continuations in the boundary graph
	for (auto& c : regularity.continuations())
	{
		force_edge(c.v0_prev, c.v0);
		force_edge(c.v1, c.v1_next);

        // For the edges after the continuation we want to add the edges to the boundary but 
        // not force them necessarily.
        const int v0_next = PathUtils::opposite(V.size(), c.v0_prev, c.v0);
        const int v1_prev = PathUtils::opposite(V.size(), c.v1_next, c.v1);
        
        vec2i e_0_next, e_1_prev;
        e_0_next = (Circular(V, (Eigen::Index)c.v0 + 1) == v0_next) ? vec2i(V(c.v0), V(v0_next)) : vec2i(V(v0_next), V(c.v0));
        e_1_prev = (Circular(V, (Eigen::Index)c.v1 - 1) == v1_prev) ? vec2i(V(v1_prev), V(c.v1)) : vec2i(V(c.v1), V(v1_prev));

        bool next_0_exists = false, prev_1_exists = false;
        for (size_t i = 0; i < boundary_graph.size(); ++i) {
            next_0_exists = next_0_exists || (boundary_graph[i].v0 == e_0_next(0) && boundary_graph[i].v1 == e_0_next(1));
            prev_1_exists = prev_1_exists || (boundary_graph[i].v0 == e_1_prev(0) && boundary_graph[i].v1 == e_1_prev(1));
        }

        if (!next_0_exists) {
            BoundaryGraph::Edge e;
            e.v0 = e_0_next(0);
            e.v1 = e_0_next(1);
            const vec2 d_error = GeomRaster::distance_bounds_from_points_with_slack(B, B.col(e.v0), B.col(e.v1), e.v0, e.v1);
            e.d_min = d_error(0);
            e.d_max = d_error(1);
            e.flags = 0x0;
            boundary_graph.emplace_back(e);
        }

        if (!prev_1_exists) {
            BoundaryGraph::Edge e;
            e.v0 = e_1_prev(0);
            e.v1 = e_1_prev(1);
            const vec2 d_error = GeomRaster::distance_bounds_from_points_with_slack(B, B.col(e.v0), B.col(e.v1), e.v0, e.v1);
            e.d_min = d_error(0);
            e.d_max = d_error(1);
            e.flags = 0x0;
            boundary_graph.emplace_back(e);
        }
	}	
}

// ----------------------------------------------------------------------------
std::vector<Continuation> Regularity::find_continuation_candidates(
	const mat2x& B,
	const vecXi& V,
	const RegularityInformation& reg,
	const bool circular,
	const bool allow_move
) {
	std::vector<Continuation> result;
	mat2x P;
	BoundaryGraph::trace_to_points(B, V, P);
	double resolution_scaling_factor = GeomRaster::get_resolution_scaling_factor_32x32(B);	

	std::vector<int> C;
	PathUtils::compute_convexities(P, C);

	for (int i = 0; i < V.size(); ++i) 
	{
		if (C[i] != -1)
			continue; //continuations only at concave corners
		auto moves_i = get_vertex_move_candidates(B, V, i, circular, reg, allow_move);
		for (int move_i : moves_i)
		{			
			auto i_boundary_index = Circular(B, V(i) + move_i);
			if (i_boundary_index == CircularAt(V, i + 1) || i_boundary_index == CircularAt(V, i - 1))
				continue; // do not allow collapsed vertices for now
			P.col(i) = B.col(i_boundary_index);
			for (int j = i + PV_REGULARITY_CONTINUATION_MIN_SEPARATION; j < std::min(V.size(), j - PV_REGULARITY_CONTINUATION_MIN_SEPARATION + V.size()); ++j) {
				if (C[j] != -1)
					continue; //continuations only at concave corners
				auto moves_j = get_vertex_move_candidates(B, V, j, circular, reg, allow_move);
				for (int move_j : moves_j)
				{
					auto j_boundary_index = Circular(B, V(j) + move_j);
					if (j_boundary_index == CircularAt(V, j + 1) || j_boundary_index == CircularAt(V, j - 1))
						continue; // do not allow collapsed vertices for now
					P.col(j) = B.col(j_boundary_index);

					// Given that the polyline is ordered, there are only two possible valid continuations 
					// which respect the ordering of the points					
					Continuation r1((int)Circular(V, i - 1), i, move_i, j, move_j, (int)Circular(V, j + 1));
					Continuation r2((int)Circular(V, i + 1), i, move_i, j, move_j, (int)Circular(V, j - 1));

#if USE_CONTINUATION_LOGIC_2018
                    bool r1_valid = is_continuation_valid_2018(P, V, r1, circular, resolution_scaling_factor, reg);
                    bool r2_valid = is_continuation_valid_2018(P, V, r2, circular, resolution_scaling_factor, reg);
#else
					bool r1_valid = is_continuation_valid(P, V, r1, circular, resolution_scaling_factor, reg);
					bool r2_valid = is_continuation_valid(P, V, r2, circular, resolution_scaling_factor, reg);
#endif

					if (r1_valid || r2_valid)
					{
						bool are_points_opposite = are_points_facing_inside(P, i, j);
						if (!are_points_opposite)
							continue;

						if (r1_valid)
							result.emplace_back(r1);
						if (r2_valid)
							result.emplace_back(r2);
					}					
				}
				//reset j to original position
				P.col(j) = B.col(V(j));
			}			
		}
		//reset i to original position
		P.col(i) = B.col(V(i));
	}
	return result;
}