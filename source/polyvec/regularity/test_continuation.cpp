// Regularity
#include <polyvec/regularity/test_continuations.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/options.hpp>

// libc++
#include <algorithm>

using namespace polyvec;
using namespace polyfit;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

namespace {
    bool test_degenerate_edges(
        const mat2x& PP_0,
        const mat2x& PP_1,
        const Continuation& r
    ) {
        return GeomRaster::is_edge_degenerate(PP_0.col(r.v0), PP_0.col(r.v0_prev)) ||
            GeomRaster::is_edge_degenerate(PP_1.col(r.v1), PP_1.col(r.v1_next));
    }
}

bool test_continuation(
    const Eigen::Matrix2Xd& PP_0,
    const Eigen::Matrix2Xd& PP_1,
    Continuation& r
) {
    // Should we even match this?
    if (test_degenerate_edges(PP_0, PP_1, r)) {
        return false;
    }

    compute_continuation_accuracy_metrics(PP_0, PP_1, r);
    compute_continuation_curvature_metrics(PP_0, PP_1, r);
    bool fail = false;

    double max_distance_midpoints_sq = 2 * Options::get()->regularity_continuation_distance_max_32x32;
    max_distance_midpoints_sq *= max(GeomRaster::get_resolution_scaling_factor_32x32(PP_0), GeomRaster::get_resolution_scaling_factor_32x32(PP_1));
    max_distance_midpoints_sq *= max_distance_midpoints_sq;

    // Test distance between midpoints projected on the tangents
    if (r.distance_midpoints_sq > max_distance_midpoints_sq) {
        fail = true;
    }

    // Test separation angle
    //PF_DEV_F("r.angle_separation_0 %f < %f", PF_DEG(r.angle_separation_0), PF_DEG(Options::get()->regularity_continuation_angle_intersection_max;
    if (min(r.angle_separation_0, r.angle_separation_1) < Options::get()->regularity_continuation_angle_separation_max - PF_EPS) {
        fail = true;
    }

    // Testing the continuation angle
    if (r.angle_continuation_difference_0 < -PF_RAD(5) || r.angle_continuation_difference_1 < -PF_RAD(5)) {
        fail = true;
    }

    if ((r.angle_continuation_0 > PF_EPS&& r.angle_continuation_0 < Options::get()->regularity_continuation_angle_max) ||
        (r.angle_continuation_1 > PF_EPS&& r.angle_continuation_1 < Options::get()->regularity_continuation_angle_max)) {
        fail = true;
    }

    // Testing the the two curves are in different quadrants
    // ____/    /___ GOOD
    // ____/    \___ BAD
    if (r.angle_polygon_0 < M_PI_4 - PF_EPS && r.angle_polygon_1 < M_PI_4 - PF_EPS) {
        fail = true;
    }

    // Handling angles at perfectly 45 degrees
    const bool angle_tangent_0_PI_4 = abs(r.angle_polygon_0 - M_PI_4) < PF_EPS;
    const bool angle_tangent_1_PI_4 = abs(r.angle_polygon_1 - M_PI_4) < PF_EPS;
    if ((angle_tangent_0_PI_4 && r.angle_polygon_1 < M_PI_4 - PF_EPS) ||
        (angle_tangent_1_PI_4 && r.angle_polygon_0 < M_PI_4 - PF_EPS)) {
        fail = true;
    }

    // Testing the flat polygon angle
    if (min(r.angle_polygon_0, r.angle_polygon_1) < Options::get()->regularity_continuation_angle_polygon_max) {
        fail = true;
    }

    // Testing the intersection angle
    if (r.intersection_angle < Options::get()->regularity_continuation_angle_intersection_max) {
        fail = true;
    }

    return !fail;
}

void compute_continuation_accuracy_metrics(
    const Eigen::Matrix2Xd& PP_0, const Eigen::Matrix2Xd& PP_1,
    Continuation& r
) {
    const vec2 p0 = PP_0.col(r.v0);
    const vec2 p1 = PP_1.col(r.v1);
    const vec2 p0_prev = PP_0.col(r.v0_prev);
    const vec2 p1_next = PP_1.col(r.v1_next);
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

void compute_continuation_curvature_metrics(
    const Eigen::Matrix2Xd& PP_0, const Eigen::Matrix2Xd& PP_1,
    Continuation& r
) {
    // Continuation angle
    const vec2 p0 = PP_0.col(r.v0);
    const vec2 p1 = PP_1.col(r.v1);

    const int v0_next_expected = PathUtils::opposite(PP_0.cols(), r.v0_prev, r.v0);
    const int v1_prev_expected = PathUtils::opposite(PP_1.cols(), r.v1_next, r.v1);
    const int v0_next_guessed = PathUtils::find_next_non_degenerate_corner(PP_0, r.v0, v0_next_expected, true);
    const int v1_prev_guessed = PathUtils::find_next_non_degenerate_corner(PP_1, r.v1, v1_prev_expected, true);
    const vec2 p0_next = PP_0.col(v0_next_guessed);
    const vec2 p1_prev = PP_1.col(v1_prev_guessed);
    r.v0_next_guessed = v0_next_guessed;
    r.v1_prev_guessed = v1_prev_guessed;

    // Angle between the next polygon edge and the continuation edge
    r.angle_separation_0 = AngleUtils::spanned_shortest(p1, p0, p0_next);
    r.angle_separation_1 = AngleUtils::spanned_shortest(p0, p1, p1_prev);

    // Polygon angle
    r.angle_polygon_0 = M_PI - AngleUtils::spanned_shortest(PP_0.col(r.v0_prev), p0, p0_next);
    r.angle_polygon_1 = M_PI - AngleUtils::spanned_shortest(PP_1.col(r.v1_next), p1, p1_prev);

    // Angle between the curve tangent and the continuation edge
    r.angle_continuation_0 = AngleUtils::spanned_shortest(PP_0.col(r.v0_prev), p0, p1);
    r.angle_continuation_1 = AngleUtils::spanned_shortest(PP_1.col(r.v1_next), p1, p0);

    double angle_continuation_difference_00 = M_PI - AngleUtils::spanned_shortest(p0_next, p0, PP_0.col(r.v0_prev));
    double angle_continuation_difference_01 = M_PI - AngleUtils::spanned_shortest(p1, p0, PP_0.col(r.v0_prev));
    double angle_continuation_difference_10 = M_PI - AngleUtils::spanned_shortest(PP_1.col(r.v1_next), p1, p1_prev);
    double angle_continuation_difference_11 = M_PI - AngleUtils::spanned_shortest(p0, p1, PP_1.col(r.v1_next));

    // This should not really be done here
    angle_continuation_difference_00 = abs(angle_continuation_difference_00 - M_PI) < PF_EPS_MEDIUM ? 0. : angle_continuation_difference_00;
    angle_continuation_difference_01 = abs(angle_continuation_difference_01 - M_PI) < PF_EPS_MEDIUM ? 0. : angle_continuation_difference_01;
    angle_continuation_difference_10 = abs(angle_continuation_difference_10 - M_PI) < PF_EPS_MEDIUM ? 0. : angle_continuation_difference_10;
    angle_continuation_difference_11 = abs(angle_continuation_difference_11 - M_PI) < PF_EPS_MEDIUM ? 0. : angle_continuation_difference_11;

    r.angle_continuation_difference_0 = angle_continuation_difference_00 - angle_continuation_difference_01;
    r.angle_continuation_difference_1 = angle_continuation_difference_10 - angle_continuation_difference_11;

    const vec2 dir0 = (PP_1.col(r.v1_next) - PP_1.col(r.v1)).normalized();
    const vec2 dir1 = (PP_0.col(r.v0_prev) - PP_0.col(r.v0)).normalized();
    r.intersection_angle = acos(dir0.dot(dir1));
}

void filter_same_side_continuations_greedy(
    const std::vector<const Eigen::Matrix2Xd*>& polygons,
    std::vector<Continuation>& R,
    const bool circular
) {
    sort(R.begin(), R.end(), [&polygons](const Continuation& lhs, const Continuation& rhs) {
        if (abs(lhs.distance_midpoints_sq - rhs.distance_midpoints_sq) < PF_EPS) {
            const double d_sq_opposite_lhs = (polygons[lhs.polygon_0]->col(lhs.v0) - polygons[lhs.polygon_1]->col(lhs.v1)).squaredNorm();
            const double d_sq_opposite_rhs = (polygons[rhs.polygon_0]->col(rhs.v0) - polygons[rhs.polygon_1]->col(rhs.v1)).squaredNorm();
            return d_sq_opposite_lhs < d_sq_opposite_rhs;
        }

        return lhs.distance_midpoints_sq < rhs.distance_midpoints_sq;
        });

    std::vector<Continuation> R_unique;
    for (size_t i = 0; i < R.size(); ++i) {
        bool discard = false;

        for (size_t j = 0; !discard && j < R_unique.size(); ++j) {
            const bool same_polygon_0 = R[i].polygon_0 == R_unique[j].polygon_0;
            const bool same_polygon_1 = R[i].polygon_1 == R_unique[j].polygon_1;
            const bool same_side_00 = R[i].v0 == R_unique[j].v0 && R[i].v0_prev == R_unique[j].v0_prev;
            const bool same_side_01 = R[i].v0 == R_unique[j].v1 && R[i].v0_prev == R_unique[j].v1_next;
            const bool same_side_10 = R[i].v1 == R_unique[j].v0 && R[i].v1_next == R_unique[j].v0_prev;
            const bool same_side_11 = R[i].v1 == R_unique[j].v1 && R[i].v1_next == R_unique[j].v1_next;
            discard = ((same_polygon_0 && (same_side_00 || same_side_01)) ||
                      (same_polygon_1 && (same_side_10 || same_side_11)));
        }

        if (discard) {
            continue;
        }

        R_unique.emplace_back(R[i]);
    }

    R = move(R_unique);
}

NAMESPACE_END(polyfit)
NAMESPACE_END(Regularity)