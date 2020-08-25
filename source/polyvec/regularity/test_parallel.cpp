// Polyvec
#include <polyvec/regularity/test_parallel.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/core/constants.hpp>

using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

void test_and_set_parallel_alignment(
    const Eigen::Matrix2Xd& B,
    const std::vector<int>& C,
    const Eigen::VectorXi& V,
    Parallel& r
) {
    const vec2 p00 = B.col(V(r.v00));
    const vec2 p01 = B.col(V(r.v01));
    const vec2 p10 = B.col(V(r.v10));
    const vec2 p11 = B.col(V(r.v11));

    const vec2 d_axis_0 = (p00 - p01).normalized();
    const vec2 d_axis_1 = (p11 - p10).normalized();
    const double a_axis_0 = atan2(d_axis_0.y(), d_axis_0.x());
    const double a_axis_1 = atan2(d_axis_1.y(), d_axis_1.x());
    const double a_axis_lerp = .5 * (a_axis_0 + a_axis_1);
    vec2 d_axis_lerp(cos(a_axis_lerp), sin(a_axis_lerp));
    d_axis_lerp.normalize(); // ehm...

    const double t_project_00 = abs(LineUtils::project_t(p00, p11, p11 + d_axis_lerp));
    const double t_project_11 = abs(LineUtils::project_t(p11, p00, p00 + d_axis_lerp));
    const double t_project_01 = abs(LineUtils::project_t(p01, p10, p10 + d_axis_lerp));
    const double t_project_10 = abs(LineUtils::project_t(p10, p01, p01 + d_axis_lerp));

    double t_project_scaled = PV_REGULARITY_MAX_ALIGNMENT_DISTANCE;
    t_project_scaled *= GeomRaster::get_resolution_scaling_factor_32x32(B);

    r.aligned_00_11 = t_project_00 < t_project_scaled && t_project_11 < t_project_scaled;
    r.aligned_01_10 = t_project_01 < t_project_scaled && t_project_10 < t_project_scaled;

    r.aligned_00_11 = r.aligned_00_11 || GeomRaster::are_points_axis_aligned(p00, p11);
    r.aligned_01_10 = r.aligned_01_10 || GeomRaster::are_points_axis_aligned(p01, p10);

    // The need to be both convex
    r.aligned_00_11 = r.aligned_00_11 && C[r.v00] == -1 && C[r.v11] == -1;
    r.aligned_01_10 = r.aligned_01_10 && C[r.v01] == -1 && C[r.v10] == -1;
}

bool test_parallel(
    const Eigen::Matrix2Xd& B,
    const Eigen::VectorXi& V,
    const int vi0,
    const int vi1,
    const int vj0,
    const int vj1,
    polyfit::Regularity::Parallel& r
) {
    const vec2 pi0 = B.col(V(vi0));
    const vec2 pi1 = B.col(V(vi1));
    const vec2 pj0 = B.col(V(vj0));
    const vec2 pj1 = B.col(V(vj1));

    auto d = (pj0 - pi1); //the vector connecting the parallel vertex pair	

    // Connecting line must be (almost) axis-aligned or at a 45° angle
    if (d.cwiseAbs().minCoeff() > 1. + PF_EPS && std::abs(std::abs(d.x()) - std::abs(d.y())) > 1.0 + PF_EPS) {
        return false;
    }

    // They should be roughly aligned
    const vec2 di = (pi1 - pi0);
    const vec2 dj = (pj1 - pj0);
    const double deviation = acos(-di.dot(dj) / (d.norm() * dj.norm()));
    if (deviation > PV_REGULARITY_PARALLEL_MAX_DEVIATION_ANGLE) {
        return false;
    }

    if ((di.cwiseAbs().minCoeff() < PF_EPS && dj.cwiseAbs().minCoeff() >= PF_EPS) || (dj.cwiseAbs().minCoeff() < PF_EPS && di.cwiseAbs().minCoeff() >= PF_EPS))
        return false; //if one of them is axis-aligned and the other is not

    // ... and the thickness ratio of both ends should be close to 1
    // one of the thicknesses is len(d)
    // find the other by projecting the endpoints onto the other line in the direction of d
    double endThickness = d.norm();
    double farThickness;
    double t1, t2;
    LineUtils::intersect(pi0, pi0 + d, pj0, pj1, t1, t2);
    if (t1 > 0 && t2 >= 0 && t2 <= 1)
        farThickness = t1 * endThickness;
    LineUtils::intersect(pj1, pj1 - d, pi0, pi1, t1, t2);
    if (t1 > 0 && t2 >= 0 - PF_EPS && t2 <= 1 + PF_EPS)
        farThickness = t1 * endThickness;
    if (std::max(endThickness, farThickness) / std::min(endThickness, farThickness) >= 2.0 - PF_EPS_MEDIUM)
        return false;

    // Skipping tiny matches
    const double len0_sq = (pi0 - pi1).squaredNorm();
    const double len1_sq = (pj0 - pj1).squaredNorm();

    if (len0_sq < 1. + PF_EPS_MEDIUM || len1_sq < 1. + PF_EPS_MEDIUM) {
        return false;
    }

    const double t00 = LineUtils::project_t(pi0, pj0, pj1);
    const double t01 = LineUtils::project_t(pi1, pj0, pj1);
    const double t10 = LineUtils::project_t(pj0, pi0, pi1);
    const double t11 = LineUtils::project_t(pj1, pi0, pi1);

    double overlap0, overlap1;
    if (!Num::test_and_calculate_interval_overlap(0., 1., min(t00, t01), max(t00, t01), overlap0) ||
        !Num::test_and_calculate_interval_overlap(0., 1., min(t10, t11), max(t10, t11), overlap1)) {
        return false;
    }

    const double length0 = std::sqrt(len0_sq);
    const double length1 = std::sqrt(len1_sq);
    const double overlap0_relative = overlap0 * length1 / min(length1, length0);
    const double overlap1_relative = overlap1 * length0 / min(length1, length0);

    const double min_overlap_thr = .75;
    if (overlap0_relative < min_overlap_thr || overlap1_relative < min_overlap_thr) {
        return false;
    }

    // checking if the match is visible in the local neighborhood
    const double distance = (((pi0 + pi1) * .5) - (pj0 + pj1) * .5).norm();
    if (distance > max(length0, length1) * .75) {
        return false;
    }

    // check if the raster supports this parallelity
    // find the primary direction of the two edges (the axis with maximum extent)
    auto get_primary_direction = [](const vec2& d) {
        if (std::abs(std::abs(d.x()) - std::abs(d.y())) < PF_EPS)
            return -1; //45°
        else if (std::abs(d.x()) > std::abs(d.y()))
            return 0;
        else
            return 1;
    };

    int primary_i = get_primary_direction(di);
    int primary_j = get_primary_direction(dj);
    if (primary_i == -1)
        primary_i = primary_j;
    if (primary_j == -1)
        primary_j = primary_i;
    if (primary_i != primary_j)
        return false; //edges lie in different octants
    int primary = primary_i == -1 ? 0 : primary_i;

    // now walk along the raster in the primary direction and count the 
    // steps

    // counts the maximum and minimum step in the primary direction
    auto count_max_step = [&](int v_from, int v_to)
    {
        vec2 maxSteps(0, 0);
        int last_step = v_from;
        int v = v_from;
        vec2 minmax(0, 0);
        while (v != v_to)
        {
            auto last_v = v;
            v = Circular(B, v + 1);
            if (std::abs((B.col(last_v) - B.col(v))(1 - primary)) > PF_EPS || v == v_to)
            {
                double current_steps = (B.col(v) - B.col(last_step))(primary);
                last_step = v;
                if (current_steps < minmax(0))
                    minmax(0) = current_steps;
                if (current_steps > minmax(1))
                    minmax(1) = current_steps;
            }
        }
        return minmax;
    };

    auto minmaxi = count_max_step(V(vi0), V(vi1));
    auto minmaxj = count_max_step(V(vj0), V(vj1));

    // one of the lines is in the opposite direction
    minmaxj = -minmaxj;
    std::swap(minmaxj(0), minmaxj(1));

    if ((minmaxi - minmaxj).cwiseAbs().maxCoeff() > 1. + PF_EPS)
        return false;

    r.v00 = vi0;
    r.v01 = vi1;
    r.v10 = vj0;
    r.v11 = vj1;
    r.connected_00_11 = false;//CircularDist(V, vj1, vi0) <= 1;
    r.connected_01_10 = true;// CircularDist(V, vi1, vj0) <= 1;
    return true;
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)