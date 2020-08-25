// Polyvec
#include <polyvec/regularity/are_points_facing_inside.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/num.hpp>

using namespace Eigen;
using namespace polyfit;
using namespace std;

NAMESPACE_BEGIN(polyvec)

bool are_points_facing_inside(
    const Eigen::Matrix2Xd& P,
    const int v0,
    const int v1
) {
    const vec2 p0 = P.col(v0);
    const vec2 p1 = P.col(v1);

    double p_t = INFINITY, s_t = INFINITY;
    for (size_t i = 0; i < P.cols(); ++i) {
        // Ignore immediately neighboring segments
		if (Circular(P, i + 1) == v0 || i == v0)
			continue;
		if (Circular(P, i + 1) == v1 || i == v1)
			continue;
        
        const vec2 s0 = P.col(i);
        const vec2 s1 = CircularAt(P, i + 1);

        //printf("p0 %f %f\n", p0(0), p0(1));
        //printf("p1 %f %f\n", p1(0), p1(1));
        //printf("s0 %f %f\n", s0(0), s0(1));
        //printf("s1 %f %f\n", s1(0), s1(1));
        if (LineUtils::intersect(p0, p1, s0, s1, p_t, s_t) &&
            p_t > -PF_EPS && p_t < 1. + PF_EPS && s_t > -PF_EPS && s_t < 1. + PF_EPS) {
            return false;
        }
    }

    return true;
}

bool are_points_facing_inside(
    const std::vector<const Eigen::Matrix2Xd*>& P,
    const int polygon_0,
    const int v0,
    const int polygon_1,
    const int v1
) {
    const vec2 p0 = P[polygon_0]->col(v0);
    const vec2 p1 = P[polygon_1]->col(v1);

    double p_t = INFINITY, s_t = INFINITY;
    for (size_t polygon = 0; polygon < P.size(); ++polygon) {
        for (Index i = 0; i < P[polygon]->cols(); ++i) {

            // Ignore immediately neighboring segments in the same polygon
            const int d0 = min(CircularDist(*P[polygon], v0, i), CircularDist(*P[polygon], i, v0));
            const int d1 = min(CircularDist(*P[polygon], v1, i), CircularDist(*P[polygon], i, v1));
            if ((polygon_0 == polygon && d0 <= 1) || (polygon_1 == polygon && d1 <= 1)) {
                continue;
            }

            const vec2 s0 = P[polygon]->col(i);
            const vec2 s1 = CircularAt(*P[polygon], i + 1);

            // There can be other polygons with the same points (we should not check for them in the first place...)
            if ((s0 - p0).squaredNorm() < PF_EPS || (s1 - p0).squaredNorm() < PF_EPS ||
                (s0 - p1).squaredNorm() < PF_EPS || (s1 - p1).squaredNorm() < PF_EPS) {
                continue;
            }

            //printf("p0 %f %f\n", p0(0), p0(1));
            //printf("p1 %f %f\n", p1(0), p1(1));
            //printf("s0 %f %f\n", s0(0), s0(1));
            //printf("s1 %f %f\n", s1(0), s1(1));
            if (LineUtils::intersect(p0, p1, s0, s1, p_t, s_t) &&
                p_t > -PF_EPS && p_t < 1. + PF_EPS && s_t > -PF_EPS && s_t < 1. + PF_EPS) {
                return false;
            }
        }
    }

    return true;
}

NAMESPACE_END(polyvec)