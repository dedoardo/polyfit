// polyvec
#include <polyvec/curve-tracer/curvature_variation_optimize_coordinate_descent.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>

using namespace std;
using namespace polyfit;
using namespace Eigen;

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(CurveTracer)

namespace {
    inline vec2 find_distances_jaklic (
        const double p1x, const double p1y,
        const double t1x, const double t1y,
        const double p2x, const double p2y,
        const double t2x, const double t2y) {

        //Jaklic et al., Curvature variation minimizing cubic Hermite interpolants
        //http://valjhun.fmf.uni-lj.si/~emil/clanki/curvatureAMC.pdf
        const double c01 = t1x * t2y - t1y * t2x;
        const double dpx = p2x - p1x;
        const double dpy = p2y - p1y;
        const double c0 = t1x * dpy - t1y * dpx;
        const double c1 = t2x * dpy - t2y * dpx;

        if (c01 == 0)
            return vec2::Zero();

        return vec2(-2 * c1 / c01 / 3, 2 * c0 / c01 / 3);
    }

    inline Eigen::Matrix<double, 8, 1> evaluate_derivatives (
        const double t, 
        const double p1x, const double p1y, 
        const double t1x, const double t1y, 
        const double p2x, const double p2y, 
        const double t2x, const double t2y) {

        const double t2 = t * t;
        const double t3 = t2 * t;
        const double invT = 1 - t;
        const double invT2 = invT * invT;
        const double invT3 = invT2 * invT;

        const double px = invT3 * p1x + 3 * invT2 * t * (p1x + t1x) + 3 * invT * t2 * (p2x - t2x) + t3 * p2x;
        const double py = invT3 * p1y + 3 * invT2 * t * (p1y + t1y) + 3 * invT * t2 * (p2y - t2y) + t3 * p2y;

        const double diff1x = 3 * (invT2 * t1x + 2 * invT * t * (-t2x - t1x + p2x - p1x) + t2 * t2x);
        const double diff1y = 3 * (invT2 * t1y + 2 * invT * t * (-t2y - t1y + p2y - p1y) + t2 * t2y);

        const double diff2x = 6 * (invT * (p2x - p1x - 2 * t1x - t2x) + t * (2 * t2x + t1x - p2x + p1x));
        const double diff2y = 6 * (invT * (p2y - p1y - 2 * t1y - t2y) + t * (2 * t2y + t1y - p2y + p1y));

        const double diff3x = 6 * (3 * t2x + 3 * t1x - 2 * p2x + 2 * p1x);
        const double diff3y = 6 * (3 * t2y + 3 * t1y - 2 * p2y + 2 * p1y);

        Matrix<double, 8, 1> derivatives;
        derivatives << px, py, diff1x, diff1y, diff2x, diff2y, diff3x, diff3y;
        return derivatives;
    }

    inline double curvature_variation_sq(
        const double t,
        const double p1x, const double p1y,
        const double t1x, const double t1y,
        const double p2x, const double p2y,
        const double t2x, const double t2y) {

        Matrix<double, 8, 1> diffs = evaluate_derivatives(t, p1x, p1y, t1x, t1y, p2x, p2y, t2x, t2y);

        double diff1x = diffs(2), diff1y = diffs(3);
        double diff2x = diffs(4), diff2y = diffs(5);
        double diff3x = diffs(6), diff3y = diffs(7);

        double diff1NormSq = diff1x * diff1x + diff1y * diff1y;
        double diff1Norm4 = diff1NormSq * diff1NormSq;
        double diff1Norm10 = diff1Norm4 * diff1Norm4 * diff1NormSq;

        double diff1CrossDiff3 = diff1x * diff3y - diff1y * diff3x;
        double diff1DotDiff2 = diff1x * diff2x + diff1y * diff2y;
        double diff1CrossDiff2 = diff1x * diff2y - diff1y * diff2x;
        double psi = diff1NormSq * diff1CrossDiff3 - 3 * diff1DotDiff2 * diff1CrossDiff2;

        return psi * psi / diff1Norm10;
    }

    inline double curvature_variation(
        const double p1x, const double p1y,
        const double t1x, const double t1y,
        const double p2x, const double p2y,
        const double t2x, const double t2y) {

        const double M = 20;
        const double h = 1 / (2 * M);

        double sum = 0;
        sum += h / 3 * (curvature_variation_sq(0, p1x, p1y, t1x, t1y, p2x, p2y, t2x, t2y) + curvature_variation_sq(1, p1x, p1y, t1x, t1y, p2x, p2y, t2x, t2y));
        for (int i = 1; i <= M - 1; ++i)
        {
            double t = 2 * i * h;
            sum += 2 * h / 3 * curvature_variation_sq(t, p1x, p1y, t1x, t1y, p2x, p2y, t2x, t2y);
        }
        for (int i = 1; i <= M - 1; ++i)
        {
            double t = 2 * (i - 1) * h;
            sum += 4 * h / 3 * curvature_variation_sq(t, p1x, p1y, t1x, t1y, p2x, p2y, t2x, t2y);
        }
        return sum;
    }
}

void minimize_curvature_variation_coordinate_descent(
    const double p1x, const double p1y,
    const double t1x, const double t1y,
    const double p2x, const double p2y,
    const double t2x, const double t2y,
    double& t1_scale, double& t2_scale
) {
    //Planar cubic G 1 and quintic G 2 Hermite interpolations via curvature variation minimization
    vec2 distances = find_distances_jaklic(p1x, p1y, t1x, t1y, p2x, p2y, t2x, t2y);

    // do coordinate descent
    double dpx = p2x - p1x;
    double dpy = p2y - p1y;
    double dp = sqrt(dpx * dpx + dpy * dpy);

    // mmmmmagic
    distances(0) = dp / 2;
    distances(1) = dp / 2;

    double lowerbound = 0.2 * dp;
    double upperbound = 5 * dp;

    for (int i = 0; i < 2; ++i) {
        if (distances(i) < lowerbound)
            distances(i) = lowerbound;
        if (distances(i) > upperbound)
            distances(i) = upperbound;
    }

    int descentDimension = 0;
    bool running = true;
    auto energy = [=] (const double d0, const double d1){ return curvature_variation(p1x, p1y, d0 * t1x, d0 * t1y, p2x, p2y, d1 * t2x, d1 * t2y); };
    double currentEnergy = energy(distances(0), distances(1));

    int iterations = 0;
    while (running) {
        ++iterations;
        double t = distances(descentDimension);
        double lastT = t;

        //line search in two directions
        for (int iSign = 0; iSign < 2; ++iSign) {
            int sign = iSign == 0 ? 1 : -1;
            double step = 1 * sign;
            while (abs(step) > 0.01) {
                double nextT = t + step;
        
                double nextE;
                if (descentDimension == 0) {
                    nextE = energy(nextT, distances(1));
                } else {
                    nextE = energy(distances(0), nextT);
                }

                //printf("next energy before %f < %f\n", nextE, currentEnergy);
                if (nextT >= lowerbound && nextT <= upperbound && nextE < currentEnergy) {
                    currentEnergy = nextE;
                    t = nextT;
                    //printf("next energy %f\n", nextE);
                }
                else {
                    step *= 0.5;
                }
            }
            if (t != lastT)
                break;
        }

        distances(descentDimension) = t;
        descentDimension = (descentDimension + 1) % 2;

        if (t < lowerbound)
            t = lowerbound;
        if (t > upperbound)
            t = upperbound;

        running = t != lastT || iterations == 1;
    }

    t1_scale = distances(0);
    t2_scale = distances(1);
}

void minimize_curvature_variation_coordinate_descent(
    polyvec::BezierCurve& bezier
) {
    mat24 control_points = bezier.get_control_points();
    vec2 p0 = control_points.col(0);
    vec2 p3 = control_points.col(3);
    vec2 t0 = bezier.dposdt(0.).normalized();
    vec2 t3 = bezier.dposdt(1.).normalized();

    double t0_scale, t3_scale;
    minimize_curvature_variation_coordinate_descent (
        p0(0), p0(1), t0(0), t0(1), p3(0), p3(1), t3(0), t3(1), t0_scale, t3_scale
    );

    // For t=0 and t=1 solve for P2 and P3 to obtain the second and third control points
    control_points.col(1) = p0 + t0 * t0_scale;
    control_points.col(2) = p3 - t3 * t3_scale;
    bezier.set_control_points(control_points);
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvec)