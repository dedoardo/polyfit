#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curvature_variation_optimize_coordinate_descent.hpp>
#include <polyvec/io/pdf.hpp>
#include <polyvec/curve-tracer/curve_objectives.hpp>
#include "drawing.hpp"

#include <memory>

using namespace std;
using namespace polyvec;
using namespace polyfit;

namespace {
    void print_curve(const BezierCurve& curve) {
        const mat24& cp = curve.get_control_points();
        printf("control points\n");
        printf("P0: %.3f %.3f\n", cp(0, 0), cp(1, 0));
        printf("P1: %.3f %.3f\n", cp(0, 1), cp(1, 1));
        printf("P2: %.3f %.3f\n", cp(0, 2), cp(1, 2));
        printf("P3: %.3f %.3f\n", cp(0, 3), cp(1, 3));

        const vec2 t0 = curve.dposdt(0.);
        printf("T0: %.3f %.3f\n", t0(0), t0(1));

        const vec2 t3 = curve.dposdt(1.);
        printf("T3: %.3f %.3f\n", t3(0), t3(1));
    }

    void write_curve(BezierCurve& curve, const char* uri) {
        DevicePDF* pdf = new DevicePDF(uri, 1, 1);
        draw::curve(&curve, 2.5, colors::talking_orange);

        const mat24& cp = curve.get_control_points();
        mat2x B(2, 4);
        B.col(0) = vec2(0., 0.);
        B.col(1) = vec2(10., 0.);
        B.col(2) = vec2(10., 10.);
        B.col(3) = vec2(0., 10.);
        draw_raster_background(B, Style::outline(colors::gray * 1.75, 2.5));
        draw::line(cp.col(0), cp.col(1), Style::outline(colors::red, 1.5));
        draw::line(cp.col(1), cp.col(2), Style::outline(colors::red, 1.5));
        draw::line(cp.col(2), cp.col(3), Style::outline(colors::red, 1.5));
        pdf->draw(0, 0);
        delete pdf;
    }

    double print_total_curvature_variation(BezierCurve& curve) {
        std::shared_ptr<BezierCurve> curve_ptr = make_shared<BezierCurve>(curve);
        GlobFitBezierAngleBasedParametrization curve_param(curve_ptr);
        GlobFitObjective_CurvatureVariation opt_obj;
        opt_obj.set_params(&curve_param);

        vecXd obj;
        matXd jac;
        opt_obj.compute_objective_and_jacobian(obj, jac);
        printf("total curvature variation %f\n", obj.norm());
        return obj.norm();
    }
}

int main() {
    BezierCurve curve;
    mat24 control_points;
    control_points.col(0) = vec2(0., 0.);
    control_points.col(1) = vec2(0., 5.);
    control_points.col(2) = vec2(1., 5.);
    control_points.col(3) = vec2(1., 0.);
    curve.set_control_points(control_points);

    print_curve(curve);
    const double k_variation_before = print_total_curvature_variation(curve);
    write_curve(curve, "D:/data/polyvec/out/before_optimization.svg");

    polyvec::CurveTracer::minimize_curvature_variation_coordinate_descent(curve);
    printf("\nafter optimization\n");

    print_curve(curve);
    const double k_variation_after = print_total_curvature_variation(curve);
    write_curve(curve, "D:/data/polyvec/out/after_optimization.svg");

    if (k_variation_before >= k_variation_after - PF_EPS) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}