/*
    Simple app to debug/draw individual curves
*/
// Polyvec
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curvature_variation_optimize_coordinate_descent.hpp>
#include <polyvec/io/pdf.hpp>
#include "../dev/drawing.hpp"

// libc++
#include <cstdlib>

using namespace std;
using namespace polyvec;
using namespace polyfit;

#define CANVAS_SIZE 32
#define CURVE_STYLE 2.5, colors::turtle_purple
#define CONTROL_POLYGON_STYLE Style::outline(colors::talking_orange, 2.5)
#define WRITE_FILE(uri) "D:/data/polyvec/out/" uri ".svg"

double evaluate_curvature_variation(BezierCurve& curve) {
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

void draw_curve(const char* uri, GlobFitCurve& c) {
    DevicePDF* pdf = new DevicePDF(uri, 1, 1);
    
    mat24 B;
    B.col(0) = vec2(0., 0.);
    B.col(1) = vec2(CANVAS_SIZE, 0.);
    B.col(2) = vec2(CANVAS_SIZE, CANVAS_SIZE);
    B.col(3) = vec2(0., CANVAS_SIZE);

    draw_raster_background(B, Style::outline(colors::gray * 1.75, 2.5));
    if (c.get_type() == GLOBFIT_CURVE_BEZIER) {
        BezierCurve& bez = dynamic_cast<BezierCurve&>(c);
        draw::line(bez.get_control_points().col(0), bez.get_control_points().col(1), CONTROL_POLYGON_STYLE);
        draw::line(bez.get_control_points().col(1), bez.get_control_points().col(2), CONTROL_POLYGON_STYLE);
        draw::line(bez.get_control_points().col(2), bez.get_control_points().col(3), CONTROL_POLYGON_STYLE);
    }

    const double k_variation = evaluate_curvature_variation(dynamic_cast<BezierCurve&>(c));
    draw::text(c.pos(.5), to_string(k_variation), draw::font_pdf, Style::text());

    draw::curve(&c, CURVE_STYLE);
    pdf->draw(0, 0);
    delete pdf;
}

int main(int argc, char* argv) {
    BezierCurve bez;
    mat24 cp;

    {
        cp.col(0) << 27.448683, 9.816228;
        cp.col(1) << 24.602633, 8.867544;
        cp.col(2) << 24.288854,  10.394427;
        cp.col(3) << 21.605573, 9.052786;

        bez.set_control_points(cp);
        draw_curve(WRITE_FILE("bez00"), bez);
        polyvec::CurveTracer::minimize_curvature_variation_coordinate_descent(bez);
        draw_curve(WRITE_FILE("bez01"), bez);
    }

    {
        cp.col(0) << 21.605573, 18.947214;
        cp.col(1) << 24.288854, 17.605573;
        cp.col(2) << 24.602633, 19.132456;
        cp.col(3) << 27.448683, 18.183772;

        bez.set_control_points(cp);
        draw_curve(WRITE_FILE("bez10"), bez);
        polyvec::CurveTracer::minimize_curvature_variation_coordinate_descent(bez);
        draw_curve(WRITE_FILE("bez11"), bez);
    }


    return EXIT_SUCCESS;
}