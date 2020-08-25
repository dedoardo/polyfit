#pragma once

#include <Eigen/Core>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/utils/matrix.hpp>

namespace draw_paper_figures {
    struct Data {
        const Eigen::Matrix2Xd* B;
        const Eigen::Matrix2Xd* PP;
        const Eigen::VectorXi*  P;
        const std::vector<polyfit::BoundaryGraph::Edge>* E;
        const std::vector<polyvec::CurvePrimitive>* curves;
        const std::vector<polyvec::TangentFitType>* tangent_fits;
        const std::vector<polyvec::AttemptInfo>* curves_attempts;
		Eigen::Vector4d color = Eigen::Vector4d(0, 0, 0, 1);
    };

    inline void grid(const Data& data) {
        const auto& raster = *data.B;
        const auto& style = Style::outline(colors::gray * 1.75, .05);

        Eigen::Vector2d min, max;
        min.setConstant(std::numeric_limits<double>::infinity());
        max = -min;

        for (int i = 0; i < raster.cols(); ++i) {
            for (int j = 0; j < 2; ++j) {
                if (raster.coeff(j, i) < min(j)) {
                    min(j) = raster.coeff(j, i);
                }

                if (raster.coeff(j, i) > max(j)) {
                    max(j) = raster.coeff(j, i);
                }
            }
        }

        min.x() = std::floor(min.x() + 0.5 - 1);
        min.y() = std::floor(min.y() + 0.5 - 1);
        max.x() = std::ceil(max.x() + 0.5 + 1);
        max.y() = std::ceil(max.y() + 0.5 + 1);

        for (int i = (int)min(0); i <= (int)max(0); ++i) {
            draw::line(Eigen::Vector2d((double)i - 0.5, min(1) - 0.5), Eigen::Vector2d((double)i - 0.5, max(1) - 0.5), style);
        }

        for (int i = (int)min(1); i <= (int)max(1); ++i) {
            draw::line(Eigen::Vector2d(min(0) - 0.5, (double)i - 0.5), Eigen::Vector2d(max(0) - 0.5, (double)i - 0.5), style);
        }
    }

    inline void raster_boundary(const Data& data) {
        const auto& raster = *data.B;
        std::vector<real2> points;

        for (int i = 0; i < raster.cols(); ++i) {
            points.emplace_back(raster.col(i));
        }

        draw::polygon(points, Style::fill(real3(61.f/255, 104.f/255, 174.f/255)));
    }

    inline void graph(const Data& data) {
        const auto& B = *data.B;
        const auto& E = *data.E;

        for (size_t i = 0; i < E.size(); ++i) {
            const auto p0 = B.col(E[i].v0);
            const auto p1 = B.col(E[i].v1);
            draw::line(p0, p1, Style::outline(colors::forest_green, .03));
        }
    }

    inline void polygon(const Data& data) {
        draw::polyline(*data.PP, Style::outline(colors::black, 0.075));
    }

    inline void polygon_curves(const Data& data) {
        const auto& PP = *data.PP;
        for (Eigen::Index i = 0; i < PP.cols(); ++i) {
            const auto pp = PP.col(i);
            const auto pn = polyfit::CircularAt(PP, i + 1);
            const auto t = pn - pp;
            Eigen::Vector2d n(-t.y(), t.x());
            n.normalize();
            const auto mid = .5 * (pp + pn);
            draw::line(mid, mid + n * .2, Style::outline(colors::black, .1));
            draw::line(mid, mid - n * .2, Style::outline(colors::black, .1));
        }

        polygon(data);
    }

    inline void curve_primitives_classification(const Data& data) {
        using namespace polyfit;
        const auto& PP = *data.PP;
        const auto& B = *data.B;

        for (Eigen::Index i = 0; i < PP.cols(); ++i) {
            const int tangent_fit = data.tangent_fits->at(i);

            //if (tangent_fit == TANGENT_FIT_CONSTANT) {
            //    const auto pp = .5 * (CircularAt(PP, i - 1) + CircularAt(PP, i));
            //    const auto pn = .5 * (CircularAt(PP, i) + CircularAt(PP, i + 1));
            //    draw::line(CircularAt(PP, i), pp, Style::outline(colors::red, .1));
            //    draw::line(CircularAt(PP, i), pn, Style::outline(colors::red, .1));
            //    draw::point(pp, .1, Style::fill(colors::black));
            //    draw::point(pn, .1, Style::fill(colors::black));
            //    draw::point(CircularAt(PP, i), .1, Style::fill(colors::black));
            //}else {
                const auto& curves = data.curves_attempts->at(tangent_fit).attempt[0].primitive_seq.primitives;
                for (size_t j = 0; j < curves.size(); ++j) {
                    if (curves[j].corner == i) {
                        switch (curves[j].curve->get_curve()->get_type()) {
                        case GLOBFIT_CURVE_BEZIER:
                            draw::curve(curves[j].curve->get_curve().get(), .1, colors::calm_blue);
                            break;
                        case GLOBFIT_CURVE_LINE:
                            draw::curve(curves[j].curve->get_curve().get(), .1, colors::red);
                            break;
                        }

                        draw::point(curves[j].curve->get_curve()->pos(0.), .1, Style::fill(colors::black));
                        draw::point(curves[j].curve->get_curve()->pos(1.), .1, Style::fill(colors::black));
                    }
                }
            //}
        }

        //draw_curve_primitives(curves, PrimitiveTypeColorFunctor());

        //for (size_t i = 0; i < curves.size(); ++i) {
            //draw::point(curves[i].curve->get_curve()->pos(1.), .1, Style::fill(colors::black));
        //}
    }

    inline void curve_primitives(const Data& data) {
        draw_curve_primitives(*data.curves, PrimitiveTypeColorFunctor());

        for (size_t i = 0; i < data.curves->size(); ++i) {
            draw::point(data.curves->at(i).curve->get_curve()->pos(1.), .1, Style::fill(colors::black));
        }
    }

    inline void curve_fill_classification(const Data& data) {
        draw_curve_primitives_closed(data.curves_attempts->at(0).attempt[0].primitive_seq.primitives, colors::black);
    }

    inline void curve_fill(const Data& data) {
        draw_curve_primitives_closed(*data.curves, data.color);
    }
}