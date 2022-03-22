#ifndef POLYVEC_USE_POTRACE
#define POLYVEC_USE_POTRACE
#endif

#include "common_includes.hpp"

#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>
#include <polyvec/io/pdf.hpp>

#include <algorithm>

#define RASTER_STYLE Style::fill ( colors::black )
#define BOUNDARY_STYLE Style::outline ( colors::gray, 5. )
#define PATH_STYLE Style::outline(colors::black, 1. )
#define TANGENT_STYLE Style::outline(colors::turtle_purple, .1 )
#define TANGENT_PROJECT_STYLE Style::outline(colors::turtle_purple, .05 )
#define CURVATURE_STYLE Style::outline(colors::calm_blue, 3. )
#define CURVATURE_STYLE_CLAMPED Style::outline(colors::royal_blue, 3. )

#define ERROR_MIN_COLOR colors::forest_green
#define ERROR_MAX_COLOR colors::red

#define ACCURACY_LEGEND_LEN 20

namespace polyvec {
    inline void draw_triplet ( const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const real3& color, const double thickness ) {
        draw::point ( p0, thickness * 2, Style::fill ( colors::forest_green ) );
        draw::point ( p1, thickness * 2, Style::fill ( colors::forest_green ) );
        draw::point ( p2, thickness * 2, Style::fill ( colors::forest_green ) );
        draw::line ( p0, p1, Style::outline ( color, thickness ) );
        draw::line ( p1, p2, Style::outline ( color, thickness ) );
    }

    inline void draw_curve_primitives_accuracy ( const Eigen::MatrixXd& B, const std::vector<polyfit::BoundaryGraph::Edge>& E,
            const std::vector<polyfit::Vertex>& VE, std::vector<CurvePrimitive>& primitives ) {
        using namespace polyfit;
        using namespace std;

        for ( int i = 0; i < primitives.size(); ++i ) {
            auto& c = primitives[i].curve;

            if ( primitives[i].corner == -1 ) {
                PF_ASSERT ( c->get_curve()->get_type() == GLOBFIT_CURVE_LINE );
                draw::line ( c->get_curve()->pos ( 0. ), c->get_curve()->pos ( 1. ), Style::outline ( colors::forest_green, 1. ) );
                continue;
            }

            double e_edge = E[VE[primitives[i].corner]].d_max;
            double e_curve = primitives[i].error.accuracy.max_error();
            double e_max = .75;

            double delta = min ( e_max, abs ( e_curve - e_edge ) );

            PF_LOGF ( "corner %d error-edge %f error-curve %f", primitives[i].corner, e_edge, e_curve );

            vec3 color_min = colors::forest_green;
            vec3 color_max = colors::red;

            const string text = misc::sfmt ( "(%d) poly-d-max %.3f curve-d-max %.3f curve-d-acc %.3f", primitives[i].corner, e_edge, e_curve,
                                             primitives[i].error.accuracy );

            vec3 color = misc::lerp ( color_min, color_max, delta / e_max );
            draw::curve ( c->get_curve().get(), 5., color );
            draw::text ( c->get_curve()->pos ( .35 ), text, draw::font_pdf / 3, Style::text() );
        }
    }

    inline void draw_curve_primitives_curvature ( std::vector<CurvePrimitive>& prims ) {
        using namespace Eigen;

        const double min_curvature_radius = VectorOptions::get()->spline_fit_curvature_g0_thr;
        const double max_curvature_vis = 2.;

        for ( int i = 0; i < prims.size(); ++i ) {
            auto& curve = prims[i].curve;

            draw::curve ( curve->get_curve().get(), 5., colors::black );

            if ( curve->get_curve()->get_type() == GLOBFIT_CURVE_LINE ) {
                continue;
            }

            const double t_step = .025;
            double t = t_step;
            double max_curvature_radius_vis = 10.;

            double k_min = INFINITY;
            double k_max = -INFINITY;

            while ( t < Curve::END + constants::Eps ) {
                t = misc::saturate ( t );
                const Vector2d pos = curve->get_curve()->pos ( t );
                Vector2d normal = curve->get_curve()->dposdt ( t ).normalized();
                normal = Vector2d ( -normal.y(), normal.x() );
                double k;
                SmoothCurveUtil::curvature_eval ( curve->get_curve()->dposdt ( t ), curve->get_curve()->dposdtdt ( t ), k );
                k = std::min ( abs ( k ), max_curvature_radius_vis );

                const double r = 1. / k;

                real3 color;

                if ( r < min_curvature_radius ) {
                    color = colors::red;
                } else {
                    color = misc::lerp ( ERROR_MAX_COLOR, ERROR_MIN_COLOR, r / max_curvature_vis );
                }

                k_min = std::min ( r, k_min );
                k_max = std::max ( r, k_max );

                draw::line ( pos, pos + normal * std::min ( r, 5. ), Style::outline ( color, .75 ) );
                t += t_step;
            }

            const str text = misc::sfmt ( "max %.3f min %.3f", k_max, k_min );
            draw::text ( curve->get_curve()->pos ( .35 ), text, draw::font_pdf / 3, Style::text() );
        }
    }

    inline void draw_curve_primitives_closed ( const std::vector<CurvePrimitive>& primitives, const Eigen::Vector4d& color) {
        std::vector<GlobFitCurve*> curves;

        for ( const CurvePrimitive& prim : primitives ) {
            curves.emplace_back ( prim.curve->get_curve().get() );
        }

        draw::curves_fill ( curves, color );
    }

    inline void draw_curve_primitives_closed(const std::vector<CurvePrimitive>& primitives, const Eigen::Vector3d& color = colors::black) {
        draw_curve_primitives_closed(primitives, Eigen::Vector4d(color(0), color(1), color(2), 1.0));
    }

    inline void draw_pdf ( const str& uri, const std::function<void() >& fn ) {
        DevicePDF* pdf = new DevicePDF ( uri.c_str(), 1, 1 );
        fn();
        pdf->draw ( 0, 0 );
        delete pdf;
    }

    inline void draw_curve_primitives_corner_error ( const Eigen::Matrix2Xd& P, const std::vector<CurvePrimitive>& C ) {
        using namespace polyfit;
        using namespace Eigen;

        const double max_error_accuracy = VectorOptions::get()->spline_fit_distance_point_g0_thr;
        ( void ) max_error_accuracy;
        const double min_curvature_radius = VectorOptions::get()->spline_fit_curvature_g0_thr;
        const double max_curvature_vis = 3.;

        for ( int i = 0; i < ( int ) C.size(); ++i ) {
            draw::curve ( C[i].curve->get_curve().get(), 1., i % 2 ? colors::forest_green : colors::talking_orange );
            //draw::text(C[i].curve->pos(.5), std::to_string(i), draw::font_pdf, Style::text());

            PF_LOGF ( "Curve %d corner %d", i, C[i].corner );
            PF_LOGF ( "accuracy pos %f neg %f tot %f", C[i].error.accuracy.e_pos, C[i].error.accuracy.e_neg, C[i].error.accuracy );

            if ( C[i].corner != -1 ) {
                int c = C[i].corner;
                vec2 pp = CircularAt ( P, c - 1 );
                vec2 p = CircularAt ( P, c );
                vec2 pn = CircularAt ( P, c + 1 );

                //draw::line(p, p + misc::lerp(p, pp, .1), Style::outline(colors::red, 1.));
                //draw::line(p, p + misc::lerp(p, pn, .1), Style::outline(colors::red, 1.));
                double a = acos ( ( pp - p ).normalized().dot ( ( pn - p ).normalized() ) );
                PF_LOGF ( "points %f %f %f %f %f %f", pp.x(), pp.y(), p.x(), p.y(), pn.x(), pn.y() );
                PF_LOGF ( "angle %f / %f ", a, PF_RAD ( 135 ) );

                if ( a < PF_RAD ( 135 ) ) {
                    GlobFitCurve* curve = C[i].curve->get_curve().get();

                    const double t_step = .025;
                    double t = t_step;
                    double max_curvature_radius_vis = 10.;

                    while ( t < Curve::END + constants::Eps ) {
                        t = misc::saturate ( t );
                        const Vector2d pos = curve->pos ( t );
                        Vector2d normal = curve->dposdt ( t ).normalized();
                        normal = Vector2d ( -normal.y(), normal.x() );
                        double k;
                        SmoothCurveUtil::curvature_eval ( curve->dposdt ( t ), curve->dposdtdt ( t ), k );
                        k = std::min ( abs ( k ), max_curvature_radius_vis );

                        const double r = 1. / k;

                        real3 color;

                        if ( r < min_curvature_radius ) {
                            color = colors::red;
                        } else {
                            color = misc::lerp ( ERROR_MAX_COLOR, ERROR_MIN_COLOR, r / max_curvature_vis );
                        }

                        draw::line ( pos, pos + normal * std::min ( r, 5. ), Style::outline ( color, .45 ) );
                        t += t_step;
                    }

                }
            }

            draw::text ( C[i].curve->get_curve()->pos ( .5 ), misc::sfmt ( "accuracy pos %.3f neg %.3f", C[i].error.accuracy.e_pos, C[i].error.accuracy.e_neg ),
                         draw::font_pdf / 2, Style::text() );
        }
    }

    inline void draw_curve_indices ( const std::vector<CurvePrimitive>& primitives ) {
        for ( int i = 0; i < primitives.size(); ++i ) {
            const auto& c = primitives[i].curve.get();
            draw::text ( ( c->get_curve()->pos ( 0. ) + c->get_curve()->pos ( 1. ) ) / 2, std::to_string ( i ), draw::font_pdf, Style::text() );
        }
    }

    struct AlternatingColorFunctor {
        real3 operator() ( int primitiveId, const CurvePrimitive& prim ) {
            return primitiveId % 2 ? colors::forest_green : colors::talking_orange;
        }
    };

    struct AlternatingCornerColorFunctor {
        real3 operator() ( int primitiveId, const CurvePrimitive& prim ) {
            return prim.corner % 2 ? colors::forest_green : colors::talking_orange;
        }
    };

    struct ErrorHighlightFunctor {
        real3 operator() ( int primitiveId, const CurvePrimitive& prim ) {
            const double max_error_accuracy = VectorOptions::get()->spline_fit_distance_point_g0_thr;
            const double min_curvature_radius = VectorOptions::get()->spline_fit_curvature_g0_thr;

            bool accuracyValid = prim.error.accuracy.max_error() < max_error_accuracy;
            bool curvatureValid = prim.error.curvature.r_min > min_curvature_radius;

            if ( !accuracyValid && !curvatureValid ) {
                return colors::red;
            } else if ( !accuracyValid ) {
                return colors::fuchsia;
            } else if ( !curvatureValid ) {
                return colors::blue;
            } else {
                return colors::forest_green;
            }
        }
    };

    struct TangentFitColorFunctor {
		const real3 colors[TANGENT_FIT_SAMPLES_COUNT] = {
			colors::forest_green,
			colors::talking_orange,
			colors::red
		};

		const std::vector<TangentFitType>& tangent_fits;

		TangentFitColorFunctor(const std::vector<TangentFitType>& tangent_fits) :
			tangent_fits(tangent_fits) { }

        real3 operator() ( int primitiveId, const CurvePrimitive& prim ) {
			PF_ASSERT(prim.corner >= 0 && prim.corner < (int)tangent_fits.size());
			return colors[(int)tangent_fits[prim.corner]];
        }
    };

    struct ConstantColorFunctor {
        const real3 color;

        ConstantColorFunctor(const real3 color) : color(color) { }

        real3 operator() (int primitiveId, const CurvePrimitive& prim) {
            return color;
        }
    };

    struct PrimitiveTypeColorFunctor {
        PrimitiveTypeColorFunctor() { }

        real3 operator() (int primitiveId, const CurvePrimitive& prim) {
            switch (prim.curve->get_curve()->get_type()) {
            case GLOBFIT_CURVE_BEZIER:
                return real3(62.f/255, 106.f/255, 175.f/255);
            case GLOBFIT_CURVE_LINE:
                return real3(122.f/255, 181.f/255, 66.f/255);
            default:
                return colors::black;
            }
        }
    };

//TColorFunctor: real3(int primitiveId, const CurvePrimitive& prim)
    template <typename TColorFunctor>
    inline void draw_curve_primitives ( const std::vector<CurvePrimitive>& primitives, TColorFunctor&& color, bool visualize_parametrization = false) {
        using namespace Eigen;

        const double max_curvature_vis = 3.;

        for ( int i = 0; i < ( int ) primitives.size(); ++i ) {
            const auto& prim = primitives[i];

            //draw::text(primitives[i].curve->get_curve()->pos(.5), std::to_string(i), draw::font_pdf, Style::text());
            auto this_color = std::forward<TColorFunctor> ( color ) ( i, prim );
            draw::curve ( primitives[i].curve->get_curve().get(), 1.2, this_color );

			if (visualize_parametrization) {
				draw::point(prim.curve->get_curve()->pos(0.25), 0.4, Style::fill(colors::red));
				draw::point(prim.curve->get_curve()->pos(0.50), 0.4, Style::fill(colors::red));
				draw::point(prim.curve->get_curve()->pos(0.75), 0.4, Style::fill(colors::red));
			}

#if 0 // curvature
            GlobFitCurve* curve = primitives[i].curve.get();

            const double t_step = .025;
            double t = t_step;
            double max_curvature_radius_vis = 10.;

            while ( t < Curve::END + constants::Eps ) {
                t = misc::saturate ( t );
                const Vector2d pos = curve->pos ( t );
                Vector2d normal = curve->dposdt ( t ).normalized();
                normal = Vector2d ( -normal.y(), normal.x() );
                double k;
                SmoothCurveUtil::curvature_eval ( curve->dposdt ( t ), curve->dposdtdt ( t ), k );
                k = std::min ( abs ( k ), max_curvature_radius_vis );

                const double r = 1. / k;

                real3 color;

                if ( r < min_curvature_radius ) {
                    color = colors::red;
                } else {
                    color = misc::lerp ( ERROR_MAX_COLOR, ERROR_MIN_COLOR, r / max_curvature_vis );
                }

                draw::line ( pos, pos + normal * std::min ( r, 5. ), Style::outline ( color, .1 ) );
                t += t_step;
            }

#endif

#if 0 // bezier control points

            if ( primitives[i].curve->get_type() == GLOBFIT_CURVE_BEZIER ) {
                auto c0 = primitives[i].curve->as_bezier()->get_bezier().get_control_points().col ( 0 );
                auto c1 = primitives[i].curve->as_bezier()->get_bezier().get_control_points().col ( 1 );
                auto c2 = primitives[i].curve->as_bezier()->get_bezier().get_control_points().col ( 2 );
                auto c3 = primitives[i].curve->as_bezier()->get_bezier().get_control_points().col ( 3 );

                draw::line ( c0, c1, Style::outline ( colors::talking_orange, .05 ) );
                draw::line ( c1, c2, Style::outline ( colors::talking_orange, .05 ) );
                draw::line ( c2, c3, Style::outline ( colors::talking_orange, .05 ) );
            }

#endif

#if 0 // fitting data

            //if (prim.fit_midpoints && i == primitives.size() - 1) {
            if ( prim.fitting_info ) {
                for ( int j = 0; j < ( int ) prim.fitting_info->fit_midpoints.size(); ++j ) {
                    draw::point ( prim.fitting_info->fit_midpoints.at ( j ), .1, Style::fill ( i % 2 ? colors::forest_green : colors::talking_orange ) );
                }
            }

            //}

            if ( prim.fitting_info && ( i == ( int ) primitives.size() - 1 ) ) {
                for ( int j = 0; j < ( int ) prim.fitting_info->fit_tangents.size(); ++j ) {
                    //if (prim.fit_tangent_points) {
                    //draw::point(prim.fit_tangent_points->at(j), .02, Style::fill(colors::forest_green));
                    //}
                    //draw::line(prim.fit_tangent_points->at(j), prim.fit_tangent_points->at(j) + prim.fit_tangents->at(j), Style::outline(colors::calm_blue, .1));
                }
            }

#endif
        }
    }

    inline void draw_sequence_error_tangents ( CurvePrimitiveSequence& seq ) {
        using namespace Eigen;

        for ( Index i = 0; i < ( Index ) seq.primitives.size(); ++i ) {
            CurvePrimitive& primitive = seq.primitives[i];

            draw::curve ( primitive.curve->get_curve().get(), .1, i % 2 ? colors::forest_green : colors::talking_orange );

            double max_tangent_error = 0.1;

            for ( Index j = 0; j < ( int ) primitive.fitting_info.dense_tangents.fit_tangents.size(); ++j ) {
                const Vector2d tangent_target = primitive.fitting_info.dense_tangents.fit_tangents.at ( j );
                const double tangent_t = misc::saturate ( primitive.fit_tangents_t[j] );
                Vector2d tangent_curve;
                SmoothCurveUtil::tangent_eval ( primitive.curve->get_curve()->dposdt ( tangent_t ), tangent_curve );

                const Vector2d tangent_error = ( tangent_target - tangent_curve ).cwiseAbs();
                //dbg::info(FMT("Tangent error %.5f %.5f", tangent_error.x(), tangent_error.y()));

                const Vector3d error_color = misc::lerp ( ERROR_MIN_COLOR, ERROR_MAX_COLOR, misc::saturate ( tangent_error.maxCoeff() / max_tangent_error ) );
                const Vector2d curve_pt = primitive.curve->get_curve()->pos ( tangent_t );

                Vector2d dir = tangent_target.normalized();
                dir = Vector2d ( -dir.y(), dir.x() );
                draw::line ( curve_pt, curve_pt + dir, Style::outline ( error_color, .01 ) );
            }

            //draw::curve(primitive.curve.get(), .025, misc::lerp(ERROR_MIN_COLOR, ERROR_MAX_COLOR, error_t));
        }
    }

    inline void draw_text ( int row, const str& text, const real3& color = colors::black ) {
        const double h = ( draw::font_pdf ) * 16;

        //draw::box ( { 0., 0., 128., 32. }, colors::black );
        draw::text ( {5., h * row + 2.5 }, text, h, Style::text ( color ) );
    }

    inline void draw_raster_closed ( const Eigen::Matrix2Xd& raster, const real4& color ) {
        std::vector<real2> points;

        for ( int i = 0; i < raster.cols(); ++i ) {
            points.emplace_back ( raster.col ( i ) );
        }

        draw::polygon ( points, Style::fill ( color ) );
    }


    inline void draw_raster_closed ( const Eigen::Matrix2Xd& raster, const real3& color = colors::black ) {
        std::vector<real2> points;

        for ( int i = 0; i < raster.cols(); ++i ) {
            points.emplace_back ( raster.col ( i ) );
        }

        draw::polygon ( points, Style::fill ( color ) );
    }

    inline void draw_raster_dashed(const Eigen::Matrix2Xd& raster, const real4& color) {
        for (int i = 0; i < raster.cols(); ++i) {
            // dbg::info(FMT("Point %f %f", raster.col(i), raster.col((i + 1) % raster.cols())));
            draw::line(raster.col(i), raster.col((i + 1) % raster.cols()), Style::outline(color.segment(0, 3), 2., LineType::Dash));
        }
    }
    
    inline void draw_raster ( const Eigen::Matrix2Xd& raster, const real4& color ) {
        for ( int i = 0; i < raster.cols(); ++i ) {
            // dbg::info(FMT("Point %f %f", raster.col(i), raster.col((i + 1) % raster.cols())));
            draw::line ( raster.col ( i ), raster.col ( ( i + 1 ) % raster.cols() ), Style::outline ( color.segment ( 0, 3 ), 6. ) );
        }
    }
    
    inline void draw_raster(const Eigen::Matrix2Xd& raster, const real3& color) {
        real4 color4;
        color4.segment(0, 3) = color;
        color4(3) = 1.0;
        draw_raster(raster, color4);
    }

    inline void draw_raster ( const Eigen::Matrix2Xd& raster, const Style& style = BOUNDARY_STYLE ) {
		draw::polyline(raster, style);
    }

    inline void draw_raster_background ( const Eigen::Matrix2Xd& raster, const Style& style, bool dots = true ) {
        Eigen::Vector2d min, max;
        min.setConstant ( std::numeric_limits<double>::infinity() );
        max = -min;

        for ( int i = 0; i < raster.cols(); ++i ) {
            for ( int j = 0; j < 2; ++j ) {
                if ( raster.coeff ( j, i ) < min ( j ) ) {
                    min ( j ) = raster.coeff ( j, i );
                }

                if ( raster.coeff ( j, i ) > max ( j ) ) {
                    max ( j ) = raster.coeff ( j, i );
                }
            }
        }

        min.x() = std::floor ( min.x() + 0.5 - 1 );
        min.y() = std::floor ( min.y() + 0.5 - 1 );
        max.x() = std::ceil ( max.x() + 0.5 + 1 );
        max.y() = std::ceil ( max.y() + 0.5 + 1 );

        for ( int i = ( int ) min ( 0 ); i <= ( int ) max ( 0 ); ++i ) {
            draw::line ( Eigen::Vector2d ( ( double ) i - 0.5, min ( 1 ) - 0.5 ), Eigen::Vector2d ( ( double ) i - 0.5, max ( 1 ) - 0.5 ), style );
        }

        for ( int i = ( int ) min ( 1 ); i <= ( int ) max ( 1 ); ++i ) {
            draw::line ( Eigen::Vector2d ( min ( 0 ) - 0.5, ( double ) i - 0.5 ), Eigen::Vector2d ( max ( 0 ) - 0.5, ( double ) i - 0.5 ), style );
        }

		if(dots)
			for ( int i = min ( 0 ); i < max ( 0 ); ++i )
				for ( int j = min ( 1 ); j < max ( 1 ); ++j ) {
					draw::point ( Eigen::Vector2d ( i, j ), 0.03, Style::fill ( style.line_color ) );
				}
    }

    inline void draw_points ( const Eigen::Matrix2Xd& P ) {
        for ( int i = 0; i < P.cols(); ++i ) {
            draw::point ( P.col ( i ), 1., Style::fill ( colors::black ) );
        }
    }

    inline void draw_raster_path ( const Eigen::Matrix2Xd& raster, const real3& color = colors::black, const double w = .5 ) {
        for ( int i = 0; i < raster.cols(); ++i ) {
            // dbg::info(FMT("Point %f %f", raster.col(i), raster.col((i + 1) % raster.cols())));
            draw::line ( raster.col ( i ), raster.col ( ( i + 1 ) % raster.cols() ), Style::outline ( color, w ) );
        }
    }

    inline void draw_accuracy_legend() {
        const double t_step = .01;
        double t = t_step;

        while ( t <= 1. + constants::Eps ) {
            real3 color = misc::lerp ( ERROR_MIN_COLOR, ERROR_MAX_COLOR, t );
            real2 pos_prev = misc::lerp ( real2 ( 0., 0. ), real2 ( ACCURACY_LEGEND_LEN, 0. ), t - t_step );
            real2 pos = misc::lerp ( real2 ( 0., 0. ), real2 ( ACCURACY_LEGEND_LEN, 0. ), t );

            draw::line ( pos_prev, pos, Style::outline ( color, 1. ) );

            t += t_step;
        }
    }

    inline void draw_curves_curvature ( const std::vector<CurvePrimitive>& curves ) {
        using namespace Eigen;

        for ( int i = 0; i < curves.size(); ++i ) {
            auto curve = curves[i].curve->get_curve();

            const double t_step = .025;
            double t = t_step;
            double max_r = 20.;

            while ( t < Curve::END + constants::Eps ) {
                t = misc::saturate ( t );
                const Vector2d pos = curve->pos ( t );
                Vector2d normal = curve->dposdt ( t ).normalized();
                normal = Vector2d ( -normal.y(), normal.x() );
                double k;
                SmoothCurveUtil::curvature_eval ( curve->dposdt ( t ), curve->dposdtdt ( t ), k );
				if(std::abs(k) >= PF_EPS)
				{
					bool clamped = false;
					double r = 0.2 * 1.0 / k;
					if (r < -max_r)
					{
						r = -max_r;
						clamped = true;
					}
					if (r > max_r)
					{
						r = max_r;
						clamped = true;
					}

					draw::line(pos, pos + normal * r, clamped ? CURVATURE_STYLE_CLAMPED : CURVATURE_STYLE);
				}
					
                t += t_step;
            }
        }
    }

    inline void draw_curves_invalid ( const std::vector<GlobFitCurve*>& curves, const std::vector<bool>& invalid ) {
        assert_break ( ( int ) curves.size() == ( int ) invalid.size() );

        for ( int i = 0; i < curves.size(); ++i ) {
            if ( invalid[i] ) {
                draw::curve ( curves[i], .25, colors::red );
            } else {
                draw::curve ( curves[i], .25, colors::black );
            }
        }
    }

    inline void draw_path ( const Eigen::Matrix2Xd& path ) {
        for ( int i = 0; i < path.cols(); ++i ) {
            draw::line ( path.col ( i ), path.col ( ( i + 1 ) % path.cols() ), PATH_STYLE );
        }
    }

    inline void draw_path ( const Eigen::Matrix2Xd& path, const bool is_closed, const Style& style ) {
        const int end_index = ( int ) ( is_closed ? path.cols() : ( int ) path.cols() -1 );

        for ( int i = 0; i < end_index; ++i ) {
            draw::line ( path.col ( i ), path.col ( ( i + 1 ) % path.cols() ), style );
        }
    }

    inline void draw_raster_convexities(const Eigen::Matrix2Xd& raster) {
        std::vector<int> convexities;
        polyfit::PathUtils::compute_convexities(raster, convexities);

        for (size_t i = 0; i < raster.cols(); ++i) {
            draw::text(raster.col(i), std::to_string(convexities[i]), draw::font_pdf, Style::text());
        }
    }

    inline void draw_raster_indices ( const Eigen::Matrix2Xd& raster ) {
        for ( int i = 0; i < raster.cols(); ++i ) {
            // dbg::info(FMT("Point %f %f", raster.col(i), raster.col((i + 1) % raster.cols())));
            //draw::line(raster.col(i), raster.col((i + 1) % raster.cols()), BOUNDARY_STYLE);
            const real2 p0 = raster.col(i);
            const real2 p1 = raster.col((i + 1) % raster.cols());
            real2 off = (p0 - p1).normalized();
            off = real2(off(1), -off(0));
            draw::text ( raster.col ( i ), std::to_string ( i ), draw::font_pdf / 4, Style::text() );
        }
    }

	inline void draw_raster_polygon_indices(const Eigen::Matrix2Xd& raster, const Eigen::VectorXi& polygon) {
		for (int i = 0; i < polygon.size(); ++i) {
			const std::string label = misc::sfmt("%d(%d)", i, polygon(i));
			draw::text(raster.col(polygon(i)), label, draw::font_pdf, Style::text());
		}
	}

    inline void draw_control_polygon ( const std::vector<GlobFitCurve*>& curves ) {
        for ( int i = 0; i < curves.size(); ++i ) {
            auto bezier = dynamic_cast<BezierCurve*> ( curves[i] );

            if ( bezier == nullptr ) {
                continue;
            }

            const Eigen::Vector2d c0 = bezier->get_control_points().col ( 0 );
            const Eigen::Vector2d c1 = bezier->get_control_points().col ( 1 );
            const Eigen::Vector2d c2 = bezier->get_control_points().col ( 2 );
            const Eigen::Vector2d c3 = bezier->get_control_points().col ( 3 );

            draw::line ( c0, c1, Style::outline ( colors::talking_orange, .1 ) );
            draw::line ( c1, c2, Style::outline ( colors::talking_orange, .1 ) );
            draw::line ( c2, c3, Style::outline ( colors::talking_orange, .1 ) );
        }
    }

    inline void draw_polygon ( const Eigen::Matrix2Xd& raster, const std::vector<Eigen::Index>& polygon ) {
        for ( int j = 0; j < polygon.size(); ++j ) {
            int corner_j = ( int ) polygon[j];
            int corner_j_next = ( int ) polygon[j + 1];

            if ( corner_j >= raster.cols() ) {
                corner_j -= ( int ) raster.cols();
            }

            if ( corner_j_next >= raster.cols() ) {
                corner_j_next -= ( int ) raster.cols();
            }

            real2 src, dst;

            if ( corner_j >= raster.cols() ) {
                src = .5 * ( raster.col ( corner_j - raster.size() ) + raster.col ( corner_j - raster.size() + 1 ) );
            } else {
                src = raster.col ( corner_j );
            }

            if ( corner_j_next >= raster.size() ) {
                dst = .5 * ( raster.col ( corner_j_next - raster.size() ) + raster.col ( corner_j_next - raster.size() + 1 ) );
            } else {
                dst = raster.col ( corner_j_next );
            }

            draw::line ( src, dst, Style::outline ( colors::talking_orange, 1. ) );

            // Drawing the line of the polygon file to which the corner corresponds
            //draw::text ( src, std::to_string ( j + 1 ), draw::font_pdf * 3, Style::text() );
            //draw::text ( src, std::to_string ( polygon[j] ), draw::font_pdf * 1, Style::text() );
            real3 color;
            ( void ) color;
            double radius = 1.;
            ( void ) radius;

#if 0

            if ( labels[j] == 0 ) {
                color = colors::red;
                radius = 2.5;
            } else if ( labels[j] == 2 ) {
                color = colors::forest_green;
            } else {
                assert_break ( false );
            }

            draw::point ( src, radius, Style::fill ( color ) );
#endif
        }
    }


    inline void draw_curves ( const std::vector<Eigen::Matrix2Xd>& curves ) {
        auto build_bezier = [&] ( const Eigen::Matrix<double, 2, 4> ctrl ) {
            ::polyvec::BezierCurve bz;
            bz.set_control_points ( ctrl );
            return bz;
        };

        for ( int i = 0; i < ( int ) curves.size(); ++i ) {
            if ( curves[i].cols() == 4 ) {
                Eigen::Vector3d color =
                    i % 2 == 0 ? colors::talking_orange : colors::enemys_blood;
                auto bz = build_bezier ( curves[i] );
                auto tess = bz.get_tesselation2();

                for ( int j = 0; j < ( int ) tess.cols() - 1; ++j )
                    draw::line ( tess.col ( j ), tess.col ( j + 1 ),
                                 Style::outline ( color, 7 ) );
            } else {
                Eigen::Vector3d color = i % 2 == 0 ? colors::forest_green : colors::green;
                draw::line ( curves[i].col ( 0 ), curves[i].col ( 1 ), Style::outline ( color, 7 ) );
            }
        }
    }

    inline void draw_curves_fill ( const std::vector<Eigen::Matrix2Xd>& curves ) {
        auto build_bezier = [&] ( const Eigen::Matrix<double, 2, 4> ctrl ) {
            ::polyvec::BezierCurve bz;
            bz.set_control_points ( ctrl );
            return bz;
        };

        std::vector<real2> polygon;

        for ( int i = 0; i < ( int ) curves.size(); ++i ) {
            if ( curves[i].cols() == 4 ) {
                auto bz = build_bezier ( curves[i] );
                auto tess = bz.get_tesselation2();

                for ( int j = 0; j < ( int ) tess.cols() - 1; ++j ) {
                    polygon.push_back ( tess.col ( j ) );
                }
            } else {
                for ( int j = 0; j < ( int ) curves[i].cols() - 1; ++j ) {
                    polygon.push_back ( curves[i].col ( j ) );
                }
            }
        }

        draw::polygon ( polygon, Style::fill ( colors::dark_gray ) );
    }

}