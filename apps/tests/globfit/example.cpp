
#include <polyvec/io/vtk_curve_writer.hpp> // Optional
#include <polyvec/misc.hpp>
#include <polyvec/api.hpp>
#include <polyvec/curve_objectives.hpp>
#include <polyvec/curve_bezier.hpp>
#include <polyvec/curve_solver.hpp>

#include <memory>

using namespace polyvec;
using namespace std;

namespace {

    void
    example_impl () {

        //
        // Create input data.
        // Three bezier curves_outline .
        // Some points for bezier 0 to fit distance.
        // Some points for bezier 1 to fit tangent.
        // Same pos transition between curves_outline 0 and 1
        // Same tangent transition between curves_outline 1 and 2
        //
        VtkCurveWriter writer;

        //
        // Bezier Curves and a line
        //
		auto curve0 = make_shared<BezierCurve>();
		auto curve1 = make_shared<BezierCurve>();
		auto curve2 = make_shared<GlobFitCurve_Line>();

        // Curve 0
        {
            Eigen::Matrix<double, 2, 4> control_points;
            control_points.row ( 0 ) << 0, 0.5, 0.4, 1;
            control_points.row ( 1 ) << 0, 0.3, 0.7, 1;
            curve0->set_control_points ( control_points );

            for ( int i=0; i<4; ++i ) {
                writer.add_point ( control_points.col ( i ) );
            }

            writer.add_polyline ( curve0->get_tesselation2() );
            writer.add_polyline ( control_points );
        }

        // Curve 1
        {
            Eigen::Matrix<double, 2, 4> control_points;
            control_points.row ( 0 ) << 1.1, 0.7, 0.2, 0;
            control_points.row ( 1 ) << 1.1, 1.1, 1.3, 2;
            curve1->set_control_points ( control_points );

            for ( int i=0; i<4; ++i ) {
                writer.add_point ( control_points.col ( i ) );
            }

            writer.add_polyline ( curve0->get_tesselation2() );
            writer.add_polyline ( control_points );
        }

        // Curve 2
        {
            Eigen::Matrix<double, 2, 2> end_points;
            end_points.row ( 0 ) << -0.1, 0.5; // x values
            end_points.row ( 1 ) << 2.1,  3; // y values
            curve2->set_points ( end_points );

            for ( int i=0; i<2; ++i ) {
                writer.add_point ( end_points.col ( i ) );
            }

            writer.add_polyline ( curve2->get_tesselation2() );
        }

		GlobFitBezierAngleBasedParametrization param0(curve0);
		GlobFitBezierAngleBasedParametrization param1(curve1);
		GlobFitLineShiftBasedParametrization param2(curve2);

        //
        // Objectives
        //
        // cannot put the points in a block as the lambda will access out of
        // scope variables
        // Eigen::Vectorization is disabled, so it is fine to have vec<Eigen::Vec>
        GlobFitObjective_SamePosition obj_same_pos;
        GlobFitObjective_SameTangent obj_same_tangent;
        //
        GlobFitObjective_FitPointLength obj_fix_pos_curve2;
        GlobFitObjective_FitPointLength obj_fix_pos_curve0;
        GlobFitObjective_FitTangentLength obj_fix_tagent_curve0;
        //
        GlobFitObjective_PointLieOnLine obj_bezier_point_lie_on_line_curve0;

        //
        std::vector<Eigen::Vector2d> points_to_fit;
        std::vector<Eigen::Vector2d> tangents_to_fit_pos;
        std::vector<Eigen::Vector2d> tangents_to_fit_dir;
        std::vector<double> weights;
        points_to_fit = {Eigen::Vector2d ( 0.2, 0.2 ),
                         Eigen::Vector2d ( 0.4, 0.6 ),
                         Eigen::Vector2d ( 0.45, 0.55 ),
                         Eigen::Vector2d ( 0.47, 0.36 ),
                         Eigen::Vector2d ( 0.7, 0.5 ),
                         Eigen::Vector2d ( 0.99, 0.58 )
                        };
        tangents_to_fit_pos = {Eigen::Vector2d ( 0.2, 0.8 ), Eigen::Vector2d ( 0.4, 0.9 ), Eigen::Vector2d ( 0.7, 1.3 ) };
        tangents_to_fit_dir = {Eigen::Vector2d ( 1, 1 ), Eigen::Vector2d ( -1, 1 ), Eigen::Vector2d ( 0, 1 ) };
        weights = {1., 1., 1.};

        // same pos
		obj_same_pos.set_params(&param0, &param1);
        obj_same_pos.set_weight ( 100 );

        // same tangent
		obj_same_tangent.set_params(&param1, &param2);
        obj_same_tangent.set_weight ( 100 );

        // Fix stuff on curve 0
		obj_fix_pos_curve0.set_params(&param0, 0., curve0->pos(0.));
        obj_fix_pos_curve0.set_weight ( 100 );
        //
		obj_fix_tagent_curve0.set_params(&param0, 0., curve0->dposdt(0).normalized());
        obj_fix_tagent_curve0.set_weight ( 100 );


        // Fix stuff on curve 2
		obj_fix_pos_curve2.set_params(&param2, 1., curve2->pos(1.));
        obj_fix_pos_curve2.set_weight ( 100 );

        // Point lie on line curve 0
        const Eigen::Vector2d obj_bezier_point_lie_on_line_curve0_LINE_POINT ( 1, -0.3 );
        const Eigen::Vector2d obj_bezier_point_lie_on_line_curve0_LINE_DIR ( 0.5, 0.2  );
		obj_bezier_point_lie_on_line_curve0.set_params(&param0, 0.1, obj_bezier_point_lie_on_line_curve0_LINE_POINT, obj_bezier_point_lie_on_line_curve0_LINE_DIR);

        {
            writer.add_line ( obj_bezier_point_lie_on_line_curve0_LINE_POINT-3*obj_bezier_point_lie_on_line_curve0_LINE_DIR,
                              obj_bezier_point_lie_on_line_curve0_LINE_POINT+3*obj_bezier_point_lie_on_line_curve0_LINE_DIR );
            writer.add_point ( obj_bezier_point_lie_on_line_curve0_LINE_POINT );
            writer.dump ( "test_dump/globfit_example_0.vtk" );
        }

        // Run the fitter
        auto callback = [&] ( int idx ) {
            writer.clear();
            writer.add_polyline ( curve0->get_tesselation2() );
            writer.add_polyline ( curve1->get_tesselation2() );
            writer.add_polyline ( curve2->get_tesselation2() );
            //
            writer.add_line ( obj_bezier_point_lie_on_line_curve0_LINE_POINT-3*obj_bezier_point_lie_on_line_curve0_LINE_DIR,
                              obj_bezier_point_lie_on_line_curve0_LINE_POINT+3*obj_bezier_point_lie_on_line_curve0_LINE_DIR );
            writer.add_point ( obj_bezier_point_lie_on_line_curve0_LINE_POINT );
            //
            writer.dump ( polyvec_str ( "test_dump/globfit_example_" << idx+1 << ".vtk" ) );
        };

        GlobFitter globfitter;
        globfitter.set_curves ( {&param0, &param1, &param2} );
        globfitter.set_objectives ( {
            &obj_same_tangent,
            &obj_same_pos,
            &obj_fix_pos_curve2,
            &obj_fix_pos_curve0,
            &obj_fix_tagent_curve0,
            &obj_bezier_point_lie_on_line_curve0
        } );
        globfitter.setup ( 30 /*max iter*/ ); // don't forget to call
        globfitter.run_fitter ( stdout, callback );


    } //end of function
} // end of anonymus

namespace polyvectest {
    namespace globfit {
        int example ( int /*argc*/, char** /*argv*/ ) {
            example_impl();
            return 0;
        } // end of function
    } // globfitter
} // deadline code test
