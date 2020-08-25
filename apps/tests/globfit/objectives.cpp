#include <polyvec/io/vtk_curve_writer.hpp> // Optional
#include <polyvec/misc.hpp>
#include <polyvec/api.hpp>
#include <polyvec/curve_objectives.hpp>
#include <polyvec/curve_bezier.hpp>
#include <polyvec/curve_solver.hpp>
#include <polyvec/utils/deriv_test.hpp>

using namespace polyvec;
using namespace std;

namespace {



// Must have called setup() on globfit before calling this
    void
    subtest_globfit_objective ( GlobFitter& globfit, const std::string obj_name ) {

        assert_break ( globfit.is_setup() );
        printf ( "-- Testing Objective: %s \n", obj_name.c_str() );

        Eigen::VectorXd params0 = globfit.get_params();
        Eigen::VectorXd delta = misc::randn ( globfit.n_parameters(), 1 );
        double h0 = 1;
        int n_halving = 8;
        Eigen::VectorXd obj;
        //Eigen::MatrixXd jac;
		Eigen::SparseMatrix<double> jac;


        // Test change of variables for fit
        derivtest::run ( params0,
                         delta,
        [&] ( const Eigen::VectorXd & new_params ) -> Eigen::VectorXd {
            globfit.set_params ( new_params );
            globfit.compute_objective_and_jacobian ( obj, jac );
            return obj;
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::VectorXd  { //
            return jac* hdelta;
        },
        std::cout,
        n_halving,
        h0 );
    } // All done with subtest_globfit_objective

}

namespace polyvectest {
    namespace globfit {
        int test_objectives ( int /*argc*/, char** /*argv*/ ) {

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
            // Bezier Curves
            //
			auto bz0 = make_shared<BezierCurve>();
			auto bz1 = make_shared<BezierCurve>();
			auto bz2 = make_shared<BezierCurve>();

            // Curve 0
            {
                BezierCurve bzc0;
                Eigen::Matrix<double, 2, 4> control_points;
                control_points.row ( 0 ) << 0, 0.5, 0.4, 1;
                control_points.row ( 1 ) << 0, 0.3, 0.7, 1;
                bz0->set_control_points ( control_points );

                for ( int i=0; i<4; ++i ) {
                    writer.add_point ( control_points.col ( i ) );
                }

                writer.add_polyline ( bzc0.get_tesselation2() );
                writer.add_polyline ( control_points );
            }

            // Curve 1
            {
                BezierCurve bzc1;
                Eigen::Matrix<double, 2, 4> control_points;
                control_points.row ( 0 ) << 1.1, 0.7, 0.8, 0;
                control_points.row ( 1 ) << 1.1, 1.1, 1.3, 2;
                bz1->set_control_points ( control_points );

                for ( int i=0; i<4; ++i ) {
                    writer.add_point ( control_points.col ( i ) );
                }

                writer.add_polyline ( bzc1.get_tesselation2() );
                writer.add_polyline ( control_points );
            }

            // Curve 2
            {
                BezierCurve bzc2;
                Eigen::Matrix<double, 2, 4> control_points;
                control_points.row ( 0 ) << -0.1, 0.5, -0.5, 0.5;
                control_points.row ( 1 ) << 2.1, 2.2, 2.5, 3;
                bz2->set_control_points ( control_points );

                for ( int i=0; i<4; ++i ) {
                    writer.add_point ( control_points.col ( i ) );
                }

                writer.add_polyline ( bzc2.get_tesselation2() );
                writer.add_polyline ( control_points );
            }

			GlobFitBezierAngleBasedParametrization param0(bz0);
			GlobFitBezierAngleBasedParametrization param1(bz1);
			GlobFitBezierAngleBasedParametrization param2(bz2);

            //
            // Objectives
            //
            // cannot put the points in a block as the lambda will access out of
            // scope variables
            // Eigen::Vectorization is disabled, so it is fine to have vec<Eigen::Vec>
            GlobFitObjective_SamePosition obj_same_pos;
            GlobFitObjective_SameTangent obj_same_tangent;
            GlobFitObjective_FitPointLength obj_fix_pos;
            GlobFitObjective_FitTangentLength obj_fix_tangent_cross_prod_formulation;
            GlobFitObjective_FitTangentLength obj_fix_tangent_entry_diff_formulation;
            GlobFitObjective_PointLieOnLine obj_point_lie_on_line;

            //
            std::vector<Eigen::Vector2d> points_to_fit;
            std::vector<Eigen::Vector2d> tangents_to_fit_pos;
            std::vector<Eigen::Vector2d> tangents_to_fit_dir;
            std::vector<double> weights;
            points_to_fit = {Eigen::Vector2d ( 0.2, 0.2 ), Eigen::Vector2d ( 0.4, 0.6 ), Eigen::Vector2d ( 0.7, 0.5 ) };
            tangents_to_fit_pos = {Eigen::Vector2d ( 0.2, 0.8 ), Eigen::Vector2d ( 0.4, 0.9 ), Eigen::Vector2d ( 0.7, 1.3 ) };
            tangents_to_fit_dir = {Eigen::Vector2d ( 1, 1 ), Eigen::Vector2d ( -1, 1 ), Eigen::Vector2d ( 0, 1 ) };
            weights = {1., 1., 1.};

            // same pos
			obj_same_pos.set_params(&param0, &param1);
            obj_same_pos.set_weight ( 1.2 );

            // same tangent 
            obj_same_tangent.set_formulation(GlobFitObjective_SameTangent::FORMULATION_VECTOR_DIFF);
			obj_same_tangent.set_params(&param1, &param2);
            obj_same_tangent.set_weight ( 0.9 );
            
			// Fix pos
			obj_fix_pos.set_params(&param0, 0., misc::randn(2, 1));
            obj_fix_pos.set_weight ( 0.85 );


            // Fix tangent
			obj_fix_tangent_cross_prod_formulation.set_params(&param2, 0.6568, misc::randn(2, 1).normalized());
            obj_fix_tangent_cross_prod_formulation.set_weight ( 1.3 );

            // point lie on line
			obj_point_lie_on_line.set_params(&param0, 0.1, Eigen::Vector2d(1, -0.3), Eigen::Vector2d(0.5, 0.2));
            obj_point_lie_on_line.set_weight ( 10. );

            //
            // Test same pos
            //
            {
                GlobFitter globfitter;
                globfitter.set_curves ( {&param0, &param1, &param2} );
                globfitter.set_objectives ( {&obj_same_pos} );
                globfitter.setup ( 60 /* maxiter -- does not matter in this case */ );
                globfitter.set_params ( globfitter.get_params() );
                subtest_globfit_objective ( globfitter, "SAMEPOS FITTER OBJECTIVE" );
            }

            //
            // Test same tangent
            //
            {
                GlobFitter globfitter;
				globfitter.set_curves({ &param0, &param1, &param2 });
                globfitter.set_objectives ( {&obj_same_tangent} );
                globfitter.setup ( 60 /* maxiter -- does not matter in this case */ );
                globfitter.set_params ( globfitter.get_params() );
                subtest_globfit_objective ( globfitter, "SAME TANGENT FITTER OBJECTIVE" );
            }

            //
            // Test fix pos
            //
            {
                GlobFitter globfitter;
				globfitter.set_curves({ &param0, &param1, &param2 });
                globfitter.set_objectives ( {&obj_fix_pos} );
                globfitter.setup ( 60 /* maxiter -- does not matter in this case */ );
                globfitter.set_params ( globfitter.get_params() );
                subtest_globfit_objective ( globfitter, "FIX POS OBJECTIVE" );
            }

            //
            // Test fix tangent
            //
            {
                GlobFitter globfitter;
				globfitter.set_curves({ &param0, &param1, &param2 });
                globfitter.set_objectives ( {&obj_fix_tangent_cross_prod_formulation} );
                globfitter.setup ( 60 /* maxiter -- does not matter in this case */ );
                globfitter.set_params ( globfitter.get_params() );
                subtest_globfit_objective ( globfitter, "FIX TANGENT CROSS PROD  OBJECTIVE" );
            }

            //
            // Test the bezier dev from init
            //
            {
                GlobFitter globfitter;
				globfitter.set_curves({ &param0, &param1, &param2 });
                globfitter.set_objectives ( {&obj_point_lie_on_line} );
                globfitter.setup ( 60 /* maxiter -- does not matter in this case */ );
                globfitter.set_params ( globfitter.get_params() ); // why am I doing this?
                subtest_globfit_objective ( globfitter, "CURVE POINT LIE ON LINE" );
            }

            //
            // Now all together
            //
            {
                GlobFitter globfitter;
				globfitter.set_curves({ &param0, &param1, &param2 });
                globfitter.set_objectives ( {
                    &obj_same_tangent,
                    &obj_same_pos,
                    &obj_fix_pos,
                    &obj_fix_tangent_cross_prod_formulation,
                    &obj_point_lie_on_line
                } );
                globfitter.setup ( 60 /* maxiter -- does not matter in this case */ );
                globfitter.set_params ( globfitter.get_params() );
                subtest_globfit_objective ( globfitter, "ALL TOGETHER" );
            }

            writer.dump ( "test_dump/test_objectives.vtk" );

            return 0;
        }
    } // globfitter
} // deadline code test