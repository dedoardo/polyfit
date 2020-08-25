#include <polyvec/io/vtk_curve_writer.hpp> // Optional
#include <polyvec/misc.hpp>
#include <polyvec/api.hpp>
#include <polyvec/curve_objectives.hpp>
#include <polyvec/curve_bezier.hpp>
#include <polyvec/curve_solver.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/utils/deriv_test.hpp>

using namespace polyvec;

// also test arc
#define WITH_ARC

namespace {

using Vec1d = Eigen::Matrix<double, 1, 1>;
Vec1d vec1d ( const double v ) {
    Vec1d vv;
    vv << v;
    return vv;
};


void
test_param_derivatives ( GlobFitCurve& curve ) {

    Eigen::VectorXd params0 = curve.get_params();
    Eigen::VectorXd delta = misc::randn ( curve.n_params(), 1 );
    double h0 = 0.2;
    int n_halving = 5;
    const double t = 0.65684;

    printf ( "-- Testing derivatives of pos w.r.t t \n" );

    if ( 1 ) {
        double perturbed_dt;
        derivtest::run (
            vec1d ( t ),
            vec1d ( 1. ),
        [&] ( const Eigen::VectorXd & perturbed_params ) -> Eigen::Vector2d {
            perturbed_dt = perturbed_params ( 0 );
            return curve.pos ( perturbed_dt );
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
            return curve.dposdt ( perturbed_dt ) *  hdelta ;
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
            return curve.dposdtdt ( perturbed_dt ) * hdelta* hdelta;
        },
        std::cout,
        n_halving,
        h0 );
    }

    printf ( "-- Testing derivatives of pos w.r.t curve params \n" );

    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_params ) -> Eigen::Vector2d {
        curve.set_params ( perturbed_params );
        return curve.pos ( t );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
        return curve.dposdparams ( t ) * hdelta;
    },
    std::cout,
    n_halving,
    h0 );

    printf ( "-- Testing derivatives of dposdt w.r.t curve params \n" );

    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_params ) -> Eigen::Vector2d {
        curve.set_params ( perturbed_params );
        return curve.dposdt ( t );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
        return curve.dposdtdparams ( t ) * hdelta;
    },
    std::cout,
    n_halving,
    h0 );

    printf ( "-- Testing derivatives of dposdtdt w.r.t curve params \n" );

    if ( curve.get_type() != GLOBFIT_CURVE_ARC ) {
    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_params ) -> Eigen::Vector2d {
        curve.set_params ( perturbed_params );
        return curve.dposdtdt ( t );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
        return curve.dposdtdtdparams ( t ) * hdelta;
    },
    std::cout,
    n_halving,
    h0 );

    }

    if ( curve.get_type() == GLOBFIT_CURVE_BEZIER ) {

        BezierCurve& curve_bz = *static_cast<BezierCurve*> ( &curve );

        printf ( "-- Testing derivatives of LENGTH w.r.t params \n" );

        derivtest::run (
            params0,
            delta,
        [&] ( const Eigen::VectorXd & perturbed_params ) -> Vec1d {
            curve_bz.set_params ( perturbed_params );
            return vec1d ( curve_bz.length() );
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Vec1d {
            return curve_bz.dlengthdparams().transpose() * hdelta;
        },
        std::cout,
        n_halving,
        h0 );
    }
} // All done with param derivatives


void
test_projection ( GlobFitCurve& curve, const Eigen::Vector2d& project_me ) {
//test_bezier_projection ( const Eigen::Matrix2Xd& control_points, const Eigen::Vector2d& project_me ) {

    printf ( "-- Testing projection \n" );
    static int n_called = -1;
    ++n_called;

    //BezierCurve bz;
    VtkCurveWriter writer;
    //GlobFitCurve_Bezier globfit_bz;
    //bz.set_control_points ( control_points );
    //globfit_bz.set_bezier ( bz );


    {
        const double tproject = curve.project ( project_me );
        const Eigen::Vector2d projection = curve.pos ( tproject );
        writer.add_line ( project_me, projection );

        if ( curve.get_type() == GLOBFIT_CURVE_BEZIER ) {
            writer.add_polyline ( static_cast<BezierCurve*> ( &curve )->get_control_points() );
        }

        writer.add_polyline ( curve.get_tesselation2() );
        writer.add_point ( projection );
        writer.add_point ( project_me );
    }

    writer.dump ( polyvec_str ( "test_dump/globfit_curve_projection_" << n_called << ".vtk" ) );

    Eigen::VectorXd params0 = curve.get_params();
    Eigen::VectorXd delta = misc::randn ( curve.n_params(), 1 );

    double h0 = 0.2;
    int n_halving = 5;
    double t_proj = -1e10;

    using Vec1d = Eigen::Matrix<double, 1, 1>;
    auto vec1d = [] ( const double v ) {
        Vec1d vv;
        vv << v;
        return vv;
    };

    printf ( "- testing dtproj dparams \n" );
    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_params ) -> Vec1d {
        curve.set_params ( perturbed_params );
        t_proj = curve.project ( project_me );
        return vec1d ( t_proj );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Vec1d { //
        return curve.dtprojectdparams ( t_proj, project_me ).transpose() * hdelta;
    },
    std::cout,
    n_halving,
    h0 );

    printf ( "- testing dposproj dparams \n" );
    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_params ) -> Eigen::Vector2d {
        curve.set_params ( perturbed_params );
        t_proj = curve.project ( project_me );
        return curve.pos ( t_proj );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { //
        const Eigen::VectorXd dtprojdparams = curve.dtprojectdparams ( t_proj, project_me );
        return curve.dposprojectdparams ( t_proj, dtprojdparams ) * hdelta;
    },
    std::cout,
    n_halving,
    h0 );

    printf ( "- testing ddtanjproj dparams \n" );
    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_params ) -> Eigen::Vector2d {
        curve.set_params ( perturbed_params );
        t_proj = curve.project ( project_me );
        return curve.dposdt ( t_proj );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { //
        const Eigen::VectorXd dtprojdparams = curve.dtprojectdparams ( t_proj, project_me );
        return curve.dposdtprojectdparams ( t_proj, dtprojdparams ) * hdelta;
    },
    std::cout,
    n_halving,
    h0 );
}// All done with param derivatives

}// end of anonymus

namespace polyvectest {
namespace globfit {

int
test_curves ( int /*argc*/, char** /*argv*/ ) {

    printf ( "----- TESTING BEZIER - \n" );
    {
        Eigen::Matrix<double, 2, 4> control_points;
        control_points.row ( 0 ) << 0, 1, 2, 0;
        control_points.row ( 1 ) << 0, 0.1, 1, 1;
        BezierCurve bz;
        bz.set_control_points ( control_points );

        test_param_derivatives ( bz );
        test_projection ( bz, Eigen::Vector2d ( 2, 1 ) );
    }

    printf ( "----- TESTING LINE - \n" );
    {
        Eigen::Matrix<double, 2, 2> control_points;
        control_points.row ( 0 ) << 0.2, 0.9;
        control_points.row ( 1 ) << 0.1, 1.3;
        GlobFitCurve_Line globfit_line;
        globfit_line.set_points ( control_points );
        test_param_derivatives ( globfit_line );
        test_projection ( globfit_line, Eigen::Vector2d ( 0.5, 0.7 ) );
    }

    return 0;
}
} // globfitter
} // deadlinecodetest