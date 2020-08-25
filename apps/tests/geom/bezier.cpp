#include <sstream>

#include <polyvec/api.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/curve_bezier.hpp>
#include <polyvec/utils/deriv_test.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>

using namespace polyvec;

namespace {

void
test_draw() {

    printf ( "-- Drawing a bezier plus tangents and normal \n" );

    // Test drawing
    BezierCurve bz;
    VtkCurveWriter writer;

    Eigen::Matrix<double, 2, 4> control_points;

    control_points.row ( 0 ) << 0, 1, 2, 0;
    control_points.row ( 1 ) << 0, 0, 1, 1;

    bz.set_control_points ( control_points );

    //
    // Draw the bezier and the control points
    //
    writer.add_point ( control_points.col ( 0 ) );
    writer.add_point ( control_points.col ( 1 ) );
    writer.add_point ( control_points.col ( 2 ) );
    writer.add_point ( control_points.col ( 3 ) );
    //
    writer.add_polyline ( control_points );
    //
    writer.add_polyline ( bz.get_tesselation2() );
    writer.dump ( "test_dump/draw_bezier.vtk" );

    //
    // Draw the bezier tangets
    //
    std::vector<double> t_values = {0, 0.25, 0.5, 0.75, 1.};

    writer.clear();

    for ( unsigned i = 0; i < t_values.size(); ++i ) {
        const double t = t_values[i];
        const Eigen::Vector2d pos = bz.pos ( t );
        const Eigen::Vector2d tangent = bz.dposdt ( t ).normalized();

        writer.add_line ( pos, pos + 0.2 * tangent );
    }

    writer.dump ( "test_dump/draw_bezier_tangent.vtk" );


    //
    // Test bezier normals
    //
    writer.clear();

    for ( unsigned i = 0; i < t_values.size(); ++i ) {
        const double t = t_values[i];
        const Eigen::Vector2d pos = bz.pos ( t );
        const Eigen::Vector2d tangent = bz.dposdt ( t ).normalized();
        const Eigen::Vector2d der2 = bz.dposdtdt ( t );
        const Eigen::Vector2d normal = ( der2 - der2.dot ( tangent ) * tangent ).normalized();

        writer.add_line ( pos, pos + 0.2 * normal );
    }

    writer.dump ( "test_dump/draw_bezier_normal.vtk" );


} // All done test draw



void
test_t_derivatives() {
    printf ( "-- Testing derivatives w.r.t t \n" );

    BezierCurve bz;

    Eigen::Matrix<double, 2, 4> control_points;

    control_points.row ( 0 ) << 0, 1, 2, 0;
    control_points.row ( 1 ) << 0, 0, 1, 1;

    bz.set_control_points ( control_points );

    using Vec1d = Eigen::Matrix<double, 1, 1>;
    Vec1d t0, delta;
    double t_perturbed;
    double h0 = 0.2;
    int n_halving = 5;
    t0 << 0.6841;
    delta << 1;

    derivtest::run (
        t0,
        delta,
    [&] ( const Eigen::VectorXd & t_perturbed_in ) -> Eigen::Vector2d {
        assert ( t_perturbed_in.size() == 1 );
        t_perturbed = t_perturbed_in ( 0 );
        return bz.pos ( t_perturbed );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
        assert ( hdelta.size() == 1 );
        return bz.dposdt ( t_perturbed ) * hdelta ( 0 );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
        assert ( hdelta.size() == 1 );
        Vec1d ans;
        return bz.dposdtdt ( t_perturbed ) * hdelta ( 0 ) * hdelta ( 0 );
    },
    std::cout,
    n_halving,
    h0 );
} // All done with t derivatives


void
test_param_derivatives() {

    BezierCurve bz;

    Eigen::Matrix<double, 2, 4> control_points;

    control_points.row ( 0 ) << 0, 1, 2, 0;
    control_points.row ( 1 ) << 0, 0, 1, 1;

    bz.set_control_points ( control_points );

    Eigen::VectorXd params0 = misc::reshaped ( bz.get_control_points(), 8, 1 );
    Eigen::VectorXd delta = misc::randn ( 8, 1 );
    double h0 = 0.2;
    int n_halving = 5;
    const double t = 0.65684;

    printf ( "-- Testing derivatives of pos w.r.t curve params \n" );

    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Eigen::Vector2d {
        bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
        return bz.pos ( t );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { return bz.dposdparams ( t ) * hdelta; },
    std::cout,
    n_halving,
    h0 );

    printf ( "-- Testing derivatives of dposdt w.r.t curve params \n" );

    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Eigen::Vector2d {
        bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
        return bz.dposdt ( t );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { return bz.dposdtdparams ( t ) * hdelta; },
    std::cout,
    n_halving,
    h0 );

    printf ( "-- Testing derivatives of dposdtdtdt w.r.t curve params \n" );

    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Eigen::Vector2d {
        bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
        return bz.dposdtdt ( t );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { //
    return bz.dposdtdtdparams ( t ) * hdelta; 
    },
    std::cout,
    n_halving,
    h0 );

    printf ( "-- Testing derivatives of LENGTH w.r.t params \n" );

    using Vec1d = Eigen::Matrix<double, 1, 1>;
    auto vec1d = [] ( const double v ) {
        Vec1d vv;
        vv << v;
        return vv;
    };

    derivtest::run (
        params0,
        delta,
    [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Vec1d {
        bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
        return vec1d ( bz.length() );
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Vec1d { return bz.dlengthdparams().transpose() * hdelta; },
    std::cout,
    n_halving,
    h0 );

} // All done with param derivatives


void
test_projection ( const Eigen::Matrix2Xd& control_points, const Eigen::Vector2d& project_me, bool test_derivative ) {

    printf ( "-- Testing projection \n" );
    static int n_called = -1;
    ++n_called;

    BezierCurve bz;
    VtkCurveWriter writer;
    bz.set_control_points ( control_points );

    {
        const double tproject = bz.project ( project_me );
        const Eigen::Vector2d projection = bz.pos ( tproject );
        writer.add_line ( project_me, projection );
        writer.add_polyline ( control_points );
        writer.add_polyline ( bz.get_tesselation2() );
        writer.add_point ( projection );
        writer.add_point ( project_me );
    }

    writer.dump ( polyvec_str ( "test_dump/bezier_projection_" << n_called << ".vtk" ) );

    if ( test_derivative ) {
        Eigen::VectorXd params0 = misc::reshaped ( bz.get_control_points(), 8, 1 );
        Eigen::VectorXd delta = misc::randn ( 8, 1 );
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
        [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Vec1d {
            bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
            t_proj = bz.project ( project_me );
            return vec1d ( t_proj );
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Vec1d { //
            return bz.dtprojectdparams ( t_proj, project_me ).transpose() * hdelta;
        },
        std::cout,
        n_halving,
        h0 );

        printf ( "- testing dposproj dparams \n" );
        derivtest::run (
            params0,
            delta,
        [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Eigen::Vector2d {
            bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
            t_proj = bz.project ( project_me );
            return bz.pos ( t_proj );
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { //
            const Eigen::VectorXd dtprojdparams = bz.dtprojectdparams ( t_proj, project_me );
            return bz.dposprojectdparams ( t_proj, dtprojdparams ) * hdelta;
        },
        std::cout,
        n_halving,
        h0 );

        printf ( "- testing ddtanjproj dparams \n" );
        derivtest::run (
            params0,
            delta,
        [&] ( const Eigen::VectorXd & perturbed_control_points ) -> Eigen::Vector2d {
            bz.set_control_points ( misc::reshaped ( perturbed_control_points, 2, 4 ) );
            t_proj = bz.project ( project_me );
            return bz.dposdt ( t_proj );
        },
        [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d { //
            const Eigen::VectorXd dtprojdparams = bz.dtprojectdparams ( t_proj, project_me );
            return bz.dposdtprojectdparams ( t_proj, dtprojdparams ) * hdelta;
        },
        std::cout,
        n_halving,
        h0 );
    }


} // All done with param derivatives

} // end of anonymus

namespace polyvectest {
namespace geom {
int
test_bezier ( int /*argc*/, char** /*argv*/ ) {

    test_draw();

    test_t_derivatives();
    test_param_derivatives();

    {
        Eigen::Matrix<double, 2, 4> control_points;
        control_points.row ( 0 ) << 0, 1, 2, 0;
        control_points.row ( 1 ) << 0, .2, 1, 1;
        test_projection ( control_points, Eigen::Vector2d ( 2, 1 ), true );
        test_projection ( control_points, Eigen::Vector2d ( 0.5, 0.5 ), true );
        test_projection ( control_points, Eigen::Vector2d ( -1, 0 ), true );
        test_projection ( control_points, Eigen::Vector2d ( 0, 1.3 ), true );
    }

    return 0;
}
}
}
