#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/geometry/line.hpp>

using namespace polyvec;

namespace {
void
run_test ( const  Eigen::Matrix2d& line1, const Eigen::Matrix2d& line2 ) {
    static int n_called = 0;

    VtkCurveWriter writer;
    bool do_intersect;
    double t_intersection_line_1, t_intersection_line_2;

    do_intersect = LineUtils::intersect (
        line1.col ( 0 ), line1.col ( 1 ), line2.col ( 0 ), line2.col ( 1 ),
        t_intersection_line_1, t_intersection_line_2 );

    writer.add_polyline ( line1 );
    writer.add_polyline ( line2 );

    if ( do_intersect ) {
        printf ( "- Testing intersection: POSITIVE at line1 param %.10g \n", t_intersection_line_2 );
        const Eigen::Vector2d intersection_2 = LineUtils::line_at ( line2.col ( 0 ), line2.col ( 1 ),
                                               t_intersection_line_2 );
        const Eigen::Vector2d intersection_1 = LineUtils::line_at ( line1.col ( 0 ), line1.col ( 1 ),
                                               t_intersection_line_1 );
        const double diff = ( intersection_2 - intersection_1 ).norm();
        assert_break ( diff < 1e-8 );
        writer.add_point ( intersection_1 );
        writer.add_line ( line1.col ( 0 ), intersection_1 );
        writer.add_line ( line2.col ( 0 ), intersection_1 );
    } else {
        printf ( "- Testing intersection: NEGATIVE \n" );
    }

    writer.dump ( polyvec_str ( "test_dump/share_test_line_intersection_" << n_called << ".vtk" ) );

    ++n_called;
}

} // end of anonymus

namespace polyvectest {
namespace geom {
int test_line_intersection ( int /*argc*/, char** /*argv*/ ) {

    VtkCurveWriter writer;

    auto run_test_d = [] ( double v0, double v1, double v2, double v3, double v4, double v5, double v6, double v7 ) {
        Eigen::Matrix2Xd line1 ( 2, 2 );
        Eigen::Matrix2Xd line2 ( 2, 2 );
        line1.col ( 0 ) = Eigen::Vector2d ( v0, v1 );
        line1.col ( 1 ) = Eigen::Vector2d ( v2, v3 );
        line2.col ( 0 ) = Eigen::Vector2d ( v4, v5 );
        line2.col ( 1 ) = Eigen::Vector2d ( v6, v7 );
        run_test ( line1, line2 );
    };

    run_test_d ( 0, 0, 0.5, 0,
                 4, -1, 4, 0.2 );
    run_test_d ( 0, 0, 0.5, 0,
                 4, 0.2, 4, -1 );
    run_test_d ( 0.2, 0.3, 0.5, 0.1,
                 4, -1, 4.1, -0.6 );
    run_test_d ( 0.2, 0.3, 0.5, 0.1,
                 4.1, -0.6, 4, -1 );
    run_test_d ( 0.2, 0.3, 0.5, 0.1,
                 3.2, 3.3, 3.5, 3.1 );


    return 0;
} // end of function
} // end of share
} // end of deadline code test
