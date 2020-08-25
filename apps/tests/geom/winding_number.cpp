#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/winding_number.hpp>

using namespace polyvec;

namespace {
    void
    run_test_winding ( const  Eigen::Matrix2Xd& polygon, const Eigen::Vector2d& point ) {
        static int n_called = 0;

        VtkCurveWriter writer;
        bool is_trustable;
        double winding_number;

        writer.add_polyline ( polygon );
        writer.add_line ( polygon.col ( polygon.cols()-1 ), polygon.col ( 0 ) );

        writer.add_point ( point );

        WindingNumber::compute_winding ( polygon, point, winding_number, is_trustable );
        printf ( "-(%d) Testing WINDING #: value: %8.4f   is_trustable: %d \n", n_called, winding_number, ( int ) is_trustable );

        writer.dump ( polyvec_str ( "test_dump/test_winding_number" << n_called << ".vtk" ) );

        ++n_called;
    }

    void
    run_test_orientation ( const  Eigen::Matrix2Xd& polygon ) {
        static int n_called = 0;

        VtkCurveWriter writer;
        bool is_ccw;

        writer.add_polyline ( polygon );
        writer.add_line ( polygon.col ( polygon.cols()-1 ), polygon.col ( 0 ) );


        WindingNumber::compute_orientation ( polygon, is_ccw );
        printf ( "-(%d) Testing WINDING #: is_ccw: %d  \n", n_called, ( int ) is_ccw );

        writer.dump ( polyvec_str ( "test_dump/test_orientation" << n_called << ".vtk" ) );

        ++n_called;
    }

} // end of anonymus

namespace polyvectest {
    namespace geom {
        int test_winding_number ( int /*argc*/, char** /*argv*/ ) {

            const unsigned n_points = 4;
            Eigen::Matrix2Xd polygon ( 2, n_points );
            polygon.col(0) << 8, -1;
            polygon.col(1) << 8, +1;
            polygon.col(2) << -8, +1;
            polygon.col(3) << -8, -1;

            run_test_winding(polygon, Eigen::Vector2d(-8, -1));
            run_test_winding(polygon, Eigen::Vector2d(-4, -1));
            run_test_winding(polygon, Eigen::Vector2d(-12, -1));
            run_test_winding(polygon, Eigen::Vector2d(-4, -0.9999));
            run_test_winding(polygon, Eigen::Vector2d(-2.6565, -1.3));
            run_test_winding(polygon, Eigen::Vector2d(-2.6565, -0.74));

            run_test_orientation(polygon); // 1
            polygon.row(0) *= -1;
            run_test_orientation(polygon); // -1
            polygon.row(0) *= -1;
            polygon.col(0) *= -1;
            run_test_orientation(polygon); // -1
            polygon.col(0) *= -1;
            run_test_orientation(polygon); // 1


            return 0;
        } // end of function
    } // end of share
} // end of deadline code test