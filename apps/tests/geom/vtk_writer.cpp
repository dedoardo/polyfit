#include <polyvec/utils/deriv_test.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>

namespace polyvectest {
namespace geom {
int test_vtk_writer ( int /*argc*/, char** /*argv*/ ) {

    ::polyvec::VtkCurveWriter writer;

    Eigen::Matrix2Xd polygon1 ( 2, 4 );
    Eigen::Matrix2Xd polygon2 ( 2, 2 );
    Eigen::Matrix2Xd points ( 2, 3 );

    polygon1.row ( 0 ) << 0, 1, 1, 2;
    polygon1.row ( 1 ) << 0, 0, 1, 1;

    polygon2.row ( 0 ) << 1, 0;
    polygon2.row ( 1 ) << 0, 1;

    points.row ( 0 ) << 0.5, 1.5, -0.5;
    points.row ( 1 ) << 0.5, 1.5, -0.5;


    writer.add_line ( polygon2.col ( 0 ), polygon2.col ( 1 ) );
    writer.add_polyline ( polygon1 );
    writer.add_point ( points.col ( 0 ) );
    writer.add_point ( points.col ( 1 ) );
    writer.add_point ( points.col ( 2 ) );


    writer.dump ( "test_dump/test_vtk_writer.vtk" );

    return 0;
}
}
}
