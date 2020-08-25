#include <cstdio>

#include <polyvec/io/vtk_curve_writer.hpp>

namespace polyvec {
  
void
VtkCurveWriter::add_line ( const Eigen::Vector2d& p0, const Eigen::Vector2d& p1 ) {
    Eigen::Matrix2Xd polyline ( 2, 2 );
    polyline.col ( 0 ) = p0;
    polyline.col ( 1 ) = p1;

    add_polyline ( polyline );
}

void
VtkCurveWriter::add_point ( const Eigen::Vector2d& pt ) {
    _points.push_back ( pt );
}

void
VtkCurveWriter::add_polyline ( const Eigen::Matrix2Xd& polyline , const bool is_close) {
    _lines.push_back ( polyline );
    if(is_close){
        Eigen::Matrix2Xd extra(2,2) ;
        extra.col(0) =  polyline.col( (int)polyline.cols()-1);
        extra.col(1) =  polyline.col(0) ;
        _lines.push_back( extra );
    }
}

void
VtkCurveWriter::clear() {
    // Don't clear() to keep the memory
    _lines.resize ( 0 );
    _points.resize ( 0 );
}

void
VtkCurveWriter::dump ( const std::string fname ) {
    FILE* fl = fopen ( fname.c_str(), "wb" );

    if ( !fl ) {
        printf ( "File %s was not found \n", fname.c_str() );
        return;
    }

    // Count the number of points
    int n_total_points = 0;

    for ( unsigned i = 0; i < _lines.size(); ++i ) {
        n_total_points += ( int ) _lines[i].cols();
    }

    n_total_points += ( int ) _points.size();

    // write the header
    fprintf ( fl, "# vtk DataFile Version 2.0\n" );
    fprintf ( fl, "Shayan's output VTK line file \n" );
    fprintf ( fl, "ASCII\n" );
    fprintf ( fl, "DATASET POLYDATA\n" );
    fprintf ( fl, "\n" );

    // write the vertices
    fprintf ( fl, "POINTS %d float\n", n_total_points );

    for ( unsigned pid = 0; pid < _points.size(); ++pid ) {
        fprintf ( fl, "%.12g %.12g 0 \n", _points[pid].x(), _points[pid].y() );
    }

    for ( unsigned cid = 0; cid < _lines.size(); ++cid ) {
        for ( unsigned vid = 0; vid < _lines[cid].cols(); ++vid ) {
            fprintf ( fl, "%.12g %.12g 0 \n", _lines[cid].col ( vid ).x(), _lines[cid].col ( vid ).y() );
        }
    }

    fprintf ( fl, "\n" );

    //
    // write the points and lines
    //

    int point_counter = 0;

    //
    fprintf ( fl, "VERTICES %d %d \n", int ( _points.size() ), 2 * int ( _points.size() ) );

    for ( unsigned pid = 0; pid < _points.size(); ++pid ) {
        fprintf ( fl, "1 %d \n", point_counter );
        ++point_counter;
    }

    fprintf ( fl, "\n" );

    //
    fprintf ( fl, "LINES %d %d \n", int ( _lines.size() ), int ( _lines.size() ) + int ( n_total_points-_points.size() ) );

    for ( unsigned cid = 0; cid < _lines.size(); ++cid ) {
        fprintf ( fl, "%d ", int ( _lines[cid].cols() ) );

        for ( unsigned vid = 0; vid < _lines[cid].cols(); ++vid ) {
            fprintf ( fl, "%d ", point_counter );
            ++point_counter;
        }

        fprintf ( fl, "\n" );
    }

    fprintf ( fl, "\n" );



    fclose ( fl );
}

}
