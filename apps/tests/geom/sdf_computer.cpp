#include <sstream>
#include <iostream>

#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/sdf_computer.hpp>
#include <polyvec/geometry/winding_number.hpp>

#include <dev/canvas.hpp>

using namespace polyvec;

namespace {

    const std::string apple_pts = R"(19.5  0.5
20.5  0.5
21.5  0.5
22.5  0.5
22.5  1.5
22.5  2.5
22.5  3.5
22.5  4.5
21.5  4.5
21.5  5.5
20.5  5.5
20.5  6.5
19.5  6.5
18.5  6.5
18.5  7.5
17.5  7.5
17.5  8.5
17.5  9.5
18.5  9.5
19.5  9.5
19.5  8.5
20.5  8.5
21.5  8.5
22.5  8.5
23.5  8.5
24.5  8.5
25.5  8.5
25.5  9.5
26.5  9.5
27.5  9.5
27.5 10.5
28.5 10.5
28.5 11.5
29.5 11.5
29.5 12.5
30.5 12.5
30.5 13.5
30.5 14.5
30.5 15.5
30.5 16.5
30.5 17.5
30.5 18.5
30.5 19.5
30.5 20.5
30.5 21.5
30.5 22.5
29.5 22.5
29.5 23.5
29.5 24.5
29.5 25.5
28.5 25.5
28.5 26.5
27.5 26.5
27.5 27.5
27.5 28.5
26.5 28.5
26.5 29.5
25.5 29.5
25.5 30.5
24.5 30.5
23.5 30.5
23.5 31.5
22.5 31.5
22.5 32.5
21.5 32.5
20.5 32.5
19.5 32.5
18.5 32.5
17.5 32.5
17.5 31.5
16.5 31.5
15.5 31.5
15.5 32.5
14.5 32.5
13.5 32.5
12.5 32.5
11.5 32.5
10.5 32.5
10.5 31.5
 9.5 31.5
 9.5 30.5
 8.5 30.5
 7.5 30.5
 7.5 29.5
 6.5 29.5
 6.5 28.5
 5.5 28.5
 5.5 27.5
 5.5 26.5
 4.5 26.5
 4.5 25.5
 3.5 25.5
 3.5 24.5
 3.5 23.5
 3.5 22.5
 2.5 22.5
 2.5 21.5
 2.5 20.5
 2.5 19.5
 2.5 18.5
 2.5 17.5
 2.5 16.5
 2.5 15.5
 2.5 14.5
 2.5 13.5
 2.5 12.5
 3.5 12.5
 3.5 11.5
 4.5 11.5
 4.5 10.5
 5.5 10.5
5.5 9.5
6.5 9.5
6.5 8.5
7.5 8.5
8.5 8.5
9.5 8.5
10.5  8.5
11.5  8.5
12.5  8.5
12.5  9.5
13.5  9.5
14.5  9.5
14.5  8.5
14.5  7.5
13.5  7.5
13.5  6.5
12.5  6.5
11.5  6.5
11.5  5.5
10.5  5.5
10.5  4.5
9.5 4.5
9.5 3.5
9.5 2.5
10.5  2.5
11.5  2.5
12.5  2.5
12.5  3.5
13.5  3.5
13.5  4.5
14.5  4.5
15.5  4.5
16.5  4.5
16.5  3.5
16.5  2.5
17.5  2.5
17.5  1.5
18.5  1.5
18.5  0.5)";

    const std::string evil_cat_pts = R"(
  25.5 29.5
26.5 29.5
26.5 28.5
27.5 28.5
27.5 27.5
27.5 26.5
27.5 25.5
27.5 24.5
27.5 23.5
26.5 23.5
26.5 22.5
26.5 21.5
25.5 21.5
25.5 20.5
25.5 19.5
24.5 19.5
24.5 18.5
23.5 18.5
23.5 17.5
24.5 17.5
25.5 17.5
25.5 18.5
26.5 18.5
27.5 18.5
27.5 19.5
28.5 19.5
28.5 20.5
28.5 21.5
29.5 21.5
29.5 22.5
29.5 23.5
29.5 24.5
30.5 24.5
30.5 25.5
30.5 26.5
30.5 27.5
29.5 27.5
29.5 28.5
29.5 29.5
29.5 30.5
28.5 30.5
28.5 31.5
27.5 31.5
26.5 31.5
25.5 31.5
25.5 32.5
24.5 32.5
23.5 32.5
22.5 32.5
21.5 32.5
20.5 32.5
19.5 32.5
18.5 32.5
17.5 32.5
16.5 32.5
15.5 32.5
14.5 32.5
13.5 32.5
12.5 32.5
11.5 32.5
10.5 32.5
 9.5 32.5
 9.5 31.5
 9.5 30.5
10.5 30.5
10.5 29.5
10.5 28.5
10.5 27.5
10.5 26.5
10.5 25.5
 9.5 25.5
 9.5 24.5
 9.5 23.5
 8.5 23.5
 8.5 22.5
 8.5 21.5
 7.5 21.5
 7.5 20.5
 7.5 19.5
 6.5 19.5
 6.5 18.5
 6.5 17.5
 5.5 17.5
 5.5 16.5
 5.5 15.5
 5.5 14.5
 4.5 14.5
 4.5 13.5
 4.5 12.5
 4.5 11.5
 4.5 10.5
4.5 9.5
3.5 9.5
3.5 8.5
2.5 8.5
2.5 7.5
2.5 6.5
3.5 6.5
3.5 5.5
3.5 4.5
3.5 3.5
4.5 3.5
4.5 2.5
4.5 1.5
5.5 1.5
5.5 2.5
6.5 2.5
6.5 3.5
7.5 3.5
8.5 3.5
9.5 3.5
10.5  3.5
11.5  3.5
11.5  2.5
12.5  2.5
12.5  1.5
13.5  1.5
13.5  2.5
13.5  3.5
14.5  3.5
14.5  4.5
14.5  5.5
15.5  5.5
15.5  6.5
15.5  7.5
15.5  8.5
14.5  8.5
14.5  9.5
14.5 10.5
13.5 10.5
13.5 11.5
14.5 11.5
14.5 12.5
15.5 12.5
15.5 13.5
16.5 13.5
16.5 14.5
17.5 14.5
18.5 14.5
18.5 15.5
19.5 15.5
20.5 15.5
20.5 16.5
21.5 16.5
21.5 17.5
21.5 18.5
22.5 18.5
22.5 19.5
23.5 19.5
23.5 20.5
23.5 21.5
23.5 22.5
24.5 22.5
24.5 23.5
24.5 24.5
24.5 25.5
24.5 26.5
24.5 27.5
24.5 28.5
23.5 28.5
23.5 29.5
24.5 29.5)";

    const std::string apple_polypath = R"(
17.5  7.5
17.5  9.5
22.5  8.5
27.5  9.5
30.5 12.5
30.5 22.5
27.5 28.5
22.5 32.5
17.5 32.5
16.5 31.5
15.5 32.5
10.5 32.5
 5.5 28.5
 2.5 22.5
 2.5 12.5
6.5 8.5
12.5  8.5
14.5  9.5
14.5  7.5
9.5 4.5
9.5 2.5
12.5  2.5
13.5  4.5
16.5  4.5
18.5  0.5
22.5  0.5
22.5  4.5 )";


    const std::string evil_cat_polypath = R"(
  9.5 32.5
10.5   28
 4.5 14.5
4.5 9.5
2.5 7.5
4.5 1.5
6.5 3.5
11.5  3.5
13.5  1.5
15.5  5.5
15.5  8.5
13.5 11.5
15.5 13.5
20.5 16.5
22.5 19.5
24.5 25.5
23.5 29.5
26.5 29.5
27.5 28.5
27.5 23.5
23.5 17.5
25.5 17.5
28.5 19.5
30.5   26
28.5 31.5)";

    Eigen::Matrix2Xd create_polygon ( const std::string coord_str ) {
        std::istringstream ss ( coord_str, std::ios::in );
        Eigen::Matrix2Xd ans;
        unsigned n_points = 0;

        // std::cerr << coord_str;

        for ( ;; ) {
            double x, y;
            ss >> x >> y;

            if ( !ss.fail() ) {
                ++n_points;
            } else {
                break;
            }
        }

        ss.clear();
        ss.str ( coord_str );
        ans.resize ( 2, n_points );
        n_points = 0;

        for ( ;; ) {
            double x, y;
            ss >> x >> y;

            if ( !ss.fail() ) {
                ans.col ( n_points ) << x, y;
                ++n_points;
            } else {
                break;
            }
        }


        return ans;
    }

    void
    cast_one_ray ( const  Eigen::Matrix2Xd& polygon, const std::string name, const unsigned pt_id, const real2& dir ) {
        printf ( "- Testing %s, pt %d, dir %g %g \n", name.c_str(), ( int ) pt_id, dir ( 0 ), dir ( 1 ) );
        SdfComputer sdf_computer ( polygon );

        VtkCurveWriter writer;
        writer.add_polyline ( polygon );
        writer.add_line ( polygon.col ( polygon.cols()-1 ), polygon.col ( 0 ) );
        writer.add_point ( polygon.col ( pt_id ) );

        bool does_hit=false;
        real2 dest;
        sdf_computer.cast_ray_from_inside ( polygon.col ( pt_id ), dir, does_hit, dest );

        if ( does_hit ) {
            writer.add_point ( dest );
            writer.add_line ( polygon.col ( pt_id ), dest );
        } else {
            writer.add_line ( dir*4+polygon.col ( pt_id ), polygon.col ( pt_id ) );
        }

        writer.dump ( polyvec_str ( "test_dump/sdf_computer_SINGLE_ray_" << name << "_" << pt_id << ".vtk" ) );
    }

    void
    sample_cone ( const real2& dir, const double angle ) {
        static int n_called = 0;
        VtkCurveWriter writer;
        real2 zero ( 0, 0 );

        printf ( "- Sampling ray dir (%g, %g) angle %g \n", dir.x(), dir.y(), angle );

        writer.add_point ( zero );

        std::vector<real2> ray_dirs;
        SdfComputer::sample_rays_in_cone ( dir, angle, 40, ray_dirs );

        for ( unsigned i = 0 ; i < ray_dirs.size() ; ++i ) {
            writer.add_line ( zero, zero + ray_dirs[i] );
        }

        writer.dump ( polyvec_str ( "test_dump/sample_ray" << "_" << n_called << ".vtk" ) );
        ++n_called;
    }


    void
    cast_ray_from_all_corners ( const  Eigen::Matrix2Xd& polygon, const std::string name ) {
        printf ( "- Testing %s \n", name.c_str() );
        const unsigned n_rays = 20;
        SdfComputer sdf_computer ( polygon );

        bool is_ccw;
        const double cone_angle = polyvec::geom::radians ( 90 );
        polyvec::WindingNumber::compute_orientation ( polygon, is_ccw );


        for ( unsigned i = 0 ; i < polygon.cols() ; ++i ) {

            VtkCurveWriter writer;
            writer.add_polyline ( polygon );
            writer.add_line ( polygon.col ( polygon.cols()-1 ), polygon.col ( 0 ) );
            writer.add_point ( polygon.col ( i ) );

            real2 pt0 = polygon.col ( i );
            real2 ptnext = polygon.col ( ( i + 1 )  % polygon.cols() );
            real2 ptprev = polygon.col ( ( i - 1 + polygon.cols() ) % polygon.cols() );
            real2 tang = ( ( ptnext-pt0 ).normalized() + ( pt0-ptprev ).normalized() ).normalized();
            real2 normal = SdfComputer::tangent2inwardnormal ( tang, is_ccw );
            writer.add_point ( normal * 2 + polygon.col ( i ) );
            writer.add_point ( normal * 1 + polygon.col ( i ) );
            writer.add_point ( normal * 3 + polygon.col ( i ) );
            writer.add_line ( polygon.col ( i ), normal * 3 + polygon.col ( i ) );

            std::vector<real2> ray_dirs;
            SdfComputer::sample_rays_in_cone ( normal, cone_angle, n_rays, ray_dirs );

            for ( unsigned j = 0 ; j < ray_dirs.size() ; ++j ) {
                real2 dest ( -1000, -1000 );
                real2 direction = ray_dirs[j];
                bool does_hit=false;
                sdf_computer.cast_ray_from_inside ( polygon.col ( i ), direction, does_hit, dest );

                if ( does_hit ) {
                    writer.add_point ( dest );
                    writer.add_line ( polygon.col ( i ), dest );
                }
            }

            writer.dump ( polyvec_str ( "test_dump/sdf_computer_MANY_ray_" << name << "_" << i << ".vtk" ) );
        }
    }

    void find_robust_mean ( const Eigen::VectorXd& values, const double expectation ) {
        double robust_mean, mean, std_dev, median;
        SdfComputer::compute_robust_mean ( values, robust_mean, mean, std_dev, median );
        const double diff = std::abs ( robust_mean-expectation );
        printf ( "-- Robust mean,  expected: %10.5g  acheived: %10.10g  diff %10.10g \n",  expectation, robust_mean, diff );
        printf ( "      mean: %10.5g, std_dev: %10.5g, median: %10.5g \n", mean, std_dev, median );
    }

    void find_sdf ( const std::string outname, const Eigen::Matrix2Xd& polygpath ) {

#if 0
        std::vector<real2> voxels;
        voxels.reserve ( voxels_mat2d.cols() );

        for ( unsigned i = 0 ; i < voxels_mat2d.cols() ; ++i ) {
            voxels.push_back ( voxels_mat2d.col ( i ) );
        }

        // ---------------------------

        Polygon polygon;
        polygon.read_from_points ( outname, voxels.data(), ( polyvec::index ) voxels.size(), true );
        PolyGraph graph ( polygon.segment() );

        polyvec::EdgeError error;
        error.clear();
        error.add ( "angle", errors::angle, 1.0, geom::radians ( 135 ) );
        error.add ( "accuracy", errors::accuracy, 1.0, 0. );
        error.add ( "continuity", errors::continuity, 0.5, geom::radians ( 120 ) );
        error.add ( "inflection", errors::inflection, 0.25, geom::radians ( 90 ) );
        error.normalize_weights();

        Pathfinder pathfinder ( &graph, &error );
        pathfinder.find_shortest();

        PolyPath path = pathfinder.shortest();

        // ---------------------------

        Eigen::Matrix2d path_pts ( 2, path.size() );

        for ( polyvec::index cid = 0 ; cid < ( polyvec::index ) path.size(); ++cid ) {
            path_pts.col ( cid ) = path.point ( cid );
        }

#endif

        SdfComputer sdf_computer ( polygpath );

        // ---------------------------

        const int extra_plots = 2;
        const double pdf_dw = 1.;
        const double pdf_dh = 1. / double ( polygpath.cols() + extra_plots );
        PDF* pdf = vgfx::pdf_open ( polyvec_str ( "test_dump/sdf_computer_"<<outname<<".svg" ), 1./pdf_dh, pdf_dw );

        // SDF VALUES
        {
            for ( polyvec::index cid = 0 ; cid < ( polyvec::index ) polygpath.cols(); ++cid ) {
                double sdf_value;
                real2 src;
                std::vector<real2> dests;
                sdf_computer.compute_sdf ( ( int ) cid, 0.5, sdf_value, &src, &dests );

                draw::line ( polygpath.col ( cid ), polygpath.col ( ( cid+1 ) % polygpath.cols() ), Style::outline ( colors::gray, 3 ) );
                draw::text ( src, polyvec_str ( sdf_value ), draw::font_pdf, Style::text() );
            }

            vgfx::pdf_draw ( pdf, { 0, ( 0 ) *pdf_dh, pdf_dw, pdf_dh } );
        }


        // LOG OF SDF VALUES
        {
            for ( polyvec::index cid = 0 ; cid < ( polyvec::index ) polygpath.cols(); ++cid ) {
                double sdf_value;
                real2 src;
                std::vector<real2> dests;
                sdf_computer.compute_sdf ( ( int ) cid, 0.5, sdf_value, &src, &dests );

                draw::line ( polygpath.col ( cid ), polygpath.col ( ( cid+1 ) % polygpath.cols() ), Style::outline ( colors::gray, 3 ) );
                draw::text ( src, polyvec_str ( log(sdf_value) ), draw::font_pdf, Style::text() );
            }

            vgfx::pdf_draw ( pdf, { 0, ( 1 ) *pdf_dh, pdf_dw, pdf_dh } );
        }


        for ( polyvec::index cid = 0 ; cid < ( polyvec::index ) polygpath.cols(); ++cid ) {
            double sdf_value;
            real2 src;
            std::vector<real2> dests;

            sdf_computer.compute_sdf ( ( int ) cid, 0.5, sdf_value, &src, &dests );

            for ( polyvec::index cid = 0 ; cid < ( polyvec::index ) polygpath.cols(); ++cid ) {
                draw::line ( polygpath.col ( cid ), polygpath.col ( ( cid+1 ) % polygpath.cols() ), Style::outline ( colors::gray, 3 ) );
            }

            for ( polyvec::index hit_id = 0 ; hit_id < ( polyvec::index ) dests.size() ; ++hit_id ) {
                draw::text ( src, polyvec_str ( sdf_value ), draw::font_pdf, Style::text() );
                draw::line ( src, dests[hit_id], Style::outline ( colors::enemys_blood, 1 ) );
            }

            vgfx::pdf_draw ( pdf, { 0, ( cid+extra_plots ) *pdf_dh, pdf_dw, pdf_dh } );
        }


        vgfx::pdf_close ( pdf );
    }

} // end of anonymus

namespace polyvectest {
    namespace geom {
        int test_sdf_computer ( int /*argc*/, char** /*argv*/ ) {
            cast_one_ray ( create_polygon ( evil_cat_pts ), "evil_cat", 0, real2 ( 0, 1 ) );
            //
            sample_cone ( real2 ( 1, 1 ), ::polyvec::geom::radians ( 20 ) );
            sample_cone ( real2 ( 1, 1 ), ::polyvec::geom::radians ( 120 ) );
            sample_cone ( real2 ( -1, 1 ), ::polyvec::geom::radians ( 20 ) );
            sample_cone ( real2 ( -1, 1 ), ::polyvec::geom::radians ( 120 ) );
            sample_cone ( real2 ( 1, -1 ), ::polyvec::geom::radians ( 20 ) );
            sample_cone ( real2 ( 1, -1 ), ::polyvec::geom::radians ( 120 ) );
            //
            // cast rays in evil cat
            //cast_ray_from_all_corners ( create_polygon ( evil_cat_pts ), "evil_cat" );
            // check changing ori
            //cast_ray_from_all_corners ( create_polygon ( evil_cat_pts ).rowwise().reverse().eval(), "evil_cat_rev" );
            //cast_ray_from_all_corners ( create_polygon ( apple_pts ), "apple" );
            //
            {
                Eigen::VectorXd values ( 5 );
                values << 1, 2, 3, 4, 5;
                find_robust_mean ( values, 3 );
            }
            {
                Eigen::VectorXd values ( 6 );
                values << 1, 2, 3, 4, 5, 6;
                find_robust_mean ( values, 3.5 );
            }
            {
                Eigen::VectorXd values ( 5 );
                values << 4, 7, 9, 100, -2;
                find_robust_mean ( values, 4.5 );
            }
            //
            find_sdf ( "evil_cat", create_polygon ( evil_cat_polypath ) ) ;

            return 0;
        } // end of function
    } // end of share
} // end of deadline code test