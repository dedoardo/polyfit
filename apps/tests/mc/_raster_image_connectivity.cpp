#include <iostream>
//
#include <polyvec/api.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/misc.hpp>
//
#include <polyvec/mc/get_bounding_box.hpp>
#include <polyvec/mc/get_pixel_regions.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
//
#include <Eigen/Geometry>

NAMESPACE_BEGIN()
NAMESPACE_END()

NAMESPACE_BEGIN ( polyvectest )
NAMESPACE_BEGIN ( mc )

int raster_image_connectivity ( int, char** ) {
    std::vector<std::string> image_addrs;
    image_addrs.emplace_back ( POLYVEC_TEST_PATH "/mc/abc.png" );
    image_addrs.emplace_back ( POLYVEC_TEST_PATH "/mc/holes001.png" );
    // image_addrs.emplace_back ( POLYVEC_TEST_PATH "/mc/chess.png" );

    printf ( "==== Testing raster image connectivity \n" );

    polyfit::Log::open ( "test_dump/log-export-multicolor-traces.txt" );

    for ( int i = 0 ; i < (int)image_addrs.size() ; ++i ) {
        std::vector<Eigen::Matrix2Xi> boundaries;
        std::vector<Eigen::Vector4d> colors;
        polyfit::mc::get_pixel_regions_from_bmp ( image_addrs[i], boundaries, colors );

        // Get the bounding box of the whole picture
        Eigen::AlignedBox<int, 2> bbox = polyfit::mc::get_bounding_box ( boundaries );

        // Can you build a connectivity
        polyfit::mc::RasterImageConnectivity raster_conn = polyfit::mc::RasterImageConnectivity::build ( boundaries );
        const bool is_mesh_okay = raster_conn.connectivity().check_sanity_slowly ( true );
        printf ( "Mesh okay status: %d \n", ( int ) is_mesh_okay );
        raster_conn.dump_as_vtk ( ::polyvec::misc::sfmt ( "test_dump/raster_image_conn_%d_ALL.vtk", i, 0 ) ) ;

        // WRITE REGIONS
        for ( int j = 0; j < ( int ) raster_conn.n_regions(); ++j ) {
            ::polyvec::VtkCurveWriter writer;
            std::vector<int> holes;
            int outer_boundary;
            raster_conn.region_to_polygon ( j, outer_boundary, holes );
            assert_break ( outer_boundary == j );

            writer.add_polyline ( raster_conn.get_polygon_points ( outer_boundary ).cast<double>() , true);
            writer.dump ( ::polyvec::misc::sfmt ( "test_dump/raster_image_conn_%d_REGION%04d.vtk", i, j ) ) ;
        } // end of regions

        // WRITE HOLES
        for ( int j = raster_conn.n_regions(); j < ( int ) raster_conn.n_polygons(); ++j ) {
            ::polyvec::VtkCurveWriter writer;
            writer.add_polyline ( raster_conn.get_polygon_points ( j ).cast<double>() , true);
            writer.dump ( ::polyvec::misc::sfmt ( "test_dump/raster_image_conn_%d_HOLE%04d.vtk", i, j ) ) ;
        } // end of holes

        // WRITE HOLES AND REGIONS
        for ( int j = 0; j < raster_conn.n_regions(); ++j ) {
            ::polyvec::VtkCurveWriter writer;
            std::vector<int> holes;
            int outer_boundary;
            raster_conn.region_to_polygon ( j, outer_boundary, holes );
            assert_break ( outer_boundary == j );

            writer.add_polyline ( raster_conn.get_polygon_points ( outer_boundary ).cast<double>() , true);

            for ( int h : holes ) {
                Eigen::Matrix2Xi bdry = raster_conn.get_polygon_points ( h );
                bdry.rowwise().reverseInPlace();
                writer.add_polyline ( bdry.cast<double>() , true);
            }

            writer.dump ( ::polyvec::misc::sfmt ( "test_dump/raster_image_conn_%d_BOTH%04d.vtk", i, j ) ) ;
        } // end of holes  + regions

    } // end of images

    return EXIT_FAILURE;
}

NAMESPACE_END ( mc )
NAMESPACE_END ( polyvectest )
