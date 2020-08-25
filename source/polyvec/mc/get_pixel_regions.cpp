#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/image-segment/image_segment.hpp>
#include <polyvec/io/image.hpp>
#include <polyvec/mc/get_pixel_regions.hpp>
#include <polyvec/mc/prune_identical_pixel_polygons.hpp>
#include <polyvec/mc/get_bounding_box.hpp>
#include <polyvec/mc/random_color.hpp>

// eigen
#include <Eigen/Geometry>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )

void
post_process_pixel_regions (
    const std::vector<Eigen::Matrix2Xd>& boundaries_d,
    std::vector<Eigen::Matrix2Xi>& boundaries,
    std::vector<int>& out2in ) {


    // Convert all to int
    boundaries.resize ( boundaries_d.size() );

    for ( int i = 0; i < ( int ) boundaries_d.size(); ++i ) {
        // Make sure the double values are actually int
        Eigen::Matrix2Xd boundary_d_current = boundaries_d[i];

        auto is_really_int = [&boundary_d_current] {
            return  ( boundary_d_current.cast<int>().cast<double>() - boundary_d_current ).norm() < 1e-10 ;
        };

        if ( !is_really_int() ) {
            boundary_d_current.array() -= 0.5;
        }

        boundaries[i] = boundary_d_current.cast<int>();

        assert_break ( is_really_int() );
    }

    // Find the bbox
    Eigen::AlignedBox<int, 2> bbox = get_bounding_box( boundaries );

    // Set the minimum value to 0,0
    for ( int i = 0; i < ( int ) boundaries.size(); ++i ) {
        boundaries[i].colwise() -= bbox.min();
    }


    //
    // Orientation and duplicacy
    //

    //
    // Method 1 randomly kill one duplicate -- does not respect colors
    //  This is broken. Paths are not necessarily duplicates.
    //
    #if 0
        // Make sure that the orientation of all is CCW (when coordinate 0,0 is
        // assumed to be at lower left)
        for ( int i = 0; i < ( int ) boundaries.size(); ++i ) {
            bool is_ccw;
            // Just give it the double version
            polyvec::WindingNumber::compute_orientation ( boundaries_d[i], is_ccw );

            // Reverse if it is not ccw
            if ( !is_ccw ) {
                boundaries[i].rowwise().reverseInPlace();
            }
        }

        // Now kill the duplicate versions and set the indexing
        {
            std::vector<Eigen::Matrix2Xi> boudaries_bup = boundaries;
            prune_identical_pixel_polygons ( boudaries_bup, boundaries, out2in );
        }
    }
    #endif
    //
    // Method 2 -- prune based on orienation (ignore holes)
    // only keep ccw's
    // Assumes that Ed's code exports stuff as ccw, which seems to be the case.
    //
    {
        std::vector<Eigen::Matrix2Xi> boudaries_bup = boundaries;
        boundaries.resize ( 0 );
        out2in.resize ( 0 );

        for ( int i = 0; i < ( int ) boudaries_bup.size(); ++i ) {
            bool is_ccw;
            // Just give it the double version
            polyvec::WindingNumber::compute_orientation ( boudaries_bup[i].cast<double>(), is_ccw );

            // Reverse if it is not ccw
            if ( is_ccw ) {
                boundaries.push_back( boudaries_bup[i] );
                out2in.push_back( i );
            }
        }
    } // end of duplicate handling
}


void
get_pixel_regions_from_bmp (
    const std::string fname,
    std::vector<Eigen::Matrix2Xi>& boundaries,
    std::vector<Eigen::Vector4d>& colors ) {

    std::vector<Eigen::Matrix2Xd> boundaries_d;
    std::vector<Eigen::Vector4d> colors_d;

    // Read and process the image
    IO::Image image;
    IO::read_image ( fname.c_str(), image );
    ImageSegment::expand_and_cleanup ( image );
    ImageSegment::extract_closed_regions ( image, boundaries_d, colors_d );

    // Process regions
    std::vector<int> out2in;
    post_process_pixel_regions ( boundaries_d, boundaries, out2in );

    // set the colors
    // The duplicate boundaries might come with different colors
    // which can break this
    colors.resize ( out2in.size() );

    for ( int i = 0 ; i < ( int ) boundaries.size() ; ++i ) {
        colors[i] = colors_d[ out2in[ i ]  ] ;
    }

    // or just assign random colors for now
    // ...

}

NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
