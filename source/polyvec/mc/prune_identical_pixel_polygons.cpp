#include <polyvec/mc/prune_identical_pixel_polygons.hpp>
#include <polyvec/mc/are_pixel_polygons_identical.hpp>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )

void prune_identical_pixel_polygons (
    const std::vector<Eigen::Matrix2Xi>& in,
    std::vector<Eigen::Matrix2Xi>& out,
    std::vector<int>& out2in ) {

    out.resize ( 0 );
    out2in.resize ( 0 );

    for ( int i = 0 ; i < ( int ) in.size() ; ++i ) {
        bool is_duplicate = false;

        for ( int j = 0 ; j < ( int ) out.size() ; ++j ) {
            const bool are_identical = are_pixel_polygons_identical ( in[i], out[j] );

            if ( are_identical ) {
                is_duplicate = true;
                break;
            }
        }

        if ( !is_duplicate ) {
            out.push_back ( in[i] );
            out2in.push_back ( i );
        }
    }

}


NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
