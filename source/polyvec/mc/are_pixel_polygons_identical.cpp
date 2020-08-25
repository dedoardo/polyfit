// polyvec multicolor
#include <polyvec/mc/are_pixel_polygons_identical.hpp>
#include <polyvec/mc/sort_points.hpp>


NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )


bool are_pixel_polygons_identical ( const Eigen::Matrix2Xi& polygon1, const Eigen::Matrix2Xi& polygon2 ) {
    const Eigen::Matrix2Xi polygon1_sorted = sort_points ( polygon1 );
    const Eigen::Matrix2Xi polygon2_sorted = sort_points ( polygon2 );

    if ( polygon1_sorted.cols() == polygon2_sorted.cols() ) {
        const Eigen::Matrix2Xi diff = polygon1_sorted - polygon2_sorted;
        const bool are_same =   diff.norm() < 1e-12;
        return are_same;
    } else {
        return false;
    }
}

NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
