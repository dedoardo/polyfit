// polyvec multicolor
#include <polyvec/mc/sort_points.hpp>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN(mc)

Eigen::Matrix2Xi sort_points ( const Eigen::Matrix2Xi& points ) {

    // Sort points in 2D to get a canonical order
    auto is_smaller_than = [] ( const Eigen::Vector2i &p1,
    const Eigen::Vector2i &p2 ) -> bool {

        // This would not distinguish swapping x and y values
        // const Eigen::Vector2i p1s( std::min(p1.x(), p1.y()), std::max(p1.x(), p1.y()) );
        // const Eigen::Vector2i p2s( std::min(p2.x(), p2.y()), std::max(p2.x(), p2.y()) );
        // This would
        const Eigen::Vector2i p1s = p1;
        const Eigen::Vector2i p2s = p2;

        if ( p1s.x() != p2s.x() ) {
            return p1s.x() < p2s.x();
        } else if ( p1s.y() != p2s.y() ) {
            return p1s.y() < p2s.y();
        } else {
            return false;
        }
    };

    // Put the points in a vecotr with stl iterators
    std::vector<Eigen::Vector2i> points_vec(points.cols());
    for ( int i = 0; i < ( int ) points.cols(); ++i ) {
        points_vec[i] = points.col ( i );
    }

    // Sort them 
    std::sort ( points_vec.begin(), points_vec.end(), is_smaller_than);

    // Put them back in a mX2d an return them
    Eigen::Matrix2Xi ans(2, points.cols());
    for ( int i = 0; i < ( int ) points.cols(); ++i ) {
        ans.col(i) = points_vec[ i ];
    }

    return ans;
}

NAMESPACE_END (mc)
NAMESPACE_END ( polyfit )
