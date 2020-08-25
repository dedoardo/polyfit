#include <polyvec/api.hpp>

#include <polyvec/mc/get_bounding_box.hpp>

#include <Eigen/Geometry>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)


Eigen::AlignedBox<int, 2>
get_bounding_box(const std::vector<Eigen::Matrix2Xi>& regions) {
   Eigen::AlignedBox<int, 2> bbox;

    for ( int i = 0; i < ( int ) regions.size(); ++i ) {
        for ( int j = 0; j < ( int ) regions[i].cols(); ++j ) {
            bbox.extend ( regions[i].col ( j ) );
        }
    }
    return bbox;
}

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)
