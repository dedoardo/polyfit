// Computer the winding number of a polygon at a point
// Author: Shayan Hoshyari

#ifndef polyvec_winding_number_
#define polyvec_winding_number_

#include <Eigen/Core>

namespace polyvec {
    namespace WindingNumber {
        // Input:
        // polygon: polygon 2 x numpoints matrix
        // point: the point to compute the winding number at
        // Returns:
        // winding_number: value of winding number
        // is_trustable: is the number trustable, or are we too close to the boundary?
        // NOTE: assumes that the polygon points are sorted in CCW order.
        // multiply the answer by -1 if the order is CW.
        void compute_winding ( const Eigen::Matrix2Xd& polygon, const Eigen::Vector2d& point,  double& winding_number, bool& is_trustable );

        // Input:
        // polygon: polygon 2 x numpoints matrix
        // Returns:
        // winding_number: is ccw
       void compute_orientation ( const Eigen::Matrix2Xd& polygon,  bool &is_ccw, double &area );
       void compute_orientation ( const Eigen::Matrix2Xd& polygon,  bool &is_ccw );

    }
} // end of polyvec


#endif