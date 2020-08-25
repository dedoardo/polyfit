#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/geom.hpp>

void polyvec::WindingNumber::compute_winding (
    const Eigen::Matrix2Xd& polygon,
    const Eigen::Vector2d& point,
    double& winding_number,
    bool& is_trustable ) {

    is_trustable = true;
    winding_number = 0;
    const double tol = 1e-10;

    // Check that the point is not on one of the edges
    for ( unsigned i = 0 ; i < polygon.cols() ; ++i ) {
        const unsigned ip1 = ( i + 1 ) % polygon.cols();
        const Eigen::Vector2d pt0 = polygon.col ( i );
        const Eigen::Vector2d pt1 = polygon.col ( ip1 );
        const double angle = geom::angle ( ( pt0-point ).normalized(), ( pt1-point ).normalized() );
        winding_number += angle;

        // Check if we are very close to this line
        const Eigen::Vector2d dir =  ( pt1 - pt0 );
        const double boundless_t = dir.dot ( point - pt0.col ( 0 ) ) / dir.squaredNorm();
        const double t_proj =  std::min ( 1., std::max ( 0., boundless_t ) );

        if ( ( point-LineUtils::line_at ( pt0, pt1, t_proj ) ).squaredNorm() < tol ) {
            is_trustable = false;
        }

    }

    winding_number /= constants::PI2;

}

void polyvec::WindingNumber::compute_orientation (
    const Eigen::Matrix2Xd& polygon,
    bool& is_ccw,
    double& area ) {

    area = 0;

    // Check that the point is not on one of the edges
    for ( unsigned i = 0 ; i < polygon.cols() ; ++i ) {
        const unsigned ip1 = ( i + 1 ) % polygon.cols();
        const Eigen::Vector2d pt0 = polygon.col ( i );
        const Eigen::Vector2d pt1 = polygon.col ( ip1 );
        const Eigen::Vector2d tangent = ( pt1-pt0 );
        const Eigen::Vector2d normal ( -tangent.y(), tangent.x() );
        const double x_mid = ( pt0.x() + pt1.x() ) / 2.;
        area += x_mid * normal.x();
    }

    if ( area >= 0 ) {
        is_ccw = false;
    } else {
        is_ccw = true;
    }

    area = std::abs ( area );
}

void polyvec::WindingNumber::compute_orientation (
    const Eigen::Matrix2Xd& polygon,
    bool& is_ccw ) {

    double area = 0;
    polyvec::WindingNumber::compute_orientation ( polygon, is_ccw, area );
}
