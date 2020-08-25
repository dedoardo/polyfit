#include <polyvec/geometry/sdf_computer.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/winding_number.hpp>

polyvec::SdfComputer::SdfComputer ( const Eigen::Matrix2Xd& points ) :
    _points ( points ) {
    WindingNumber::compute_orientation ( _points, _is_ccw );
}

bool
polyvec::SdfComputer::is_ray_inside ( const real2& pt0, const real2& pt1 ) {
    //const std::vector<double> t_values_to_test = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    const std::vector<double> t_values_to_test = {0.5};
    const double tol = 1e-8;
    bool is_trustable;
    double winding_number;

    for ( double t: t_values_to_test ) {
        const real2 point_at_t = LineUtils::line_at ( pt0, pt1, t );
        WindingNumber::compute_winding ( _points, point_at_t, winding_number, is_trustable );

        if ( is_trustable ) {
            if ( std::abs ( winding_number ) > tol ) {
                return true;
            } else {
                return false;
            }
        }
    }

    // this line is on an edge.
    return false;
}

void
polyvec::SdfComputer::cast_ray_from_inside (
    const real2& source,
    const real2& direction,
    bool& does_hit_from_inside,
    real2& dest ) {

    const double winding_number_tol = 1e-8;
    real2 direction_norm = direction.normalized();

    // perturb source slightly
    const double perturbation_eps = 1e-3;
    real2 source_pert = source + direction_norm * perturbation_eps;

    // check if the perturbed source is inside
    {
        bool is_trustable;
        double winding_number;
        WindingNumber::compute_winding ( _points, source_pert, winding_number, is_trustable );

        // We are either outside or on the boundary
        if ( ( !is_trustable ) || ( std::abs ( winding_number ) < winding_number_tol )   ) {
            does_hit_from_inside = false;
            dest.setConstant ( std::numeric_limits<double>::max() );
            return;
        }
    }

    // if perturbed is inside, find the closest location where the ray intersects
    does_hit_from_inside = false;
    dest.setConstant ( std::numeric_limits<double>::max() );
    double dist = std::numeric_limits<double>::max() ;

    for ( unsigned i = 0 ; i < _points.cols() ; ++i ) {

        // Find intersection points
        real2 intersection_point;
        bool does_intersect;
        const unsigned ip1 = ( i + 1 ) % _points.cols();
        const real2 poly_pt0 = _points.col ( i );
        const real2 poly_pt1 = _points.col ( ip1 );
        const Eigen::Vector2d ray_pt0 = source_pert;
        double dist_on_ray;
        does_intersect = geom::ray_intersect ( ray_pt0, direction, geom::segment ( poly_pt0, poly_pt1 ), dist_on_ray );

        if ( does_intersect && ( dist_on_ray > 0 ) ) {
            does_hit_from_inside = true;
            intersection_point =dist_on_ray*direction+ray_pt0;
            double current_dist = ( intersection_point-source_pert ).norm();

            if ( is_ray_inside ( ray_pt0, intersection_point ) && ( current_dist < dist ) ) {
                dest = intersection_point;
                dist = current_dist;
            }
        }

    }

    assert_break ( does_hit_from_inside == true );
}


void polyvec::SdfComputer::sample_rays_in_cone (
    const real2& cone_center_dir,
    const double cone_angle,
    const unsigned n_rays,
    std::vector<real2>& ray_dirs ) {

    ray_dirs.clear();
    ray_dirs.reserve ( n_rays );

    const real2 cone_center_dir_norm = cone_center_dir.normalized();
    const double step = 2. * cone_angle / ( n_rays - 1 );
    const double cone_center_dir_angle = geom::angle ( cone_center_dir_norm );

    for ( unsigned i = 0 ; i < n_rays ; ++i ) {
        const double ray_angle =  cone_center_dir_angle - cone_angle + i * step;
        ray_dirs.push_back ( real2 ( cos ( ray_angle ), sin ( ray_angle ) ) );
    }

}

double polyvec::SdfComputer::compute_robust_mean ( const Eigen::VectorXd& ray_lengths ) {
    double robust_mean, mean, std_dev, median;
    compute_robust_mean ( ray_lengths, robust_mean, mean, std_dev, median );
    return robust_mean;
}

void polyvec::SdfComputer::compute_robust_mean ( const Eigen::VectorXd& ray_lengths, double& robust_mean, double& mean, double& standard_dev, double& median ) {

    assert_break ( ray_lengths.size() > 1 );

    // Find standard deviation
    mean = ray_lengths.mean();
    const double sigma2 = ( ( ray_lengths.array() - mean ).square().sum() ) / ( ray_lengths.size() - 1 );
    standard_dev = sqrt ( sigma2 );

    // Find median
    Eigen::VectorXd ray_lengths_sorted = ray_lengths;
    std::sort ( ray_lengths_sorted.data(), ray_lengths_sorted.data()+ray_lengths.size() );

    if ( ray_lengths.size() % 2 == 1 ) {
        median = ray_lengths_sorted[ ray_lengths.size() / 2];
    } else {
        median = 0.5 * ( ray_lengths_sorted[ ray_lengths.size() / 2 -1  ] + ray_lengths_sorted[ ray_lengths.size() / 2 ] );
    }

    // Find mean of values withing one standard dev of mean
    robust_mean = 0;
    unsigned n_accepted_rays= 0;

    for ( unsigned i = 0 ; i < ray_lengths.size() ; ++i ) {
        const double len = ray_lengths[i];

        if ( std::abs ( len-median ) < standard_dev ) {
            robust_mean += len;
            ++n_accepted_rays;
        }
    }

    robust_mean /= n_accepted_rays;
}


polyvec::real2 polyvec::SdfComputer::tangent2inwardnormal ( const real2& tang, const bool is_ccw ) {
    return  real2 ( -tang.y(), tang.x() ) * ( is_ccw ? 1 : -1. );
}


double polyvec::SdfComputer::compute_sdf ( const int point_idx, const double t_value ) {
    double ans;
    compute_sdf ( point_idx, t_value, ans, nullptr, nullptr );
    return ans;
}

void polyvec::SdfComputer::compute_sdf ( const int point_idx, const double t_value, double& sdf, real2* __ray_source, std::vector<real2>* __ray_dests ) {

    const double tol = 1e-10;
    const int n_rays = 40;
    const double cone_angle_div2 = geom::radians(60);

    assert_break ( t_value >= 0 );
    assert_break ( t_value < 1-tol );

    // next and prev points
    const int point_next  = ( point_idx + 1 ) % _points.cols();
    const int point_prev  = ( point_idx - 1 + _points.cols() ) % _points.cols();

    // find ray source
    const real2 ray_source  = LineUtils::line_at ( _points.col ( point_idx ), _points.col ( point_next ), t_value );

    if ( __ray_source ) {
        *__ray_source = ray_source;
    }

    // find normal
    real2 normal;

    if ( std::abs ( t_value ) < tol ) {
        const real2 tangent_1 =  _points.col ( point_next ) - _points.col ( point_idx ) ;
        const real2 tangent_2 =  _points.col ( point_idx ) - _points.col ( point_prev ) ;
        const real2 tangent = ( tangent_1+tangent_2 ).normalized();
        normal = tangent2inwardnormal ( tangent, _is_ccw );
    } else {
        const real2 tangent = ( _points.col ( point_next ) - _points.col ( point_idx ) ).normalized();
        normal = tangent2inwardnormal ( tangent, _is_ccw );
    }

    // find directions
    std::vector<real2> ray_dirs;
    sample_rays_in_cone ( normal, cone_angle_div2, n_rays, ray_dirs );

    // find intersections and lengths
    std::vector< double > ray_lengths;

    if ( __ray_dests ) {
        __ray_dests->clear();
    }

    for ( unsigned i = 0 ; i < n_rays ; ++i ) {
        real2 dest;
        bool does_hit;
        cast_ray_from_inside ( ray_source, ray_dirs[i], does_hit, dest );

        if ( does_hit ) {
            ray_lengths.push_back ( ( dest-ray_source ).norm() );

            if ( __ray_dests ) {
                __ray_dests->push_back ( dest );
            }
        }
    }

    // find robust average
    if ( ray_lengths.size() < 2 ) {
        sdf = std::numeric_limits<double>::max();
    } else {
        sdf = compute_robust_mean ( Eigen::VectorXd::ConstMapType ( ray_lengths.data(), ray_lengths.size() ) );
    }
}
