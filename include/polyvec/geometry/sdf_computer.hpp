// Computes the Shape Diameter Function (sdf) for a closed polygon
// Consistent Mesh Partitioning and Skeletonisation using the Shape Diameter Function
// Author: Shayan Hoshyari

#ifndef polyvec_sdf_computer_
#define polyvec_sdf_computer_

#include <polyvec/api.hpp>

// TODO: weight the rays computing sigma and sdf
// TODO: viz accepted and not accepted rays
// TODO: compute proper average rather than just the sdf at the midpoint. 

namespace polyvec {
    class SdfComputer {
    public:
        SdfComputer ( const Eigen::Matrix2Xd& points ) ;

        bool is_ray_inside ( const real2& pt0, const real2& pt1 );

        void cast_ray_from_inside ( const real2& source,
                                    const real2& direction,
                                    bool& does_hit_from_inside,
                                    real2& dest );

        static void sample_rays_in_cone ( const real2& cone_center_dir,
                                          const double cone_angle_div2,
                                          const unsigned n_rays,
                                          std::vector<real2>& ray_dirs );


        static double compute_robust_mean ( const Eigen::VectorXd &ray_lengths );
        static void compute_robust_mean ( const Eigen::VectorXd &ray_lengths, double &robust_mean, double &mean, double &standard_dev, double &median );

        static real2 tangent2inwardnormal ( const real2 &tangent, const bool is_ccw);

        double compute_sdf( const int point_idx, const double t_value);
        void compute_sdf( const int point_idx, const double t_value, double &sdf, real2 *ray_source, std::vector<real2> *ray_dst );

    private:
        Eigen::Matrix2Xd _points;
        bool _is_ccw;
    };
}

#endif