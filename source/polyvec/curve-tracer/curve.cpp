// Author: A tired Shayan Hoshyari.
// Header
#include <polyvec/curve-tracer/curve.hpp>

// Polyvec
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/misc.hpp>

// Eigen
#include <unsupported/Eigen/Polynomials>

namespace polyvec {
    const char* curve_type_to_string(const GlobFitCurveType type) {
        switch (type) {
        case GLOBFIT_CURVE_LINE:
            return "line";
        case GLOBFIT_CURVE_ARC:
            return "arc";
        case GLOBFIT_CURVE_BEZIER:
            return "bezier";
        default:
            return "unknown";
        }
    }


    double GlobFitCurve::t_end = 1.;

    //constexpr double GlobFitCurve::t_end;

	// SHAYAN: Why on earth are this two functions here!?!?!?!?!??!
    double Curve::distance_to ( const real2& point ) const {
        return ( point - pos ( project ( point ) ) ).squaredNorm();
    }

    double Curve::distance_to_sq ( const real2& point ) const {
        return sqrt ( distance_to ( point ) );
    }

// =============================================================
//                        GlobFitCurve
// =============================================================
    double GlobFitCurve::distance_to ( const Eigen::Vector2d& point ) const {
        return ( point - pos ( project ( point ) ) ).squaredNorm();
    }

    double GlobFitCurve::distance_to_sq ( const Eigen::Vector2d& point ) const {
        return sqrt ( distance_to ( point ) );
    }

    void GlobFitCurve::approximate_intersect (
        const real2& line_src,
        const real2& line_dst,
        std::vector<double>& line_at,
        std::vector<double>&  this_curve_at ) {

        std::vector<int> polygon_seg_no;

        const Eigen::VectorXd tesselationt = get_tesselationt();
        LineUtils::intersect_polyline ( line_src, line_dst, get_tesselation2(), line_at, polygon_seg_no, this_curve_at );

        // Fix the t value from polyline segment to curve by linear interpolation
        for ( int i = 0 ; i < int ( polygon_seg_no.size() ); ++i  ) {
            assert_break ( this_curve_at[i] >= 0. );
            assert_break ( this_curve_at[i] <= 1. );
            assert_break ( polygon_seg_no[i] >= 0 );
            assert ( polygon_seg_no[i] <  ( int ) tesselationt.size()-1 );
            this_curve_at[i] = tesselationt[polygon_seg_no[i]]* ( 1. - this_curve_at[i] ) + tesselationt[polygon_seg_no[i]+1]* ( this_curve_at[i] );
        }
    }

	int GlobFitCurve::approximate_intersect_ray_closest(
		const real2& ray_o,
		const real2& ray_d,
		double& line_at,
		double& this_curve_at
	) {
		std::vector<double> line_ats;
		std::vector<double> this_curve_ats;
		approximate_intersect(ray_o, ray_o + ray_d, line_ats, this_curve_ats);
		PF_ASSERT(line_ats.size() == this_curve_ats.size());

		line_at = INFINITY;
		this_curve_at = INFINITY;
		for (size_t i = 0; i < line_ats.size(); ++i) {
			if (abs(line_ats[i]) < abs(line_at)) {
				line_at = line_ats[i];
				this_curve_at = this_curve_ats[i];
			}
		}

		return line_ats.size();
	}    

	std::pair<double, double> GlobFitCurve::split_t(double t, double t_split)
	{
		return std::make_pair(t / t_split, 1 - (1 - t) / (1 - t_split));
	}

}  // namespace polyvec
