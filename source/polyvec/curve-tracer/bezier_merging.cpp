
#include <iostream>
#include <map>
#include <fstream>

#include <unsupported/Eigen/KroneckerProduct>
#include <Eigen/Geometry>

#include <polyvec/curve-tracer/bezier_merging.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/shortest-path/dijkstra.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>


NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( BezierMerging )

// ===========================================
// Control points parameterization
// ===========================================

void
PointParameters::is_potracable ( bool& result, Eigen::Vector2d& oo, double& alpha, double& beta ) const {

#define undoable() do{ result = false; return; }while(0)

    // const double norm_cap = 10000;

    // if ( control_points.norm() > norm_cap ) {
    //     undoable();
    // }



    bool is_inflection_okay =  !AngleUtils::have_opposite_convexity_with_tol (
                                   control_points.col ( 0 ),
                                   control_points.col ( 1 ),
                                   control_points.col ( 2 ),
                                   control_points.col ( 3 )
                               );

    if ( !is_inflection_okay ) {
        undoable();
    }

    double lhs_at, rhs_at;
    bool do_intersect = ::polyvec::LineUtils::intersect (
                            control_points.col ( 0 ),  control_points.col ( 1 ),
                            control_points.col ( 3 ),  control_points.col ( 2 ),
                            lhs_at,
                            rhs_at );

    bool is_intersection_okay = do_intersect && ( lhs_at >= 1. )  && ( rhs_at >= 1. );

    if ( !is_intersection_okay ) {
        undoable();
    }

    oo = ::polyvec::LineUtils::line_at ( control_points.col ( 0 ),  control_points.col ( 1 ), lhs_at );
    alpha = 1. / lhs_at;
    beta = 1. / rhs_at;
    result = true;

#undef undoable
}


bool
PointParameters::is_potracable () const {
    Eigen::Vector2d oy;
    double alpha, beta;
    bool potracable;
    is_potracable ( potracable, oy, alpha, beta );
    return potracable;
}

PotraceParameters
PointParameters::as_potrace_params() const {
    bool potracable;
    Eigen::Vector2d oy;
    double alpha, beta;
    PotraceParameters ans;

    is_potracable ( potracable, oy, alpha, beta );

    if ( potracable ) {
        ans.alpha = alpha;
        ans.beta = beta;
        Mat23 p0, p1;
        //
        p0.col ( 0 ) << -1, 0;
        p0.col ( 1 ) << 0, 1;
        p0.col ( 2 ) << 1, 0;
        //
        p1.col ( 0 ) << control_points.col ( 0 );
        p1.col ( 1 ) << oy;
        p1.col ( 2 ) << control_points.col ( 3 );
        //
        find_affine_transformation ( p0, p1, ans.A, ans.b );
    } else {
        assert_break ( 0 );
    }

    return ans;
}

// ===========================================
// Potrace parameterization
// ===========================================

double
PotraceParameters::get_areaz() const {
    return 3./10* ( 2*alpha + 2*beta - alpha*beta );
}


double
PotraceParameters::get_areay() const {
    return std::abs ( polyvec::Num::determinant ( A ) ) *get_areaz();
}

PointParameters
PotraceParameters::get_zz() const {
    Mat24 ans;
    ans.col ( 0 ) << -1, 0;
    ans.col ( 1 ) << -1+alpha, alpha;
    ans.col ( 2 ) << 1-beta, beta;
    ans.col ( 3 ) << 1, 0;

    return {ans};
}

Eigen::Vector2d
PotraceParameters::get_oz() {
    return Eigen::Vector2d ( 0, 1 );
}

PointParameters
PotraceParameters::get_yy() const {
    return { ( A*get_zz().control_points ).colwise() + b};
}

Eigen::Vector2d
PotraceParameters::get_oy() const {
    return A*get_oz() + b;
}

PotraceParameters
PotraceParameters::equiparamed() const {
    PotraceParameters ans;
    ans = *this;
    ans.alpha = ans.beta = 2 - std::sqrt ( ( 2-alpha ) * ( 2-beta ) );
    return ans;
}

bool
PotraceParameters::is_ccw() const {
    return polyvec::Num::determinant ( A ) < 0;
}


// stolen from potrace <3
// calculate (p1-p0)x(p3-p2)
namespace {
	double
		_cprod(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3) {
		double x1, y1, x2, y2;

		x1 = p1.x() - p0.x();
		y1 = p1.y() - p0.y();
		x2 = p3.x() - p2.x();
		y2 = p3.y() - p2.y();

		return x1 * y2 - x2 * y1;
	}

	// calculate the point t in [0..1] on the (convex) bezier curve
	//   (p0,p1,p2,p3) which is tangent to q1-q0. Return -1.0 if there is no
	//   solution in [0..1].
	double
		_tangent(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3, const Eigen::Vector2d& q0,
			const Eigen::Vector2d& q1) {
		double A, B, C;   /* (1-t)^2 A + 2(1-t)t B + t^2 C = 0 */
		double a, b, c;   /* a t^2 + b t + c = 0 */
		double d, s, r1, r2;

		A = _cprod(p0, p1, q0, q1);
		B = _cprod(p1, p2, q0, q1);
		C = _cprod(p2, p3, q0, q1);

		a = A - 2 * B + C;
		b = -2 * A + 2 * B;
		c = A;

		d = b * b - 4 * a * c;

		if ((std::abs(a) < 1e-12)
			&& (std::abs(b) > 1e-12)
			&& (-c / b >= 0)
			&& (-c / b <= 1.)) {
			return -c / b;
		}

		if ((std::abs(a) < 1e-12) || (d < 1e-12)) {
			return -1.0;
		}

		s = sqrt(d);

		r1 = (-b + s) / (2 * a);
		r2 = (-b - s) / (2 * a);

		if (r1 >= 0 && r1 <= 1) {
			return r1;
		}
		else if (r2 >= 0 && r2 <= 1) {
			return r2;
		}
		else {
			return -1.0;
		}
	}
}

void
PotraceParameters::t_for_tangent ( const Eigen::Vector2d& tangent_y, double& t, bool& success ) const {
    const Mat24 yy = get_yy().control_points;
    t = _tangent ( yy.col ( 0 ), yy.col ( 1 ), yy.col ( 2 ), yy.col ( 3 ), Eigen::Vector2d::Zero(), tangent_y.normalized() );

    if ( t < 0 ) {
        success = false;
    } else {
        success = true;
    }

# if 0 // my own attempt. Potrace seems nicer.
    const double pi = ::polyvec::constants::PI;
    // Find the tangent in the reference coordiantes
    Eigen::Vector2d tangent_z = Num::solve_linear_system ( A, tangent_y );
    const Eigen::Vector2d zz = get_zz().control_points;
    assert_break ( tangent_z.norm() > 1e-10 );
    // Now find the angle of the tangent with the x axis
    // also find the begin and end angle of the bezier
    const double angle = std::atan2 ( tangent_z.y(), tangent_z.x() );
    const double begin_angle = std::atan2 ( zz.col ( 1 ).y() - zz.col ( 0 ).y(), zz.col ( 1 ).x() - zz.col ( 0 ).x() );
    const double end_angle = std::atan2 ( zz.col ( 3 ).y() - zz.col ( 2 ).y(), zz.col ( 3 ).x() - zz.col ( 2 ).x() );
    assert_break ( std::abs ( begin_angle-pi/4 ) < 1e-5 );
    assert_break ( std::abs ( end_angle+pi/4 ) < 1e-5 );
#endif
}

// ===========================================
// Functions
// ===========================================

void
find_affine_transformation ( const Mat23& points0, const Mat23& points1, Eigen::Matrix2d& A, Eigen::Vector2d& b ) {
    Eigen::Matrix<double, 6, 6> LHS;
    Eigen::Matrix<double, 6, 1> RHS;
    Eigen::Matrix<double, 6, 1> unknowns;

    Eigen::Matrix2d eye =  Eigen::Matrix2d::Identity();

    LHS.block<2, 4> ( 0, 0 ) = Eigen::kroneckerProduct ( points0.col ( 0 ).transpose(), eye );
    LHS.block<2, 4> ( 2, 0 ) = Eigen::kroneckerProduct ( points0.col ( 1 ).transpose(), eye );
    LHS.block<2, 4> ( 4, 0 ) = Eigen::kroneckerProduct ( points0.col ( 2 ).transpose(), eye );
    //
    LHS.block<2, 2> ( 0, 4 ) = eye;
    LHS.block<2, 2> ( 2, 4 ) = eye;
    LHS.block<2, 2> ( 4, 4 ) = eye;
    //
    RHS.segment<2> ( 0 ) = points1.col ( 0 );
    RHS.segment<2> ( 2 ) = points1.col ( 1 );
    RHS.segment<2> ( 4 ) = points1.col ( 2 );

    // Debugging
    //std::cout << LHS << std::endl;
    //std::cout << Num::determinant(LHS) << std::endl;

    unknowns = polyvec::Num::solve_linear_system ( LHS, RHS );

    A.col ( 0 ) = unknowns.segment<2> ( 0 );
    A.col ( 1 ) = unknowns.segment<2> ( 2 );
    b = unknowns.segment<2> ( 4 );

    // Make sure it is working
    const double error = ( ( A*points0 ).colwise() + b - points1 ).norm() ;
    assert_break ( error < 1e-10 );
}

namespace {
bool
skip_comments ( FILE* fl ) {
    char tmp[4096];
    bool is_there_more_to_read = false;
    assert_break ( fl );

    for ( ;; ) {
        const int c = fgetc ( fl );

        if ( c ==EOF ) {
            is_there_more_to_read = false;
            break;
        } else if ( c == '#' ) {
            fgets ( tmp, 4096, fl );
            continue;
        } else if ( ( c != ' ' ) && ( c != '\n' ) && ( c != '\r' ) ) {
            int success = fseek ( fl, -1, SEEK_CUR );
            assert ( success == 0 );
            is_there_more_to_read = true;
            break;
        }
    }

    return is_there_more_to_read;
}
}

void
dump_curve_sequence ( const std::vector<Eigen::Matrix2Xd>& curves, FILE* file ) {
    fprintf ( file, "# N curves \n %d \n", ( int ) curves.size() );

    for ( int i = 0 ; i < curves.size() ; ++ i ) {
        fprintf ( file, "# Curve %d, n points \n %d \n", i, ( int ) curves[i].cols() );

        for ( int j = 0 ; j < curves[i].cols() ; ++ j ) {
            fprintf ( file, "%.12g ", curves[i] ( 0, j ) );
        }

        fprintf ( file, "\n" );

        for ( int j = 0 ; j < curves[i].cols() ; ++ j ) {
            fprintf ( file, "%.12g ", curves[i] ( 1, j ) );
        }

        fprintf ( file, "\n" );
    }
}

std::vector<Eigen::Matrix2Xd>
read_curve_sequence ( FILE* file ) {
#define ENSURE(X) if(!(X)) {assert_break(0);}

    int n_curves;
    ENSURE ( skip_comments ( file ) );
    ENSURE ( fscanf ( file, "%d", &n_curves ) == 1 );

    std::vector<Eigen::Matrix2Xd> curves ( n_curves );

    for ( int i = 0 ; i < curves.size() ; ++ i ) {
        int n_points;
        ENSURE ( skip_comments ( file ) );
        ENSURE ( fscanf ( file, "%d", &n_points ) == 1 );
        curves[i].resize ( 2, n_points );

        ENSURE ( skip_comments ( file ) );

        for ( int j = 0 ; j < curves[i].cols() ; ++ j ) {
            ENSURE ( fscanf ( file, "%lf", &curves[i] ( 0, j ) ) == 1 );
        }

        ENSURE ( skip_comments ( file ) );

        for ( int j = 0 ; j < curves[i].cols() ; ++ j ) {
            ENSURE ( fscanf ( file, "%lf", &curves[i] ( 1, j ) ) == 1 );
        }
    }

    return curves;


#undef ENSURE
}


void merged_curve (
    const std::vector<Eigen::Matrix2d>& line_begin,
    const std::vector<Mat24>& beziers,
    const std::vector<Eigen::Matrix2d>& line_end,
    FILE* verbose,
    bool& doable,
    Mat24& ans_mat ) {

    PotraceParameters ans;


    auto eprintf = [verbose] ( const std::string strr ) {
        if ( verbose ) {
            fprintf ( verbose, "%s", strr.c_str() );
        }
    };

#define undoable() do{doable =false; eprintf(""); return;}while(0)

    assert_break ( line_begin.size() <= 1 );
    assert_break ( line_end.size() <= 1 );
    assert_break ( beziers.size() >= 1 );

    // ============
    // First check that each bezier should be potracable
    // ============
    for ( int i = 0 ; i < ( int ) beziers.size() ; ++i ) {
        if ( !PointParameters ( {beziers[i]} ).is_potracable() ) {
            eprintf ( StringUtils::fmt ( "curve[%d] is not potracable \n", i ) );
            undoable();
        }
    }

    // ============
    // Now check that transitions from each curve also does not
    // introduce an inflection
    // The assumption is that the tangents of all the transitions match
    // ============

    // line begin
    if ( line_begin.size() ) {
        if ( AngleUtils::have_opposite_convexity (
                    line_begin[0].col ( 0 ),
                    beziers[0].col ( 1 ),
                    beziers[0].col ( 2 ),
                    beziers[0].col ( 3 ) ) ) {
            eprintf ( StringUtils::fmt ( "line0 to bezier 1 inflected \n" ) );
            undoable();
        }
    }

    // line end
    if ( line_end.size() ) {
        if ( AngleUtils::have_opposite_convexity (
                    beziers.back().col ( 0 ),
                    beziers.back().col ( 1 ),
                    beziers.back().col ( 2 ),
                    line_end[0].col ( 1 ) ) ) {
            eprintf ( StringUtils::fmt ( "line_end to bezier end inflected \n" ) );
            undoable();
        }
    }

    // Beziers
    for ( int i = 0 ; i < ( int ) beziers.size()-1 ; ++i ) {
        if ( AngleUtils::have_opposite_convexity (
                    beziers[i].col ( 1 ),
                    beziers[i].col ( 2 ),
                    beziers[i+1].col ( 1 ),
                    beziers[i+1].col ( 2 ) ) ) {
            eprintf ( StringUtils::fmt ( "bezier[%d] to bezier[%d] inflected \n", i, i +1 ) );
            undoable();
        }
    }

    // ============
    //   Check that the first and last line intersection each other
    // ============
    Eigen::Vector2d y0, y3, oy;
    {
        Eigen::Matrix2d ray_begin;
        Eigen::Matrix2d ray_end;

        if ( line_begin.size() ) {
            ray_begin = line_begin.front();
        } else {
            ray_begin = beziers.front().leftCols<2>();
        }

        if ( line_end.size() ) {
            ray_end = line_end.front();
        } else {
            ray_end = beziers.back().rightCols<2>();
        }

        ray_end.rowwise().reverseInPlace();

        double lhs_at, rhs_at;
        const bool do_intersect = ::polyvec::LineUtils::intersect ( ray_begin.col ( 0 ), ray_begin.col ( 1 ), ray_end.col ( 0 ), ray_end.col ( 1 ), lhs_at, rhs_at );

        if ( ( !do_intersect )
                || ( lhs_at < 1. )
                || ( rhs_at < 1. ) ) {
            eprintf ( StringUtils::fmt ( "Intersection Problem \n" ) );
            undoable();
        }

        y0 = ray_begin.col ( 0 );
        y3 = ray_end.col ( 0 );
        oy = ::polyvec::LineUtils::line_at ( ray_begin.col ( 0 ), ray_begin.col ( 1 ), lhs_at );

        // This is intersecting on the other side of the universe
        if ( oy.norm() > 10000 ) {
            //polyvec::VtkCurveWriter writer;
            //writer.add_polyline( ray_begin );
            //writer.add_polyline( ray_end   );
            //writer.dump("DAMN.vtk");
            undoable();
        }
    }

    // ============
    //  Now find the area of the bezier sequence
    // ============
    double area_shape = 0;
    {
        auto find_carpet  = [&line_begin, &line_end, &beziers] {

            Eigen::Matrix2Xd carpet ( 2, beziers.size() + line_begin.size() + line_end.size() + 1 );
            int colid = 0;

            if ( line_begin.size() ) {
                carpet.col ( colid ) = line_begin.front().col ( 0 );
                ++colid;
            }

            for ( int i = 0 ; i < ( int ) beziers.size() ; ++i ) {
                carpet.col ( colid ) = beziers[i].col ( 0 );
                ++colid;
            }

            if ( line_end.size() ) {
                carpet.col ( colid ) = line_end.front().col ( 0 );
                ++colid;
                carpet.col ( colid ) = line_end.front().col ( 1 );
                ++colid;
            } else {
                carpet.col ( colid ) = beziers.back().col ( 3 );
                ++colid;
            }

            assert_break ( colid == ( int ) carpet.cols() );

            return carpet;
        }; // find carpet

        bool is_carpet_ccw;
        double area_carpet;
        Eigen::Matrix2Xd carpet = find_carpet ();
        ::polyvec::WindingNumber::compute_orientation ( carpet, is_carpet_ccw, area_carpet );

        double area_beziers = 0;

        for ( int i = 0 ; i < ( int ) beziers.size() ; ++i ) {
            area_beziers += PointParameters ( {beziers[i]} ).as_potrace_params().get_areay();
        }

        area_shape = area_carpet + area_beziers;
    }

    //
    // If we are here, we can finally create the merged bezier
    //
    {
        Mat23 p0, p1;
        //
        p0.col ( 0 ) << -1, 0;
        p0.col ( 1 ) << 0, 1;
        p0.col ( 2 ) << 1, 0;
        //
        p1.col ( 0 ) << y0;
        p1.col ( 1 ) << oy;
        p1.col ( 2 ) << y3;
        //
        find_affine_transformation ( p0, p1, ans.A, ans.b );

        const double area_ref  =  area_shape / std::abs ( polyvec::Num::determinant ( ans.A ) );

        assert_break( !std::isnan( area_ref ) );
        if ( ( area_ref > 1.19 ) ) {
            eprintf ( StringUtils::fmt ( "Area Problem \n" ) );
            undoable();
        }

        ans.alpha =   ans.beta  = 2. - sqrt ( 4-10*area_ref/3. ) ;
        assert_break (  ans.beta > 0 );

        if ( ans.beta > 1 ) {
            // DEBUGGING - TEMP, REMOVE
            //for (auto bz : beziers ) {
            //    std::cout << bz.row(0) << " " << bz.row(1) << std::endl;
            // }
            undoable();
        }

        assert_break ( area_ref > 0 );
        assert_break ( std::abs ( area_ref - ans.get_areaz() ) < 1e-4 );
        ans_mat = ans.get_yy().control_points;
        doable = true;
    }

#undef undoable
}

void
is_merge_acceptable (
    const std::vector<Mat24>& beziers,
    const Mat24& merged_bezier,
    const double tolerance,
    const double alphamax,
    double& penalty,
    bool& is_acceptable ) {

    //
    // HELPERS
    //

    auto get_distance = [] ( const PotraceParameters &bzp, ::polyvec::BezierCurve &bz, const Eigen::Vector2d& p0, const Eigen::Vector2d& p1,
    bool &proj_success, double &dist ) {
        double t;
        bzp.t_for_tangent ( p1-p0, t, proj_success );

        if ( !proj_success ) {
            return;
        }

        dist = ::polyvec::LineUtils::distance_from_point ( p0, p1, bz.pos ( t ) );
    };

    auto get_distance_signed = [] ( const PotraceParameters &bzp,  ::polyvec::BezierCurve &bz, const Eigen::Vector2d& p0, const Eigen::Vector2d& p1,
    const Eigen::Vector2d pforsign,  bool &proj_success, double &dist ) {
        double t;
        bzp.t_for_tangent ( p1-p0, t, proj_success );

        if ( !proj_success ) {
            return;
        }

        const double dist_pref =  ::polyvec::LineUtils::signed_distance_from_point ( p0, p1, pforsign );
        dist =  ::polyvec::LineUtils::signed_distance_from_point ( p0, p1, bz.pos ( t ) );

        if ( dist*dist_pref < 0 ) {
            dist = -std::abs ( dist );
        } else {
            dist = std::abs ( dist );
        }
    };

    auto check_distance_a = [tolerance] ( const double dist_a ) {
        assert_break ( dist_a >= 0 );

        if ( dist_a > tolerance ) {
            return false;
        }

        return true;
    };

    auto check_distance_b = [tolerance] ( const double dist_b ) {
        if ( dist_b < -tolerance ) {
            return false;
        }

        return true;
    };



    //
    // Loop over all beziers. Calculate the dist_a and dist_b distances
    // from the potrace paper.
    //
    penalty = 0;
    is_acceptable = true;
    ::polyvec::BezierCurve eval_ctx;
    eval_ctx.set_control_points ( merged_bezier );
    PotraceParameters project_ctx;
    project_ctx = PointParameters ( {merged_bezier} ).as_potrace_params();

    //
    // Check alpha
    //
    if ( project_ctx.alpha > alphamax ) {
        is_acceptable = false;
        return;
    }

    //
    // Check distace
    //
    for ( int i = 0 ; i < ( int ) beziers.size() ; ++i ) {

        // dist a
        const Eigen::Vector2d a = PointParameters ( {beziers[i]} ).as_potrace_params().get_oy();

        if ( i < ( int ) beziers.size() -1 ) {
            const Eigen::Vector2d anext = PointParameters ( {beziers[i+1]} ).as_potrace_params().get_oy();
            bool can_project;
            double dist_a;
            get_distance ( project_ctx, eval_ctx, a, anext, can_project, dist_a );

            if ( ( !can_project ) || ( !check_distance_a ( dist_a ) ) ) {
                is_acceptable = false;
                return;
            }

            penalty += dist_a * dist_a;
        } // end of dist a

        // dist b
        {
            const Eigen::Vector2d b = beziers[i].col ( 0 );
            const Eigen::Vector2d bnext =  beziers[i].col ( 3 );
            bool can_project;
            double dist_b;
            get_distance_signed ( project_ctx, eval_ctx, b, bnext, a, can_project, dist_b );

            if ( ( !can_project ) || ( !check_distance_b ( dist_b ) ) ) {
                is_acceptable = false;
                return;
            }

            penalty += dist_b * dist_b;
        } // end of dist b
    } // end of beziers

} //  merge_penalty()

void merge_penalty (
    const std::vector<Mat24>& beziers,
    const Mat24& merged_bezier,
    double& penalty ) {

    // Do a one sided projection of the points
    constexpr int n_samples = 10;

    penalty = 0;
    ::polyvec::BezierCurve big_eval_ctx;
    big_eval_ctx.set_control_points ( merged_bezier );

    for ( int i = 0 ; i < ( int ) beziers.size() ; ++i ) {
        ::polyvec::BezierCurve sub_eval_ctx;
        sub_eval_ctx.set_control_points ( beziers[i] );
        // This is a question . Should we multiply by length?
        const double weight = 1. / n_samples *  sub_eval_ctx.length();
        // Perhaps not
        // const double weight = 1. / n_samples;


        for ( int j = 1 ; j < ( n_samples+1 ) ; ++j ) {
            const double t_sub = 1. / ( n_samples+1 ) * j ;
            const Eigen::Vector2d point_sub = sub_eval_ctx.pos ( t_sub );
            const double t_big = big_eval_ctx.project ( point_sub );
            const Eigen::Vector2d point_big = big_eval_ctx.pos ( t_big );
            penalty += weight * ( point_big-point_sub ).norm();
        } // end of points
    }  // end of bezeirs
}

void
merge_curves (
    const std::vector<Eigen::Matrix2Xd>& curves_in, // lines of bezier control points
    const bool is_circular,
    const MergeRecursivelyOptions options,
	const std::set<Regularity::SymmetryPair>& symmetric_curves,
    std::vector<Eigen::Matrix2Xd>& curves_out,
    std::vector<std::vector<int>>& out2in  ) {

    const int n_curves = ( int ) curves_in.size();

    //
    // Check if a curve is right after a corner
    //
    auto is_line_after_corner = [&n_curves, &curves_in, &is_circular] ( const int i ) {
        if ( ( i == 0 ) && ( !is_circular ) ) {
            return false;
        } else {
            const int iprev = ( ( i-1 )+n_curves ) %n_curves;

            if ( ( curves_in[iprev].cols() == 2 )
                    && ( curves_in[i].cols() == 2 ) ) {
                return true;
            } else {
                return false;
            }
        }
    };

    //
    // Prepare the data for a merge curve
    //
    auto get_merge_data = [&n_curves, &curves_in, &is_circular]
                          ( const int i,
                            const int n_merge,
                            std::vector<Mat22>& line_begin,
                            std::vector<Mat24>& beziers,
                            std::vector<Mat22>& line_end,
    std::vector<int>& parent_ids ) {

        assert_break ( is_circular || ( n_merge+i<=n_curves ) );
        line_begin.resize ( 0 );
        line_end.resize ( 0 );
        beziers.resize ( 0 );
        parent_ids.resize ( 0 );

        for ( int offset = 0 ; offset < n_merge ; ++offset ) {
            const int curve_id = ( i + offset ) %n_curves;

            if ( ( curves_in[curve_id].cols() == 2 ) && ( offset==0 ) ) {
                line_begin.push_back ( curves_in[curve_id] );
            } else if ( ( curves_in[curve_id].cols() == 2 ) && ( offset==n_merge-1 ) ) {
                line_end.push_back ( curves_in[curve_id] );
            } 
            else {
                assert_break ( curves_in[curve_id].cols() == 4 );
                beziers.push_back ( curves_in[curve_id] );
            }

            parent_ids.push_back ( curve_id );
        }

    };

    //
    // Convert a curves end or begin point into node ids
    //
    auto n_graph_nodes = [&is_circular, &curves_in, &n_curves]() {
        if ( is_circular  ) {
            return n_curves;
        } else {
            return n_curves+1;
        }
    };
    auto curve_end_node = [&is_circular, &curves_in, &n_curves] ( const int i ) {
        if ( is_circular && ( i == n_curves-1 ) ) {
            return 0;
        } else {
            return i+1;
        }
    };
    auto curve_begin_node = [&is_circular, &curves_in, &n_curves] ( const int i ) {
        return i;
    };

    //
    // All the possible merge situations
    //
    std::vector<ShortestPath::Node> graph_nodes ( n_graph_nodes() );
    std::map<std::pair<int, int>, int> mc_map;
    std::vector<double> mc_penalty;
    std::vector<std::vector<int>> mc_parent_ids;
    std::vector<Eigen::Matrix2Xd> mc_pts;

    //
    // Create all the merge candidates
    //
    {
        std::vector<Mat22> line_begin;
        std::vector<Mat24> beziers;
        std::vector<Mat22> line_end;
        ShortestPath::Node node;
        std::vector<int> parent_ids;
        Mat24 merged_bezier;
        bool is_merge_doable;
        double potrace_penalty;
        double distace_penalty;		

        for ( int cbeginid = 0 ; cbeginid < n_curves ; ++cbeginid ) {
            const int max_n_merge = is_circular ? n_curves : n_curves - cbeginid;

			// find the next symmetry pair after the current curve
			auto next_symmetry = std::lower_bound(symmetric_curves.begin(), symmetric_curves.end(), Regularity::SymmetryPair(cbeginid, cbeginid, -1));
			if (next_symmetry == symmetric_curves.end())
				next_symmetry = symmetric_curves.begin(); //wrap around

			// We employ the following strategy to promote symmetries:
			// A merge is only valid if all its symmetries are completely within the
			// merge or none is (i.e., it is ok if only one curve of symmetry pairs
			// is within the merge). This is evaluated independently for each 
			// participating symmetric region

			struct SymmetryInfo
			{
				// stores the number of symmetric pairs that are completely
				// contained
				int n_complete_symmetries = 0;

				// stores the number of symmetric paris where only one part is
				// contained
				int n_partial_symmetries = 0;
			};

			std::map<int, SymmetryInfo> symmetry_info;

            for ( int n_merge = 1 ; n_merge <= max_n_merge ; ++n_merge ) {

				const int cendid = (cbeginid + n_merge - 1) % n_curves;
				
				// check if we have new symmetries
				while (next_symmetry != symmetric_curves.end() && next_symmetry->first == cendid)
				{
					auto& info = symmetry_info[next_symmetry->region];
					// check if the other curve is already in the merge
					if (PathUtils::contains_closed(curves_in.size(), cbeginid, cendid, next_symmetry->second, is_circular))
					{
						++info.n_complete_symmetries;
						--info.n_partial_symmetries;
					}
					else
						++info.n_partial_symmetries;

					++next_symmetry;
					if (next_symmetry == symmetric_curves.end())
						next_symmetry = symmetric_curves.begin(); //wrap around
				}

                // If n_merge is just 1 then we don't have to do much
                // The fit is accepted and the error is 0
                if ( n_merge==1 ) {
                    mc_pts.push_back ( curves_in[cbeginid] );
                    mc_penalty.push_back ( 0 );
                    mc_parent_ids.push_back ( {cbeginid} );
                } else {
					
                    // Check that we do not pass a corner                    
                    const int cendidm1 = ( cbeginid+n_merge-2 ) %n_curves;

                    // Handle lines
                    if ( options.allow_merging_lines ) {
                        // If were were fitting to potrace data this should be enough
                        if (  is_line_after_corner ( cendid ) ) {
                            break;
                        }
                        // For ed's assymetric case, we also need this
                        if ( ( curves_in[cendidm1].cols() == 2 ) && ( cendidm1 != cbeginid) ) {
                            break;
                        }
                    }

                    // New strategy, should work for Ed too
                    // Simply don't let cbeginid or cendid to become line
                    // basically we are not attempting to merge lines with anything
                    if ( !options.allow_merging_lines ) {
                        if ( ( curves_in[cbeginid].cols() == 2 ) || ( curves_in[cendid].cols() == 2 ) ) {
                            break;
                        }
                    }

					bool symmetry_valid = true;
					for(auto& info : symmetry_info)
						if (info.second.n_complete_symmetries > 0 && info.second.n_partial_symmetries > 0)
							symmetry_valid = false; // violates symmetry
					if (!symmetry_valid)
						continue;

                    // Get info needed for merging
                    get_merge_data ( cbeginid, n_merge, line_begin, beziers, line_end, parent_ids );

                    // Merge the curve
                    merged_curve ( line_begin, beziers, line_end, nullptr, is_merge_doable, merged_bezier );

                    if ( !is_merge_doable ) {
                        continue;
                    }

                    // Check if acceptable
                    is_merge_acceptable ( beziers, merged_bezier, options.distance_tolerance, options.alpha_max, potrace_penalty, is_merge_doable );

                    if ( !is_merge_doable ) {
                        continue;
                    }

                    // Find the penalty
                    merge_penalty ( beziers, merged_bezier, distace_penalty );

                    // Add to the curve to nodes
                    mc_pts.push_back ( merged_bezier );
                    mc_penalty.push_back ( distace_penalty );
                    mc_parent_ids.push_back ( parent_ids );

                } // end of more than one curve

            } // end of n_merge
        } // end of begin id
    } // end of creating the merge candidates


    //
    // Now create a graph
    //
    for ( int node_id = 0 ; node_id < n_graph_nodes() ; ++node_id  ) {
        graph_nodes[node_id].v = node_id;
    }

    for ( int mcid = 0 ; mcid < ( int ) mc_parent_ids.size() ; ++mcid ) {
        const int node_from = curve_begin_node ( mc_parent_ids[mcid].front() );
        const int node_to = curve_end_node ( mc_parent_ids[mcid].back() );
        graph_nodes[node_from].add_neighbor(node_to, mc_penalty[mcid] + options.per_bezier_const );
        mc_map.insert ( std::make_pair ( std::make_pair ( node_from, node_to ), mcid ) );
    }


    // Find the shortest path
    {
        std::vector<Eigen::Index> P;
        ShortestPath::State state;

        if ( is_circular ) {
            ShortestPath::find_cycle ( state, graph_nodes, {0}, P );
            // Debugging
            // ShortestPath::print_graph ( graph_nodes, stderr );
        } else {
            ShortestPath::find ( state, graph_nodes, 0, n_graph_nodes()-1, P );
        }

        const int max_node_offset = is_circular ? ( int ) P.size() : ( int ) P.size()-1;

        curves_out.resize ( 0 );
        out2in.resize ( 0 );

        for ( int nodeoffset = 0 ; nodeoffset <  max_node_offset ; ++nodeoffset ) {
            const int node_from = (int)P[nodeoffset];
            const int node_to = (int)P[ ( nodeoffset+1 ) %P.size()];
            auto mc_id_it = mc_map.find ( std::make_pair ( node_from, node_to ) );
            assert_break ( mc_id_it != mc_map.end() );
            const int mc_id = (int)mc_id_it->second;

            curves_out.push_back ( mc_pts[ mc_id ] );
            out2in.push_back ( mc_parent_ids[ mc_id ] );
        }

    }

} // merge_recursively()

NAMESPACE_END ( BezierMerging )
NAMESPACE_END ( polyfit )

