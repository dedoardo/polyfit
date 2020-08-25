#include <Eigen/Core>

#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/utils/potrace.hpp>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( PotraceUtil )

std::vector<Eigen::Matrix2Xd>
decode_potrace_curve ( potrace_curve_t* potrace_curve ) {
#if defined(POLYVEC_USE_POTRACE)

    int n = potrace_curve->n;
    int* tag = potrace_curve->tag;
    std::vector<Eigen::Matrix2Xd> ans;
    ans.reserve ( n );
    potrace_dpoint_t ( *c ) [3] = potrace_curve->c;
    Eigen::Matrix2Xd curve_points;

    for ( int i = 0; i < n; i++ ) {
        switch ( tag[i] ) {
        case POTRACE_CORNER: {
            curve_points.resize ( 2, 2 );
            {
                if ( i == 0 ) {
                    curve_points.col ( 0 ) << c[n - 1][2].x, c[n - 1][2].y;
                } else {
                    curve_points.col ( 0 ) << c[i - 1][2].x, c[i - 1][2].y;
                }

                curve_points.col ( 1 ) << c[i][1].x, c[i][1].y;
                ans.push_back ( curve_points );
            }
            {
                curve_points.col ( 0 ) << c[i][1].x, c[i][1].y;
                curve_points.col ( 1 ) << c[i][2].x, c[i][2].y;
                ans.push_back ( curve_points );
            }
            break;
        }

        case POTRACE_CURVETO: {
            curve_points.resize ( 2, 4 );

            if ( i == 0 ) {
                curve_points.col ( 0 ) << c[n - 1][2].x, c[n - 1][2].y;
            } else {
                curve_points.col ( 0 ) << c[i - 1][2].x, c[i - 1][2].y;
            }

            curve_points.col ( 1 ) << c[i][0].x, c[i][0].y;
            curve_points.col ( 2 ) << c[i][1].x, c[i][1].y;
            curve_points.col ( 3 ) << c[i][2].x, c[i][2].y;

            ans.push_back ( curve_points );
            break;
        }
        }
    }  // end of for over curves

    return ans;
#else
    assert_break ( 0 );
    return std::vector<Eigen::Matrix2Xd>();
#endif  //  defined(POLYVEC_USE_POTRACE)
}  // decode_potrace_curve


// Using the potrace function seemed to be harder than copy pasting the
// code for this case.
// copied from trace.c in potrace.
std::vector<Eigen::Matrix2Xd>
fit_polygon ( const Eigen::Matrix2Xd& polygon, const double alphamax ) {
    assert_break ( 0 && "I have not tested this yet" );
#if 0
#define sign(x) ((x) > 0 ? 1 : (x) < 0 ? -1 : 0)

    std::vector<Eigen::Matrix2Xd> curves_out ;

    auto mod = [] ( int i, int n ) {
        return ( i + n ) % n;
    };
    auto interval = [] ( double lambda, const Eigen::Vector2d& a,
    const Eigen::Vector2d& b ) {
        return a + lambda * ( b - a );
    };

    auto dorth_infty = [] ( const Eigen::Vector2d& p0, const Eigen::Vector2d& p2 ) {
        Eigen::Vector2d r;

        r.y() = sign ( p2.x() - p0.x() );
        r.x() = -sign ( p2.y() - p0.y() );

        return r;
    };

    auto ddenom = [dorth_infty] ( const Eigen::Vector2d& p0,
    const Eigen::Vector2d& p2 ) {
        Eigen::Vector2d r = dorth_infty ( p0, p2 );

        return r.y() * ( p2.x() - p0.x() ) - r.x() * ( p2.y() - p0.y() );
    };

    auto dpara = [] ( const Eigen::Vector2d& p0, const Eigen::Vector2d& p1,
    const Eigen::Vector2d& p2 ) {
        double x1, y1, x2, y2;

        x1 = p1.x() - p0.x();
        y1 = p1.y() - p0.y();
        x2 = p2.x() - p0.x();
        y2 = p2.y() - p0.y();

        return x1 * y2 - x2 * y1;
    };

    int m = ( int ) polygon.cols();

    int i, j, k;
    double dd, denom, alpha;
    Eigen::Vector2d p2, p3, p4;
    int make_corner, make_smooth;

    /* examine each vertex and find its best fit */
    for ( i = 0; i < m; i++ ) {
        j = mod ( i + 1, m );
        k = mod ( i + 2, m );
        p4 = interval ( 1 / 2.0, polygon.col ( k ), polygon.col ( j ) );

        denom = ddenom ( polygon.col ( i ), polygon.col ( k ) );

        if ( denom != 0.0 ) {
            dd = dpara ( polygon.col ( i ), polygon.col ( j ), polygon.col ( k ) ) / denom;
            dd = fabs ( dd );
            alpha = dd > 1 ? ( 1 - 1.0 / dd ) : 0;
            alpha = alpha / 0.75;
        } else {
            alpha = 4 / 3.0;
        }

        // curve->alpha0[j] = alpha; /* remember "original" value of alpha */

        make_corner = ( alpha >= alphamax );
        make_smooth = ( alpha < alphamax );

        // Shayan: this is the part that I have not yet figured out how to do.
        if ( make_corner ) { /* pointed corner */
            // curve->tag[j] = POTRACE_CORNER;
            // curve->c[j][1] = curve->vertex[j];
            // curve->c[j][2] = p4;
            curves_out.push_back ( Eigen::Matrix2Xd ( 2, 2 ) ) ;
            curves_out.push_back ( Eigen::Matrix2Xd ( 2, 2 ) ) ;
            curves_out.back().push_back (
        } else if ( make_smooth ) {
            if ( alpha < 0.55 ) {
                alpha = 0.55;
            } else if ( alpha > 1 ) {
                alpha = 1;
            }

            p2 = interval ( .5 + .5 * alpha, polygon.col ( i ), polygon.col ( j ) );
            p3 = interval ( .5 + .5 * alpha, polygon.col ( k ), polygon.col ( j ) );
            curve->tag[j] = POTRACE_CURVETO;
            curve->c[j][0] = p2;
            curve->c[j][1] = p3;
            curve->c[j][2] = p4;
        } else {
            assert ( 0 );
        }

        // curve->alpha[j] = alpha; /* store the "cropped" value of alpha */
        // curve->beta[j] = 0.5;
    }

    return;

#undef sign
#endif
}

std::vector<Eigen::Matrix2Xd>
run_pipeline ( const Eigen::Matrix2Xd& path,
               bool should_adjust_path,
               bool optimize_curve ) {
#if defined(POLYVEC_USE_POTRACE)

    //
    // Trace with potrace
    //

    // create the parameters
    param_unique_ptr potrace_params;
    potrace_params.reset ( potrace_param_default() );
    potrace_params->alphamax = 1.0;
    potrace_params->opticurve = optimize_curve;
    potrace_params->opttolerance = 1.0;

    // create the path
    path_unique_ptr potrace_path;
    potrace_path.reset ( potrace_path_create() );

    //
    // set the pixels
    //
    Eigen::Matrix2Xi path_int = path.cast<int>();
    potrace_path_manual ( potrace_path.get(), ( int ) path_int.cols(),
                          path_int.data() );

    //
    // Fit the polygon
    //
    potrace_ipolygon_automatic ( potrace_path.get(), potrace_params.get() );

    //
    // adjust the path
    //
    if ( should_adjust_path ) {
        potrace_dpolygon_automatic ( potrace_path.get(), potrace_params.get() );
    } else {
        potrace_dpolygon_naive ( potrace_path.get() );
    }

    // Now smooth out
    Eigen::VectorXi path_transitions ( path.cols() );
    path_transitions.setConstant ( POTRACE_TRANSITION_GUESS );
    potrace_smooth ( potrace_path.get(), potrace_params.get(),
                     ( int ) path_transitions.size(), path_transitions.data() );

    // Now get the curves
    return decode_potrace_curve ( &potrace_path->curve );

#else
    assert_break ( 0 );
    return std::vector<Eigen::Matrix2Xd>();
#endif

}  // end of run_pipeline()

NAMESPACE_END ( PotraceUtil )
NAMESPACE_END ( polyfit )
