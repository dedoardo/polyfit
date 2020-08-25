namespace polyvec {
    namespace geom {
        //
        // Axis-aligned bounding box
        //
        polyvec_inline aabb::aabb ( double tl_x, double tl_y, double w, double h ) {
            min = real2 ( tl_x, tl_y );
            max = min + real2 ( w, h );
        }

        polyvec_inline aabb::aabb ( const real2& size ) {
            min = real2::Zero();
            max = size;
        }

        polyvec_inline aabb::aabb ( const real2& topleft, const real2& size ) {
            min = topleft;
            max = topleft + size;
        }

        polyvec_inline aabb::aabb ( const real2* pts, index n_pts ) {
            if ( pts == nullptr ) {
                return;
            }

            min = real2_hi;
            max = real2_lo;

            for ( index i = 0; i < n_pts; ++i ) {
                real2 pt = pts[i];
                min.x() = ::std::min ( min.x(), pt.x() );
                min.y() = ::std::min ( min.y(), pt.y() );
                max.x() = ::std::max ( max.x(), pt.x() );
                max.y() = ::std::max ( max.y(), pt.y() );
            }
        }

        polyvec_inline aabb& aabb::add ( const real2& p ) {
            min.x() = ::std::min ( min.x(), p.x() );
            min.y() = ::std::min ( min.y(), p.y() );
            max.x() = ::std::max ( max.x(), p.x() );
            max.y() = ::std::max ( max.y(), p.y() );
            return *this;
        }

        polyvec_inline double aabb::width() const {
            return max.x() - min.x();
        }

        polyvec_inline double aabb::height() const {
            return max.y() - min.y();
        }

        polyvec_inline double aabb::resolution() const {
            return ::std::max ( width(), height() );
        }

        polyvec_inline double aabb::left() const {
            return min.x();
        }

        polyvec_inline double aabb::right() const {
            return max.x();
        }

        polyvec_inline double aabb::top() const {
            return min.y();
        }

        polyvec_inline double aabb::bottom() const {
            return max.y();
        }

        polyvec_inline real2 aabb::tl() const {
            return min;
        }

        polyvec_inline real2 aabb::tr() const {
            return real2 ( max.x(), min.y() );
        }

        polyvec_inline real2 aabb::br() const {
            return max;
        }

        polyvec_inline real2 aabb::bl() const {
            return real2 ( min.x(), max.y() );
        }

        polyvec_inline real2 aabb::dims() const {
            return real2 ( width(), height() );
        }

        polyvec_inline bool aabb::valid() const {
            return width() > 0 && height() > 0;
        }

        //
        // Segment
        //
        polyvec_inline segment::segment ( const real2& src, const real2& dst ) :
            src ( src ),
            dst ( dst ) {
        }

        polyvec_inline segment& segment::scale ( double s ) {
            src += s * ( src - dst );
            dst += s * ( dst - src ); // i know i know.. doesn't matter
            return *this;
        }

        polyvec_inline double segment::length() const {
            return ( dst - src ).norm();
        }

        polyvec_inline double segment::orientation() const {
            real2 dir = dst - src;
            return atan2 ( dir.y(), dir.x() );
        }

        polyvec_inline real2 segment::pos ( double t ) const {
            assert_break ( t >= 0 && t <= 1. );
            return src + ( dst - src ) * t;
        }

        polyvec_inline real2 segment::dd() const {
            return ( dst - src );
        }

        polyvec_inline real2 segment::dir() const {
            return dd().normalized();
        }

        polyvec_inline real2 segment::midpoint() const {
            return .5 * ( src + dst );
        }

        polyvec_inline aabb segment::bb() const {
            return aabb().add ( src ).add ( dst );
        }

        //
        // Intersection/Evaluation
        //
        // Returns true if the edge connecting the two points lies inside the polygon
        // the points in `polygon` are assumed to be circular
        polyvec_inline bool is_edge_inside_polygon ( const Eigen::Matrix2Xd& polygon, const real2& edge_src, const real2& edge_dst ) {
            using namespace Eigen;

            // moving the points a little inside the shape
            const real2 src = edge_src + ( edge_dst - edge_src ) * constants::Eps;
            const real2 dst = edge_dst + ( edge_src - edge_dst ) * constants::Eps;
            (void)(src);
            (void)(dst);

            // Checking if the line intersects any of the edges
            for ( Eigen::Index i = 0; i < polygon.cols(); ++i ) {
                const real2 polygon_edge_src = polygon.col ( i );
                const real2 polygon_edge_dst = polygon.col ( i + 1 );

                double polygon_iset_t, edge_iset_t;

                if ( !LineUtils::intersect ( polygon_edge_src, polygon_edge_dst, edge_src, edge_dst, polygon_iset_t, edge_iset_t ) ) {

                }
            }

            break_here;
            return true;
        }

        //
        // Angle routines
        //
        polyvec_inline double radians ( double degrees ) {
            return ( degrees * constants::PI ) / 180.;
        }
        polyvec_inline double degrees ( double radians ) {
            return ( radians * 180 ) / constants::PI;
        }

        polyvec_inline double angle ( const real2& v ) {
            return atan2 ( v[1], v[0] );
        }
        polyvec_inline double angle ( const real2& v1, const real2& v2 ) {
            return atan2 ( v1[0] * v2[1] - v1[1] * v2[0], v1.dot ( v2 ) );
        }

        //These functions bring an angle into the range [0,2Pi] or [rangeStart, rangeStart+2Pi]
        //They assume we're not too far on the negative side of the range
        polyvec_inline double to_range ( double angle ) {
            return fmod ( angle + 8 * constants::PI, constants::PI2 );
        }
        polyvec_inline double to_range ( double angle, double rangeStart ) {
            return fmod ( angle + 16 * constants::PI - rangeStart, constants::PI2 ) + rangeStart;
        }

        polyvec_inline double axis_deviation ( const real2& d ) {
            return std::fmod ( ::std::atan2 ( d.y(), d.x() ) + constants::PI2, constants::PI_half );
        }

		// clockwise angle from (1, 0)
        polyvec_inline double orientation_2PI (const real2& d) {
            double angle = atan2(d.y(), d.x());
            return d.y() < 0 ? (2. * M_PI - std::abs(angle)) : angle;
        }

		polyvec_inline Eigen::Vector2d direction_from_angle(const double radians) {
			return Eigen::Vector2d(cos(radians), sin(radians));
		}

		polyvec_inline double smallest_angle_between_points(const Eigen::Vector2d& pt0, const Eigen::Vector2d& pt1, const Eigen::Vector2d& pt2) {
			const Eigen::Vector2d d0 = (pt0 - pt1).normalized();
			const Eigen::Vector2d d1 = (pt2 - pt1).normalized();
			return acos(d0.dot(d1));
		}

        polyvec_inline double angle_between_vectors ( const real2& lhs, const real2& rhs ) {
            return acos ( std::max ( -1.0, std::min ( lhs.dot ( rhs ), 1.0 ) ) );
        }

        polyvec_inline double closest_horizontal_angle ( double angle ) {
            angle = to_range ( angle );

            if ( ( angle < constants::PI_half ) || ( angle > 3 * constants::PI_half ) ) {
                return 0;
            } else {
                return constants::PI;
            }
        }

        polyvec_inline double closest_vertical_angle ( double angle ) {
            angle = to_range ( angle );

            if ( angle < constants::PI ) {
                return constants::PI_half;
            } else {
                return 3 * constants::PI_half;
            }
        }
    }

//
// CurveParams
//
    inline const CurveParams::ParamVec& CurveParams::params() const {
        return _curve;
    }

    inline void CurveParams::set_params ( const ParamVec& curve ) {
        _curve = curve;
        on_params_changed();
    }

//
// Curve
//
    inline Line& Curve::line() {
        return * ( ( Line* ) this );
    }

    inline Arc& Curve::arc() {
        return * ( ( Arc* ) this );
    }

    inline Clothoid& Curve::clothoid() {
        return * ( ( Clothoid* ) this );
    }

    inline const Line& Curve::line() const {
        return * ( ( Line* ) this );
    }

    inline const Arc& Curve::arc() const {
        return * ( ( Arc* ) this );
    }

    inline const Clothoid& Curve::clothoid() const {
        return * ( ( Clothoid* ) this );
    }


//
// Line
//
// Constructors
    inline Line::Line ( const real2& src, const real2& dst ) {
        _curve.resize ( n_params() );
        _curve.head<2>() = src;
        real2 dir = dst - src;
        double len = dir.norm();

        if ( fabs ( len ) < Eigen::NumTraits<double>::dummy_precision() ) { //fail gracefully on zero-length line
            _curve[LENGTH] = 0;
            _curve[ANGLE] = 0;
            _der = real2 ( 1., 0. );
        } else {
            _curve[LENGTH] = len;
            _curve[ANGLE] = geom::angle ( dir );
            _der = dir * ( 1. / len );
        }
    }

// Interfaces
    inline int Line::n_params() const {
        return 4;
    }

    inline void Line::on_params_changed() {
        _der = real2 ( cos ( _src_angle() ), sin ( _src_angle() ) );
    }

    inline const char* Line::type_str() const {
        return "line";
    }

    inline int Line::type() const {
        return LINE;
    }

    inline Curve* Line::clone() const {
        Line* copy = new Line;
        copy->set_params ( params() );
        return copy;
    }

    inline real2 Line::src() const {
        return _src();
    }

    inline real2 Line::dst() const {
        return pos ( length() );
    }

    inline real2 Line::pos ( double s ) const {
        return src() + s * _der;
    }

    inline real2 Line::der ( double s ) const {
        param_unused ( s );
        return _der;
    }

    inline real2 Line::der2 ( double s ) const {
        param_unused ( s );
        return real2::Zero();
    }

    inline double Line::project ( const real2& point ) const {
        return std::min ( _length(), std::max ( 0., _der.dot ( point - _src() ) ) );
    }

    inline bool Line::trace_ray ( const real2& rayorig, const real2& raydir, double& sray, double& scurve, double from,
                                  double to ) const {
        const double tol = 1e-6;

        const Eigen::Vector2d normal = Eigen::Vector2d ( _der ( 1 ), -_der ( 0 ) );
        const double denom = normal.dot ( raydir );

        if ( std::abs ( denom ) < tol ) {
            return false;
        }

        sray = ( _src() - rayorig ).dot ( normal ) / denom;

        const Eigen::Vector2d point = rayorig + sray * raydir;
        scurve = _der.dot ( point - _src() );

        if ( ( scurve < to + tol ) && ( scurve > from - tol ) ) {
            return true;
        } else {
            return false;
        }
    }

    inline double Line::length() const {
        return _length();
    }

    inline double Line::angle ( double s ) const {
        param_unused ( s );
        return _src_angle();
    }

    inline double Line::curvature ( double s ) const {
        param_unused ( s );
        return 0;
    }

    inline double Line::src_angle() const {
        return _src_angle();
    }

    inline double Line::src_curvature() const {
        return 0.;
    }

    inline double Line::dst_angle() const {
        return _src_angle();
    }

    inline double Line::dst_curvature() const {
        return 0.;
    }

    inline void Line::derivative_at ( double s, ParamDer2& out, ParamDer2& out_tan ) const {
        out_tan = out = ParamDer2::Zero ( 2, 4 );
        out ( 0, X ) = 1;
        out ( 1, Y ) = 1;
        out ( 0, ANGLE ) = -s * _der ( 1 );
        out ( 1, ANGLE ) = s * _der ( 0 );
        out_tan ( 0, ANGLE ) = -_der ( 1 );
        out_tan ( 1, ANGLE ) = _der ( 0 );
    }

    inline void Line::derivative_at_end ( int continuity, EndDer& out ) const {
        out = EndDer::Zero ( 2 + continuity, 4 );
        out ( 0, X ) = 1;
        out ( 1, Y ) = 1;
        out ( 0, LENGTH ) = _der ( 0 );
        out ( 1, LENGTH ) = _der ( 1 );
        out ( 0, ANGLE ) = -_length() * _der ( 1 );
        out ( 1, ANGLE ) = _length() * _der ( 0 );

        if ( continuity >= 1 ) {
            out ( 2, ANGLE ) = 1.;
        }
    }

    inline real3 Line::implicit_equation_coeffs() const {
        Eigen::Vector2d normal ( -_der ( 1 ), _der ( 0 ) );
        Eigen::Vector2d x0 ( _curve[X], _curve[Y] );

        Eigen::Vector3d abc;

        abc ( 0 ) = normal ( 0 );
        abc ( 1 ) = normal ( 1 );
        abc ( 2 ) = -normal.dot ( x0 );

        return abc;
    }

    inline Curve::ParamDer3 Line::implicit_equation_coeffs_derivative() const {
        ParamDer3 dabc_dparams ( 3, n_params() );
        dabc_dparams.setZero();

        // a = -sin angle
        // b = cos angle
        // c = sin * x - cos * y

        const double sin_angle = sin ( _curve[ANGLE] );
        const double cos_angle = cos ( _curve[ANGLE] );
        const double x = _curve[X];
        const double y = _curve[Y];

        //
        dabc_dparams ( 0, ANGLE ) = -cos_angle;

        //
        dabc_dparams ( 1, ANGLE ) = -sin_angle;

        //
        dabc_dparams ( 2, X ) = sin_angle;
        dabc_dparams ( 2, Y ) = -cos_angle;
        dabc_dparams ( 2, ANGLE ) = x * cos_angle + y * sin_angle;

        return dabc_dparams;
    }

//
// Arc
//
    inline Arc::Arc ( const real2& src, double angle, double length, double curvature ) {
        _curve.resize ( n_params() );
        _curve.head<2>() = src;
        _curve[ANGLE] = geom::to_range ( angle );
        _curve[LENGTH] = length;
        _curve[CURVATURE] = curvature;

        on_params_changed();
    }

    inline Arc::Arc ( const real2& start, const real2& mid, const real2& end, bool& success ) {
        success = false;
        _curve.resize ( n_params() );

        for ( int i = 0; i < n_params(); ++i ) {
            _curve[i] = 0;
        }

        real2 mid1 = mid - start;
        real2 end1 = end - start;

        double twiceSignedArea = ( mid1[0] * end1[1] - mid1[1] * end1[0] );
        double abc = sqrt ( mid1.squaredNorm() * end1.squaredNorm() * ( mid - end ).squaredNorm() );

        // There seems to be a bug here.
        // When this is combined with the clothoid projector
        // Maybe use _flat instead.
        if ( fabs ( twiceSignedArea ) < 1e-16 || abc < 1e-16 ) {
            return;    //degenerate arc
        }

        _curve.head<2>() = start;
        _curve[CURVATURE] = 2. * twiceSignedArea / abc;
        double halfArcAngle = asin ( fabs ( 0.5 * end1.norm() * _curve[CURVATURE] ) );

        if ( mid1.dot ( end - mid ) < 0. ) {
            halfArcAngle = constants::PI - halfArcAngle;
        }

        _curve[LENGTH] = fabs ( 2. * halfArcAngle / _curve[CURVATURE] );

        if ( twiceSignedArea < 0. ) {
            _curve[ANGLE] = geom::angle ( end1 ) + halfArcAngle;
        } else {
            _curve[ANGLE] = geom::angle ( end1 ) - halfArcAngle;
        }

        success = true;
        on_params_changed();
    }

    inline int Arc::n_params() const {
        return 5;
    }

    inline void Arc::on_params_changed() {
        _tangent = real2 ( cos ( _src_angle() ), sin ( _src_angle() ) );
        _angle_diff = _length() * _curvature();
        _flat = fabs ( _curvature() ) < 1e-3;

        if ( !_flat ) {
            real2 to_center ( -_tangent.y(), _tangent.x() );
            _radius = 1. / _curvature();
            _center = _src() + _radius * to_center;
        } else {
            //   dbg::warning ( FMT ( "Flat arc" ) );
        }
    }

    inline const char* Arc::type_str() const {
        return "arc";
    }

    inline int Arc::type() const {
        return ARC;
    }

    inline Curve* Arc::clone() const {
        Arc* copy = new Arc;
        copy->set_params ( params() );
        return copy;
    }

    inline real2 Arc::src() const {
        return _src();
    }

    inline real2 Arc::dst() const {
        return pos ( length() );
    }

    inline real2 Arc::pos ( double s ) const {
        const double angle = _src_angle() + s * _curvature();

        if ( _flat ) {
            return _src() + _tangent * s;
        } else {
            return _center + _radius * real2 ( sin ( angle ), -cos ( angle ) );
        }
    }

    inline real2 Arc::der ( double s ) const {
        const double angle = _src_angle() + s * _curvature();
        return real2 ( cos ( angle ), sin ( angle ) );
    }

    inline real2 Arc::der2 ( double s ) const {
        const double angle = _src_angle() + s * _curvature();
        return real2 ( -sin ( angle ), cos ( angle ) ) * _curvature();
    }

    inline double Arc::project ( const real2& point ) const {
        double t;

        if ( _flat ) {
            t = ( point - _src() ).dot ( _tangent );
        } else {
            double projAngle = atan2 ( point[1] - _center[1], point[0] - _center[0] );

            //To compute the projection, get the angle difference into the range centered on the midpoint of
            //the arc.
            if ( _curve[CURVATURE] > 0. ) {
                double projAngleDiff = constants::PI_half + projAngle - _src_angle();
                t = geom::to_range ( projAngleDiff, 0.5 * _angle_diff - constants::PI ) * _radius;
            } else {
                double projAngleDiff = constants::PI_half - projAngle + _src_angle();
                t = -geom::to_range ( projAngleDiff, -0.5 * _angle_diff - constants::PI ) * _radius;
            }
        }


        return std::max ( 0., std::min ( _length(), t ) );
    }

    inline bool Arc::trace_ray ( const real2& rayorig, const real2& raydir, double& sray, double& scurve, double from,
                                 double to ) const {
        if ( _flat ) {
            return false;
        }

        const double tol = 1e-8;
        sray = 1e10, scurve = -1e10;
        bool success = false;
        double alpha[2];

        const Eigen::Vector2d omc = rayorig - _center;
        const double aa = raydir.squaredNorm();
        const double bb = raydir.dot ( omc );
        const double cc = omc.squaredNorm() - _radius * _radius;
        const double delta = bb * bb - aa * cc;

        if ( delta < 0 ) {
            return false;
        } else if ( delta < tol ) {
            alpha[0] = alpha[1] = -bb / aa;
        } else {

            const double sqrt_delta = sqrt ( delta );
            alpha[0] = ( -bb + sqrt_delta ) / aa;
            alpha[1] = ( -bb - sqrt_delta ) / aa;
        }


        for ( int i = 0; i < 2; i++ ) {
            const Eigen::Vector2d point = rayorig + alpha[i] * raydir;
            double candid_scurve;

            const double projAngle = atan2 ( point[1] - _center[1], point[0] - _center[0] );

            if ( _curve[CURVATURE] > 0. ) {
                const double projAngleDiff = constants::PI_half + projAngle - _src_angle();
                candid_scurve = geom::to_range ( projAngleDiff, 0.5 * _angle_diff - constants::PI ) * _radius;
            } else {
                const double projAngleDiff = constants::PI_half - projAngle + _src_angle();
                candid_scurve = -geom::to_range ( projAngleDiff, -0.5 * _angle_diff - constants::PI ) * _radius;
            }

            if ( ( candid_scurve < to + tol ) &&
                    ( candid_scurve > from - tol ) &&
                    ( std::abs ( alpha[i] ) < std::abs ( sray ) ) ) {
                sray = alpha[i];
                scurve = candid_scurve;
                success = true;
            }
        }

        return success;
    }

    inline double Arc::length() const {
        return _length();
    }

    inline real2 Arc::center() const {
        return _center;
    }

    inline double Arc::radius() const {
        return _radius;
    }

    inline double Arc::radius_abs() const {
        return std::abs ( _radius );
    }

    inline double Arc::angle ( double s ) const {
        return _src_angle() + s * _curvature();
    }

    inline double Arc::curvature ( double s ) const {
        param_unused ( s );
        return _curvature();
    }

    inline double Arc::src_angle() const {
        return _src_angle();
    }

    inline double Arc::src_curvature() const {
        return _curvature();
    }

    inline double Arc::dst_angle() const {
        return _src_angle() + _angle_diff;
    }

    inline double Arc::dst_curvature() const {
        return _curvature();
    }

    inline void Arc::derivative_center ( ParamDer2& out ) const {
        const double eps = 1e-6;
        assert_break ( !_flat );

        Arc* aa = ( Arc* ) this->clone();
        out.resize ( 2, n_params() );

        ParamVec mod = _curve;
        real2 plus, minus;

        for ( int i = 0; i < n_params(); i++ ) {
            mod = _curve;
            mod ( i ) += eps;
            aa->set_params ( mod );
            plus = aa->center();

            mod ( i ) -= 2.*eps;
            aa->set_params ( mod );
            minus = aa->center();

            out.col ( i ) = ( plus - minus ) / 2. / eps;
        }

        delete aa;
    }

    inline void Arc::derivative_radius ( ParamDer1& out ) const {
        const double eps = 1e-6;
        assert_break ( !_flat );

        Arc* aa = ( Arc* ) this->clone();
        out.resize ( 1, n_params() );

        ParamVec mod = _curve;
        double plus, minus;

        for ( int i = 0; i < n_params(); i++ ) {
            mod = _curve;
            mod ( i ) += eps;
            aa->set_params ( mod );
            plus = aa->radius();

            mod ( i ) -= 2.*eps;
            aa->set_params ( mod );
            minus = aa->radius();

            out ( 0, i ) = ( plus - minus ) / 2. / eps;
        }

        delete aa;
    }

    inline void Arc::derivative_at ( double s, ParamDer2& out, ParamDer2& out_tan ) const {
        out_tan = out = ParamDer2::Zero ( 2, 5 );
        out ( 0, X ) = 1;
        out ( 1, Y ) = 1;

        real2 diff = pos ( s ) - _src();
        out ( 0, ANGLE ) = -diff[1];
        out ( 1, ANGLE ) = diff[0];

        if ( _flat ) {
            out_tan.col ( ANGLE ) = real2 ( -_tangent[1], _tangent[0] );
            out.col ( CURVATURE ) = ( 0.5 * s * s ) * real2 ( -_tangent[1], _tangent[0] );
            out_tan.col ( CURVATURE ) = s * real2 ( -_tangent[1], _tangent[0] );
        } else {
            double angle = _src_angle() + s * _curve[CURVATURE];

            double cosa, sina;
            cosa = cos ( angle );
            sina = sin ( angle );

            out_tan.col ( ANGLE ) = real2 ( -sina, cosa );
            out ( 0, CURVATURE ) = ( s * cosa + ( _tangent[1] - sina ) * _radius ) * _radius;
            out ( 1, CURVATURE ) = ( s * sina - ( _tangent[0] - cosa ) * _radius ) * _radius;
            out_tan ( 0, CURVATURE ) = -s * sina;
            out_tan ( 1, CURVATURE ) = s * cosa;
        }
    }

    inline void Arc::derivative_at_end ( int continuity, EndDer& out ) const {
        out = EndDer::Zero ( 2 + continuity, 5 );

        ParamDer2 pDer, dummy;
        derivative_at ( _length(), pDer, dummy );
        real2 endTan = der ( _length() );

        // Set \partial (x_end,y_end) / \partial (params)
        out.block ( 0, 0, 2, 5 ) = pDer;


        // Set \partial (ALL) / \partial (LENGTH)
        out.col ( LENGTH ).head<2>() = endTan;

        // Set \partial (angle_end) / \partial (params)
        if ( continuity >= 1 ) {
            out ( 2, ANGLE ) = 1.;
            out ( 2, LENGTH ) = _curve[CURVATURE];
            out ( 2, CURVATURE ) = _curve[LENGTH];
        }

        // Set \partial (\kappa_end) / \partial (params)
        if ( continuity == 2 ) {
            out ( 3, CURVATURE ) = 1.;
        }
    }
}