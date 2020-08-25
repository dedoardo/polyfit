#ifndef polyvec_geom_h_
#define polyvec_geom_h_

// C++ STL
#include <deque>

// polyvec
#include <polyvec/api.hpp>
#include <polyvec/core/macros.hpp>
#include <polyvec/geometry/line.hpp>

namespace polyvec {

    namespace geom {
        // ---
        // Geometric primitive axis-aligned bounding box
        struct aabb {
            aabb ( const double tl_x, const double tl_y, const double w, const double h );
            aabb ( const real2& size );
            aabb ( const real2& topleft, const real2& size );
            aabb ( const real2* pts = nullptr, index n_pts = 0 );

            aabb& add ( const real2& p );

            double width() const;
            double height() const;
            double resolution() const;

            double left() const;
            double right() const;
            double top() const;
            double bottom() const;

            real2 tl() const;
            real2 tr() const;
            real2 br() const;
            real2 bl() const;
            real2 dims() const;

            bool valid() const;

            real2 min = real2::Constant(std::numeric_limits<double>::infinity());
            real2 max = real2::Constant(-std::numeric_limits<double>::infinity());
        };

        // ---
        // Geometric primitive segment between two points
        struct segment {
            segment ( const real2& from = real2_0, const real2& to = real2_0 );

            real2 src = real2_0;
            real2 dst = real2_0;

            segment& scale ( double s );
            double   length() const;
            double   orientation() const;
            real2    pos ( double t ) const ;
            real2    dd() const;
            aabb     bb() const;
            real2    dir() const;
            real2    midpoint() const;
        };


        // ---
        // Misc Intersection/Evaluation
        polyvec_inline bool is_edge_inside_polygon ( const Eigen::Matrix2Xd& polygon,
                const real2& point_src, const real2& point_dst );

        // moved to line utils.hpp
        // polyvec_inline double project_t ( const real2& point, const segment& line );
        // polyvec_inline real2 line_at ( const real2& src, const real2& dst, const double t );
        // polyvec_inline bool intersect ( const real2& lhs_src, const real2& lhs_dst,
        //                                 const real2& rhs_src, const real2& rhs_dst,
        //                                double& lhs_at, double& rhs_at );

        // shayan: moved to misc.hpp
        // polyvec_inline Eigen::MatrixXd reshaped ( const Eigen::MatrixXd& in, int m, int n );

        // sign of cross
        polyvec_inline int side_of_line(const Eigen::Vector2d& src, const Eigen::Vector2d& dst, const Eigen::Vector2d& pt) {
            const double x = (pt.x() - src.x()) * (dst.y() - src.y()) - (pt.y() - src.y()) * (dst.x() - src.x());
            return abs(x) < constants::Eps ? 0 : ((0 < x) - (x < 0));
        }

        // discards dot sign
        polyvec_inline bool lines_are_aligned(const Eigen::Vector2d& d0, const Eigen::Vector2d& d1) {
            return (1. - abs(d0.dot(d1))) < constants::Eps;
        }        

        // do they share one of the two coordinates?
        polyvec_inline bool points_are_axis_aligned(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) {
            return abs((p0.x() - p1.x())) < constants::Eps || 
                   abs((p0.y() - p1.y())) < constants::Eps;
        }
        
        polyvec_inline double angle_deviation_horz (const double rad) {
            return std::min(abs(rad), M_PI - abs(rad));
        }

        polyvec_inline double points_match(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1) {
            return (p0 - p1).norm() < 1e-1;
        }

            //polyvec_inline Eigen::VectorXd projection_derivatives(const int n_params,
            //                                                      const double time,
            //                                                      const double tend,
            //                                                      const Eigen::Vector2d &outside_point,
            //                                                      const Eigen::Vector2d &rr,
            //                                                      const Eigen::Vector2d &drdt,
            //                                                      const Eigen::Vector2d &drdtdt,
            //                                                      const Eigen::Matrix2Xd &drdparams,
            //                                                      const Eigen::Matrix2Xd &drdtdparams);


        // ---
        // _Highly optimized_ angle utilities
        polyvec_inline double radians ( double degrees );
        polyvec_inline double degrees ( double radians );
        polyvec_inline double angle ( const real2& v );
        polyvec_inline double angle ( const real2& v1, const real2& v2 );
        polyvec_inline double to_range ( double angle );
        polyvec_inline double to_range ( double angle, double rangeStart );
        polyvec_inline double axis_deviation ( const real2& d );
        polyvec_inline double angle_between_vectors ( const real2& lhs, const real2& rhs );
        polyvec_inline double closest_horizontal_angle ( double angle );
        polyvec_inline double closest_vertical_angle ( double angle );

        // ---
        // Misc routines
        real2   project ( const real2& pt, const segment& s );
        bool    ray_intersect ( const real2& orig, const real2& dir, const segment& segment, double& ray_d );
        bool    point_on_line ( const segment& segment, const real2& pt );

        template <typename T>
        T sq ( const T& in ) {
            return in * in;
        }

        // ---
        // Fresnel functions from Cornucopia
        //almost full double-precision accuracy, using rational approximations
        void fresnel ( double xxa, double* ssa, double* cca );
        void fresnel ( const Eigen::VectorXd& t, Eigen::VectorXd* s, Eigen::VectorXd* c );

        //roughly single-precision accuracy, using polynomial approximations
        void fresnel_approx ( double xxa, double* ssa, double* cca );
        void fresnel_approx ( const Eigen::VectorXd& t, Eigen::VectorXd* s, Eigen::VectorXd* c ); //sse vectorized
    }

// ---
// OLD CODE, REMOVE
// ---
// Optimization parameterization of a curve
    struct CurveParams {
        typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::AutoAlign, 6, 1> ParamVec;

        //derivative of x and y w.r.t. parameters
        typedef Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor, 1, 6> ParamDer1;
        typedef Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::AutoAlign, 2, 6> ParamDer2;
        typedef Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::AutoAlign, 3, 6> ParamDer3;

        //derivative of x, y, angle, curvature w.r.t. parameters
        typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign, 4, 6> EndDer;

        enum Param {
            X = 0,     // start x
            Y,         // start y
            ANGLE,     // start angle
            LENGTH,    // arc length
            CURVATURE, // curvature
            DCURVATURE // curvature'
        };

        const ParamVec& params() const;
        void            set_params ( const ParamVec& params );

        virtual int     n_params() const = 0;
        virtual void    on_params_changed() = 0;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    protected:
        inline double _length() const {
            return _curve[LENGTH];
        }
        inline real2  _src() const {
            return _curve.head<2>();
        }
        inline double _src_angle() const {
            return _curve[ANGLE];
        }
        inline double _curvature() const {
            return _curve[CURVATURE];
        }
        inline double _curvature_der() const {
            return _curve[DCURVATURE];
        }

        ParamVec _curve;
    };

    struct Line;
    struct Arc;
    struct Clothoid;

    struct Curve : public CurveParams {
        enum Type {
            NO_CURVE,
            LINE = 1 << 0,
            ARC = 1 << 1,
            CLOTHOID = 1 << 2,
            BEZIER = 1 << 3,
            ALL = ( LINE | ARC | CLOTHOID | BEZIER )
        };

        static constexpr double START = 0.;
        static constexpr double END = 1.;

        static const Type types[3];

        // must have virtual destructor
        virtual ~Curve() = default;

        Line&     line();
        Arc&      arc();
        Clothoid& clothoid();

        const Line&     line() const;
        const Arc&      arc() const;
        const Clothoid& clothoid() const;

        virtual const char* type_str() const = 0;
        virtual int         type() const = 0;
        virtual Curve*      clone() const = 0;

        virtual real2 src() const = 0;
        virtual real2 dst() const = 0;

        virtual real2 pos ( double s ) const = 0;
        virtual real2 der ( double s ) const = 0;
        virtual real2 der2 ( double s ) const = 0;

        virtual double distance_to_sq ( const real2& point ) const;
        virtual double distance_to ( const real2& point ) const;
        virtual double project ( const real2& point ) const = 0;
        virtual bool trace_ray ( const real2& rayorig, const real2& raydir, double& sray, double& scurve, double from,
                                 double to ) const = 0;

        virtual double length() const = 0;
        virtual double angle ( double s ) const = 0;
        virtual double curvature ( double s ) const = 0;
        virtual double src_angle() const = 0;
        virtual double src_curvature() const = 0;
        virtual double dst_angle() const = 0;
        virtual double dst_curvature() const = 0;

        virtual void derivative_at ( double s, ParamDer2& out, ParamDer2& out_tan ) const = 0;
        virtual void derivative_at_end ( int continuity, EndDer& out ) const = 0;
        virtual void to_end_curvature_der ( Eigen::MatrixXd& der ) const = 0;
    };

    struct Line : public Curve {
        Line() = default; // clone()
        Line ( const real2& src, const real2& dst );

        // CurveParams interface
        int n_params() const override;
        void on_params_changed() override;

        // Curve interface
        const char* type_str() const override;
        int         type() const override;
        Curve*      clone() const override;

        real2 src() const override;
        real2 dst() const override;

        real2 pos ( double s ) const override;
        real2 der ( double s ) const override;
        real2 der2 ( double s ) const override;

        double project ( const real2& point ) const override;
        bool   trace_ray ( const real2& rayorig, const real2& raydir, double& sray, double& scurve, double from,
                           double to ) const override;

        double length() const override;
        double angle ( double s ) const override;
        double curvature ( double s ) const override;
        double src_angle() const override;
        double src_curvature() const override;
        double dst_angle() const override;
        double dst_curvature() const override;

        void derivative_at ( double s, ParamDer2& out, ParamDer2& out_tan ) const override;
        void derivative_at_end ( int continuity, EndDer& out ) const override;
        void to_end_curvature_der ( Eigen::MatrixXd& der ) const override {
            param_unused ( der );
        };

        // Line-specific
        real3     implicit_equation_coeffs() const;
        ParamDer3 implicit_equation_coeffs_derivative() const;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        real2 _der;
    };

    struct Arc : public Curve {
        Arc() = default;
        Arc ( const real2& src, double angle, double length, double curvature );
        Arc ( const real2& start, const real2& mid, const real2& end, bool& success );

        // CurveParams interface
        int n_params() const override;
        void on_params_changed() override;

        // Curve interface
        const char* type_str() const override;
        int         type() const override;
        Curve*      clone() const override;

        real2 src() const override;
        real2 dst() const override;

        real2 pos ( double s ) const override;
        real2 der ( double s ) const override;
        real2 der2 ( double s ) const override;

        double project ( const real2& point ) const override;
        bool   trace_ray ( const real2& rayorig, const real2& raydir, double& sray, double& scurve, double from,
                           double to ) const override;

        double length() const override;
        real2  center() const;
        double radius() const;
        double radius_abs() const;

        double angle ( double s ) const override;
        double curvature ( double s ) const override;
        double src_angle() const override;
        double src_curvature() const override;
        double dst_angle() const override;
        double dst_curvature() const override;

        void derivative_at ( double s, ParamDer2& out, ParamDer2& out_tan ) const override;
        void derivative_at_end ( int continuity, EndDer& out ) const override;
        void to_end_curvature_der ( Eigen::MatrixXd& der ) const override {
            param_unused ( der );
        };

        // Arc-specific
        void derivative_center ( ParamDer2& out ) const;
        void derivative_radius ( ParamDer1& out ) const;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        real2  _tangent; // start
        real2  _center; // if arc is not flat
        double _radius; // if arc is not flat, 1 / curvature
        double _angle_diff; // length * curvature
        bool   _flat; // if true arc is almost flat and we should use an approximation
    };

}

#include "geom.inl"


#endif // polyvec_geom_h_
