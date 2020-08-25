// Header
#include <polyvec/curve-tracer/curve_bezier.hpp>

// C++ stl
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>

// Eigen
#include <Eigen/Core>

// Polyvec
#include <polyvec/geom.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/geometry/line.hpp>

namespace polyvec {

	// Taken from https://github.com/inkscape/lib2geom: src/2geom/bezier.h
	// Compute the value of a Bernstein-Bezier polynomial.
	// This method uses a Horner-like fast evaluation scheme.
	// param t Time value
	// param control_points control points (one per column)
	//
	template <typename Derived>
	Eigen::Matrix<double, Derived::RowsAtCompileTime,1>
		_bezier_value_at(const Eigen::MatrixBase<Derived>& control_points, const double t) {

		assert_break(t <= 1.);
		assert_break(t >= 0);

		double u = 1.0 - t;
		int degree = control_points.cols() - 1;

		Eigen::Matrix<double, Derived::RowsAtCompileTime, 1> ans;
		
		double bc = 1;
		double tn = 1;

		ans = control_points.col(0) * u;

		for (unsigned i = 1; i < degree; i++) {
			tn = tn * t;
			bc = bc * (degree - i + 1) / i;
			ans = (ans + tn * bc * control_points.col(i)) * u;
		}
		ans = (ans + tn * t * control_points.col(degree));

		return ans;
	}

	// From https://github.com/inkscape/lib2geom/: src/2geom/bezier.cpp
	Eigen::Matrix2Xd
		_derivative_bezier_curve(const Eigen::Matrix2Xd& control_points) {
		assert_break(control_points.cols() >= 2);
		Eigen::Matrix2Xd ans(2, control_points.cols() - 1);

		const int order = (int)control_points.cols() - 1;

		for (unsigned i = 0; i < ans.cols(); ++i) {
			ans.col(i) = order * (control_points.col(i + 1) - control_points.col(i));
		}

		return ans;
	}

	Eigen::Matrix3Xd
		_tesselate(const Eigen::Matrix2Xd& control_points, const int n_sampling) {
		Eigen::Matrix3Xd ans(3, n_sampling);

		for (int i = 0; i < n_sampling; ++i) {
			const double t = 1. / (n_sampling - 1) * i;
			ans.col(i) (0) = t;
			ans.col(i).tail<2>() = _bezier_value_at(control_points, t);
		}

		return ans;
	}

	const std::vector<double>&
		_get_gauss_quad_locs() {
		static std::vector<double> ans;
		static bool is_init = false;

		if (!is_init) {
			ans.resize(6);
			ans[0] = -9.3246951420315202781230155449399e-01L;
			ans[1] = -6.6120938646626451366139959501991e-01L;
			ans[2] = -2.3861918608319690863050172168071e-01L;
			ans[3] = -ans[2];
			ans[4] = -ans[1];
			ans[5] = -ans[0];

			for (double& a : ans) {
				a = a / 2. + 0.5;
			}

			is_init = true;
		}

		return ans;
	}

	const std::vector<double>&
		_get_gauss_quad_weights() {
		static std::vector<double> ans;
		static bool is_init = false;

		if (!is_init) {
			ans.resize(6);
			ans[0] = 1.7132449237917034504029614217273e-01L;
			ans[1] = 3.6076157304813860756983351383772e-01L;
			ans[2] = 4.6791393457269104738987034398955e-01L;
			ans[3] = ans[2];
			ans[4] = ans[1];
			ans[5] = ans[0];

			for (double& a : ans) {
				a = a / 2.;
			}

			is_init = true;
		}

		return ans;
	}

// ===================================================================
//                            Bezier Curve
// ===================================================================

    Eigen::Vector2d
    BezierCurve::pos ( const double t ) const {
        return _bezier_value_at ( _control_points_d0, t );
    }

    Eigen::Vector2d
    BezierCurve::dposdt ( const double t ) const {
        return _bezier_value_at ( _control_points_d1, t );
    }

    Eigen::Vector2d
    BezierCurve::dposdtdt ( const double t ) const {
        return _bezier_value_at ( _control_points_d2, t );
    }

	Eigen::Vector2d BezierCurve::dposdtdtdt(const double t) const
	{
		return _bezier_value_at(_control_points_d3, t);
	}

    Eigen::Matrix2Xd
    BezierCurve::dposdparams ( const double t ) const {
        const int n_points = n_control_points();

        // Eigen::Matrix2Xd ans(2, n_points);

        // Coefficient of bernstein monomials in a bernstein polynomial
        const double b0[] = {1., 0., 0., 0.};
        const double b1[] = {0., 1., 0., 0.};
        const double b2[] = {0., 0., 1., 0.};
        const double b3[] = {0., 0., 0., 1.};

        // Value of derivatives
		Eigen::MatrixXd control(4, 4);
		control <<
			1., 0., 0., 0.,
			0., 1., 0., 0.,
			0., 0., 1., 0.,
			0., 0., 0., 1.;

		auto dxdp = _bezier_value_at(control, t);

        Eigen::Matrix2Xd ans = Eigen::Matrix2Xd::Zero ( 2, n_points * 2 );
        ans.row ( 0 ) << dxdp(0), 0, dxdp(1), 0, dxdp(2), 0, dxdp(3), 0;
        ans.row ( 1 ) << 0, dxdp(0), 0, dxdp(1), 0, dxdp(2), 0, dxdp(3);

        return ans;
    }

    Eigen::Matrix2Xd
    BezierCurve::dposdtdparams ( const double t ) const {
        const int order = 3;
        const double drder = order;
        const int n_points = n_control_points();

        // Eigen::Matrix2Xd ans(2, n_points);

        // Coefficient of bernstein monomials in a bernstein polynomial
        const double b0_deriv[] = {-drder, 0., 0.};
        const double b1_deriv[] = {drder, -drder, 0.};
        const double b2_deriv[] = {0., drder, -drder};
        const double b3_deriv[] = {0., 0., drder};

		Eigen::MatrixXd control(4, 3);
		control <<
			-drder, 0., 0.,
			drder, -drder, 0.,
			0., drder, -drder,
			0., 0., drder;

        // Value of derivatives
		auto dxdtdp = _bezier_value_at(control, t);        

        Eigen::Matrix2Xd ans = Eigen::Matrix2Xd::Zero ( 2, n_points * 2 );
        ans.row ( 0 ) << dxdtdp(0), 0, dxdtdp(1), 0, dxdtdp(2), 0, dxdtdp(3), 0;
        ans.row ( 1 ) << 0, dxdtdp(0), 0, dxdtdp(1), 0, dxdtdp(2), 0, dxdtdp(3);

        return ans;
    }

    Eigen::Matrix2Xd
    BezierCurve::dposdtdtdparams ( const double t ) const {
        const int order = 3;
        const double drder = order*(order-1);
        const int n_points = n_control_points();

        // Eigen::Matrix2Xd ans(2, n_points);

        // Coefficient of bernstein monomials in a bernstein polynomial
		Eigen::MatrixXd control(4, 2);
		control <<
			drder, 0.,
			-2 * drder, drder,
			drder, -2 * drder,
			0., drder;

        // Value of derivatives
		auto dxdtdtdp = _bezier_value_at(control, t);

        Eigen::Matrix2Xd ans = Eigen::Matrix2Xd::Zero ( 2, n_points * 2 );
        ans.row ( 0 ) << dxdtdtdp(0), 0, dxdtdtdp(1), 0, dxdtdtdp(2), 0, dxdtdtdp(3), 0;
        ans.row ( 1 ) << 0, dxdtdtdp(0), 0, dxdtdtdp(1), 0, dxdtdtdp(2), 0, dxdtdtdp(3);

        return ans;
    }

	Eigen::Matrix2Xd BezierCurve::dposdtdtdtdparams(const double t) const
	{
		const int order = 3;
		const double factor = 6;
		const int n_points = n_control_points();

		// Coefficient of bernstein monomials in a bernstein polynomial
		Eigen::MatrixXd control(4, 1);
		control <<
			-factor,
			3 * factor,
			-3 * factor,
			factor;

		// Value of derivatives
		auto dxdtdp = _bezier_value_at(control, t);		

		Eigen::Matrix2Xd ans = Eigen::Matrix2Xd::Zero(2, n_points * 2);
		ans.row(0) << dxdtdp(0), 0, dxdtdp(1), 0, dxdtdp(2), 0, dxdtdp(3), 0;
		ans.row(1) << 0, dxdtdp(0), 0, dxdtdp(1), 0, dxdtdp(2), 0, dxdtdp(3);

		return ans;
	}

    double
    BezierCurve::project ( const Eigen::Vector2d& point ) const {

        // int closest_idx = -1;
        double closest_t = -1;
        double closest_dist2 = std::numeric_limits<double>::max();

        // First find an initial guess for the closest point
        for ( int i = 0; i < _tesselation.cols(); ++i ) {
            double dist2 = ( point - _tesselation.col ( i ).tail<2>() ).squaredNorm();

            if ( dist2 < closest_dist2 ) {
                // closest_idx = i;
                closest_t = _tesselation.col ( i ) ( 0 );
                closest_dist2 = dist2;
            }
        }

        // Now refine your guess using newton iterations
        auto newton_iteration = [&] ( double guess ) -> double {
            Eigen::Vector2d post, dert, der2t;
            post = this->pos ( guess );
            dert = this->dposdt ( guess );
            der2t = this->dposdtdt ( guess );

            double dot = dert.dot ( point - post );
            double dotDer = der2t.dot ( point - post ) - dert.squaredNorm();

            // Make sure the iteration does not make things worse
            if ( dotDer >= -1e-30 ) {
                return guess;
            }

            // Make sure the iteration does not shoot us out of the range
            return std::max ( 0., std::min ( 1., guess - dot / dotDer ) );
        };
        closest_t = newton_iteration ( closest_t );
        closest_t = newton_iteration ( closest_t );
        closest_t = newton_iteration ( closest_t );
        closest_t = newton_iteration ( closest_t );

        return closest_t;
    }

    Eigen::VectorXd
    BezierCurve::dtprojectdparams ( const double time, const Eigen::Vector2d& point ) const {
        const int n_params = n_control_points() * 2;
        const double tend = 1.;

        return SmoothCurveUtil::projection_derivatives (
                   n_params,
                   time,
                   tend,
                   point,
                   pos ( time ),
                   dposdt ( time ),
                   dposdtdt ( time ),
                   dposdparams ( time ),
                   dposdtdparams ( time ) );
    }

    Eigen::Matrix2Xd
    BezierCurve::dposprojectdparams ( const double t, const Eigen::VectorXd& dtprojectdparams ) const {
        return dposdparams ( t ) + dposdt ( t ) * dtprojectdparams.transpose();
    }

    Eigen::Matrix2Xd
    BezierCurve::dposdtprojectdparams ( const double t, const Eigen::VectorXd& dtprojectdparams ) const {
        return dposdtdparams ( t ) + dposdtdt ( t ) * dtprojectdparams.transpose();
    }


    void
    BezierCurve::set_control_points ( const Eigen::Matrix2Xd& control_points_d0_in ) {
        assert_break ( control_points_d0_in.cols() == n_control_points() );
        _control_points_d0 = control_points_d0_in;
        _control_points_d1 = _derivative_bezier_curve ( _control_points_d0 );
        _control_points_d2 = _derivative_bezier_curve ( _control_points_d1 );
		_control_points_d3 = _derivative_bezier_curve ( _control_points_d2 );
        _tesselation = _tesselate ( _control_points_d0, _n_tesselation );
    }

    const Eigen::Matrix2Xd&
    BezierCurve::get_control_points() const {
        return _control_points_d0;
    }

	void BezierCurve::set_params(const Eigen::VectorXd & params)
	{
		Eigen::Matrix2Xd controlPoints(2, 4);
		controlPoints << params[0], params[2], params[4], params[6], params[1], params[3], params[5], params[7];
		set_control_points(controlPoints);
	}

	Eigen::VectorXd BezierCurve::get_params() const
	{
		Eigen::VectorXd params(n_params());
		params << 
			_control_points_d0.coeff(0, 0), _control_points_d0.coeff(1, 0),
			_control_points_d0.coeff(0, 1), _control_points_d0.coeff(1, 1),
			_control_points_d0.coeff(0, 2), _control_points_d0.coeff(1, 2),
			_control_points_d0.coeff(0, 3), _control_points_d0.coeff(1, 3);
		return params;
	}

	//const Eigen::Matrix3Xd&
    //BezierCurve::get_tesselation3() {
    //    return _tesselation;
    //}

    Eigen::VectorXd
    BezierCurve::get_tesselationt() const {
        return _tesselation.row(0).transpose();
    }

    Eigen::Matrix2Xd
    BezierCurve::get_tesselation2() const {
        return _tesselation.bottomRows<2>();
    }


    double
    BezierCurve::length() const {
        const std::vector<double>& ww = _get_gauss_quad_weights();
        const std::vector<double>& tt = _get_gauss_quad_locs();
        double len = 0;

        for ( unsigned i = 0; i < tt.size(); ++i ) {
            len += ww[i] * dposdt ( tt[i] ).norm();
        }

        return len;
    }


    Eigen::VectorXd
    BezierCurve::dlengthdparams() {
        const int n_points = n_control_points();
        const double tol = 1e-10;
        const std::vector<double>& ww = _get_gauss_quad_weights();
        const std::vector<double>& tt = _get_gauss_quad_locs();

        Eigen::VectorXd dlendprams = Eigen::VectorXd ( n_points * 2 );
        dlendprams.setZero();

        for ( unsigned i = 0; i < tt.size(); ++i ) {
            Eigen::Vector2d tang = dposdt ( tt[i] );
            double tang_norm = std::max ( tang.norm(), tol );
            dlendprams += ww[i] * 1. / tang_norm * dposdtdparams ( tt[i] ).transpose() * tang;
        }

        return dlendprams;
    }    

	geom::aabb BezierCurve::get_bounding_box() const
	{
		geom::aabb b;
		for (int i = 0; i < 4; ++i)
			b.add(_control_points_d0.col(i));
		return b;
	}

	BezierCurve::BezierCurve(const Eigen::Matrix2Xd& C) {
		set_control_points(C);
	}

	std::pair<GlobFitCurve*, GlobFitCurve*> BezierCurve::split(double t) const
	{
		Eigen::Matrix2Xd cLeft(2, 4), cRight(2, 4);
		
		auto& p = _control_points_d0;

		// left curve
		cLeft.col(0) = p.col(0);
		cLeft.col(1) = _bezier_value_at(p.leftCols<2>(), t);
		cLeft.col(2) = _bezier_value_at(p.leftCols<3>(), t);
		cLeft.col(3) = _bezier_value_at(p.leftCols<4>(), t);
		
		// right curve
		cRight.col(0) = _bezier_value_at(p.rightCols<4>(), t);
		cRight.col(1) = _bezier_value_at(p.rightCols<3>(), t);
		cRight.col(2) = _bezier_value_at(p.rightCols<2>(), t);
		cRight.col(3) = p.col(3);

		return std::make_pair(new BezierCurve(cLeft), new BezierCurve(cRight));
	}
}
