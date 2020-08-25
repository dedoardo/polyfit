#ifndef _polyvec_curve_bezier_h_
#define _polyvec_curve_bezier_h_

// Eigen
#include <Eigen/Core>

#include "curve.hpp"

namespace polyvec {

//The intrinsic parameters of a Bezier are the control points in the order 
//  c1x, c1y, c2x, c2y, c3x, c3y, c4x, c4y
class BezierCurve : public GlobFitCurve{
public:
    BezierCurve() = default;

	GlobFitCurveType get_type() const { return GLOBFIT_CURVE_BEZIER; }
	int n_params() const { return 8; }

    Eigen::Vector2d pos ( const double t ) const override;
    Eigen::Vector2d dposdt ( const double t ) const override;
    Eigen::Vector2d dposdtdt ( const double t ) const override;
	Eigen::Vector2d dposdtdtdt ( const double t ) const override;

    Eigen::Matrix2Xd dposdparams ( const double t ) const override;
    Eigen::Matrix2Xd dposdtdparams ( const double t ) const override;
    Eigen::Matrix2Xd dposdtdtdparams ( const double t ) const override;
	Eigen::Matrix2Xd dposdtdtdtdparams ( const double t ) const override;

    double project ( const Eigen::Vector2d& point ) const override;

    Eigen::VectorXd dtprojectdparams ( const double t, const Eigen::Vector2d& point ) const override;
    Eigen::Matrix2Xd dposprojectdparams ( const double t, const Eigen::VectorXd& dtprojectdparams ) const override;
    Eigen::Matrix2Xd dposdtprojectdparams ( const double t, const Eigen::VectorXd& dtprojectdparams ) const override;

    void set_control_points ( const Eigen::Matrix2Xd& control_points_d0_in );
    const Eigen::Matrix2Xd& get_control_points() const;
    constexpr static int n_control_points() {
        return 4;
    }

	void set_params(const Eigen::VectorXd& params) override;
	Eigen::VectorXd get_params() const override;

    Eigen::Matrix2Xd get_tesselation2() const override;
    Eigen::VectorXd  get_tesselationt() const override;

    double length() const;
    Eigen::VectorXd dlengthdparams();

	GlobFitCurve* clone() const { return new BezierCurve(*this); }

	geom::aabb get_bounding_box() const;

	BezierCurve(const Eigen::Matrix2Xd& C);

	std::pair<GlobFitCurve*, GlobFitCurve*> split(double t) const;

private:
    constexpr static unsigned _n_tesselation = 50;    

    Eigen::Matrix2Xd _control_points_d0;
    Eigen::Matrix2Xd _control_points_d1;
    Eigen::Matrix2Xd _control_points_d2;
	Eigen::Matrix2Xd _control_points_d3;
    Eigen::Matrix3Xd _tesselation;
};

}

#endif // _polyvec_curve_bezier_h_