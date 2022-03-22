#ifndef _polyvec_curve_line_h_
#define _polyvec_curve_line_h_

// C++ STL
#include <functional>
#include <vector>

// Eigen
#include <Eigen/Core>

#include <polyvec/curve-tracer/curve.hpp>

namespace polyvec {
	// Rather than wrapping curnucopia's line, it seems to be easier to
	// just rewrite this one.
	// The better practice .
	class GlobFitCurve_Line : public  GlobFitCurve {
	public:
		// GlobFitCurve Implementation
		GlobFitCurveType get_type() const override;
		int n_params() const override;

		Eigen::Vector2d pos(const double t) const override;
		Eigen::Vector2d dposdt(const double t) const override;
		Eigen::Vector2d dposdtdt(const double t) const override;
		Eigen::Vector2d dposdtdtdt(const double t) const override;

		Eigen::Matrix2Xd dposdparams(const double t) const override;
		Eigen::Matrix2Xd dposdtdparams(const double t) const override;
		Eigen::Matrix2Xd dposdtdtdparams(const double t) const override;
		Eigen::Matrix2Xd dposdtdtdtdparams(const double t) const override;

		double project(const Eigen::Vector2d& point) const override;

		Eigen::VectorXd dtprojectdparams(const double t, const Eigen::Vector2d& point) const override;
		Eigen::Matrix2Xd dposprojectdparams(const double t, const Eigen::VectorXd& dtprojectdparams) const override;
		Eigen::Matrix2Xd dposdtprojectdparams(const double t, const Eigen::VectorXd& dtprojectdparams) const override;

		void set_params(const Eigen::VectorXd& params) override;
		Eigen::VectorXd get_params() const override;

		Eigen::Matrix2Xd get_tesselation2() const override;
		Eigen::VectorXd get_tesselationt() const override;
		
		GlobFitCurve* clone() const override;

		geom::aabb get_bounding_box() const override;

		// Line functions
		void set_points(const Eigen::Matrix2d& points);
		void set_points(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1);
		Eigen::Matrix2d get_points();

		// Constructors
		GlobFitCurve_Line() = default;
		GlobFitCurve_Line(const Eigen::Matrix2d& points);
		GlobFitCurve_Line(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1);

		std::pair<GlobFitCurve*, GlobFitCurve*> split(double t) const override;

	private:
		static constexpr int _n_params = 4;
		bool _are_points_set = false;
		// First column is first point, second point is the second point
		Eigen::Matrix2d _points = Eigen::Matrix2d::Constant(1e10);
	};
}

#endif