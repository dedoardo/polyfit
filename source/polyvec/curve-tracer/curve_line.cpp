#include <polyvec/curve-tracer/curve_line.hpp>
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/misc.hpp>

using namespace polyvec;

GlobFitCurveType GlobFitCurve_Line::get_type() const {
	return GLOBFIT_CURVE_LINE;
}

int GlobFitCurve_Line::n_params() const {
	return _n_params;
}

Eigen::Vector2d GlobFitCurve_Line::pos(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	return _points.col(0) + (_points.col(1) - _points.col(0)) * t;
}

Eigen::Vector2d GlobFitCurve_Line::dposdt(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	return (_points.col(1) - _points.col(0));
}

Eigen::Vector2d GlobFitCurve_Line::dposdtdt(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	return Eigen::Vector2d::Zero();
}

Eigen::Vector2d GlobFitCurve_Line::dposdtdtdt(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	return Eigen::Vector2d::Zero();
}


Eigen::Matrix2Xd GlobFitCurve_Line::dposdparams(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	Eigen::Matrix<double, 2, _n_params> ans;
	ans.row(0) << 1 - t, 0, t, 0;
	ans.row(1) << 0, 1 - t, 0, t;
	return ans;
}

Eigen::Matrix2Xd GlobFitCurve_Line::dposdtdparams(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	Eigen::Matrix<double, 2, _n_params> ans;
	ans.row(0) << -1, 0, 1, 0;
	ans.row(1) << 0, -1, 0, 1;
	return ans;
}

Eigen::Matrix2Xd GlobFitCurve_Line::dposdtdtdparams(const double t) const {
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	Eigen::Matrix<double, 2, _n_params> ans;
	ans.setZero();
	return ans;
}

Eigen::Matrix2Xd polyvec::GlobFitCurve_Line::dposdtdtdtdparams(const double t) const
{
	assert_break(_are_points_set);
	assert_break(t <= GlobFitCurve::t_end);
	assert_break(t >= 0);
	Eigen::Matrix<double, 2, _n_params> ans;
	ans.setZero();
	return ans;
}

double GlobFitCurve_Line::project(const Eigen::Vector2d& point) const {
	assert_break(_are_points_set);
	const Eigen::Vector2d dir = (_points.col(1) - _points.col(0));
	const double boundless_t =
		dir.dot(point - _points.col(0)) / dir.squaredNorm();
	return std::min(t_end, std::max(0., boundless_t));
}

Eigen::VectorXd GlobFitCurve_Line::dtprojectdparams(
	const double time, const Eigen::Vector2d& point) const {
	return SmoothCurveUtil::projection_derivatives(
		_n_params, time, GlobFitCurve::t_end, point, pos(time), dposdt(time),
		dposdtdt(time), dposdparams(time), dposdtdparams(time));
}

Eigen::Matrix2Xd GlobFitCurve_Line::dposprojectdparams(
	const double t, const Eigen::VectorXd& dtprojectdparams) const {
	return dposdparams(t) + dposdt(t) * dtprojectdparams.transpose();
}

Eigen::Matrix2Xd GlobFitCurve_Line::dposdtprojectdparams(
	const double t, const Eigen::VectorXd& dtprojectdparams) const {
	return dposdtdparams(t) + dposdtdt(t) * dtprojectdparams.transpose();
}

void GlobFitCurve_Line::set_params(const Eigen::VectorXd& params) {
	_are_points_set = true;
	_points = misc::reshaped(params, 2, 2);
}

Eigen::VectorXd GlobFitCurve_Line::get_params() const {
	return misc::reshaped(_points, 4, 1);
}

Eigen::Matrix2Xd GlobFitCurve_Line::get_tesselation2() const {
	return _points;
}

Eigen::VectorXd GlobFitCurve_Line::get_tesselationt() const {

	return Eigen::Vector2d(0, 1);
}

void GlobFitCurve_Line::set_points(const Eigen::Matrix2d& points) {
	_are_points_set = true;
	_points = points;
}

void GlobFitCurve_Line::set_points(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) {
	_are_points_set = true;
	_points.col(0) = p0;
	_points.col(1) = p1;
}

Eigen::Matrix2d GlobFitCurve_Line::get_points() {
	return _points;
}

GlobFitCurve* GlobFitCurve_Line::clone() const {
	return new GlobFitCurve_Line(*this);
}

geom::aabb polyvec::GlobFitCurve_Line::get_bounding_box() const
{
	geom::aabb b;
	for (int i = 0; i < 2; ++i)
		b.add(_points.col(i));
	return b;
}

GlobFitCurve_Line::GlobFitCurve_Line(const Eigen::Matrix2d& points) {
	_are_points_set = true;
	_points = points;
}

GlobFitCurve_Line::GlobFitCurve_Line(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) {
	_are_points_set = true;
	_points.col(0) = p0;
	_points.col(1) = p1;
}

std::pair<GlobFitCurve*, GlobFitCurve*> polyvec::GlobFitCurve_Line::split(double t) const
{
	throw std::runtime_error("Not implemented");
}
