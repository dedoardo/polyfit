// polyvec
#include <polyvec/curve-tracer/curve_parametrizations.hpp>

// libc++
#include <algorithm>

using namespace polyvec;

#define ENABLE_DERIVATIVE_CHECK 0

// =====================================================
//         GlobFitBezierAngleBasedParametrization
// =====================================================

polyvec::GlobFitBezierAngleBasedParametrization::GlobFitBezierAngleBasedParametrization(
    std::shared_ptr<BezierCurve> curve, 
    const bool invert_front_coordinate_system,
    const bool invert_back_coordinate_system)
	: GlobFitCurveParametrization(curve),
    front_system_inverted(invert_front_coordinate_system),
    back_system_inverted(invert_back_coordinate_system)
{
	auto& pts = curve->get_control_points();
	Eigen::Vector2d beg_tangent = pts.col(1) - pts.col(0);
	Eigen::Vector2d end_tangent = pts.col(3) - pts.col(2);

	auto get_angle = [](const Eigen::Vector2d in) -> double {
		const double len = in.squaredNorm();
		const double sn = in.y() / len;
		const double cs = in.x() / len;
		return std::atan2(sn, cs);
	};

	startP = pts.col(0);
	endP = pts.col(3);
	front_system.set_secondary(beg_tangent.normalized() * (front_system_inverted ? -1. : +1.));
	back_system.set_secondary(end_tangent.normalized() * (back_system_inverted ? -1. : +1.));
	params.beg_pt.setZero();
	params.end_pt.setZero();
	params.beg_tangent_len = beg_tangent.norm();
	params.end_tangent_len = end_tangent.norm();
	params.beg_tangent_angle = get_angle(beg_tangent);
	params.end_tangent_angle = get_angle(end_tangent);

	calculate_jacobian();

	initial_params = get_params();
}

void GlobFitBezierAngleBasedParametrization::set_params(const Eigen::VectorXd & p)
{
	params.beg_tangent_len = p[0];
	params.end_tangent_len = p[1];
	params.beg_pt.x() = p[2];
	params.beg_pt.y() = p[3];
	params.end_pt.x() = p[4];
	params.end_pt.y() = p[5];
	params.beg_tangent_angle = p[6];
	params.end_tangent_angle = p[7];

	set_curve_parameters();
	calculate_jacobian();

#if ENABLE_DERIVATIVE_CHECK
	static bool checking = false;
	if (!checking)
	{
		checking = true;
		check_derivatives();
		checking = false;
	}
#endif
}

Eigen::VectorXd GlobFitBezierAngleBasedParametrization::get_params() const
{
	Eigen::VectorXd p(n_params());
	p << params.beg_tangent_len, 
        params.end_tangent_len, 
        params.beg_pt.x(),
        params.beg_pt.y(),
		params.end_pt.x(), 
        params.end_pt.y(),
        params.beg_tangent_angle, 
        params.end_tangent_angle;
	return p;
}

void GlobFitBezierAngleBasedParametrization::reverse() {
    std::swap(params.beg_pt, params.end_pt);
    std::swap(params.beg_tangent_angle, params.end_tangent_angle);
    std::swap(params.beg_tangent_len, params.end_tangent_len);
    std::swap(startP, endP);
    std::swap(front_system, back_system);
    std::swap(front_system_inverted, back_system_inverted);
    params.beg_tangent_angle += M_PI;
    params.end_tangent_angle += M_PI;
    set_curve_parameters();
    calculate_jacobian();
}

GlobFitCurveParametrization * polyvec::GlobFitBezierAngleBasedParametrization::clone()
{
	auto copy = new GlobFitBezierAngleBasedParametrization(std::shared_ptr<BezierCurve>(static_cast<BezierCurve*>(get_curve()->clone())));
	copy->params = params;
	copy->startP = startP;
	copy->endP = endP;
	copy->front_system = front_system;
	copy->back_system = back_system;
    copy->front_system_inverted = front_system_inverted;
    copy->back_system_inverted = back_system_inverted;
	copy->copy_fixed_parameter_from(*this);
	return copy;
}

GlobFitCurveParametrization * GlobFitBezierAngleBasedParametrization::create_for_curve(std::shared_ptr<GlobFitCurve> curve) const
{
	return new GlobFitBezierAngleBasedParametrization(std::static_pointer_cast<BezierCurve>(curve));
}

void GlobFitBezierAngleBasedParametrization::reduce_degrees_of_freedom(DofOptions::Type options)
{
	if(options & DofOptions::FIX_FRONT_TANGENT)
		fix_parameter(6);

	if (options & DofOptions::KEEP_FRONT_ON_BISECTOR)
		fix_parameter(3, 0.0);

	if(options & DofOptions::FIX_BACK_TANGENT)
		fix_parameter(7);

	if(options & DofOptions::KEEP_BACK_ON_BISECTOR)
		fix_parameter(5, 0.0);
}

void GlobFitBezierAngleBasedParametrization::set_curve_parameters()
{
	const double cs_beg = cos(params.beg_tangent_angle);
	const double sn_beg = sin(params.beg_tangent_angle);
	const double cs_end = cos(params.end_tangent_angle);
	const double sn_end = sin(params.end_tangent_angle);

	Eigen::VectorXd p(get_curve()->n_params());

	auto p0 = startP + front_system.matrix() * params.beg_pt;
	auto p3 = endP + back_system.matrix() * params.end_pt;
	p.block<2, 1>(0, 0) = p0;
	p.block<2, 1>(2, 0) = p0 + params.beg_tangent_len * Eigen::Vector2d(cs_beg, sn_beg);
	p.block<2, 1>(4, 0) = p3 - params.end_tangent_len * Eigen::Vector2d(cs_end, sn_end);
	p.block<2, 1>(6, 0) = p3;

	get_curve()->set_params(p);
}

void GlobFitBezierAngleBasedParametrization::calculate_jacobian()
{
	const int xid = 0;
	const int yid = 1;
	const int n_dim = 2;
	dCurveParamsdParams.resize(get_curve()->n_params(), n_params());

	const double cs_beg = cos(params.beg_tangent_angle);
	const double sn_beg = sin(params.beg_tangent_angle);
	const double cs_end = cos(params.end_tangent_angle);
	const double sn_end = sin(params.end_tangent_angle);

	const double lbeg = params.beg_tangent_len;
	const double lend = params.end_tangent_len;

	dCurveParamsdParams.row(0 * n_dim + xid) << 0, 0, front_system.matrix().coeff(0, 0), front_system.matrix().coeff(0, 1), 0, 0, 0, 0;
	dCurveParamsdParams.row(0 * n_dim + yid) << 0, 0, front_system.matrix().coeff(1, 0), front_system.matrix().coeff(1, 1), 0, 0, 0, 0;
	//
	dCurveParamsdParams.row(1 * n_dim + xid) << cs_beg, 0, front_system.matrix().coeff(0, 0), front_system.matrix().coeff(0, 1), 0, 0, -lbeg * sn_beg, 0;
	dCurveParamsdParams.row(1 * n_dim + yid) << sn_beg, 0, front_system.matrix().coeff(1, 0), front_system.matrix().coeff(1, 1), 0, 0, lbeg* cs_beg, 0;
	//
	dCurveParamsdParams.row(2 * n_dim + xid) << 0, -cs_end, 0, 0, back_system.matrix().coeff(0, 0), back_system.matrix().coeff(0, 1), 0, lend* sn_end;
	dCurveParamsdParams.row(2 * n_dim + yid) << 0, -sn_end, 0, 0, back_system.matrix().coeff(1, 0), back_system.matrix().coeff(1, 1), 0, -lend * cs_end;
	//
	dCurveParamsdParams.row(3 * n_dim + xid) << 0, 0, 0, 0, back_system.matrix().coeff(0, 0), back_system.matrix().coeff(0, 1), 0, 0;
	dCurveParamsdParams.row(3 * n_dim + yid) << 0, 0, 0, 0, back_system.matrix().coeff(1, 0), back_system.matrix().coeff(1, 1), 0, 0;
}

// =====================================================
//         GlobFitBezierShiftBasedParametrization
// =====================================================

polyvec::GlobFitBezierShiftBasedParametrization::GlobFitBezierShiftBasedParametrization(std::shared_ptr<BezierCurve> curve)
	: GlobFitCurveParametrization(curve)
{
	auto controlPoints = curve->get_control_points();
	startP = controlPoints.col(0);
	endP = controlPoints.col(3);
	startTang = (controlPoints.col(1) - controlPoints.col(0)); //will be normalized shortly
	endTang = (controlPoints.col(3) - controlPoints.col(2)); //will be normalized shortly

	params.beg_tangent_len = startTang.norm();
	startTang /= params.beg_tangent_len;
	params.end_tangent_len = endTang.norm();
	endTang /= params.end_tangent_len;

	params.beg_shift = 0;
	params.end_shift = 0;

	calculate_jacobian();

	initial_params = get_params();
}

void polyvec::GlobFitBezierShiftBasedParametrization::set_params(const Eigen::VectorXd & p)
{
	params.beg_tangent_len = p(0);
	params.end_tangent_len = p(1);
	params.beg_shift = p(2);
	params.end_shift = p(3);

	set_curve_parameters();
	calculate_jacobian();

#if ENABLE_DERIVATIVE_CHECK
	static bool checking = false;
	if (!checking)
	{
		checking = true;
		check_derivatives();
		checking = false;
	}
#endif
}

Eigen::VectorXd polyvec::GlobFitBezierShiftBasedParametrization::get_params() const
{
	Eigen::VectorXd p(4);
	p << params.beg_tangent_len, params.end_tangent_len, params.beg_shift, params.end_shift;
	return p;
}

GlobFitCurveParametrization * polyvec::GlobFitBezierShiftBasedParametrization::clone()
{
	auto copy = new GlobFitBezierShiftBasedParametrization(std::shared_ptr<BezierCurve>(static_cast<BezierCurve*>(get_curve()->clone())));
	copy->params = params;
	copy->copy_fixed_parameter_from(*this);
	return copy;
}

GlobFitCurveParametrization * GlobFitBezierShiftBasedParametrization::create_for_curve(std::shared_ptr<GlobFitCurve> curve) const
{
	return new GlobFitBezierShiftBasedParametrization(std::static_pointer_cast<BezierCurve>(curve));
}

void polyvec::GlobFitBezierShiftBasedParametrization::set_curve_parameters()
{
	Eigen::VectorXd p(get_curve()->n_params());
	p.block<2, 1>(0, 0) = startP + params.beg_shift * Eigen::Vector2d(startTang.y(), -startTang.x());
	p.block<2, 1>(6, 0) = endP + params.end_shift * Eigen::Vector2d(endTang.y(), -endTang.x());
	p.block<2, 1>(2, 0) = p.block<2, 1>(0, 0) + params.beg_tangent_len * startTang;
	p.block<2, 1>(4, 0) = p.block<2, 1>(6, 0) - params.end_tangent_len * endTang;
	get_curve()->set_params(p);
}

void polyvec::GlobFitBezierShiftBasedParametrization::calculate_jacobian()
{
	dCurveParamsdParams.resize(get_curve()->n_params(), n_params());

	dCurveParamsdParams.row(0) << 0, 0, startTang.y(), 0;
	dCurveParamsdParams.row(1) << 0, 0, -startTang.x(), 0;
	dCurveParamsdParams.row(6) << 0, 0, 0, endTang.y();
	dCurveParamsdParams.row(7) << 0, 0, 0, -endTang.x();

	dCurveParamsdParams.row(2) << startTang.x(), 0, startTang.y(), 0;
	dCurveParamsdParams.row(3) << startTang.y(), 0, -startTang.x(), 0;
	dCurveParamsdParams.row(4) << 0, -endTang.x(), 0, endTang.y();
	dCurveParamsdParams.row(5) << 0, -endTang.y(), 0, -endTang.x();
}

// =====================================================
//         GlobFitLineParametrization
// =====================================================

polyvec::GlobFitLineParametrization::GlobFitLineParametrization(
    std::shared_ptr<GlobFitCurve_Line> curve,
    const bool invert_front_coordinate_system,
    const bool invert_back_coordinate_system)
	: GlobFitCurveParametrization(curve),
    front_system_inverted(invert_front_coordinate_system),
    back_system_inverted(invert_back_coordinate_system)
{
	const auto& pts = curve->get_points();
	startP = pts.col(0);
	endP = pts.col(1);

	auto dir_norm = (endP - startP).normalized();
	front_system.set_secondary(dir_norm * (invert_front_coordinate_system ? -1. : +1.));
	back_system.set_secondary(dir_norm * (invert_back_coordinate_system ? -1. : +1.));

	original_orthogonal_dir = front_system.primary();

	params.front.setZero();
	params.back.setZero();

	calculate_jacobian();

	initial_params = get_params();
}

void polyvec::GlobFitLineParametrization::set_params(const Eigen::VectorXd & p)
{
	params.front = p.segment<2>(0);
	params.back = p.segment<2>(2);

	set_curve_parameters();
	calculate_jacobian();

#if ENABLE_DERIVATIVE_CHECK
	static bool checking = false;
	if (!checking)
	{
		checking = true;
		check_derivatives();
		checking = false;
	}
#endif
}

Eigen::VectorXd polyvec::GlobFitLineParametrization::get_params() const
{
	Eigen::VectorXd p(4);
	p << params.front, params.back;
	return p;
}

GlobFitCurveParametrization * polyvec::GlobFitLineParametrization::clone()
{
	auto copy = new GlobFitLineParametrization(std::shared_ptr<GlobFitCurve_Line>(static_cast<GlobFitCurve_Line*>(get_curve()->clone())));
	copy->params = params;
	copy->copy_fixed_parameter_from(*this);
	copy->original_orthogonal_dir = original_orthogonal_dir;
	copy->front_system = front_system;
	copy->back_system = back_system;
    copy->back_system_inverted = back_system_inverted;
    copy->front_system_inverted = front_system_inverted;
	copy->startP = startP;
	copy->endP = endP;
	return copy;
}

GlobFitCurveParametrization * GlobFitLineParametrization::create_for_curve(std::shared_ptr<GlobFitCurve> curve) const
{
	return new GlobFitLineParametrization(std::static_pointer_cast<GlobFitCurve_Line>(curve));
}

void polyvec::GlobFitLineParametrization::set_front_primary_axis(const Eigen::Vector2d & d, bool rescale)
{
	front_system.set_primary(d * (front_system_inverted ? -1.0 : 1.0));
}

void polyvec::GlobFitLineParametrization::set_back_primary_axis(const Eigen::Vector2d & d, bool rescale)
{	
	back_system.set_primary(d * (back_system_inverted ? -1.0 : 1.0));
}

void GlobFitLineParametrization::merge_with(GlobFitLineParametrization * other)
{
	this->back_system = other->back_system;
	this->params.back = other->params.back;
	this->endP = other->endP;
	
	if (!other->parameter_info.empty())
	{
		parameter_info.resize(n_params());
		this->parameter_info[2] = other->parameter_info[2];
		this->parameter_info[3] = other->parameter_info[3];
	}

	set_curve_parameters();
}

void GlobFitLineParametrization::reverse() {
    std::swap(params.front, params.back);
    std::swap(front_system, back_system);
    std::swap(front_system_inverted, back_system_inverted);
    std::swap(startP, endP);
    set_curve_parameters();
    calculate_jacobian();
}

void GlobFitLineParametrization::reduce_degrees_of_freedom(DofOptions::Type options)
{
	if (options & DofOptions::FIX_FRONT_TANGENT || options & DofOptions::FIX_BACK_TANGENT)
	{
		couple_parameter(0, this, 2);

		// if the secondary axes have an impact on orthogonal movement, fix them to zero
		if (std::abs(original_orthogonal_dir.dot(front_system.secondary())) > PF_EPS)
			fix_parameter(1, 0.0);
		if (std::abs(original_orthogonal_dir.dot(back_system.secondary())) > PF_EPS)
			fix_parameter(3, 0.0);
	}

	if (options & DofOptions::KEEP_FRONT_ON_BISECTOR)
		fix_parameter(1, 0.0);

	if (options & DofOptions::KEEP_BACK_ON_BISECTOR)
		fix_parameter(3, 0.0);
}

void polyvec::GlobFitLineParametrization::set_curve_parameters()
{
	Eigen::VectorXd p(get_curve()->n_params());
	p <<
		startP.x() + front_system.primary().x() * params.front.x() + front_system.secondary().x() * params.front.y(),
		startP.y() + front_system.primary().y() * params.front.x() + front_system.secondary().y() * params.front.y(),
		endP.x() + back_system.primary().x() * params.back.x() + back_system.secondary().x() * params.back.y(),
		endP.y() + back_system.primary().y() * params.back.x() + back_system.secondary().y() * params.back.y();
	get_curve()->set_params(p);
}

void polyvec::GlobFitLineParametrization::calculate_jacobian()
{
	dCurveParamsdParams.resize(get_curve()->n_params(), n_params());

	dCurveParamsdParams.row(0) << front_system.primary().x(), front_system.secondary().x(), 0, 0;
	dCurveParamsdParams.row(1) << front_system.primary().y(), front_system.secondary().y(), 0, 0;
	dCurveParamsdParams.row(2) << 0, 0, back_system.primary().x(), back_system.secondary().x();
	dCurveParamsdParams.row(3) << 0, 0, back_system.primary().y(), back_system.secondary().y();
}

void CoordinateSystem::set_primary(const Eigen::Vector2d &p)
{
	_matrix.col(0) = p;
	_matrix.col(1) = Eigen::Vector2d(-p.y(), p.x());
}

void CoordinateSystem::set_secondary(const Eigen::Vector2d &s)
{
	_matrix.col(1) = s;
	_matrix.col(0) = Eigen::Vector2d(s.y(), -s.x());
}
