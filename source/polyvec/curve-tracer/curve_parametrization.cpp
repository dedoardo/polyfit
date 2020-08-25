#include <polyvec/curve-tracer/curve_parametrization.hpp>

#include <iostream>

using namespace polyvec;

GlobFitCurveParametrization::GlobFitCurveParametrization(std::shared_ptr<GlobFitCurve> curve)
	: curve(curve)
{
	initial_params = curve->get_params();
}

int polyvec::GlobFitCurveParametrization::n_params() const
{
	return curve->n_params();
}

Eigen::Matrix2Xd polyvec::GlobFitCurveParametrization::dposdparams(const double t) const
{
	if (dCurveParamsdParams.size() == 0)
		return curve->dposdparams(t);
	else
		return curve->dposdparams(t) * dCurveParamsdParams;
}

Eigen::Matrix2Xd polyvec::GlobFitCurveParametrization::dposdtdparams(const double t) const
{
	if (dCurveParamsdParams.size() == 0)
		return curve->dposdtdparams(t);
	else
		return curve->dposdtdparams(t) * dCurveParamsdParams;
}

Eigen::Matrix2Xd polyvec::GlobFitCurveParametrization::dposdtdtdparams(const double t) const
{
	if (dCurveParamsdParams.size() == 0)
		return curve->dposdtdtdparams(t);
	else
		return curve->dposdtdtdparams(t) * dCurveParamsdParams;
}

Eigen::Matrix2Xd polyvec::GlobFitCurveParametrization::dposdtdtdtdparams(const double t) const
{
	if (dCurveParamsdParams.size() == 0)
		return curve->dposdtdtdtdparams(t);
	else
		return curve->dposdtdtdtdparams(t) * dCurveParamsdParams;
}

void polyvec::GlobFitCurveParametrization::set_params(const Eigen::VectorXd & params)
{
	curve->set_params(params);
}

Eigen::VectorXd polyvec::GlobFitCurveParametrization::get_params() const
{
	return curve->get_params();
}

std::vector<polyvec::GlobFitCurveParametrization::ParameterInfo>& polyvec::GlobFitCurveParametrization::get_parameter_info()
{
	parameter_info.resize(n_params());
	return parameter_info;
}

GlobFitCurveParametrization * polyvec::GlobFitCurveParametrization::clone()
{
	auto copy = new GlobFitCurveParametrization(std::shared_ptr<GlobFitCurve>(get_curve()->clone()));	
	copy->copy_fixed_parameter_from(*this);
	return copy;
}

GlobFitCurveParametrization * GlobFitCurveParametrization::create_for_curve(std::shared_ptr<GlobFitCurve> curve) const
{
	return new GlobFitCurveParametrization(curve);
}

void polyvec::GlobFitCurveParametrization::couple_parameter(int thisParameter, GlobFitCurveParametrization * curve, int otherParameter)
{
	parameter_info.resize(n_params());
	auto& info = parameter_info[thisParameter];
	info.coupled_parameters.emplace_back(curve, otherParameter);
}

void polyvec::GlobFitCurveParametrization::uncouple_all_parameters()
{
	parameter_info.resize(n_params());
	for (auto& p : parameter_info)
		p.coupled_parameters.clear();
}

void polyvec::GlobFitCurveParametrization::fix_parameter(int parameter, double value)
{
	parameter_info.resize(n_params());
	if (std::isnan(value))
		value = this->get_params()(parameter);
	parameter_info[parameter].fixedValue = value;
}

void polyvec::GlobFitCurveParametrization::copy_fixed_parameter_from(const GlobFitCurveParametrization & other)
{
	if (other.parameter_info.empty())
		return;
	parameter_info.resize(n_params());
	for (int i = 0; i < n_params(); ++i)
	{
		parameter_info[i].fixedValue = other.parameter_info[i].fixedValue;
	}
}

void polyvec::GlobFitCurveParametrization::check_derivatives() {
	auto x = get_params();
	auto curve_params = curve->get_params();

	auto jac = dCurveParamsdParams;
		
	const double eps = 0.001;

	Eigen::VectorXd curve_paramsPlus, curve_paramsMinus;

	for (int ip = 0; ip < x.size(); ++ip)
	{
		auto xCopy = x;
		xCopy(ip) += eps;
		set_params(xCopy);
		curve_paramsPlus = curve->get_params();

		xCopy(ip) -= 2 * eps;
		set_params(xCopy);
		curve_paramsMinus = curve->get_params();

		Eigen::VectorXd expectedDeriv = (curve_paramsPlus - curve_paramsMinus) / (2 * eps);
		for (int i = 0; i < curve_params.size(); ++i)
		{
			auto error = jac.coeff(i, ip) - expectedDeriv(i);
			auto relativeError = error / std::max(std::abs(curve_params(i)), std::abs(jac.coeff(i, ip)));

			if (std::abs(relativeError) > 0.001 && std::abs(error) > 1e-10)
			{
				std::cout << "Wrong parametrization derivative for " << typeid(*this).name() << " at parameter " << ip << ", curve parameter " << i << ". Numerical derivative: " << expectedDeriv(i) <<
					", computed derivative: " << jac.coeff(i, ip) << ", curve parameter value: " << curve_params(i) << ", parameter value: " << x(ip) << std::endl;				
			}
		}
	}
	set_params(x);
}