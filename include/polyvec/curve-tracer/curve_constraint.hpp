#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>

NAMESPACE_BEGIN(polyvec)

//Represents a constraint of the form
//  curve parameter = f(other curve parameters)
class GlobFitConstraint
{
public:
	//Calculates the target parameter and its gradient from the current source parameters
	virtual std::pair<double, Eigen::VectorXd> f() = 0;
	
	const GlobFitCurveParametrization::ParameterAddress& get_target_param() const { return target_param; }
	const std::vector<GlobFitCurveParametrization::ParameterAddress>& get_source_params() const { return source_params; }
protected:
	GlobFitCurveParametrization::ParameterAddress target_param;
	std::vector<GlobFitCurveParametrization::ParameterAddress> source_params;
};

NAMESPACE_END(polyvec)