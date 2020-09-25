#pragma once

#include <polyvec/curve-tracer/curve_constraint.hpp>
#include <polyvec/curve-tracer/curve_parametrizations.hpp>

NAMESPACE_BEGIN(polyvec)

//Calculates the tangent direction of a line. The output of this constraint is the tangent angle w.r.t. the x-axis.
class GlobFitConstraint_LineDirection : public GlobFitConstraint
{
public:
	GlobFitConstraint_LineDirection(GlobFitLineParametrization* line, const GlobFitCurveParametrization::ParameterAddress& target);

	std::pair<double, Eigen::VectorXd> f();
	
private:	
	GlobFitLineParametrization* line;
};

NAMESPACE_END(polyvec)