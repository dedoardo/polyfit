#include <polyvec/curve-tracer/curve_constraints.hpp>

using namespace polyvec;

GlobFitConstraint_LineDirection::GlobFitConstraint_LineDirection(GlobFitLineParametrization* line, GlobFitCurveParametrization::ParameterAddress& target)
	: line(line)
{	
	this->source_params.push_back(GlobFitCurveParametrization::ParameterAddress(line, 0));
	this->source_params.push_back(GlobFitCurveParametrization::ParameterAddress(line, 1));
	this->source_params.push_back(GlobFitCurveParametrization::ParameterAddress(line, 2));
	this->source_params.push_back(GlobFitCurveParametrization::ParameterAddress(line, 3));

	this->target_param = target;
}

std::pair<double, Eigen::VectorXd> GlobFitConstraint_LineDirection::f()
{
	Eigen::Vector2d diff = line->get_curve()->pos(1.0) - line->get_curve()->pos(0.0);	
	auto dDiffdParams = line->dposdparams(1.0) - line->dposdparams(0.0);
	
	double angle = std::atan2(diff.y(), diff.x());
	
	double xSq = diff.x() * diff.x();
	double ySq = diff.y() * diff.y();
	double denom = xSq + ySq;

	Eigen::Vector2d dAngledDiff(-diff.y() / denom, diff.x() / denom);
	Eigen::VectorXd dAngledParams = (dAngledDiff.transpose() * dDiffdParams).transpose();

	return std::make_pair(angle, dAngledParams);
}