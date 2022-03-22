/*
	Abstracts the storage of the curve parameters and the computation of the jacobian
	for all the curves which are used in the curve fitting.

	It is constructed from an existing curve object and transforms the parameter in the
	representation more useful for optimization with a possibly different number of degrees
	of freedom. 
*/
#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curve_line.hpp>

NAMESPACE_BEGIN (polyvec)

namespace DofOptions
{
	using Type = unsigned int;
	static const Type FIX_FRONT_TANGENT = 1;
	static const Type FIX_BACK_TANGENT = 2;
	static const Type KEEP_FRONT_ON_BISECTOR = 4;
	static const Type KEEP_BACK_ON_BISECTOR = 8;
};

class GlobFitCurveParametrization {
public:


	//we identify parameters by the curve pointer and the parameter id internal to that curve
	struct ParameterAddress
	{
		GlobFitCurveParametrization* curve = nullptr;
		int internal_parameter = -1;

		ParameterAddress() { }
		ParameterAddress(GlobFitCurveParametrization* curve, int p)
			: curve(curve), internal_parameter(p)
		{ }
	};

	struct ParameterInfo
	{
		//the id of the variable in the optimization problem
		int variable_id = -1;

		//the id of the hard constraint that produces this parameter
		int constraint_id = -1;		
					
		std::vector<ParameterAddress> coupled_parameters;

		//specifies the variable if the variable is fixed
		double fixedValue = std::numeric_limits<double>::quiet_NaN();

		double fixedValuePropagated = std::numeric_limits<double>::quiet_NaN();
	};

	virtual ~GlobFitCurveParametrization() = default;
	GlobFitCurveParametrization(std::shared_ptr<GlobFitCurve> curve);

	virtual int n_params() const;

	Eigen::Matrix2Xd dposdparams(const double t) const;
	Eigen::Matrix2Xd dposdtdparams(const double t) const;
	Eigen::Matrix2Xd dposdtdtdparams(const double t) const;
	Eigen::Matrix2Xd dposdtdtdtdparams(const double t) const;

	virtual void set_params(const Eigen::VectorXd& params);
	virtual Eigen::VectorXd get_params() const;

	std::shared_ptr<GlobFitCurve> get_curve() const { return curve; }

	const Eigen::VectorXd get_initial_parameters() const { return initial_params; }

	std::vector<ParameterInfo>& get_parameter_info();

	//Creates an exact deep copy of the current parametrization.
	virtual GlobFitCurveParametrization* clone();

	//Creates a new parametrization of the same type for a new curve.
	virtual GlobFitCurveParametrization* create_for_curve(std::shared_ptr<GlobFitCurve> curve) const;

    virtual void reverse() { };

	//Couples a parameter of this curve with a parameter of another curve. Therefore, no new variable for this parameter
	//will be created, but the variable for the original curve's parameter will be used.
	void couple_parameter(int thisParameter, GlobFitCurveParametrization* curve, int otherParameter);		

	void uncouple_all_parameters();

	//Fixes a given parameter to a given value. If no value is given, the current value will be used.
	void fix_parameter(int parameter, double value = std::numeric_limits<double>::quiet_NaN());
	void copy_fixed_parameter_from(const GlobFitCurveParametrization& other);	

	virtual void reduce_degrees_of_freedom(DofOptions::Type /*options*/) { throw std::runtime_error("Tangent fixing must be implemented in a subclass"); }

protected:
	//Jacobian of the curve parameters w.r.t. the parameters of this parametrization.
	//If this is empty, the curve's parametrization is used.
	Eigen::MatrixXd dCurveParamsdParams;

	Eigen::VectorXd initial_params;

	void check_derivatives();

	std::shared_ptr<GlobFitCurve> curve;
	
	std::vector<ParameterInfo> parameter_info;
};

NAMESPACE_END(polyvec)
