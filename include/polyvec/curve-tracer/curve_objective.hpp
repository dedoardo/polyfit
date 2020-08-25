#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>

NAMESPACE_BEGIN(polyvec)

enum GlobFitObjectiveType {
	// A curve tries to fit a series of points as good as possible
	GLOBFIT_OBJECTIVE_FIT_POINTS = 0,
	// A curve tries to fit  tangents at a series of points as good as possible
	GLOBFIT_OBJECTIVE_FIT_TANGENTS,
	// A curve tries to have a certain position at a certain t value
	GLOBFIT_OBJECTIVE_FIX_POSITION,
	// A curve tries to have a certain tangent at a certain t value
	GLOBFIT_OBJECTIVE_FIX_TANGENT,
	// two curves_outline try to attain the same position at certain t values
	GLOBFIT_OBJECTIVE_SAME_POSITION,
	// two curves_outline try to attain the same tangents at certain t values
	GLOBFIT_OBJECTIVE_SAME_TANGENT,
	// two curves_outline try to attain the same  curvature at a certain t value
	GLOBFIT_OBJECTIVE_SAME_CURVATURE,
	// Objecitves to control bezier fairness
	GLOBFIT_OBJECTIVE_BEZIER_FAIRNESS,
	// A position on curve should lie on a certain line
	GLOBFIT_OBJECTIVE_POSITION_ON_CURVE_LIE_ON_LINE,
	// Regularize the scale and transition only beziers  or lines
	GLOBFIT_OBJECTIVE_SCALE_AND_TRANSITION_ONLY_REGULARIZER,
	// Minimize curvature variation
	GLOBFIT_OBJECTIVE_CURVATURE_VARIATION,
	GLOBFIT_OBJECTIVE_PARAMETER_BOUND,
	GLOBFIT_OBJECTIVE_LINE_REGULARIZATION,
	GLOBFIT_OBJECTIVE_LINE_ANGLE,

	GLOBFIT_OBJECTIVE_COUNT,
};

std::string globfitobjectivetype_as_string(GlobFitObjectiveType);

// =============================================================
//                           BASE CLASS
// =============================================================


class GlobFitObjective {
public:
	void set_weight(const double);
	double get_weight() const;
	double get_sqrt_weight() const;

	virtual int n_equations() const = 0;
	virtual GlobFitObjectiveType get_type() const = 0;

	virtual void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams) const = 0;

	GlobFitObjective() = default;
	virtual ~GlobFitObjective() = default;

	const std::vector<GlobFitCurveParametrization*>& get_curves() const { return _curves; }

protected:
	double _weight = 1;
	double _sqrt_weight = 1;

	// does not own them
	std::vector<GlobFitCurveParametrization*> _curves;

	int _n_sum_curve_params() const;
	virtual bool _is_all_data_provided() const;
};

NAMESPACE_END(polyvec)