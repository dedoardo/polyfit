#pragma once

#include <polyvec/api.hpp>

#include <polyvec/curve-tracer/spline_types.hpp>
#include <polyvec/curve-tracer/curve_constraints.hpp>
#include <polyvec/curve-tracer/curve_objectives.hpp>
#include <polyvec/curve-tracer/curve_parametrizations.hpp>

NAMESPACE_BEGIN(polyvec)

struct FitOptions;

// Responsible for optimizing the parameters of a sequence of primitives.
// Owns least squares objectives and wraps the solver invocation.
class CurveFitter {
public:
	//Adds a curve to the fitter. The fitter will only keep a reference
	//to the curve and does not assume ownership.
	int  add_curve(CurvePrimitive* primitive);

	// Restricts the tangential movement of one of the endpoints of a line to not degenerate or flip directions
	void limit_endpoint_movement(GlobFitLineParametrization* line, bool limit_back);

	void make_g0(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next);
	void make_g1(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next);
	void make_g2(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next);
	void prescribe_angle(GlobFitLineParametrization* line, double angle);
	double get_angle_residuals() const;
	void increase_angle_weights(double factor);
    void make_overlap_start_to_end(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next);
    void make_overlap_end_to_start(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next);
	void solve();
    void solve_with_stiffening();
	void eval_error();
    void reset();

	CurveFitter(double raster_aabb_diagonal, bool final_fit);

	const polyvec::FitOptions& get_fit_options() const;

private:
	//Pointers to primitives to fit - the fitter does not own them
	std::vector<CurvePrimitive*> primitives;
	std::vector<GlobFitObjective*> objs;
	//std::vector<GlobFitObjective_BezierFairness*> obj_fairness;
	//std::vector<GlobFitObjective_FitPointsProject> obj_points_project;
	std::vector<GlobFitObjective_FitPointLength> obj_points_t;
	//std::vector<GlobFitObjective_FitTangentsProject> obj_tangents;
	std::vector<GlobFitObjective_FitTangentLength> obj_tangents_t;
	//std::vector<GlobFitObjective_PointLieOnLine> obj_t_on_line;
	std::vector<GlobFitObjective_SamePosition> obj_g0;
	std::vector<GlobFitObjective_SameTangent> obj_g1;
	std::vector<GlobFitObjective_CurvatureVariation> obj_curvVar;
	std::vector<GlobFitObjective_ParameterBound> obj_parameter_bound;
	std::vector<GlobFitObjective_SameCurvature> obj_g2;
	std::vector<GlobFitObjective_LineAngle> obj_lineAngle;
	std::vector<GlobFitObjective_LineRegularization> obj_lineRegularization;
	//std::vector<GlobFitObjective_ScaleAndTranslationOnlyRegularizer> obj_linreg;

	std::vector<GlobFitConstraint_LineDirection> line_g1_constraints;

	bool final_fit;
};

NAMESPACE_END(polyvec)