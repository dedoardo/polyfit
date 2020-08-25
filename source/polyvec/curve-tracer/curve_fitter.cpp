#include <polyvec/curve-tracer/curve_fitter.hpp>

#include <polyvec/curve-tracer/curve_solver.hpp>
#include <polyvec/curve-tracer/fit_options.hpp>

#include <polyvec/core/log.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/utils/curve_sampling.hpp>

using namespace polyvec;

CurveFitter::CurveFitter(double raster_aabb_diagonal, bool final_fit)
	: final_fit(final_fit)
{
}

const polyvec::FitOptions& CurveFitter::get_fit_options() const
{
	if (final_fit)
		return polyvec::final_fit;
	else
		return polyvec::classification_fit;
}

// -------------------------------------------------------------------------
int CurveFitter::add_curve(CurvePrimitive* primitive) {
	assert_break(primitive && primitive->curve);

	if (primitive->curve->get_curve()->get_type() == GLOBFIT_CURVE_BEZIER) {
		assert_break(primitive->fitting_info.fit_midpoints.size() && primitive->fitting_info.dense_tangents.fit_tangents.size());
		//assert_break(primitive->fit_tangents->size() == primitive->fit_tangent_ts->size());
	}

	GlobFitCurveParametrization* curve = primitive->curve.get();
	const int curve_id = (int)primitives.size();
	primitives.emplace_back(primitive);

	curve->uncouple_all_parameters();
	curve->set_params(curve->get_initial_parameters());

	auto& fit_options = get_fit_options();

	if (fit_options.fix_endtangents)
	{
		curve->reduce_degrees_of_freedom(
			DofOptions::FIX_FRONT_TANGENT |
			DofOptions::FIX_BACK_TANGENT);
	}

	if(fit_options.keep_endpoints_on_bisectors)
	{
		curve->reduce_degrees_of_freedom(
			DofOptions::KEEP_FRONT_ON_BISECTOR |
			DofOptions::KEEP_BACK_ON_BISECTOR);
	}

	FittingInfo::Tangents* fit_tangents = &primitive->fitting_info.dense_tangents;	
	double fit_tangent_weight = fit_options.tangents_weight / fit_tangents->fit_tangents.size();
	double fit_point_weight = polyvec::classification_fit.points_weight / primitive->fitting_info.fit_midpoints.size();

	if (primitive->hold_front_tangent)
		curve->reduce_degrees_of_freedom(DofOptions::FIX_FRONT_TANGENT);
	if (primitive->hold_end_tangent)
		curve->reduce_degrees_of_freedom(DofOptions::FIX_BACK_TANGENT);

	// fit points
	for (size_t i = 0; i < primitive->fitting_info.fit_midpoints.size(); ++i) {
		const double t = curve->get_curve()->project(primitive->fitting_info.fit_midpoints.at(i));

		obj_points_t.emplace_back();
		auto& obj = obj_points_t.back();
		obj.set_params(primitive->curve.get(), t, primitive->fitting_info.fit_midpoints.at(i));
		obj.set_weight(primitive->point_weight_multiplier * fit_point_weight);
	}

	// fit tangents
#if POLYVEC_FITTER_EXPORT_TANGENTS_T
	primitive->fit_tangents_t.clear();
#endif
	if (!fit_tangents->fit_tangents.empty())
	{
		obj_tangents_t.emplace_back();
		auto& obj = obj_tangents_t.back();
		obj.set_params(primitive->curve.get());
		obj.set_weight(fit_tangent_weight);
		for (size_t i = 0; i < fit_tangents->fit_tangents.size(); ++i) {
			const double t = fit_tangents->fit_tangent_ts.at(i);
			obj.add_tangent(t, fit_tangents->fit_tangents.at(i));
			obj.set_using_cross_product_formulation(false);

#if POLYVEC_FITTER_EXPORT_TANGENTS_T
			primitive->fit_tangents_t.emplace_back(t);
#endif
		}
	}

	if(primitive->curve->get_curve()->get_type() == GLOBFIT_CURVE_BEZIER) {

		if (fit_options.fairness_formulation == polyvec::FitOptions::CURVATURE_VARIATION)
		{
			obj_curvVar.emplace_back();
			auto& obj = obj_curvVar.back();
			obj.set_params(curve);
			obj.set_weight(fit_options.curvature_variation_weight);
		}

#if REGULARIZE_BEZIER
		auto angleBasedBezier = dynamic_cast<GlobFitBezierAngleBasedParametrization*> (curve);

		if (angleBasedBezier != nullptr) {
			auto params = angleBasedBezier->get_params();

			obj_parameter_bound.emplace_back();
			auto& o1 = obj_parameter_bound.back();
			o1.set_params(angleBasedBezier, 0, 0.4, 0.0);			
			if (params(0) < 0.4)
				params(0) = 0.4;

			obj_parameter_bound.emplace_back();
			auto& o2 = obj_parameter_bound.back();
			o2.set_params(angleBasedBezier, 1, 0.4, 0.0);
			if (params(1) < 0.4)
				params(1) = 0.4;
			angleBasedBezier->set_params(params);
		}

#endif
	}
	else if (primitive->curve->get_curve()->get_type() == GLOBFIT_CURVE_LINE)
	{
		obj_lineRegularization.emplace_back();
		auto& o = obj_lineRegularization.back();
		o.set_params(primitive->curve.get());
	}

	return curve_id;
}

void CurveFitter::limit_endpoint_movement(GlobFitLineParametrization* line, bool limit_back)
{
	auto points = std::static_pointer_cast<GlobFitCurve_Line>(line->get_curve())->get_points();
	auto& coordinate_system = (limit_back ? line->get_back_coordinate_system() : line->get_front_coordinate_system());
	auto hard_limit = (limit_back ? 1 : -1) * (points.col(0) - points.col(1)).dot(coordinate_system.secondary());
	auto soft_limit = 0.1 * hard_limit;

	obj_parameter_bound.emplace_back();
	auto& o = obj_parameter_bound.back();
	o.set_params(line, (limit_back ? 3 : 1), soft_limit, hard_limit);
}

// -------------------------------------------------------------------------
void CurveFitter::make_g0(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next) {
	assert_break(curve_prev);
	assert_break(curve_next);

	auto prev_line = dynamic_cast<GlobFitLineParametrization*> (curve_prev);
	auto prev_bezier = dynamic_cast<GlobFitBezierAngleBasedParametrization*> (curve_prev);

	auto next_line = dynamic_cast<GlobFitLineParametrization*> (curve_next);
	auto next_bezier = dynamic_cast<GlobFitBezierAngleBasedParametrization*> (curve_next);

	if ((prev_line != nullptr || prev_bezier != nullptr) && (next_line != nullptr || next_bezier != nullptr)) {
		//find the internal parameters that we want to couple

		int prevParam1, prevParam2, nextParam1, nextParam2;

		if (prev_line != nullptr) {
			prevParam1 = 2;
			prevParam2 = 3;
		}
		else if (prev_bezier != nullptr) {
			prevParam1 = 4;
			prevParam2 = 5;
		}

		if (next_line != nullptr) {
			nextParam1 = 0;
			nextParam2 = 1;
		}
		else if (next_bezier != nullptr) {
			nextParam1 = 2;
			nextParam2 = 3;
		}

		//if (prev_line && next_line)
		//	return;

		curve_next->couple_parameter(nextParam1, curve_prev, prevParam1);
		curve_next->couple_parameter(nextParam2, curve_prev, prevParam2);

        if (prev_line && next_line && prev_line->is_back_reversed() != next_line->is_front_reversed()) {
            prev_line->fix_parameter(prevParam1, 0.0);
            prev_line->fix_parameter(prevParam2, 0.0);
        }
	}
	else {
		throw std::runtime_error("Cannot make primitives G0");
		///dbg::warning(STR("Using soft G0 constraint."));
		///obj_g0.emplace_back();
		///auto& g0 = obj_g0.back();
		///g0.set_params(curve_prev, curve_next);
		///g0.set_weight(polyvec::spline_fit_g0_weight);
	}
}

// -------------------------------------------------------------------------
void CurveFitter::make_g1(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next) {

	if (!get_fit_options().fix_endtangents) // only do something if end tangents aren't fixed anyway
	{
		auto prev_line = dynamic_cast<GlobFitLineParametrization*> (curve_prev);
		auto prev_bezier = dynamic_cast<GlobFitBezierAngleBasedParametrization*> (curve_prev);

		auto next_line = dynamic_cast<GlobFitLineParametrization*> (curve_next);
		auto next_bezier = dynamic_cast<GlobFitBezierAngleBasedParametrization*> (curve_next);

		if (prev_line && next_line) {
			//We can't do anything here. It is either a corner or the lines should have been merged.
		}
		else if (prev_line && next_bezier) {
			line_g1_constraints.emplace_back(prev_line, GlobFitCurveParametrization::ParameterAddress(next_bezier, 6));
		}
		else if (prev_bezier && next_line) {
			line_g1_constraints.emplace_back(next_line, GlobFitCurveParametrization::ParameterAddress(prev_bezier, 7));
		}
		else if (prev_bezier && next_bezier) {
			curve_prev->couple_parameter(7, curve_next, 6);
		}
		else
		{
			throw std::runtime_error("Cannot make primitives G0");
			/*dbg::warning(STR("Using soft G1 constraint."));
			obj_g1.emplace_back();
			auto& g1 = obj_g1.back();
			g1.set_params(curve_prev, curve_next);
			g1.set_weight(VectorOptions::get()->spline_fit_g1_weight);*/
		}
	}
}

void CurveFitter::make_g2(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next) {
	auto& options = get_fit_options();
	if(!options.add_soft_g2_constraint)
		return;
	
#if !G2_FOR_LINES
	if (std::dynamic_pointer_cast<GlobFitCurve_Line>(curve_prev->get_curve()) || std::dynamic_pointer_cast<GlobFitCurve_Line>(curve_next->get_curve()))
		return;
#endif

	obj_g2.emplace_back();
	auto& g2 = obj_g2.back();
    g2.set_params(curve_prev, 1.0, curve_next, 0.0);
    g2.set_weight(options.g2_weight);
}

void CurveFitter::prescribe_angle(GlobFitLineParametrization * line, double angle)
{
	obj_lineAngle.emplace_back();
	auto& o = obj_lineAngle.back();
	o.set_params(line, angle);
	o.set_weight(1);
}

double CurveFitter::get_angle_residuals() const
{
	double sum = 0;
	for (auto& o : obj_lineAngle)
	{
		Eigen::VectorXd obj;
		Eigen::MatrixXd jac;
		o.compute_objective_and_jacobian(obj, jac);
		sum += obj.squaredNorm();
	}
	return sum;
}

void CurveFitter::increase_angle_weights(double factor)
{
	for (auto& o : obj_lineAngle)
		o.set_weight(o.get_weight() * factor);
}

void CurveFitter::make_overlap_start_to_end(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next) {
    auto bezier_prev = dynamic_cast<GlobFitBezierAngleBasedParametrization*>(curve_prev);
    auto bezier_next = dynamic_cast<GlobFitBezierAngleBasedParametrization*>(curve_next);
    auto line_prev = dynamic_cast<GlobFitLineParametrization*>(curve_prev);
    auto line_next = dynamic_cast<GlobFitLineParametrization*>(curve_next);

    PF_ASSERT(bezier_prev || line_prev);
    PF_ASSERT(bezier_next || line_next);
    if (bezier_prev && bezier_next) {
        curve_prev->couple_parameter(2, curve_next, 4);
        curve_prev->couple_parameter(3, curve_next, 5);
    }
    else if (bezier_prev && line_next) {
        curve_prev->couple_parameter(2, curve_next, 2);
        curve_prev->couple_parameter(3, curve_next, 3);
    }
    else if (line_prev && bezier_next) {
        curve_prev->couple_parameter(0, curve_next, 4);
        curve_prev->couple_parameter(1, curve_next, 5);
    }
    else if (line_prev && line_next) {
        curve_prev->couple_parameter(0, curve_next, 2);
        curve_prev->couple_parameter(1, curve_next, 3);
    }
    else {
        PF_ABORT;
    }
}

void CurveFitter::make_overlap_end_to_start(GlobFitCurveParametrization* curve_prev, GlobFitCurveParametrization* curve_next) {
    auto bezier_prev = dynamic_cast<GlobFitBezierAngleBasedParametrization*>(curve_prev);
    auto bezier_next = dynamic_cast<GlobFitBezierAngleBasedParametrization*>(curve_next);
    auto line_prev = dynamic_cast<GlobFitLineParametrization*>(curve_prev);
    auto line_next = dynamic_cast<GlobFitLineParametrization*>(curve_next);

    PF_ASSERT(bezier_prev || line_prev);
    PF_ASSERT(bezier_next || line_next);
    if (bezier_prev && bezier_next) {
        curve_prev->couple_parameter(4, curve_next, 2);
        curve_prev->couple_parameter(5, curve_next, 3);
    }
    else if (bezier_prev && line_next) {
        curve_prev->couple_parameter(4, curve_next, 0);
        curve_prev->couple_parameter(5, curve_next, 1);
    }
    else if (line_prev && bezier_next) {
        curve_prev->couple_parameter(2, curve_next, 2);
        curve_prev->couple_parameter(3, curve_next, 3);
    }
    else if (line_prev && line_next) {
        curve_prev->couple_parameter(2, curve_next, 0);
        curve_prev->couple_parameter(3, curve_next, 1);
    }
    else {
        PF_ABORT;
    }
}

// -------------------------------------------------------------------------
void CurveFitter::solve() {
	objs.clear();	

	for (size_t i = 0; i < obj_points_t.size(); ++i) 
		objs.emplace_back(&obj_points_t[i]);
	for (index i = 0; i < obj_g0.size(); ++i) 
		objs.emplace_back(&obj_g0[i]);
	for (size_t i = 0; i < obj_g1.size(); ++i) 
		objs.emplace_back(&obj_g1[i]);
	for (size_t i = 0; i < obj_tangents_t.size(); ++i) 
		objs.emplace_back(&obj_tangents_t[i]);
	for (auto& o : obj_curvVar) 
		objs.push_back(&o);
	for (auto& o : obj_parameter_bound)
		objs.push_back(&o);
	for (index i = 0; i < obj_g2.size(); ++i)
		objs.emplace_back(&obj_g2[i]);	
	for (auto& o : obj_lineAngle)
		objs.emplace_back(&o);
	for (auto& o : obj_lineRegularization)
		objs.emplace_back(&o);

	std::vector<GlobFitConstraint*> constraints;
	for (auto& c : line_g1_constraints)
		constraints.push_back(&c);

	//for (size_t i = 0; i < obj_t_on_line.size(); ++i) {
	//	objs.emplace_back(&obj_t_on_line[i]);
	//}
	//for (size_t i = 0; i < obj_fairness.size(); ++i) {
	//	objs.emplace_back(obj_fairness[i]);
	//}
	//for (index i = 0; i < obj_points_project.size(); ++i)
	//	objs.emplace_back(&obj_points_project[i]);	
	//for (index i = 0; i < obj_tangents.size(); ++i)
	//	objs.emplace_back(&obj_tangents[i]);

	const int max_iterations = 50;

	GlobFitter solver;
	std::vector<GlobFitCurveParametrization*> curves;

	for (auto& p : primitives) {
		curves.push_back(p->curve.get());
	}

	solver.set_curves(curves);
	solver.set_objectives(objs);
	solver.set_constraints(constraints);
	solver.setup(max_iterations);

#if POLYVEC_LOG_CURVE_FITTER
	FILE* log = fopen("curve-fitter-log.txt", "w");
	assert_break(log);

	fprintf(log, "Fitting %d curves and %d objectives\n", (int)curves.size(), (int)objs.size());
	fprintf(log, "----------\n");
#else
	FILE* log = nullptr;
#endif

	auto iterations = solver.run_fitter(log, nullptr);
	PF_VERBOSE_F("Solver finished after %i iterations", iterations);

#if POLYVEC_LOG_CURVE_FITTER
	fclose(log);
#endif

	//FILE* fitter_error = fopen("least-squares-objectives.txt", "a");
	//fprintf(fitter_error, "\n\n\n");
	//solver.report_errors(stdout);
	//fclose(fitter_error);

	eval_error();
}

// -------------------------------------------------------------------------
void CurveFitter::solve_with_stiffening() {
    solve();

    // stiffening for line angles
    if (get_fit_options().consider_parallel)
    {
        for (int i = 0; i < 4; ++i)
        {
            if (get_angle_residuals() < PF_EPS)
                break; // there is nothing left to optimize
            increase_angle_weights(10);
            solve();
        }
    }

}

// -------------------------------------------------------------------------
void CurveFitter::eval_error() {
	// reset
	for (size_t i = 0; i < primitives.size(); ++i) {
		CurvePrimitive* prim = primitives[i];
		new (&prim->error) CurveFitError;

		prim->objective_energies.clear();
		prim->objective_energies.resize(GLOBFIT_OBJECTIVE_COUNT);
	}

	for (auto o : objs) {
		auto curve = o->get_curves().front();
		CurvePrimitive* prim = nullptr;

		for (auto& p : primitives)
			if (p->curve.get() == curve) {
				prim = p;
			}

		Eigen::VectorXd obj;
		Eigen::MatrixXd jac;
		o->compute_objective_and_jacobian(obj, jac);
		prim->objective_energies[o->get_type()] += obj.squaredNorm();
	}

	std::vector<double> iset_line_t;
	std::vector<double> iset_curve_t;

	for (size_t i = 0; i < primitives.size(); ++i) {
		CurvePrimitive* prim = primitives[i];
		GlobFitCurveParametrization* curve = prim->curve.get();

		if (!prim->fitting_info.fit_midpoints.empty()) {
			// holy moly please just fix this
			polyfit::mat2x P_fit(2, prim->fitting_info.fit_midpoints.size()), N_fit(2, prim->fitting_info.fit_midpoint_normals.size());

			for (size_t j = 0; j < prim->fitting_info.fit_midpoints.size(); ++j) {
				P_fit.col(j) = prim->fitting_info.fit_midpoints.at(j);
				N_fit.col(j) = prim->fitting_info.fit_midpoint_normals.at(j);
			}

			prim->error.accuracy = polyfit::CurveTracer::measure_accuracy(
				*prim->curve->get_curve().get(),
				P_fit,
				N_fit);
		}

		if (!prim->fitting_info.dense_tangents.fit_tangents.empty()) {
			polyfit::CurveTracer::measure_curvature_radius(*prim, prim->error.curvature);
		}
	}
}

void CurveFitter::reset() {
    primitives.clear();
    objs.clear();
}