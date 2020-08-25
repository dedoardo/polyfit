#include <polyvec/curve-tracer/curve_objective.hpp>

using namespace polyvec;

std::string polyvec::globfitobjectivetype_as_string(GlobFitObjectiveType in) {
#define PROCESS(XX) case XX: return #XX

	switch (in) {
		PROCESS(GLOBFIT_OBJECTIVE_FIT_POINTS);
		PROCESS(GLOBFIT_OBJECTIVE_FIT_TANGENTS);
		PROCESS(GLOBFIT_OBJECTIVE_FIX_POSITION);
		PROCESS(GLOBFIT_OBJECTIVE_FIX_TANGENT);
		PROCESS(GLOBFIT_OBJECTIVE_SAME_POSITION);
		PROCESS(GLOBFIT_OBJECTIVE_SAME_TANGENT);
		PROCESS(GLOBFIT_OBJECTIVE_SAME_CURVATURE);
		PROCESS(GLOBFIT_OBJECTIVE_BEZIER_FAIRNESS);
		PROCESS(GLOBFIT_OBJECTIVE_POSITION_ON_CURVE_LIE_ON_LINE);
		PROCESS(GLOBFIT_OBJECTIVE_SCALE_AND_TRANSITION_ONLY_REGULARIZER);
		PROCESS(GLOBFIT_OBJECTIVE_CURVATURE_VARIATION);
		PROCESS(GLOBFIT_OBJECTIVE_PARAMETER_BOUND);
		PROCESS(GLOBFIT_OBJECTIVE_LINE_REGULARIZATION);
		PROCESS(GLOBFIT_OBJECTIVE_LINE_ANGLE);
	default:
		assert_break(0);
		return "";
	}

#undef PROCESS
}

// =============================================================
//                           BASE CLASS
// =============================================================    

void GlobFitObjective::set_weight(const double weight) {
	assert_break(weight >= 0);
	_weight = weight;
	_sqrt_weight = sqrt(weight);
}

double GlobFitObjective::get_weight() const {
	return _weight;
}

double GlobFitObjective::get_sqrt_weight() const {
	return _sqrt_weight;
}

int GlobFitObjective::_n_sum_curve_params() const {
	int ans = 0;

	for (GlobFitCurveParametrization* c : _curves) {
		ans += c->n_params();
	}

	return ans;
}


bool GlobFitObjective::_is_all_data_provided() const {
	bool answer = true;

	answer = answer && (_curves.size() > 0);
	answer = answer && (_weight >= 0);

	return answer;
}