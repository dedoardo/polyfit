// Header
#include <polyvec/curve-tracer/curve_solver.hpp>
#include <polyvec/core/log.hpp>

// Eigen
#include <Eigen/Core>
#include <Eigen/Jacobi>

// Polyvec
#include <polyvec/api.hpp>
#include <polyvec/utils/union-find.hpp>

#include <nse/MathematicaFormatter.h>

#include <fstream>
#include <iostream>
#include <iomanip>

#define PRINT_JACOBIAN 0
#define PRINT_ITERATION_STATS 0

#define ENABLE_DERIVATIVE_CHECK 0

namespace polyvec {
// ============= EIGEN INTERFACE FOR MARQUART_LEVENBERG
int
GlobFitter::inputs() {
    return n_parameters();
}

int
GlobFitter::values() {
    return n_equations();
}

int
GlobFitter::operator() ( const Eigen::VectorXd& x, Eigen::VectorXd& fvec ) {
    const double tol = 1e-10;

    if ( ( _cached_params.size() ==0 ) || ( ( _cached_params - x ).norm() > tol ) ) {
        set_variables ( x );
        _cached_params = x;
        compute_objective_and_jacobian ( _cached_objective, _cached_jacobian );
    }

    fvec = _cached_objective;
    return 0;
}

int
GlobFitter::df ( const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& fjac ) {
    const double tol = 1e-10;

    if ( ( _cached_params.size() ==0 ) || ( ( _cached_params - x ).norm() > tol ) ) {
        set_variables ( x );
        _cached_params = x;
        compute_objective_and_jacobian ( _cached_objective, _cached_jacobian );
    }

    fjac = _cached_jacobian;
    return 0;
}

// Separate file as it takes ages to compile.
int
GlobFitter::run_fitter ( FILE* log_file, std::function<void ( int ) > callback ) {
    assert_break ( _is_setup ); //

    Eigen::VectorXd x0 = get_variables();
	set_variables(x0); //accomodate for fixed parameters
    Eigen::VectorXd x = x0;

	if (n_parameters() == 0 || n_equations() == 0)
		return 0;

    // Init Eigen
    Eigen::LevenbergMarquardt<GlobFitter> lm ( *this );
	lm.setGtol(1e-5);
	lm.setXtol(1e-5);
	lm.setFtol(1e-5);
    using Eigen::LevenbergMarquardtSpace::Running;
    using Eigen::LevenbergMarquardtSpace::Status;	

    Status status = lm.minimizeInit ( x );

    int iter_id = 0;

    if ( log_file ) {
        fprintf ( log_file, " ==== MARQURT-LEVENBERG for Bezier Fitting ===== \n" );
		fprintf(log_file, " %i equations, %i degrees of freedom\n", n_equations(), n_parameters());
        fprintf ( log_file, " %10s %20s %20s \n", "ITER", "OBJ", "GRAD" );
        fflush ( log_file );
    }

#if PRINT_ITERATION_STATS
	printf("------------------------\n");
#endif

    bool keep_running = true;

	while (keep_running && iter_id < _max_iterations) {
        status = lm.minimizeOneStep ( x );

#if ENABLE_DERIVATIVE_CHECK
		check_derivatives();
#endif

#if PRINT_ITERATION_STATS
		std::cout << lm.gnorm() << "  " << lm.fnorm() * lm.fnorm() << " " << status << std::endl;
#endif
        if ( log_file ) {
            if ( ( iter_id % 5 == 0 ) || ( status != Running ) ) {
                fprintf ( log_file, " %10d %20.8e %20.8e \n",
                         iter_id,
                         _cached_objective.norm(),
                         ( _cached_jacobian.transpose() * _cached_objective ).norm() );
                fflush ( log_file );
            }
        }

        ++iter_id;

        if ( callback ) {
            callback ( iter_id );
        }

        if ( status != Running ) {
            keep_running = false;
        } else if ( iter_id >= _max_iterations ) {
			if (log_file) {
				fprintf(log_file, " %10d %20.8e %20.8e \n",
					iter_id,
					_cached_objective.norm(),
					(_cached_jacobian.transpose() * _cached_objective).norm());
				fflush(log_file);
			}
        }

    }

    if ( log_file ) {
        switch ( status ) {
#define TAKECAREOF(XX) \
  case(Eigen::LevenbergMarquardtSpace::XX): fprintf(log_file, "FINISHED, Status: %s \n", #XX); break

            TAKECAREOF ( NotStarted );
            TAKECAREOF ( Running );
            TAKECAREOF ( ImproperInputParameters );
            TAKECAREOF ( RelativeReductionTooSmall );
            TAKECAREOF ( RelativeErrorTooSmall );
            TAKECAREOF ( RelativeErrorAndReductionTooSmall );
            TAKECAREOF ( CosinusTooSmall );
            TAKECAREOF ( TooManyFunctionEvaluation );
            TAKECAREOF ( FtolTooSmall );
            TAKECAREOF ( XtolTooSmall );
            TAKECAREOF ( GtolTooSmall );
            TAKECAREOF ( UserAsked );

#undef TAKECAREOF
        }

        // Print the objectives
        for(int i = 0 ; i < (int)_objectives.size() ; ++i)
        {
            fprintf ( log_file, "%s: %g \n",
             globfitobjectivetype_as_string(_objectives[i]->get_type()).c_str(),
             _cached_objective(i) );
        }

        fflush ( log_file );
    }

	set_variables(x);

#if PRINT_ITERATION_STATS
	std::cout << "Energy per objective type:" << std::endl;
	std::map<GlobFitObjectiveType, double> energy_per_type;
	for (int i = 0; i < _objectives.size(); ++i)
	{		
		energy_per_type[_objectives[i]->get_type()] += _cached_objective.segment(_xadj_equations[i], _xadj_equations[i + 1] - _xadj_equations[i]).squaredNorm();
	}
	for (auto& entry : energy_per_type)
	{
		std::cout << globfitobjectivetype_as_string(entry.first) << ": " << entry.second << std::endl;
	}


#endif

	return iter_id;
}

void GlobFitter::set_curves(const std::vector<GlobFitCurveParametrization*>& curves_in) {
	_curves = curves_in;
	_is_setup = false;
}

void GlobFitter::set_objectives(const std::vector<GlobFitObjective*>& objectives_in) {
	_objectives = objectives_in;
	_is_setup = false;
}

void GlobFitter::set_constraints(const std::vector<GlobFitConstraint*>& c) {
	_constraints = c;
	constraint_gradients.resize(c.size());
	_is_setup = false;
}

void GlobFitter::setup(const int max_iterations) {
	_max_iterations = max_iterations;
	

	curve_to_index.clear();
	for (int i = 0; i < _curves.size(); ++i)
		curve_to_index[_curves[i]] = i;

	//Assign a temporary index to each parameter
	nextParameter = 0;
	for (int i = 0; i < (int)_curves.size(); ++i) {
		auto curve = _curves[i];
		auto& info = curve->get_parameter_info();
		for (auto& entry : info) {
			entry.variable_id = nextParameter++;
			entry.constraint_id = -1;
			entry.fixedValuePropagated = std::numeric_limits<double>::quiet_NaN();
		}
	}

	//Assign constraint information
	for(int i = 0; i < _constraints.size(); ++i)
	{
		auto& target = _constraints[i]->get_target_param();
		auto& param = target.curve->get_parameter_info()[target.internal_parameter];		
		assert_break_msg(param.constraint_id == -1, "Only one constraint is permitted per parameter.");
		assert_break_msg(std::isnan(param.fixedValue), "A fixed parameter cannot have an additional constraint.");
		param.constraint_id = i;		
	}

	//Check if there are chained constraints
	for (auto& c : _constraints)
	{
		for (auto& source : c->get_source_params())
		{
			assert_break_msg(source.curve->get_parameter_info()[source.internal_parameter].constraint_id == -1, "Chained constraints are not supported.");			
		}
	}

	//Find the coupled parameters with a union-find
	struct EntryData
	{
		int variable_id;
		double fixedValue = std::numeric_limits<double>::quiet_NaN();
		bool is_constrained = false;
	};
	
	UnionFind<EntryData> uf(nextParameter);

	for (int i = 0; i < (int)_curves.size(); ++i) {
		auto curve = _curves[i];
		auto& info = curve->get_parameter_info();
		for (auto& entry : info) {
			uf[entry.variable_id].is_constrained = entry.constraint_id != -1;
			for (auto& c : entry.coupled_parameters) {
				auto& param2 = c.curve->get_parameter_info()[c.internal_parameter];
				auto coupledParamId = param2.variable_id;
				uf.merge(entry.variable_id, coupledParamId);
				assert_break_msg(entry.constraint_id == -1 && param2.constraint_id == -1, "Constraining coupled parameters is not supported.");
			}
		}
	}

	//set fixed values to parents
	for (int i = 0; i < (int)_curves.size(); ++i) {
		auto curve = _curves[i];
		auto& info = curve->get_parameter_info();
		for (auto& entry : info) {
			if (!std::isnan(entry.fixedValue)) {
				auto& parentFixed = uf[uf.getRepresentative(entry.variable_id)].fixedValue;
				assert_break(std::isnan(parentFixed) || std::abs(parentFixed - entry.fixedValue) < 0.0001);
				parentFixed = entry.fixedValue;
			}
		}
	}

	//assign actual parameter indices to parents
	nextParameter = 0;
	for (int i = 0; i < uf.size(); ++i)
		if (uf[i].parent == i && std::isnan(uf[i].fixedValue) && !uf[i].is_constrained)
			uf[i].variable_id = nextParameter++;

	//assign parameter indices to children
	for (int i = 0; i < (int)_curves.size(); ++i) {
		auto curve = _curves[i];
		auto& info = curve->get_parameter_info();
		for (auto& entry : info) {
			auto& parent = uf[uf.getRepresentative(entry.variable_id)];
			// reset the variable id
			entry.variable_id = -1;
			if (!std::isnan(parent.fixedValue))
				// this parameter is fixed
				entry.fixedValuePropagated = parent.fixedValue;
			else if (entry.constraint_id != -1)
				// this parameter is constrained
				;
			else
				// this parameter has a variable
				entry.variable_id = parent.variable_id;				
		}
	}
		

	//
	// Equation bookkeeping
	//
	_xadj_equations.resize(_objectives.size() + 1);
	_xadj_equations[0] = 0;
	variables_in_objective.clear();
	variables_in_objective.resize(_objectives.size());
	for (int i = 0; i < (int)_objectives.size(); ++i) {
		_xadj_equations[i + 1] = _xadj_equations[i] + _objectives[i]->n_equations();

		auto& o = _objectives[i];
		for (auto c : o->get_curves()) {
			for (auto& p : c->get_parameter_info()) {
				if (p.variable_id != -1)
					variables_in_objective[i].push_back(p.variable_id);
				else if (p.constraint_id != -1)
					for (auto& constraint_source : _constraints[p.constraint_id]->get_source_params())
					{
						auto& source_param = constraint_source.curve->get_parameter_info()[constraint_source.internal_parameter];
						if (source_param.variable_id != -1)
							variables_in_objective[i].push_back(source_param.variable_id);
					}
			}
		}
	}


	//
	// Initial params
	//
	_is_setup = true;

	//check_derivatives();
}

int GlobFitter::n_parameters() {
	assert_break(_is_setup);
	return nextParameter;
}

int GlobFitter::n_equations() {
	assert_break(_is_setup);
	return _xadj_equations.back();
}

bool GlobFitter::is_setup() {
	return _is_setup;
}


std::vector<GlobFitCurveParametrization*> GlobFitter::get_curves() {//
	return _curves;
}


std::vector<GlobFitObjective*> GlobFitter::get_objectives() {//
	return _objectives;
}

void GlobFitter::compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::SparseMatrix<double>& jac) {//
	assert_break(_is_setup);

	obj.setZero(n_equations());	
	if (jac.isCompressed()) {
		jac.resize(n_equations(), n_parameters());
		//reserve space for jacobian, this is only done once
		Eigen::VectorXi entriesPerColumn(n_parameters());
		entriesPerColumn.setZero();
		for (int i = 0; i < _objectives.size(); ++i)
		{
			for (auto variable : variables_in_objective[i])
				entriesPerColumn(variable) += _objectives[i]->n_equations();
		}		
		_cached_jacobian.reserve(entriesPerColumn);
	}

	Eigen::VectorXd sub_obj;
	Eigen::MatrixXd sub_jac;

	for (int objid = 0; objid < (int)_objectives.size(); ++objid) {
		const int eq_beg = _xadj_equations[objid];
		const int eq_len = _xadj_equations[objid + 1] - _xadj_equations[objid];

		// Eval objective
		_objectives[objid]->compute_objective_and_jacobian(sub_obj, sub_jac);

		// Set obj
		obj.segment(eq_beg, eq_len) = _objectives[objid]->get_sqrt_weight() * sub_obj;

		// Set Jacobian
		// Initialize to zero
		for (auto variable : variables_in_objective[objid])
		{
			for (int i = _xadj_equations[objid]; i < _xadj_equations[objid + 1]; ++i)
				jac.coeffRef(i, variable) = 0;
		}
		// Add the partial gradient contributions
		int subParamId = 0;		
		for (auto curve : _objectives[objid]->get_curves()) {
			auto& params = curve->get_parameter_info();
			for (auto& p : params) {
				if (p.variable_id != -1) {
					//the parameter is a variable
					for(int eq = eq_beg; eq != _xadj_equations[objid + 1]; ++eq)
						jac.coeffRef(eq, p.variable_id) += _objectives[objid]->get_sqrt_weight() * sub_jac.coeff(eq - eq_beg, subParamId);
				}
				else if (p.constraint_id != -1) {
					//the parameter is derived from other variables
					auto& sources = _constraints[p.constraint_id]->get_source_params();
					auto& dpdsources = constraint_gradients[p.constraint_id];
					for (int iSource = 0; iSource < sources.size(); ++iSource)
					{
						auto& source = sources[iSource];
						auto& psource = source.curve->get_parameter_info()[source.internal_parameter];
						if (psource.variable_id != -1)
						{
							for (int eq = eq_beg; eq != _xadj_equations[objid + 1]; ++eq)
								//apply chain rule  d obj / d source = d obj / d target * d target / d source
								jac.coeffRef(eq, psource.variable_id) += _objectives[objid]->get_sqrt_weight() * sub_jac.coeff(eq - eq_beg, subParamId) * dpdsources(iSource);
						}
					}
				}
				++subParamId;
			}
		} // end of curves_outline
	} // end of objectives

#if PRINT_JACOBIAN
	{
		std::ofstream file("jacobian.txt", std::ios::trunc | std::ios::out);
		file << "Variables: \n";
		file << get_variables().transpose() << "\n\n";
		file << "Obj                                                  | Jacobian: \n";
		for (int i_obj = 0; i_obj < _objectives.size(); ++i_obj)
		{
			auto row_start = _xadj_equations[i_obj];
			auto row_end = _xadj_equations[i_obj + 1];
			for (int row = row_start; row < row_end; ++row)
			{
				std::string type = (row == row_start ? globfitobjectivetype_as_string(_objectives[i_obj]->get_type()) : "");
				file << std::setw(40) << type << std::setw(12) << obj(row) << " | ";
				for (int col = 0; col < jac.cols(); ++col)
					file << std::setw(14) << jac.coeff(row, col);
				file << '\n';
			}			
		}
		file << '\n';

		file << "Curve parameters (variable id, fixed value, fixed value propagated): \n";
		for (auto& c : get_curves()) {
			for (auto& p : c->get_parameter_info()) 
				file << "(" << p.variable_id << ", " << p.fixedValue << ", " << p.fixedValuePropagated << ") ";
			file << '\n';
		}
		file << '\n';

		//file << "Sparse Jacobian: \n";
		//file << nse::util::FormatMathematica(jac);
		file.close();
	}
#endif
}

void GlobFitter::check_derivatives() {
	auto x = get_variables();
	Eigen::VectorXd obj;
	Eigen::SparseMatrix<double> jac;
	jac.resize(n_equations(), n_parameters());
	_cached_jacobian.resize(n_equations(), n_parameters());
	const double eps = 0.001;

	Eigen::VectorXd objPlus, objMinus;

	bool error = false;
	compute_objective_and_jacobian(obj, jac);	
	for (int ip = 0; ip < n_parameters(); ++ip)
	{
		auto xCopy = x;
		xCopy(ip) += eps;
		set_variables(xCopy);
		compute_objective_and_jacobian(objPlus, _cached_jacobian);

		xCopy(ip) -= 2 * eps;
		set_variables(xCopy);
		compute_objective_and_jacobian(objMinus, _cached_jacobian);

		Eigen::VectorXd expectedDeriv = (objPlus - objMinus) / (2 * eps);
		for (int i = 0; i < n_equations(); ++i)
		{
			auto error = jac.coeff(i, ip) - expectedDeriv(i);
			auto relativeError = error / std::max(std::abs(obj(i)), std::abs(jac.coeff(i, ip)));

			if (std::abs(relativeError) > 0.001 && std::abs(error) > 1e-10)
			{
				error = true;
				std::cout << "Wrong derivative at parameter " << ip << ", equation " << i << ". Numerical derivative: " << expectedDeriv(i) << 
					", computed derivative: " << jac.coeff(i, ip) << ", objective value: " << obj(i) << ", parameter value: " << x(ip) << std::endl;
				int obj_id = 0;
				while (_xadj_equations[obj_id] <= i)
					++obj_id;
				--obj_id;
				std::cout << "Equation belongs to objective " << obj_id << " (type " << globfitobjectivetype_as_string(_objectives[obj_id]->get_type()) << ", local equation " << i - _xadj_equations[obj_id] << ")" << std::endl;
				std::cout << "Parameter belongs to the following curves: ";
				for (int ic = 0; ic < _curves.size(); ++ic)
				{
					auto curve = _curves[ic];
					auto& curve_params = curve->get_parameter_info();
					for (int icp = 0; icp < curve_params.size(); ++icp)
					{
						if (curve_params[icp].variable_id == ip)
							std::cout << "curve " << ic << " (type " << typeid(*curve).name() << ", local parameter " << icp << ") ";
					}
				}
				std::cout << std::endl;
			}
		}
	}
	if (!error)
		std::cout << "Derivatives seem good." << std::endl;
	set_variables(x);
}

Eigen::VectorXd GlobFitter::get_variables() {//
	assert_break(_is_setup);
	Eigen::VectorXd ans(n_parameters());

	for(auto curve : _curves) {
		auto params = curve->get_params();
		auto& pInfo = curve->get_parameter_info();
		for (int i = 0; i < pInfo.size(); ++i) {
			if (pInfo[i].variable_id != -1)
				ans(pInfo[i].variable_id) = params(i);
		}
	}

	return ans;
}

void GlobFitter::set_variables(const Eigen::VectorXd& params_in) { //
	std::vector<Eigen::VectorXd> curve_params(_curves.size());

	//Fill the variables into curve parameters, consider fixed and coupled parameters
	for (int iCurve = 0; iCurve < _curves.size(); ++iCurve) {
		auto curve = _curves[iCurve];
		auto& pInfo = curve->get_parameter_info();
		auto& params = curve_params[iCurve];
		params.resize(pInfo.size());
		for (int i = 0; i < pInfo.size(); ++i) {
			if (pInfo[i].variable_id != -1)
				params(i) = params_in(pInfo[i].variable_id);
			else if (!std::isnan(pInfo[i].fixedValuePropagated))
				params(i) = pInfo[i].fixedValuePropagated;				
		}		
		//at this point, the constrained parameters are not yet filled in
		//set them anyway, so that the calculation of constraints can access them
		curve->set_params(params);
	}
	
	std::vector<bool> curve_has_constrained_parameters(_curves.size(), false);
	//Derive parameters from constraints
	for (int iConstraint = 0; iConstraint < _constraints.size(); ++iConstraint) {
		auto c = _constraints[iConstraint];
		auto& target = c->get_target_param();	
		auto curve_id = curve_to_index.at(target.curve);
		//We calculate the value of this constraint and assume the following:
		//  * all source parameters are set on the underlying parametrizations (no chained constraints)
		std::tie(curve_params[curve_id](target.internal_parameter), constraint_gradients[iConstraint]) = c->f();
		curve_has_constrained_parameters[curve_id] = true;
	}

	//Set the new curve parameters that are changed by constraints
	for (int iCurve = 0; iCurve < _curves.size(); ++iCurve) {
		if (!curve_has_constrained_parameters[iCurve])
			continue;
		auto curve = _curves[iCurve];
		curve->set_params(curve_params[iCurve]);
	}
}

void GlobFitter::report_errors(FILE* error_file) {
	for (size_t i = 0; i < _curves.size(); ++i) {
		fprintf(error_file, "curve id %d\n", (int)i);
		fprintf(error_file, "----------------------------------------------\n");

		for (size_t j = 0; j < _objectives.size(); ++j) {
			GlobFitObjective* obj = _objectives[j];					

			Eigen::VectorXd res;
			Eigen::MatrixXd jac;
			obj->compute_objective_and_jacobian(res, jac);
			auto type = globfitobjectivetype_as_string(obj->get_type());
			fprintf(error_file, "objective: %s\n", type);
			
			fprintf(error_file, "values: ");
			for (Eigen::Index k = 0; k < res.size(); ++k) {
				fprintf(error_file, "%f ", res(k));
			}
			fprintf(error_file, "\n");
		}
	}
}
}