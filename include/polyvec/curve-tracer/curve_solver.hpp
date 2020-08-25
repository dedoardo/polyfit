/*
    The Fitter is responsible for managing bookkeeping of the curves, their objectives
    and correctly composing the sparse jacobian and residual vector used by the 
    non-linear least squares solver.
*/
#pragma once

// polyvec
#include <polyvec/curve-tracer/curve_objective.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>
#include <polyvec/curve-tracer/curve_constraint.hpp>

// libc++
#include <vector>
#include <memory>

// Eigen
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/LevenbergMarquardt> // SparseFunctor

NAMESPACE_BEGIN(polyvec)

// Forward declarations
class GlobFitObjective;
class GlobFitCurve;

struct GlobFitter : public Eigen::SparseFunctor<double, int> {
public:
	GlobFitter() : Eigen::SparseFunctor<double, int>(0, 0) { }
    void set_curves ( const std::vector<GlobFitCurveParametrization*>& );
    void set_objectives ( const std::vector<GlobFitObjective*>& );

	//Adds hard constraints to the optimization problem. The following requirements must be
	//met for constraints  target = f(source):
	//  * target cannot be coupled to other parameters
	//  * target cannot be fixed
	//  * none of source can be a target of a constraint	
	void set_constraints(const std::vector<GlobFitConstraint*>&);

	//Prepares the optimization process
    void setup ( const int max_iterations );

	//returns the number of iterations
    int run_fitter ( FILE* log_file, std::function<void ( int ) > callback );
	
	void report_errors( FILE* error_file );

    int n_parameters();
    int n_equations();
    bool is_setup();

    std::vector<GlobFitCurveParametrization*> get_curves();
    std::vector<GlobFitObjective*> get_objectives();    

private:

	void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::SparseMatrix<double>& dobj_dcurveparams); // + matrix 

	//extracts the current set of variables from the current state of the curves
	Eigen::VectorXd get_variables();
	
	//updates the curves with the given set of variables
	void set_variables(const Eigen::VectorXd&);

	int nextParameter = 0;
    std::vector<GlobFitCurveParametrization*> _curves;
	std::map< GlobFitCurveParametrization*, size_t> curve_to_index;
	std::vector<GlobFitObjective*> _objectives;
	std::vector<GlobFitConstraint*> _constraints;
	std::vector<Eigen::VectorXd> constraint_gradients;

    bool _is_setup = false;

	//id of the first equation of the i-th objective
    std::vector<int> _xadj_equations = std::vector<int>();        
	std::vector<std::vector<int>> variables_in_objective;
    int _max_iterations;

    Eigen::VectorXd _cached_params;
    Eigen::VectorXd _cached_objective;
    Eigen::SparseMatrix<double> _cached_jacobian;

	void check_derivatives();

    // ============== API needed by Eigen Levenberg-Marquart
public:

    int inputs();
    int values();
    int operator() ( const Eigen::VectorXd& x, Eigen::VectorXd& fvec );
    int df ( const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& fjac );
};

NAMESPACE_END (polyvec)