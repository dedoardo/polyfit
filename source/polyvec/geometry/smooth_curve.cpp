// polyvec
#include <polyvec/geometry/smooth_curve.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_objectives.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curve_line.hpp>

// libc++
#include <memory>

using namespace polyvec;
using namespace polyfit;
using namespace std;

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(SmoothCurveUtil)

// Evaluates the norm curvature variation using the solver's objective object, 
// this is ugly but gets the job done.
double evaluate_total_curvature_variation(
    BezierCurve& curve
) {
    std::shared_ptr<BezierCurve> curve_ptr = make_shared<BezierCurve>(curve);
    GlobFitBezierAngleBasedParametrization curve_params(curve_ptr);
    GlobFitObjective_CurvatureVariation k_variation;
    k_variation.set_params(&curve_params);

    vecXd k_variation_objective;
    matXd k_variation_jacobian;
    k_variation.compute_objective_and_jacobian(k_variation_objective, k_variation_jacobian);
    return k_variation_objective.norm();
}

Eigen::VectorXd projection_derivatives (
    const int n_params,
    const double time,
    const double tend,
    const Eigen::Vector2d& outside_point,
    const Eigen::Vector2d& rr,
    const Eigen::Vector2d& drdt,
    const Eigen::Vector2d& drdtdt,
    const Eigen::Matrix2Xd& drdparams,
    const Eigen::Matrix2Xd& drdtdparams
) {
    const double tol = 1e-10;

    //
    // find z.drdt
    //
    const Eigen::Vector2d zz = rr - outside_point;
    const double zdotdrdt = zz.dot ( drdt );
    const bool is_bdry_minimizer = std::abs ( zdotdrdt ) > tol;

    //
    // Find the dt/dparamstilde
    //
    Eigen::MatrixXd dtdparamstilde ( 1, n_params );

    if ( is_bdry_minimizer && ( ( time <= tol ) || ( time + tol >= tend ) ) ) {
        dtdparamstilde.setZero(); // Assuming dtend and dtbeg / dparams = 0
    } else {
        double denom = drdt.dot ( drdt ) + zz.dot ( drdtdt );

        if ( std::abs ( denom ) < tol ) {
            denom = ( denom < 0. ? -tol : tol );
        }

        dtdparamstilde = - ( drdt.transpose() * drdparams + zz.transpose() * drdtdparams ) / denom;
    }

    // theoreticlly it is a rowvector, but all single row should become col
    return dtdparamstilde.transpose();
}

#define SMOOTH_ABS_EPSILON 1e-3
#define SMOOTH_ABS_ALPHA 1
#define NO_DIV_BY_ZERO_EPS 1e-8

void tangent_eval (
    const Eigen::Vector2d& drdt,
    Eigen::Vector2d& tangent
) {
    // const double tol = 1e-15;
    const double eps = SMOOTH_ABS_EPSILON;

    // double dposdt_norm2 = std::max ( drdt.squaredNorm() , eps );
    double dposdt_norm2 = drdt.squaredNorm() + eps ;
    double dposdt_norm = sqrt ( dposdt_norm2 );
    tangent = drdt / dposdt_norm;
}

void tangent_derivatives (
    const Eigen::Vector2d& drdt,
    Eigen::Vector2d& tangent,
    Eigen::Matrix2d& dtangent_ddrdt
) {
    // const double tol = 1e-15;
    const double eps = NO_DIV_BY_ZERO_EPS;

    // double dposdt_norm2 = std::max ( drdt.squaredNorm() , eps );
    double dposdt_norm2 = drdt.squaredNorm() + eps ;
    double dposdt_norm = sqrt ( dposdt_norm2 );
    tangent = drdt / dposdt_norm;
    dtangent_ddrdt =
        Eigen::Matrix2d::Identity() / dposdt_norm
        - 1. / ( dposdt_norm * dposdt_norm2 ) * drdt * drdt.transpose();
}

void smoothabs_eval ( const double in, double& abs ) {
    // VERSION 1
    const double eps = SMOOTH_ABS_EPSILON;
    abs = sqrt ( in*in+eps );
    // VERSION 2 -- overflows
    // const double alpha = SMOOTH_ABS_ALPHA;
    //abs = 1. / alpha * ( log(1. + exp(alpha*in)) + log(1. + exp(-alpha*in)) ) ; 
    // VERSION 3
    //abs = in*in;
}

void smoothabs_derivatives ( const double in, double& abs, double& dabsdin ) {
    // VERSION 1
    smoothabs_eval ( in, abs );
    dabsdin = in / abs;
    // VERSION 2 -- overflows
    // const double alpha = SMOOTH_ABS_ALPHA;
    // smoothabs_eval(in, abs);
    //dabsdin = 1. / (1. + exp( - alpha * in) ) - 1. / (1. + exp( alpha * in));
    // VERSION 3
    //abs = in * in;
    //dabsdin = 2*in;
}

void curvature_eval (
    const Eigen::Vector2d& drdt,
    const Eigen::Vector2d& drdtdt,
    double& curvature
) {
    const double eps = NO_DIV_BY_ZERO_EPS;

    double dposdt_norm2 = drdt.squaredNorm() + eps;
    double dposdt_norm = sqrt ( dposdt_norm2 );
    const double denomi =  ( dposdt_norm*dposdt_norm2 );
    const double nomi =  drdt.x() *drdtdt.y()-drdt.y() *drdtdt.x() ;
    curvature = nomi / denomi;
}

void curvature_derivatives (
    const Eigen::Vector2d& drdt,
    const Eigen::Vector2d& drdtdt,
    double& curvature,
    Eigen::RowVector2d& dcurvaure_ddrdt,
    Eigen::RowVector2d& dcurvaure_ddrdtdt
) {
    const double eps = NO_DIV_BY_ZERO_EPS;

    double dposdt_norm2 = drdt.squaredNorm() + eps;
    double dposdt_norm = sqrt ( dposdt_norm2 );
    //
    const double denomi =  ( dposdt_norm*dposdt_norm2 );
    const double nomi =  drdt.x() *drdtdt.y()-drdt.y() *drdtdt.x() ;
    //
    Eigen::RowVector2d dnomi_ddrdt ( drdtdt.y(), - drdtdt.x() );
    Eigen::RowVector2d dnomi_ddrdtdt ( -drdt.y(),  drdt.x() );
    //
    Eigen::RowVector2d dinvdenomi_ddrdt = -3. / ( dposdt_norm2*dposdt_norm2*dposdt_norm ) * drdt.transpose() ;

    curvature = nomi / denomi;
    dcurvaure_ddrdt = dnomi_ddrdt / denomi + nomi * dinvdenomi_ddrdt;
    dcurvaure_ddrdtdt = dnomi_ddrdtdt / denomi;
}

#undef SMOOTH_ABS_EPSILON
#undef SMOOTH_ABS_ALPHA
#undef NO_DIV_BY_ZERO_EPS

NAMESPACE_END(SmoothCurveUtil)
NAMESPACE_END(polyvec)