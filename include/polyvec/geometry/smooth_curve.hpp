// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>

// libc++
#include <memory>

NAMESPACE_BEGIN(polyvec)

// Forward declarations
class BezierCurve;
class GlobFitCurve;

NAMESPACE_BEGIN(SmoothCurveUtil)
void reverse_control_points(
    GlobFitCurve* curve
);

std::shared_ptr<GlobFitCurveParametrization> create_parametrization (
    std::shared_ptr<GlobFitCurve>& curve
);

// Returns norm(\int(|dk/dt|^2))
double evaluate_total_curvature_variation (
    BezierCurve& curve
);

// 
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
);

//
void tangent_eval (
    const Eigen::Vector2d& drdt,
    Eigen::Vector2d& tangent
);

//
void tangent_derivatives (
    const Eigen::Vector2d& drdt,
    Eigen::Vector2d& tangent,
    Eigen::Matrix2d& dtangent_ddrdt
);

//
void smoothabs_eval( const double in, double &abs);

//
void smoothabs_derivatives( const double in, double &abs, double &dabsdin);

//
void curvature_eval (
    const Eigen::Vector2d& drdt,
    const Eigen::Vector2d& drdtdt,
    double& curvaure
);

//
void curvature_derivatives (
    const Eigen::Vector2d& drdt,
    const Eigen::Vector2d& drdtdt,
    double& curvaure,
    Eigen::RowVector2d& dcurvaure_ddrdt,
    Eigen::RowVector2d& dcurvaure_ddrdtdt
);

NAMESPACE_END(SmoothCurveUtil)
NAMESPACE_END(polyvec)