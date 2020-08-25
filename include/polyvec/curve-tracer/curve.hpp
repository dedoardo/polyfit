/*
    Abstract interface for all curve primitives which permits to evaluate the functional
    and its derivatives. And a couple of other miscellaneous helper functions, some
    of which shouldn't be here.

    Implementations of the interface are not used directly in the optimization
    but rather wrapped in a CurveParameterization view which rearranges the curve parameters
*/
#pragma once

// polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/geom.hpp>

NAMESPACE_BEGIN (polyvec)

// Existing curve primitive implementation
enum GlobFitCurveType {
    GLOBFIT_CURVE_LINE = 0,
    GLOBFIT_CURVE_ARC,
    GLOBFIT_CURVE_BEZIER,
};
const char* curve_type_to_string(const GlobFitCurveType type);

// =============================================================
//                           BASE CLASS
// =============================================================

//
// This is the curve interface used by the globfitter.
// Each curve that is used by the glob fitter
// has to be wrapped inside one of these objects.
//
class GlobFitCurve_Line;
class GlobFitCurve_Arc;
class GlobFitCurve_Bezier;
class GlobFitCurve_BezierScaleAndTranslationOnly;

class GlobFitCurve {
  public:
    static double t_end;

    virtual GlobFitCurveType get_type() const = 0;	
	virtual int n_params() const = 0;

    virtual Eigen::Vector2d pos ( const double t ) const = 0;
    virtual Eigen::Vector2d dposdt ( const double t ) const = 0;
    virtual Eigen::Vector2d dposdtdt ( const double t ) const = 0;
	virtual Eigen::Vector2d dposdtdtdt ( const double t ) const = 0;

    virtual Eigen::Matrix2Xd dposdparams ( const double t ) const = 0;
    virtual Eigen::Matrix2Xd dposdtdparams ( const double t ) const = 0;
    virtual Eigen::Matrix2Xd dposdtdtdparams ( const double t ) const = 0;
	virtual Eigen::Matrix2Xd dposdtdtdtdparams ( const double t ) const = 0;

    virtual double project ( const Eigen::Vector2d& point ) const = 0;

    // This has a default implementation
    virtual Eigen::VectorXd dtprojectdparams ( const double t, const Eigen::Vector2d& point ) const = 0;
    virtual Eigen::Matrix2Xd dposprojectdparams ( const double t, const Eigen::VectorXd& dtprojectdparams ) const = 0;
    virtual Eigen::Matrix2Xd dposdtprojectdparams ( const double t, const Eigen::VectorXd& dtprojectdparams ) const = 0;

    virtual void set_params ( const Eigen::VectorXd& params ) = 0;
    virtual Eigen::VectorXd get_params() const = 0;

    virtual Eigen::Matrix2Xd get_tesselation2() const = 0;
    virtual Eigen::VectorXd  get_tesselationt() const = 0;

    virtual ~GlobFitCurve() {}

    virtual GlobFitCurve* clone() const = 0;

	virtual geom::aabb get_bounding_box() const = 0;

    // Non virtual	
    double distance_to ( const Eigen::Vector2d& point ) const;
    double distance_to_sq ( const Eigen::Vector2d& point ) const;

    // Find where a line intersects the curve
    // line_src, line_dst, are the positions on the line
    // line_at: t values of the line where it intersects the curve (no bounds check, can be between -inf to inf)
    // this_curve_at: values of the t for the curve, where it intersects the line (has bound check , so between 0 and 1 ) 
    void approximate_intersect(
       const real2& line_src,
       const real2& line_dst,
       std::vector<double>& line_at,
       std::vector<double>&  this_curve_at );

	// This is a more practical yet stricter version of the one above.
	// It returns the number of intersection points, but fills in values only for
	// the one closest to the source of the line.
	int approximate_intersect_ray_closest(
		const real2& ray_o,
		const real2& ray_d,
		double& ray_at,
		double& this_curve_at
	);

	// Splits this curve at position t and returns the two new curve parts.
	// The mapping from old parameter values to new parameter values is
	// a simple re-scaling:
	//   t_new_left = t_old / t
	//   t_new_right = 1 - (1 - t_old) / (1 - t)
	// See split_t()
	virtual std::pair<GlobFitCurve*, GlobFitCurve*> split(double t) const = 0;

	// Returns the new curve parameters for a curve split at t_split. The return values 
	// are the parameters for the left and right part, respectively, that produce the
	// same point as the original curve at t.
	static std::pair<double, double> split_t(double t, double t_split);
};

NAMESPACE_END(polyvec)