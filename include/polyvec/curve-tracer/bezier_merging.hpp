//
// bezier_merging.hpp
// To see the usage see the tests in apps/apps_test/bezier_merging/*
// and the app SH_export_traces.cpp
//

#pragma once

#include <polyvec/api.hpp>
#include <polyvec/regularity/find_symmetric_primitives.hpp>
#include <set>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( BezierMerging )

// =========================================================
// Public API
// =========================================================
using  Mat22 =  Eigen::Matrix<double, 2, 2>;
using  Mat23 =  Eigen::Matrix<double, 2, 3>;
using  Mat24 =  Eigen::Matrix<double, 2, 4>;

// Given an sequence of Beziers (and optionally a prepending and an appending line)
// fits a Bezier to all of them
// Input:
//   line_begin: either with size 0 or 1, the line at the beginning
//   beziers: with a size of at least 1, the beziers to fit to
//   line_end: either with size 0 or 1, the line at the end
//   verbose: a file to write the logs to (can be nullptr)
// Returns:
//   doable: was the operation successful
//   ans:    the control points of the merged Bezier
void merged_curve (
    const std::vector<Eigen::Matrix2d>& line_begin,
    const std::vector<Mat24>& beziers,
    const std::vector<Eigen::Matrix2d>& line_end,
    FILE* verbose,
    bool& doable,
    Mat24& ans );


// Checks if the merged bezier is an acceptable fit
// Input:
//  beziers: the small Beziers
//  merged_bezier: the merged Bezier
//  tolerance: slack for potrace envelope check
//  alpha_max: maximum acceptable alpha value for a Bezier
// Returns
//  potrace_penalty: a fit penalty defined by potrace. This is not really used in the
//    rest of the pipeline.
//  is_acceptable: whether the fit acceptable
void is_merge_acceptable (
    const std::vector<Mat24>& beziers,
    const Mat24& merged_bezier,
    const double tolerance,
    const double alpha_max,
    double& potrace_penalty,
    bool& is_acceptable  );

// The penalty of the fit, that is used in the shortest path formulation
// TODO: There should also be an option to pass line_begin and line_end to this function.
void merge_penalty (
    const std::vector<Mat24>& beziers,
    const Mat24& merged_bezier,
    double& penalty );

// Option for the merging function
struct MergeRecursivelyOptions {
    // Tolerance for the distance envelope check.
    double distance_tolerance = 0.5;
    // Increasing this will encourage merging more beziers.
    // A value of 0 will result in no merging.
    double per_bezier_const = 0.15;
    // Maximum acceptable alpha parameter for a bezier
    double alpha_max = 0.9;
    // Can lines also be merged into Beziers
    bool allow_merging_lines = false;
};

// The main merging function.
// Input:
//   curves_in: an array of lines and beziers
//   is_circular: is the path circular?
//   options: the options
//   curves_out: the merged curves
//   out2in: the index of the curves that were merged
void merge_curves (
    const std::vector<Eigen::Matrix2Xd>& curves_in,
    const bool is_circular,
    const MergeRecursivelyOptions options,
	const std::set<Regularity::SymmetryPair>& symmetric_curves,
    std::vector<Eigen::Matrix2Xd>& curves_out,
    std::vector<std::vector<int>>& out2in  );


// =========================================================
// Private stuff. Exposed for testing
// =========================================================

struct PotraceParameters;

struct PointParameters {
    // The control points
    Mat24 control_points;

    // Checks if the control points satisfy the potrace Bezier Parameterization.
    // ALso returns the point o and the alpha and beta parameters for this parameterization
    // (see formulas/bezier-merging.tex)
    void is_potracable ( bool& result, Eigen::Vector2d& oo, double& alpha, double& beta ) const;

    // Checks if the control points satisfy the potrace Bezier Parameterization.
    // ALso returns the point o and the alpha and beta parameters for this parameterization
    // (see formulas/bezier-merging.tex)
    bool is_potracable () const;

    // Return the potrace parameterization of the Bezier
    PotraceParameters as_potrace_params() const;
};

// Operations on the Potrace parameterization of a Bezier
struct PotraceParameters {

    // Alpha, beta parameters
    // (see formulas/bezier-merging.tex)
    double alpha;
    double beta;
    // The affine transformation
    // (see formulas/bezier-merging.tex)
    Eigen::Matrix2d A;
    Eigen::Vector2d b;

    // Get the area under the Bezier in the reference coordinates
    double get_areaz() const;
    // Get the area under the Bezier in the image coordinates
    double get_areay() const;
    // Get the control points in the reference coordinates
    PointParameters get_zz() const;
    // Get the o point in the reference coordinates
    static Eigen::Vector2d get_oz();
    // Get the control points in the image coordinates
    PointParameters get_yy() const;
    // Get the o point in the image coordinates
    Eigen::Vector2d get_oy() const;
    // Return the closest equiparameter (alpha=beta) Bezier
    PotraceParameters equiparamed() const;
    // Are the control point in a CCW order?
    bool is_ccw() const;

    // Input : tangent vector t
    // Output:
    //  Success: does the bezier acheive such tangent
    //  t: the t value at which the bezier acheives this tangent
    void t_for_tangent ( const Eigen::Vector2d& tangent, double& t, bool& success ) const;
};

// Finds an affine transformation that takes ponits0 to points1.
void find_affine_transformation ( const Mat23& points0, const Mat23& points1, Eigen::Matrix2d& A, Eigen::Vector2d& b );


// Saves and loads a sequence of curves to a file
void dump_curve_sequence ( const std::vector<Eigen::Matrix2Xd>& curves, FILE* file );
std::vector<Eigen::Matrix2Xd>  read_curve_sequence ( FILE* file );

NAMESPACE_END ( BezierMerging )
NAMESPACE_END ( polyfit )

