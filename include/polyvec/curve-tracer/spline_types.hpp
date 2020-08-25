#pragma once

#include <polyvec/api.hpp>

#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_fit_error.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyvec)

// Describes the primitive configuration for a single polygon corner
enum TangentFitType
{
	// Asymmetric fit with a single Bezier curve
	TANGENT_FIT_LERP = 0,

	// Symmetric fit with a Bezier curve and a line segment
	TANGENT_FIT_LERP_SYM,

	// Two line segments joined by a C0 corner
	TANGENT_FIT_CONSTANT,


	TANGENT_FIT_SAMPLES_COUNT
};
const char* tangent_fit_type_to_string(const TangentFitType v);

// Describes pieces of information to which a curve segment is fit (e.g. points, tangent directions...)
struct FittingInfo
{
    // Which edge does the source point of the primitive lie on
    polyfit::BoundaryGraph::EdgeID edge_src;

    // Which edge does the destination point of the primitive lie on
    polyfit::BoundaryGraph::EdgeID edge_dst;

	// Points to approximate (these are the midpoints of pixel boundaries)
	std::vector<Eigen::Vector2d> fit_midpoints;
	// Boundary normals for every point in fit_midpoints
	std::vector<Eigen::Vector2d> fit_midpoint_normals;

	struct Tangents
	{
		// Tangents to approximate at given parameter positions (see fit_tangent_ts)
		std::vector<Eigen::Vector2d> fit_tangents;
		// The parameter positions at which fit_tangents should be achieved
		std::vector<double>  fit_tangent_ts;
	};

	Tangents dense_tangents;
	Tangents sparse_tangents;

	// Adds tangent samples for given Bezier curve, such that tangent directions interpolate
	// linearly (in polar coordinates). tMiddle is the parameter position of the polygon corner,
	// i.e. where the midpoint of the interpolation is.
	void add_tangent_samples(GlobFitCurveParametrization* curve, double tMiddle);

	// Adds tangent samples for a given line segment.
	void add_tangent_samples(GlobFitCurve_Line* curve);

	// Adds point samples for a given raster boundary segment.
	void add_midpoint_samples(const Eigen::Matrix2Xd& points, const Eigen::VectorXi& polygonV, const size_t firstEdge, const double firstEdgeT, const size_t lastEdge, const double lastEdgeT, bool circular);
};

// Everything in this struct is only for debugging and not used in the fitting process
struct CurvePrimitiveDebugInfo
{	
	Eigen::Vector2d corner_pt;

	// debug (remove)
	std::vector<double> fit_tangents_t;

	// original endpoints
	Eigen::Vector2d endpoint_src;
	Eigen::Vector2d endpoint_dst;

	// line stuff
	bool hold_endpoint_dst = false;
};

// Describes a segment of the spline with additional data.
struct CurvePrimitive : public CurvePrimitiveDebugInfo {
    std::shared_ptr<GlobFitCurveParametrization> curve;
    FittingInfo fitting_info;
    
    CurveFitError       error;
	std::vector<double> objective_energies;	

	bool hold_front_tangent = false;
	bool hold_end_tangent = false;

	// the corner index of the polygon that this curve is fit to
	int corner;
	// determines if this primitive also fits the next polygon corner
	bool fits_next_corner = false;

	// Splits the underlying curve into two parts and re-distributes
	// the corresponding fitting terms.
	std::pair<CurvePrimitive, CurvePrimitive> split(double t) const;

	double point_weight_multiplier = 1.0;

	bool fits_same_corner(const CurvePrimitive& other, int polygon_corners) const;
};

struct CurvePrimitiveSequence {
	std::vector<CurvePrimitive> primitives;
};

NAMESPACE_END(polyvec)