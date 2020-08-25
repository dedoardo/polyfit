/*
	Calculates various error measurements for the curve and the raster
	boundary being approximated.

	The _with_slack function fallback to the distance from the primitive
	to the center of the in/out pixel.
*/
#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

namespace polyvec {
	struct CurvePrimitive;
	class GlobFitCurve;
}

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

// if the measurement is unsigned, the in/out values will coincide
struct AccuracyMeasurement {
	double e_pos = 0.;    // maximum positive distance violation
	double e_neg = 0.;	 // maximum negative distance violation
	double e_l1 = 0.;
	double e_l2_sq = 0.;
	int    count = 0;

	void combine(const AccuracyMeasurement& m);

	// prepares the measurement for add_*
	void reset();

	// adds new signed measurements
	void add_pos(const double d);
	void add_neg(const double d);

	// adds an unsigned measurement
	void add(const double d);

	// finishes calculating the norm and resets the errors to 0. if no measurements were provided
	void finish();

	// returns the min/max unsigned distance from any of the points
	double max_error() const;
	double min_error() const;
};

AccuracyMeasurement measure_accuracy(
	polyvec::GlobFitCurve& curve,  // primitive tested
	const mat2x& P,                      // fitting points
	const mat2x& N                      // fitting point normals
);

/*
	First the primitive is test for intersection against the ray with origin at the midpoint
	and direction along the normal. If the resulting distance exceeds 0.5, it is recalculated
	as the distance from the pixel center to the curve. The smallest between the inner/outer
	pixel centers is taken and added to 0.5
*/
//void measure_accuracy_signed_extended(
//	polyvec::GlobFitCurve& curve,  // primitive tested
//	const mat2x& P,                      // fitting points
//	const mat2x& N,                      // fitting point normals
//	AccuracyMeasurement& m               // results
//);

/*
	The distance equals the distance to the midpoints from the closest point on the curve
*/
void measure_accuracy_unsigned(const polyvec::CurvePrimitive* prim, AccuracyMeasurement& m);

/*
	If the distance from the midpoint to the closest point on the curve exceeds 0.5, it is recaulated
	as the distance from the pixel center to the curve. The smallest between the inner/outer pixel
	centers is chosen and added to 0.5
*/
void measure_accuracy_unsigned_extended(const polyvec::CurvePrimitive* prim, AccuracyMeasurement& m);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)
