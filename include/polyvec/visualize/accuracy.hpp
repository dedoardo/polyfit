/*
	Routines that draw accuracy error visualizations in the scene
*/
#pragma once

#include <polyvec/api.hpp>
#include <polyvec/curve-tracer/spline.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Visualize)

#if 0

// raster -> curve l1 measurements
void accuracy_measure_l1_signed_midpoints_with_slack(
	const std::vector<polyvec::CurvePrimitive>& prims
);

// raster -> curve l2 measurements
void accuracy_measure_l2_unsigned_midpoints(
	const std::vector<polyvec::CurvePrimitive>& prims
);

void accuracy_measure_unsigned_midpoints_with_slack(
	const std::vector<polyvec::CurvePrimitive>& prims
);

// raster -> polygon(curve) l1 measurements
void accuracy_measure_l1_signed_midpoints_with_slack_polygon(
	const mat2x& P, // polygon
	const std::vector<polyvec::CurvePrimitive>& prims
);

// raster -> polygon(curve) l2 measurements
void accuracy_measure_l2_unsigned_midpoints_polygon(
	const mat2x& P, // polygon
	const std::vector<polyvec::CurvePrimitive>& prims
);

// colors in red all the curves for which the polygon error exceeds the curves one
void accuracy_error_polygon_greater_l1(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
);

void accuracy_error_polygon_greater_l2(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
);

// colors in red all the curves for which the polygon error 
void accuracy_error_polygon_greater_l1(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
);

void accuracy_error_polygon_greater_l2(
	const mat2x& P,
	const std::vector<polyvec::CurvePrimitive>& prims
);

#endif

NAMESPACE_END(Visualize)
NAMESPACE_END(polyfit)