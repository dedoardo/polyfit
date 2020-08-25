#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	Given a list of points, it detects all edges which are
	visually grouped together. 

	A pair of straight edges is grouped together if the length of their overlap
	exceeds the distance between them. The length of the overlap is essentially
	the distance to other elements, following the principle of grouping from
	Gestalt psychology.

	The points do not necessarily need to face each other.

	Returns the number of pairs that have been detected
*/
int detect_parallel_edges_relative_proximity_axis_aligned (
	const mat2x& P,          // A list of points
	mat4xi& PL,              // pair of parallel edges, each column stores (src0, dst0, src1, dst1)
	const double min_feature_ratio,
	const double proximity_scale,
	double min_length_ratio, // Minimum required ratio of length between the two segments
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)