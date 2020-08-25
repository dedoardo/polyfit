/*
	Passing quality assurance.
	Should this be a set of free functions?
*/
#pragma once

// Polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/curve-tracer/tangent-fits.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(CurveTracer)

int collapse_asymmetric_constant_steps(
	const polyfit::mat2x& B,
	polyfit::vecXi& P,
	const polyfit::Regularity::RegularityInformation& RE,
	std::vector<TangentFitType>& tangent_fits,
	const bool circular
);

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyvec)