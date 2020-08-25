#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

/*
	For every pair of opposite corners which are directly facing each other
	(the closest ones are considered). They are moved within the single step.
*/
void align_opposite_corners (
	const mat2x& B,
	std::vector<BoundaryGraph::Edge>& E,
	polyfit::Regularity::RegularityInformation&    RE,
	vecXi& V,
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)