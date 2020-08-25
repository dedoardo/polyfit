#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

void pick_continuations(
	const mat2x& B, vecXi& V,
	std::vector<Continuation>& candidates,
	polyfit::Regularity::RegularityInformation& regularity,
	const bool circular,
	std::vector<polyfit::BoundaryGraph::Edge>& boundary_graph
);

std::vector<Continuation> find_continuation_candidates(
	const mat2x& B, const vecXi& V,
	const polyfit::Regularity::RegularityInformation& reg,
	const bool circular,
	const bool allow_move
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)