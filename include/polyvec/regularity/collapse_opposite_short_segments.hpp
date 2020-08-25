#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(Regularity)

void collapse_opposite_short_segments(
	const polyfit::mat2x& B,
	polyfit::vecXi& P,
	const polyfit::Regularity::RegularityInformation& RE,
	const bool circular
);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyvec)