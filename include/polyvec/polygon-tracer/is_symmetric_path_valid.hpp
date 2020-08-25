/*
	Compares the new candidate path against the old one with respect
	to a symmetric region.
	If the new path doesn't introduce new inflections the path is accepted.
	If an inflection is introduced we expect to be in regions where the polygonal
	approximation could have failed to capture an inflection.
*/
#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

bool is_symmetric_path_valid(
	const mat2x& B,                   // raster boundary
	const vecXi& P_old,               // current polygon
	const vecXi& P_new,               // candidate polygon
	const std::vector<Vertex>& S_ord, // corners in the symmetric region ordered w.r.t. B (hence P)
									  // the sequence doesn't have to be continuous as in the case of
									  // of disjoint symmetries
	const bool circular
);

NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)