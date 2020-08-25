#pragma once

#include <polyvec/api.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/curve-tracer/spline_types.hpp>

#include <set>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

struct SymmetryPair
{
	int first, second;
	int region;

	SymmetryPair();

	SymmetryPair(int first, int second, int region);

	bool operator<(const SymmetryPair& other) const;
};

std::set<SymmetryPair> find_symmetric_primitives(const std::vector<polyfit::Regularity::Symmetry>&edge_symmetries, const polyvec::CurvePrimitiveSequence& curves, int n_edges);

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)