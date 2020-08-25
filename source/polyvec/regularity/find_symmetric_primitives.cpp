#include <polyvec/regularity/find_symmetric_primitives.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/core/log.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

SymmetryPair::SymmetryPair()
	: first(-1), second(-1), region(-1)
{ }

SymmetryPair::SymmetryPair(int first, int second, int region)
	: first(first), second(second), region(region)
{ }

bool SymmetryPair::operator<(const SymmetryPair& other) const
{
	if (first != other.first)
		return first < other.first;
	if (second != other.second)
		return second < other.second;
	return region < other.region;
}

std::set<SymmetryPair> find_symmetric_primitives(const std::vector<polyfit::Regularity::Symmetry>&edge_symmetries, const polyvec::CurvePrimitiveSequence& curves, int n_edges)
{
	std::set<SymmetryPair> symmetric_primitives;

	// for each polygon edge, gather its primitives
	std::vector<std::vector<int>> per_edge_primitives(n_edges);

	std::set<std::pair<int, int>> symmetric_edges;
	for (auto& s : edge_symmetries)
		symmetric_edges.emplace(std::min(s.v0, s.v1), std::max(s.v0, s.v1));

	bool passed_last_edge = false; // due to circularity, we ignore curves for the last edge in the first round
	auto add_primitive_to_edge = [&](int edge, int primitive) {
		if (edge != n_edges - 1 || passed_last_edge)
			per_edge_primitives[edge].push_back(primitive);
		if (edge != n_edges - 1)
			passed_last_edge = true;
	};
	for (int i = 0; i < curves.primitives.size(); ++i)
	{
		auto& p = curves.primitives[i];

		auto src_edge_id = polyfit::BoundaryGraph::unpack_edge_id(p.fitting_info.edge_src);
		auto dst_edge_id = polyfit::BoundaryGraph::unpack_edge_id(p.fitting_info.edge_dst);

		add_primitive_to_edge(src_edge_id(0), i);

		if (dst_edge_id(0) != src_edge_id(0))
			add_primitive_to_edge(dst_edge_id(0), i);
	}

	for (auto& s : edge_symmetries)
	{
		auto e1 = s.v0;
		auto e2 = s.v1;

		if (per_edge_primitives[e1].size() != per_edge_primitives[e2].size())
		{
			PF_VERBOSE_F("Edges %d-%d are symmetric but have different numbers of primitives.", e1, e2);
			continue;
		}

		for (int i = 0; i < per_edge_primitives[e1].size(); ++i)
		{
			int pi1 = per_edge_primitives[e1][i];
			int pi2 = per_edge_primitives[e2][per_edge_primitives[e2].size() - 1 - i];

			auto& p1 = curves.primitives[pi1];
			auto& p2 = curves.primitives[pi2];

			auto src1 = polyfit::BoundaryGraph::unpack_edge_id(p1.fitting_info.edge_src)(0);
			auto src2 = polyfit::BoundaryGraph::unpack_edge_id(p2.fitting_info.edge_src)(0);
			auto dst1 = polyfit::BoundaryGraph::unpack_edge_id(p1.fitting_info.edge_dst)(0);
			auto dst2 = polyfit::BoundaryGraph::unpack_edge_id(p2.fitting_info.edge_dst)(0);

			bool on_single_edge_p1 = src1 == e1 && dst1 == e1;
			bool on_single_edge_p2 = src2 == e2 && dst2 == e2;

			if (on_single_edge_p1 != on_single_edge_p2)
				continue; //not a symmetric primitive

			if (on_single_edge_p1 && pi1 != pi2)
			{
				symmetric_primitives.emplace(pi1, pi2, s.region);
				symmetric_primitives.emplace(pi2, pi1, s.region);
			}
			else
			{
				//need to check if the other edge pair is also symmetric
				int other_e1 = src1 == e1 ? dst1 : src1;
				int other_e2 = src2 == e2 ? dst2 : src2;

				if (other_e1 == other_e2 || symmetric_edges.find(std::make_pair(std::min(other_e1, other_e2), std::max(other_e1, other_e2))) != symmetric_edges.end() && pi1 != pi2)
				{
					symmetric_primitives.emplace(pi1, pi2, s.region);
					symmetric_primitives.emplace(pi2, pi1, s.region);
				}
			}
		}
	}

	return symmetric_primitives;
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)