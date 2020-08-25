#include <polyvec/curve-tracer/regularity_handling.hpp>

#include <polyvec/utils/matrix.hpp>

#include <numeric>
#include <iostream>

#define MAKE_CONTINUATION_C0 1
#define HANDLE_SYMMETRY 1
#define HANDLE_SYMMETRIC_CORNERS_WITH_DIFFERENT_EDGE_LENGTHS 0
#define MAKE_PARALLEL_VERTICES_EQUAL 0
#define PRESERVE_AXIS_ALIGNED 1
#define MAKE_CONNECTED_PARALLEL_EQUAL 1

struct RegularityHandlingUFEntryData {
	polyvec::TangentFitType type = polyvec::TANGENT_FIT_LERP;
	bool changed = false;

	void change(polyvec::TangentFitType new_type)
	{
		changed = true;
		if(new_type > type)
			type = new_type;
	}
};
using RegularityHandlingUF = polyvec::UnionFind<RegularityHandlingUFEntryData>;

template <typename ResultType>
typename ResultType::TResult extract_result(RegularityHandlingUF& uf);

template <>
typename polyvec::HandledRegularityNewClassifications::TResult extract_result<polyvec::HandledRegularityNewClassifications>(RegularityHandlingUF& uf)
{
	polyvec::HandledRegularityNewClassifications::TResult result(uf.size());
	for (int i = 0; i < result.size(); ++i) {
		auto& parent = uf[uf.getRepresentative(i)];
		result[i] = parent.type;
	}
	return result;
}

template <>
typename polyvec::HandledRegularityListOfChanges::TResult extract_result<polyvec::HandledRegularityListOfChanges>(RegularityHandlingUF& uf)
{
	polyvec::HandledRegularityListOfChanges::TResult result;
	for (int i = 0; i < uf.size(); ++i) {
		auto& parent = uf[uf.getRepresentative(i)];
		if (!parent.changed)
			continue;
		result.emplace_back(i, parent.type);
	}
	return result;
}

// Returns if the given polygon corner has the given tangent fit type or no types are provided
bool is_classification_or_empty(const std::vector<polyvec::TangentFitType>& initial_corner_classification, const RegularityHandlingUF& uf, size_t corner, polyvec::TangentFitType compare_type)
{
	return initial_corner_classification.empty() || uf[corner].type == compare_type;
}

template <typename ResultType>
typename ResultType::TResult polyvec::get_regularity_actions(
	const Eigen::Matrix2Xd& polygon,
	const polyfit::Regularity::RegularityInformation& regularity,
	const std::vector<polyvec::TangentFitType>& corner_classification,
	const int polygon_id
) {
	// =========  We still have the UnionFind here although we don't need  =========
	// =========  it with the current approach. I'll keep it to be able    =========
	// =========  to quickly switch to a different strategy.               =========

	RegularityHandlingUF uf(polygon.cols());
	if (!corner_classification.empty())
		for (int i = 0; i < polygon.cols(); ++i)
			uf[i].type = corner_classification[i];

	auto change = [&](size_t v, polyvec::TangentFitType type)
	{
		uf[v].change(type);
	};

	auto merge = [&](size_t v1, size_t v2)
	{
		auto& node1 = uf[v1];
		auto& node2 = uf[v2];
		if (node1.changed || node2.changed)
			return false; //at least one of the nodes is changed by other regularities

		auto new_type = std::max(node1.type, node2.type);
		node1.change(new_type);
		node2.change(new_type);
		return true;
	};

	for (auto& r : regularity.continuations()) {

		// vertices around closures must be corners
#if MAKE_CONTINUATION_C0
		change(r.v0, TANGENT_FIT_CONSTANT);
		change(r.v1, TANGENT_FIT_CONSTANT);
#endif
		merge(r.v0, r.v1);
	}

	for (auto& r : regularity.parallels()) {
#if MAKE_PARALLEL_VERTICES_EQUAL
		// make the corresponding vertices equal
		if (r.aligned_00_11) {
			merge(r.v00, r.v11);
		}

		if (r.aligned_01_10) {
			merge(r.v01, r.v10);
		}
#endif
#if MAKE_CONNECTED_PARALLEL_EQUAL
		if (r.connected_00_11)
			merge(r.v00, r.v11);
		if (r.connected_01_10)
			merge(r.v01, r.v10);
#endif
	}

	for (auto& r : regularity.vertex_symmetries()) {
		if (r.is_relaxed)
			continue;
		// make the corresponding vertices equal
#if HANDLE_SYMMETRY
		if (!merge(r.v0, r.v1))
			continue;
#if HANDLE_SYMMETRIC_CORNERS_WITH_DIFFERENT_EDGE_LENGTHS
		//if the polygon lengths are not identical, this cannot be an asymmetric fit
		auto v0P = polygon.col(r.v0);
		auto v0LenNext = (polyfit::CircularAt(polygon, r.v0 + 1) - v0P).norm();
		auto v0LenPrev = (polyfit::CircularAt(polygon, r.v0 - 1) - v0P).norm();

		auto v1P = polygon.col(r.v1);
		auto v1LenNext = (polyfit::CircularAt(polygon, r.v1 + 1) - v1P).norm();
		auto v1LenPrev = (polyfit::CircularAt(polygon, r.v1 - 1) - v1P).norm();

		if (std::abs((v0LenNext - v1LenPrev)) > std::max(2.0, 0.2 * std::max(v0LenNext, v1LenPrev)) || std::abs((v0LenPrev - v1LenNext)) > std::max(2.0, 0.2 * std::max(v0LenPrev, v1LenNext))) {
			if (is_classification_or_empty(corner_classification, uf, r.v0, TANGENT_FIT_LERP))
				change(r.v0, TANGENT_FIT_LERP_SYM);

			if (is_classification_or_empty(corner_classification, uf, r.v1, TANGENT_FIT_LERP))
				change(r.v1, TANGENT_FIT_LERP_SYM);
		}
#endif
#endif
	}

	for (auto& r : regularity.important_edges()) {
#if PRESERVE_AXIS_ALIGNED
		if (is_classification_or_empty(corner_classification, uf, r.v0, TANGENT_FIT_LERP))
			change(r.v0, TANGENT_FIT_LERP_SYM);

		if (is_classification_or_empty(corner_classification, uf, polyfit::Circular(polygon, r.v0 + 1), TANGENT_FIT_LERP))
			change(polyfit::Circular(polygon, r.v0 + 1), TANGENT_FIT_LERP_SYM);
#endif
	}


	//Find the least common type for merged groups
	for (int i = 0; i < polygon.cols(); ++i) {
		auto rep = uf.getRepresentative(i);		
		if (rep != i)
		{
			auto& parent = uf[rep];
			parent.changed = true;
			uf[i].changed = true;
			parent.type = std::max(parent.type, uf[i].type);
		}
	}

	return extract_result<ResultType>(uf);
}


//explicit instantiations
template polyvec::HandledRegularityNewClassifications::TResult polyvec::get_regularity_actions<polyvec::HandledRegularityNewClassifications>(const Eigen::Matrix2Xd& polygon,
	const polyfit::Regularity::RegularityInformation& regularity, const std::vector<polyvec::TangentFitType>& corner_classification, const int polygon_id);

template polyvec::HandledRegularityListOfChanges::TResult polyvec::get_regularity_actions<polyvec::HandledRegularityListOfChanges>(const Eigen::Matrix2Xd& polygon,
	const polyfit::Regularity::RegularityInformation& regularity, const std::vector<polyvec::TangentFitType>& corner_classification, const int polygon_id);
