#pragma once

#include <polyvec/api.hpp>
#include <polyvec/curve-tracer/spline_types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/utils/union-find.hpp>

NAMESPACE_BEGIN(polyvec)

// Takes a current classification of polygon corners and produces a list of changes that the regularities impose. 
// ResultType: The format of the result (see HandledRegularityNewClassifications, HandledRegularityListOfChanges
// polygon: the vertices of the polygon
// regularity: regularity information in the polygon
// corner_classification: the current corner classification or an empty vector if no classifications are available
template <typename ResultType>
typename ResultType::TResult get_regularity_actions(
    const Eigen::Matrix2Xd& polygon, 
	const polyfit::Regularity::RegularityInformation& regularity, 
    const std::vector<polyvec::TangentFitType>& corner_classification,
    const int polygon_id
);


// Returns the handled regularities as a list of per-polygon corner classifications.
struct HandledRegularityNewClassifications
{
	using TResult = std::vector<TangentFitType>;
};

// Returns the handled regularities as a list of changes (corner id, new fit type). The 
// list of changes contains all potentially changed vertices, even if they already have
// the correct classification.
struct HandledRegularityListOfChanges
{
	using TResult = std::vector<std::pair<size_t, TangentFitType>>;
};

NAMESPACE_END(polyvec)