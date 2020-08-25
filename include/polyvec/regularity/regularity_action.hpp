#pragma once

#include <polyvec/api.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/polygon-tracer/iteration-callback.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PolygonTracer)

//Represents a regularity and the actions performed to the boundary graph.
class RegularityAction
{
public:
	virtual std::vector<size_t> get_edges_to_delete(const mat2x& raster, const vecXi& polygon, const std::vector<BoundaryGraph::Edge>& E, bool circular) const = 0;
	
    virtual bool can_introduce_degenerate_edges() const = 0;
    virtual bool can_introduce_inflections() const { return false; }
    virtual bool is_polygon_acceptable(const mat2x& raster, const vecXi& polygon_old, const vecXi& polygon_new, bool circular) const = 0;
};

// Applies regularities by modifying the edge set of a boundary graph. If a regularity
// turned out to produce an invalid result (or no result), its actions are reverted.
void apply_regularities(
	ShortestPath::State& G_state,
	const mat2x& B,
	const std::vector<Vertex>& M,
	std::vector<BoundaryGraph::Edge>& E,
	vecXi& P,
	const bool circular,
	IterationCallback on_iter,
	const std::vector<bool>& corners_C0,
	const std::vector<RegularityAction*>& regularity_actions);


NAMESPACE_END(PolygonTracer)
NAMESPACE_END(polyfit)