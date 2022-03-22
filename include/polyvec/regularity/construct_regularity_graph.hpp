#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/core/constants.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/regularity/symmetry.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>

#include <nse/IteratorRange.h>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

struct Parallel {
	int polygon_0 = -1, polygon_1 = -1;
	int v00 = -1, v01 = -1; // (src dst) of one of the dges
	int v10 = -1, v11 = -1; // (src dst) of the other edge

	bool aligned_00_11 = false;   // True if v00 and v11 are pixel aligned
	bool aligned_01_10 = false;   // True if v01 and v10 are pixel aligned
	bool connected_00_11 = false; // True if v00 and v11 are neighbors in the polygon
	bool connected_01_10 = false; // True if v01 and v10 are neighbors in the polygon

	std::string to_string() const;

	template <typename TFun>
	void reindex(TFun&& new_index_functor)
	{
		v00 = std::forward<TFun>(new_index_functor)(v00, -1);
		v01 = std::forward<TFun>(new_index_functor)(v01, 0);
		v10 = std::forward<TFun>(new_index_functor)(v10, 0);
		v11 = std::forward<TFun>(new_index_functor)(v11, 1);
	}
};

// The polypath (v0_prev, v0, v1, v1_next) is meant to be a single curve
// v0, v1 are the points facing each other
struct Continuation {
	Continuation() = default;
	Continuation(int v0_prev, int v0, int move_v0, int v1, int move_v1, int v1_next) :
		v0_prev(v0_prev), v0(v0), v1(v1), move_v0(move_v0), v1_next(v1_next), move_v1(move_v1) { }

	// Corners
	int polygon_0 = -1;
	int v0_prev = -1;
	int v0 = -1;
	int v0_next_guessed = -1;
	int move_v0;

	int polygon_1 = -1;
	int v1_prev_guessed = -1;
	int v1 = -1;
	int v1_next = -1;
	int move_v1;

	// Errors
	double distance_midpoints_sq = INFINITY;
	double angle_continuation_0 = INFINITY;
	double angle_continuation_1 = INFINITY;

	// todo: this two values can be computed from angle_continuation_*
	double angle_continuation_difference_0 = INFINITY;
	double angle_continuation_difference_1 = INFINITY;

	double angle_separation_0 = INFINITY;
	double angle_separation_1 = INFINITY;

	double angle_polygon_0 = INFINITY;
	double angle_polygon_1 = INFINITY;

	double intersection_angle = INFINITY;

	std::string to_string() const;

	template <typename TFun>
	void reindex(TFun&& new_index_functor)
	{
		v0_prev = std::forward<TFun>(new_index_functor)(v0_prev, -1);
		v0 = std::forward<TFun>(new_index_functor)(v0, 0);
		v0_next_guessed = std::forward<TFun>(new_index_functor)(v0_next_guessed, 1);

		v1_next = std::forward<TFun>(new_index_functor)(v1_next, 1);
		v1 = std::forward<TFun>(new_index_functor)(v1, 0);
		v1_prev_guessed = std::forward<TFun>(new_index_functor)(v1_prev_guessed, -1);
	}
};

struct Symmetry {
	int axis; // see Symmetry::Axis
	int size; // Number of points spanned by the symmetric region
	int region; // Index of the symmetric region
	int v0;
	int v1;

	bool is_relaxed;

	bool matches(const int v) const {
		return v0 == v || v1 == v;
	}

	bool matches(const int _v0, const int _v1) const {
		return (_v0 == v0 && _v1 == v1) || (_v0 == v1 && _v1 == v0);
	}

	// Returns the vertex that this symmetry associates to v. Can be v for self-symmetries.
	int other(int v) const
	{
		return v == v0 ? v1 : v0;
	}

	std::string to_string() const;

	template <typename TFun>
	void reindex(TFun&& new_index_functor)
	{
		v0 = std::forward<TFun>(new_index_functor)(v0, 0);
		v1 = std::forward<TFun>(new_index_functor)(v1, 0);
	}
};

struct ImportantEdge {
	int v0;

	std::string to_string() const;

	template <typename TFun>
	void reindex(TFun&& new_index_functor)
	{
		v0 = std::forward<TFun>(new_index_functor)(v0, 0);
	}
};

// Represents all regularity information for a raster / polygon
class RegularityInformation
{
public:
	void add(Parallel&& p);
	void add(Symmetry&& s);
	void add(Continuation&& s);
	void add_important_edge(int v0);
	void add_edge_symmetry(Symmetry&& s);

	const std::vector<Parallel>&		parallels() const { return _parallels; }
	const std::vector<Continuation>&	continuations() const { return _continuations; }
	
	class SymmetryIterator
	{
	public:
		SymmetryIterator(const RegularityInformation& reg, std::vector<int>::const_iterator current)
			: reg(reg), current(current)
		{ }

		bool operator!=(const SymmetryIterator& other) const { return current != other.current; }
		SymmetryIterator& operator++() { ++current; return *this; }
		const Symmetry& operator*() const { return reg._vertex_symmetries[*current]; }

	private:
		const RegularityInformation& reg;
		std::vector<int>::const_iterator current;
	};

	const std::vector<Symmetry>&	vertex_symmetries() const { return _vertex_symmetries; }
	const std::vector<Symmetry>&	edge_symmetries() const { return _edge_symmetries; }
	nse::util::IteratorRange<SymmetryIterator> vertex_symmetries(Vertex v) const;

	const std::vector<ImportantEdge>&	important_edges() const { return _important_edges; }

	void reset();

	void add_continuations(const mat2x& boundary, vecXi& polygon, std::vector<polyfit::BoundaryGraph::Edge>& E, bool circular, bool allow_polygon_modification);
	void add_continuations(const mat2x& boundary, vecXi& polygon, bool circular); // without allowing modifications to the polygon

	void clear_parallels() { _parallels.clear(); parallels_dirty = true; }
	void clear_continuations() { _continuations.clear(); per_edge_continuations_forward.clear(); per_edge_continuations_backward.clear(); }
	void clear_symmetries() { _vertex_symmetries.clear(); _edge_symmetries.clear(); per_vertex_symmetries.clear(); symmetries_dirty = true; }
	void clear_important_edges() { _important_edges.clear(); important_edges_dirty = true; }
	size_t size() const { return _parallels.size() + _continuations.size() + _vertex_symmetries.size() + _important_edges.size(); }

	// Returns the continuation of the directed edge (edge_src, edge_dst) or nullptr if the edge does not participate in a continuation.
	const Continuation* get_edge_continuation(Vertex edge_src, Vertex edge_dst) const;

	// Modifies the indices of all regularities so as to be consistent with a new vertex sequence that is produced
	// by deleting vertices from the original sequence.
	// Symmetries and important edges are invalidated
	void reindex_after_vertex_deletion(	const std::vector<int>& deleted_vertices_ordered, int old_vertex_count);

	// Re-calculates all dirty regularities except continuations
	void update(
		const mat2x& B,        // Fitting boundary
		const vecXi& V,        // Polygonal fit on B
        const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries,
        const std::vector<polyfit::Symmetry::SubPath>& raster_symmetries_local,
		const bool circular);

private:
	std::vector<Parallel>			_parallels;
	std::vector<Continuation>		_continuations;
	std::vector<Symmetry>			_vertex_symmetries;
	std::vector<Symmetry>			_edge_symmetries;
	std::vector<ImportantEdge>		_important_edges;

	std::vector<std::vector<int>>   per_vertex_symmetries;
	std::vector<int>				per_edge_continuations_forward, per_edge_continuations_backward;

	bool parallels_dirty = true;
	bool symmetries_dirty = true;
	bool important_edges_dirty = true;

	size_t polygon_vertices = 0;
};

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)