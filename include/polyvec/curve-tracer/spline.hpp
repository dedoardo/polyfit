/*
	The fitting logic is contained here...
	This needs to be separated into appropriate submodules
*/
#pragma once

// Polyvec
#include <polyvec/api.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/curve_fit_error.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/curve-tracer/spline_types.hpp>
#include <polyvec/curve-tracer/spline_validation.hpp>

#include <polyvec/curve-tracer/curve_fitter.hpp>

#include <polyvec/curve-tracer/measure_accuracy.hpp>
#include <polyvec/curve-tracer/measure_curvature.hpp>

// c++ stl
#include <vector>
#include <memory> // unique_ptr
#include <bitset>

NAMESPACE_BEGIN(polyvec)

struct CurvePrimitive;
class FitClassifier;

std::vector<bool> get_important_or_axis_aligned_edges(const Eigen::Matrix2Xd& polygon, const polyfit::Regularity::RegularityInformation& regularity);
std::vector<double> find_edge_angles_from_parallel(const Eigen::Matrix2Xd& polygon, const polyfit::Regularity::RegularityInformation& regularity);
void addCurveSequenceToFitter(polyvec::CurveFitter& curve_fitter, CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, int primitiveTo, int polygon_corners,
	const std::vector<bool>& is_edge_important_or_axisaligned, const std::vector<double>& edge_angles, bool allow_parallel_handling);

//Responsible for finding the correct primitive configurations for a polygon and constructing
//a primitive sequence. The sequence is optimized with the CurveFitter.
class CurveSequenceFitter {
public :
    CurveSequenceFitter(
        const Eigen::Matrix2Xd& points,
        const Eigen::Matrix2Xd& polygon,
        const std::vector<Eigen::Index>& polygon_vertices,
        const polyfit::Regularity::RegularityInformation& regularity,
        const int polygon_id = -1, // Used to filter out the regularities
        const std::vector<bool>& invert_edge_coordinate_system = std::vector<bool>(),
		const bool circular = true
	);

	CurvePrimitiveSequence all_with_type(TangentFitType type, bool optimize, bool fixFront, bool fixBack);

	CurvePrimitiveSequence fit_individual_corners(std::vector<std::pair<int, TangentFitType>> corners, bool circular, bool fix_front, bool fix_back, bool optimize = true, bool allow_parallel_handling = false);

	//Represents a consecutive sequence of corners on the polygon using the
	//indices of the first (inclusive) and last (inclusive) corner. The sequence
	//is circular, i.e. the first index can be greater than the last index.
	class CornerSequence
	{
	public:
		CornerSequence();

		CornerSequence(int first_incl, int last_incl);

		const int first_incl() const;
		const int last_incl() const;

		//TFunc: void(int)
		template <typename TFunc>
		void for_each_corner(int polygon_corners, TFunc&& callback);

		int included_corners(int polygon_corners) const;

	private:
		int _first_incl, _last_incl;
	};

	//Represents the results of an attempted fitting. Fittings may differ based
	//on the boundary constraints imposed.
	struct FittingAttempt
	{
		CurvePrimitiveSequence primitive_seq;

		std::vector<CurveFitError> corner_errors;			
		std::vector<FitValidity::Type> corner_valid;
		std::vector<double> accuracy_good_probability;
		std::vector<double> curvature_good_probability;
		std::vector<double> distances_to_corner;

		void resize(size_t n)
		{
			corner_errors.resize(n);
			corner_valid.resize(n);
			accuracy_good_probability.resize(n);
			curvature_good_probability.resize(n);
			distances_to_corner.resize(n, std::numeric_limits<double>::max());
		}

		// measures metrics derived from a primitive sequence
		void measure_sequence(const Eigen::Matrix2Xd& polygon);
	};

	struct EvolutionaryFittingState
	{
		EvolutionaryFittingState() { }
		EvolutionaryFittingState(const Eigen::Matrix2Xd& rasterPoints, const Eigen::VectorXi& polygonPoints, bool circular);

		//Stores the belief about each polygon corner. The bitset stores what fit types
		//are deemed viable for a given corner.
		std::vector<std::bitset<TANGENT_FIT_SAMPLES_COUNT>> cornerBelief;

		std::vector<double> polygon_corner_angles;
		std::vector<double> polygon_raster_corner_sizes;

		std::vector<double> corner_good_probability;

		bool is_corner_certain(int corner) const;
		TangentFitType get_highest_priority_possible_type(int corner) const;

		std::vector<CornerSequence> get_connected_sequences() const;
	};

	CurvePrimitiveSequence fit_evolutionary_simple(std::vector<TangentFitType>& outCornerTypes, const std::function<void(const EvolutionaryFittingState& state, const FittingAttempt* fits)>& callback);
	CurvePrimitiveSequence fit_evolutionary_simple(FitClassifier& classifier, std::vector<TangentFitType>& outCornerTypes, const std::function<void(const EvolutionaryFittingState& state, const FittingAttempt* fits)>& callback);

    // Calculates the tangent interpolation type for each corner in the polygon
    void classify_evolutionary_simple (
        FitClassifier& classifier,
        std::vector<TangentFitType>& outCornerTypes,
        const std::function<void(const EvolutionaryFittingState& state, const FittingAttempt* fits)>& callback
    );

    CurvePrimitiveSequence fit_initial_guess (
        const std::vector<TangentFitType>& tangentFits
    );

private:
	//input
	
    //raster boundary points
	const Eigen::Matrix2Xd& points;
	
    //polygon corner points
	const Eigen::Matrix2Xd& polygon;
	
    //indices of polygon in raster boundary
	const std::vector<Eigen::Index>& polygon_vertices;
	
    //indices of polygon in raster boundary
	Eigen::VectorXi polygonV;

    //id of the polygon which is fit by the curve
    const int polygon_id;

	const bool circular;
	const polyfit::Regularity::RegularityInformation& regularity;

    std::vector<bool> invert_edge_coordinate_system;

	void add_polygon_normals();

	void add_corner_fit(CurvePrimitiveSequence& seq, const Eigen::Index corner, TangentFitType fitType);

	//Optimizes the primitives [primitiveFrom, primitiveTo) in the sequence while
	//adding G0 and G1 constraints between them. Merges consecutive parallel lines.
	//The primitiveTo index is changed accordingly.
    void solve(CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, int& primitiveTo, bool allow_parallel_handling);
	void solve(CurvePrimitiveSequence& seq, bool circular, bool allow_parallel_handling);
	void solve(CurvePrimitiveSequence& seq, bool circular, int primitiveFrom, bool allow_parallel_handling);

	CurvePrimitive* first_bezier(CurvePrimitiveSequence& seq, int firstPrimitive) const;
	CurvePrimitive* last_bezier(CurvePrimitiveSequence& seq, int firstPrimitive) const;

	void fix_front(CurvePrimitiveSequence& seq, int firstPrimitive);
	void fix_back(CurvePrimitiveSequence& seq, int firstPrimitive);

	// ----- Helper methods for evolutionary fitting -----		
	std::vector<polyfit::CurveTracer::AccuracyMeasurement> get_polygon_measures() const;		

    // ---------------------------------------------------
	//Represents a curve in the initial (unoptimized) fits
	struct CornerFitPrototypePrimitive
	{
		GlobFitCurveParametrization* curve;
		FittingInfo fitting_info;

		bool hold_front_tangent = false;
		bool hold_end_tangent = false;
	};

	// Returns if the edge has a tangent we want to preserve (i.e., axis-aligned of 45°)
	bool polyvec::CurveSequenceFitter::edge_has_important_tangent(int i_edge) const;

	//Prototypes for fits for each polygon corner. Each corner is represented by a number of individual curves.
	std::vector<std::vector<CornerFitPrototypePrimitive>> initial_fits[TANGENT_FIT_SAMPLES_COUNT];
	
	//Initializes the fits for the fitting types up to FIT - 1 for all corners.
	template <TangentFitType FIT = TANGENT_FIT_SAMPLES_COUNT>
	void init_all_fits();

	//Initializes the fits for a specific fitting type for all corners.
	template <TangentFitType FIT>
	void init_fit();

	template <TangentFitType FIT>
	void init_corner_fit(int corner);

	std::vector<Eigen::Vector2d> corner_normals;
	std::vector<bool> edge_is_important;

	//Prescribes fixed angles for specific edges. Is nan for unfixed edge angles
	std::vector<double> edge_angles;	

	//Accuracy measurement of the polygon
	std::vector<polyfit::CurveTracer::AccuracyMeasurement> polygon_measures;

	double raster_aabb_diagonal;
};

std::string fitting_suffix();

// Isolated the test of merge_consecutive_parallel_lines so that it can be called from outside
bool can_merge_consecutive_lines(
    const Eigen::Vector2d& tangent_0,
    const Eigen::Vector2d& tangent_1
);

// skip_vertices is a bitvector (length of the total polygon vertices) which prevent lines to be merged across them
void merge_consecutive_parallel_lines(
    CurvePrimitiveSequence& seq, 
    bool circular, 
    int primitiveFrom, 
    int& primitiveTo,
    const std::vector<bool>& prevent_merge = std::vector<bool>()
);


NAMESPACE_END(polyvec)
