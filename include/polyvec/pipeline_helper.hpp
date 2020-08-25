#pragma once

#include <string>
#include <vector>
#include <iostream>

#include <polyvec/core/types.hpp>
#include <polyvec/regularity/construct_regularity_graph.hpp>
#include <polyvec/polygon-tracer/boundary-graph.hpp>
#include <polyvec/polygon-tracer/postprocessor.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/regularity/symmetry.hpp>

namespace polyvec
{
	void init_pipeline();

	struct InputData
	{
		// Raster boundary
		polyfit::mat2x R;

		void serialize(std::ostream& s) const;
		static InputData deserialize(std::istream& s);
	};	

	InputData load_input(const std::string& image_path);

	struct PolygonData
	{
		// Fitting points obtained augmenting the raster boundary R
		polyfit::mat2x B;

		// polygon vertices
		polyfit::vecXi P;

		// polygon points
		polyfit::mat2x PP;

        // Concavities for each boundary point
        std::vector<int> C;

		// This should disappear 
        std::vector<polyfit::Vertex> M;
        std::vector<bool> midpoints;
        std::vector<polyfit::Vertex> PV;

		// Regularity relationships between pairs of edges and corners
		polyfit::Regularity::RegularityInformation RE;

		// Raster symmetries
		std::vector<polyfit::Symmetry::SubPath> raster_symmetries;
        std::vector<polyfit::Symmetry::SubPath> raster_symmetries_local; // this name is cryptic

        std::vector<polyfit::BoundaryGraph::Edge> E;

		void serialize(std::ostream& s) const;
		static PolygonData deserialize(std::istream& s, const InputData& input);

		void compute_derived_structures();
	};

	PolygonData extract_polygon(const InputData& input);

	struct SplineData
	{
		// The type of curve which has been fit at each polygon corner
		std::vector<TangentFitType> tangent_fits;

		// The spline curves
		CurvePrimitiveSequence curves;
	};

	SplineData extract_curves(const PolygonData& polygon);

	struct AttemptInfo
	{
		CurveSequenceFitter::FittingAttempt attempt[3];
		CurveSequenceFitter::EvolutionaryFittingState state;
	};

	SplineData extract_curves_with_classifier(
		const PolygonData& polygon, 
		const char* classifier_uri,
		std::vector<AttemptInfo>* attempts = nullptr
	);

	SplineData extract_curves_with_refinement(
		PolygonData& polygon,
		const char* classifier_uri
	);

    PostProcessor post_process(PolygonData& polygon, SplineData& spline, const char* classifier_uri);
}