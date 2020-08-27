#include <polyvec/pipeline_helper.hpp>

#include <polyvec/core/log.hpp>
#include <polyvec/core/options.hpp>
#include <polyvec/core/constants.hpp>
#include <polyvec/utils/arguments.hpp>
#include <polyvec/io/image.hpp>
#include <polyvec/image-segment/image_segment.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/polygon-tracer/regularized.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>
#include <polyvec/polygon-tracer/iterative-global.hpp>
#include <polyvec/polygon-tracer/postprocessor.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/curve-tracer/refinement.hpp>

#include <polyvec/eigen_serialization.hpp>

// todo: remove
#include <polyvec/polygon-tracer/iterative-local.hpp>
#include <polyvec/polygon-tracer/minimum.hpp>

#include <iostream>

using namespace polyvec;
using namespace polyfit;
using namespace Eigen;

void InputData::serialize(std::ostream& s) const
{
	if (!s.good())
		throw std::runtime_error("Output stream is not valid.");
	s << Serialized(R);
}

InputData InputData::deserialize(std::istream& s)
{
	if (!s.good())
		throw std::runtime_error("Input stream is not valid.");
	InputData data;
	s >> Serialized(data.R);
	return data;
}

void PolygonData::serialize(std::ostream& s) const
{
	if (!s.good())
		throw std::runtime_error("Output stream is not valid.");
	s << Serialized(B) << " " << Serialized(P);
}

PolygonData PolygonData::deserialize(std::istream& s, const InputData& input)
{
	if (!s.good())
		throw std::runtime_error("Input stream is not valid.");
	PolygonData data;
	s >> Serialized(data.B) >> Serialized(data.P);
	data.compute_derived_structures();
	return data;
}

void polyvec::init_pipeline()
{
	pvec::load_options_default();
}

InputData polyvec::load_input(const std::string& image_path)
{
	// Reads an image
	IO::Image I;

	if (!IO::read_image(image_path.c_str(), I)) {
		PF_LOGF("image not found %s", image_path.c_str());
		throw std::runtime_error("Image not found");
	}

	PF_VERBOSE_F("read image %d %d from %s", I.width, I.height, image_path.c_str());

	// extracts overlapping closed boundaries
	InputData data;
	std::vector<vec4> colors;
	std::vector<polyfit::mat2x> boundaries;

	if (!I.pixels) {
		PF_LOGF("image not found %s", image_path.c_str());
		throw std::runtime_error("Image not found");
	}

	ImageSegment::expand_and_cleanup(I);
	ImageSegment::extract_closed_regions(I, boundaries, colors);
	
	if (boundaries.empty()) {
		PF_LOGS("boundary extraction failed");
		throw std::runtime_error("Boundary extraction failed.");
	}

    printf("boundaries %d\n", (int)boundaries.size());

	// finds the closed boundary most likely to be a binary image
	data.R = boundaries.at(ImageSegment::find_binary_color_region(boundaries));
	PF_VERBOSE_F("boundary has %lld points", data.R.cols());

	return data;
}

PolygonData polyvec::extract_polygon(const InputData& input)
{
	PolygonData polygon;

	// Constructs a fitting path given the inputs boundary, placing midpoints on symmetric
	// segments on length <= 2
	BoundaryGraph::create_fitting_path(input.R, polygon.B, polygon.M, true);

	polygon.raster_symmetries = Symmetry::find_longest(polygon.B, true);
    polygon.raster_symmetries_local = Symmetry::find_shortest(polygon.B, true);

	// For each pair of corners, if the line connecting the two samples is a valid approximation
	// of the underlying boundary, then it will be a candidate edge for the polygonal approximation
	BoundaryGraph::connect_valid_edges(polygon.B, polygon.M, polygon.E, true);

	// Finding an initial solution
	ShortestPath::State G_state;

    if (!polyfit::PolygonTracer::regularized(G_state, polygon.B, polygon.M, polygon.E, polygon.RE, polygon.P, polygon.raster_symmetries, polygon.raster_symmetries_local, true)) {
		PF_LOGS("failed to trace polygon");
		throw std::runtime_error("Failed to trace polygon");
	}
    
    // todo: hey jude
    polygon.PV.resize(polygon.P.size());
    for (size_t i = 0; i < polygon.PV.size(); ++i) {
        polygon.PV[i] = polygon.P(i);
    }

    BoundaryGraph::trace_to_points(polygon.B, polygon.P, polygon.PP);

	return polygon;
}

void PolygonData::compute_derived_structures()
{
	raster_symmetries = Symmetry::find_longest(B, true);
    raster_symmetries_local = Symmetry::find_shortest(B, true);

	BoundaryGraph::trace_to_points(B, P, PP);

	PV.resize(P.size());
	for (size_t i = 0; i < PV.size(); ++i) {
		PV[i] = P(i);
	}
	
	RE = polyfit::Regularity::RegularityInformation();
	RE.update(B, P, raster_symmetries, raster_symmetries_local, true);
	RE.add_continuations(B, P, true);
}

SplineData polyvec::extract_curves(const PolygonData& polygon) {
	SplineData spline;

	CurveSequenceFitter fitter(polygon.B, polygon.PP, polygon.PV, polygon.RE);
	spline.curves = fitter.fit_evolutionary_simple(spline.tangent_fits, nullptr);
	
	return spline;
}

SplineData polyvec::extract_curves_with_classifier(
	const PolygonData& polygon, 
	const char* classifier_uri,
	std::vector<AttemptInfo>* attempts
) {
	PF_ASSERT(classifier_uri);

	SplineData spline;

	FitClassifierRandomForest classifier;
	classifier.load_from_file(classifier_uri);
	CurveSequenceFitter fitter(polygon.B, polygon.PP, polygon.PV, polygon.RE);
	spline.curves = fitter.fit_evolutionary_simple(classifier, spline.tangent_fits, 
		[&] (const CurveSequenceFitter::EvolutionaryFittingState & state, const CurveSequenceFitter::FittingAttempt * fits) {
			if (attempts) {
				attempts->emplace_back();
				attempts->back().attempt[0] = fits[0];
				attempts->back().attempt[1] = fits[1];
				attempts->back().attempt[2] = fits[2];
				attempts->back().state = state;
			}
	});

	return spline;
}

SplineData polyvec::extract_curves_with_refinement(
    PolygonData& polygon,
	const char* classifier_uri
) {
	SplineData spline = extract_curves_with_classifier(polygon, classifier_uri, NULL);

	vecXi polygon_vertices = polygon.P;
	const int collapsed = polyvec::CurveTracer::collapse_asymmetric_constant_steps (
		polygon.B, polygon_vertices, polygon.RE,
		spline.tangent_fits, true
	);

	PF_VERBOSE_F("Collapsing %d asymmetric steps", collapsed);
	BoundaryGraph::trace_to_points(polygon.B, polygon_vertices, polygon.PP);

	// just what..
	std::vector<Index> polygon_indices(polygon_vertices.size());
	for (int i = 0; i < polygon_vertices.size(); ++i) {
		polygon_indices[i] = polygon_vertices(i);
	}

	// Refitting the curves
	CurveSequenceFitter fitter(
		polygon.B,
		polygon.PP,
		polygon_indices,
		polygon.RE,
        -1,
        std::vector<bool>(),
		true
	);

	std::vector<std::pair<int, TangentFitType>> corner_fits;
	for (int i = 0; i < polygon_vertices.size(); ++i) {
		corner_fits.emplace_back(i, spline.tangent_fits[i]);
	}

	spline.curves = fitter.fit_individual_corners(corner_fits, true, false, false);
	return spline;
}

PostProcessor polyvec::post_process(PolygonData& polygon, SplineData& spline, const char* classifier_uri) {
    polyvec::PostProcessor post_processor(
        polygon.B, polygon.M, polygon.E, polygon.P, polygon.PP, polygon.raster_symmetries, polygon.raster_symmetries_local, polygon.RE, spline.tangent_fits, spline.curves, classifier_uri
    );

    post_processor.collapse_short_segments();
	//post_processor.merge_curves();

    return post_processor;
}