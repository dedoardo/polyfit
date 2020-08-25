#pragma once

// polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/curve_fit_error.hpp>
#include <polyvec/pipeline_helper.hpp>

// Random forests
#include <andres/ml/decision-trees.hxx>

#define POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD 2

NAMESPACE_BEGIN(polyvec)

// Takes a string of the form "key1-value1_key2-value2_key3-value3" and extracts the numeric value for a given key.
// The delimiter characters can be arbitrary. The delimiter between key and value must be a single character.
int extract_number(const std::string& str, const std::string& key);

// Takes the metrics of a corner fit an classifies it into good/bad.
class FitClassifier
{
public:

	struct FitInfo
	{
		CurveFitError error;
		double distance_to_corner;

		FitInfo() { }
		FitInfo(const CurveFitError& error, double distance_to_corner)
			: error(error), distance_to_corner(distance_to_corner)
		{ }

		void write_to_file(const std::string& path) const;
		static FitInfo read_from_file(const std::string& path);

		bool operator<(const FitInfo& other) const;
	};

	struct PolygonCornerInfo
	{
		polyfit::CurveTracer::AccuracyMeasurement accuracy;
		double dist_prev;
		double dist_next;
		double angle;
		double raster_corner_size;
		int neighbors_with_different_convexity;

		int convexity[1 + 2 * POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD];
		
		PolygonCornerInfo() { }
		PolygonCornerInfo(const polyfit::CurveTracer::AccuracyMeasurement& accuracy, double dist_prev, double dist_next, double angle, double raster_corner_size, int neighbors_with_different_convexity);

		static void write(const polyfit::mat2x& B, const polyfit::vecXi& P, const polyfit::mat2x& PP, const std::string& path);
		static std::vector<PolygonCornerInfo> read(const std::string& path);
	};

	//Estimates if the fit characterized by the given parameters is good and returns the according probability.
	//fits: an array of at least four FitInfo objects of the fit with different boundary constraints
	//corner: information of the polygon corner underlying the fit
	virtual double evaluate_fit_probability(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner, const polyfit::mat2x& polygon, size_t polygon_corner) const = 0;
};

class FitClassifierGeometric
	: public FitClassifier
{
public:
	double evaluate_fit_probability(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner, const polyfit::mat2x& polygon, size_t polygon_corner) const;
};

class FitClassifierRandomForest
	: public FitClassifier
{
public:
	// Default filename for the serialized model
	static const std::string DEFAULT_FOREST_FILE;

	struct Sample
	{
		FitClassifier::FitInfo fits[4];
		FitClassifier::PolygonCornerInfo corner;
		unsigned char label;

		std::string model;

		size_t polygon; // the index of the corresponding polygon in the sample set
		size_t corner_idx; // the index of the corner in the corresponding polygon

		Sample() { }
		Sample(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner);

		static Sample read_from_data(const std::string& data_dir, const std::string& model, unsigned char label);
	};

	struct SampleSet
	{
		std::vector<Sample> samples;
		std::vector<PolygonData> polygons;
	};

	static size_t read_samples (
		const std::string& data_dir,           // Directory containing fitting information for each model
		const std::string& classification_dir, // Directory containing a Classification-{int}.txt file with labels 
		SampleSet& samples					   // Where the samples are appended
	);

	// Pre-computes derived data for the provided samples.
	static void prepare_samples(
		SampleSet& sample_set,					// The sample set whose samples to process
		std::vector<Sample>::iterator begin,	// The first sample to process
		bool remove_by_regularities,			// Specify whether to remove samples that are affected by regularities
		const std::string& data_dir				// The data directory of the samples, from which additional information are loaded
	);			

	void train (
		SampleSet& samples,
		const std::string& model_uri = DEFAULT_FOREST_FILE
	);

	// Train the random forests from input models and models
	SampleSet train (
		const std::string& data_dir,                            // Directory with the corner fitting data
		const std::string& classification_dir,                  // Directory with the class
		const std::string& model_uri = DEFAULT_FOREST_FILE // File where the trained model will be saved
	);

	// Load the serialized model
	void load_from_file(const std::string& file);
	void load_from_default_file();

	// Infer corner fitting type probability
	double evaluate_fit_probability(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner, const polyfit::mat2x& polygon, size_t polygon_corner) const;

	const andres::ml::DecisionForest<>& get_forest() const { return forest; }

private:
	andres::ml::DecisionForest<> forest;
};

class FitClassifierConstant
	: public FitClassifier
{
public:
	FitClassifierConstant(double constant_probability);
	double evaluate_fit_probability(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner, const polyfit::mat2x& polygon, size_t polygon_corner) const;
private:
	double constant_probability;
};

NAMESPACE_END(polyvec)