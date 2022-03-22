#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/curve-tracer/measure_accuracy_polygon.hpp>
#include <polyvec/utils/system.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/curve-tracer/regularity_handling.hpp>

#include <andres/marray.hxx>

#include <fstream>
#include <map>

#define REORDER_FITS 1
#define DUPLICATE_SAMPLES_REORDERED 0

#define OUTPUT_TRAINING_DATA 1

#define REMOVE_SAMPLES_AFFECTED_BY_REGULARITIES 1

using namespace polyvec;

const std::string FitClassifierRandomForest::DEFAULT_FOREST_FILE = "fitting_forest.txt";

namespace {
	double read_key_value_pair(std::ifstream& f, const std::string& key)
	{
		std::vector<char> read_key(key.length() + 1);
		f.get(read_key.data(), key.length() + 1);
		if (strcmp(read_key.data(), key.c_str()) != 0)
		{
			std::cerr << "Expecting key " << key << " but encountered " << read_key.data() << std::endl;
			throw std::runtime_error("Encountered a different key");
		}
		f.ignore(std::numeric_limits<std::streamsize>::max(), ':');

		double value;
		f >> value;
		f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		return value;
	}

	void read_convexity_values (std::ifstream& f, const std::string& key, int* values) {
		std::vector<char> read_key(key.length() + 1);
		f.get(read_key.data(), key.length() + 1);
		if (strcmp(read_key.data(), key.c_str()) != 0)
		{
			std::cerr << "Expecting key " << key << " but encountered " << read_key.data() << std::endl;
			throw std::runtime_error("Encountered a different key");
		}
		f.ignore(std::numeric_limits<std::streamsize>::max(), ':');

		for (int nhb = -POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD; nhb <= POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD; ++nhb) {
			const int nhb_idx = nhb + POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD;
			f >> values[nhb_idx];

			if (values[nhb_idx] != 1 && values[nhb_idx] != -1) {
				throw std::runtime_error("Convexity values should be either +1 or -1");
			}
		}
		f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
}

FitClassifier::PolygonCornerInfo::PolygonCornerInfo(const polyfit::CurveTracer::AccuracyMeasurement& accuracy, double dist_prev, double dist_next, double angle, double raster_corner_size, int neighbors_with_different_convexity)
	:accuracy(accuracy), dist_prev(dist_prev), dist_next(dist_next), angle(angle), raster_corner_size(raster_corner_size), neighbors_with_different_convexity(neighbors_with_different_convexity)
{ }

void FitClassifier::PolygonCornerInfo::write(const polyfit::mat2x & B, const polyfit::vecXi & P, const polyfit::mat2x & PP, const std::string & path)
{
	int convexity_nhbs[1 + 2 * POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD];

	std::ofstream polygon_file(path);
	for (int corner = 0; corner < PP.cols(); ++corner)
	{
		auto to_next = PP.col((corner + 1) % PP.cols()) - PP.col(corner);
		auto to_prev = PP.col((corner - 1 + PP.cols()) % PP.cols()) - PP.col(corner);

		polyfit::CurveTracer::AccuracyMeasurement error;
		polyfit::CurveTracer::measure_accuracy_signed_extended_polygon_fit_asymmetric(B, P, corner, error, true);
		polygon_file << "acc max: " << error.max_error() << "\n";
		polygon_file << "acc pos: " << error.e_pos << "\n";
		polygon_file << "acc neg: " << error.e_neg << "\n";
		polygon_file << "dist prev: " << to_prev.norm() << "\n";
		polygon_file << "dist next: " << to_next.norm() << "\n";
		polygon_file << "angle: " << std::acos(std::min(1.0, std::max(-1.0, to_next.dot(to_prev) / (to_next.norm() * to_prev.norm())))) << "\n";
		polygon_file << "corner size: " << polyfit::GeomRaster::raster_corner_size(B, P(corner)) << "\n";
		polygon_file << "convex diff: " << polyfit::AngleUtils::number_of_neighbors_with_different_convexity(PP, corner) << "\n";
		
		polyfit::AngleUtils::neighborhood_convexity(PP, corner, POLYGON_CORNER_INFO_CONVEXITY_NEIGHBORHOOD, convexity_nhbs);
	}
	polygon_file.close();
}

void FitClassifier::FitInfo::write_to_file(const std::string & path) const
{
	std::ofstream info_file(path);
	info_file << "acc max: " << error.accuracy.max_error() << "\n";
	info_file << "acc pos: " << error.accuracy.e_pos << "\n";
	info_file << "acc neg: " << error.accuracy.e_neg << "\n";
	info_file << "curv: " << error.curvature.r_min << "\n";
	info_file << "dist2corner: " << distance_to_corner << "\n";
	info_file.close();
}

std::vector<FitClassifier::PolygonCornerInfo> FitClassifier::PolygonCornerInfo::read(const std::string& path)
{
	std::ifstream f(path);
	if (!f.good())
		throw std::runtime_error("Cannot open file");

	std::vector<PolygonCornerInfo> result;

	while (!f.eof())
	{
		if (f.peek() < 0)
			break;
		PolygonCornerInfo info;
		double max_acc = read_key_value_pair(f, "acc max");
		info.accuracy.e_pos = read_key_value_pair(f, "acc pos");
		info.accuracy.e_neg = read_key_value_pair(f, "acc neg");		
		assert(std::abs(max_acc - info.accuracy.max_error()) < 0.01);
		info.dist_prev = read_key_value_pair(f, "dist prev");
		info.dist_next = read_key_value_pair(f, "dist next");
		info.angle = read_key_value_pair(f, "angle");
		info.raster_corner_size = read_key_value_pair(f, "corner size");
		info.neighbors_with_different_convexity = (int)read_key_value_pair(f, "convex diff");

		result.push_back(info);
	}

	return result;
}

FitClassifier::FitInfo FitClassifier::FitInfo::read_from_file(const std::string& _path)
{
	std::string path = _path;
	for (size_t i = 0; i < path.size(); ++i) {
		if (path[i] == '\\') {
			path[i] = '/';
		}
	}

	std::ifstream f(path);
	if (!f.good())
	{
		std::cerr << "Cannot open fit info file: " << path << std::endl;
		throw std::runtime_error("Cannot open file");
	}

	FitInfo info;
	double max_accuracy = read_key_value_pair(f, "acc max");
	info.error.accuracy.e_pos = read_key_value_pair(f, "acc pos");
	info.error.accuracy.e_neg = read_key_value_pair(f, "acc neg");
	info.error.curvature.r_min = read_key_value_pair(f, "curv");
	info.distance_to_corner = read_key_value_pair(f, "dist2corner");
	return info;
}

#define PARTIAL_COMPARE(attribute) if(attribute != other.attribute) return attribute < other.attribute;

bool FitClassifier::FitInfo::operator<(const FitInfo& other) const
{
	PARTIAL_COMPARE(error.accuracy.max_error());
	PARTIAL_COMPARE(error.curvature.r_min);
	PARTIAL_COMPARE(error.accuracy.e_pos);
	PARTIAL_COMPARE(error.accuracy.e_neg);
	return false;
}


int polyvec::extract_number(const std::string& str, const std::string& key)
{
	auto loc = str.find(key);
	std::stringstream ss(str.substr(loc + key.size() + 1));
	int number;
	ss >> number;
	return number;
}

FitClassifierRandomForest::Sample::Sample(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner)
	: corner(corner)
{
	for (int i = 0; i < 3; ++i)
		this->fits[i] = fits[i];
}
	
FitClassifierRandomForest::Sample FitClassifierRandomForest::Sample::read_from_data(const std::string& data_dir, const std::string& model, unsigned char label)
{
	FitClassifierRandomForest::Sample sample;
	for (int i = 0; i < 3; ++i)
	{
		bool fix0 = i & 1 == 1;
		bool fix1 = (i >> 1) & 1 == 1;
		auto info_file = data_dir + model + "_fix-" + std::to_string(fix0) + std::to_string(fix1) + "_info.txt";
		sample.fits[i] = FitClassifier::FitInfo::read_from_file(info_file);
	}

	sample.corner_idx = extract_number(model, "corner");	
	sample.label = label;
	sample.model = model;

	return sample;
}

struct FeatureDescription
{
	std::vector<std::function<double(const FitClassifierRandomForest::Sample&, const polyfit::mat2x& polygon, size_t polygon_corner)>> features;

	FeatureDescription()
	{
		for (int j = 0; j < 3; ++j)
		{
			features.push_back([j](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return nth_fit(s, j).error.accuracy.max_error(); });
			//features.push_back([j](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return nth_fit(s, j).error.accuracy.e_pos; });
			//features.push_back([j](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return nth_fit(s, j).error.accuracy.e_neg; });
			//features.push_back([j](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return nth_fit(s, j).distance_to_corner; });
			features.push_back([j](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return nth_fit(s, j).error.curvature.r_min; });
			features.push_back([j](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return nth_fit(s, j).error.accuracy.max_error() - s.corner.accuracy.max_error(); });
		}

		features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.angle; });
		
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.raster_corner_size; });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.neighbors_with_different_convexity; });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return std::max(s.corner.dist_next, s.corner.dist_prev) / std::min(s.corner.dist_next, s.corner.dist_prev); });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.accuracy.max_error(); });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.accuracy.e_pos; });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.accuracy.e_neg; });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.dist_next; });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return s.corner.dist_prev; });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return std::min(s.corner.dist_prev, s.corner.dist_next); });
		//features.push_back([](const FitClassifierRandomForest::Sample& s, const polyfit::mat2x& polygon, size_t polygon_corner) { return std::max(s.corner.dist_prev, s.corner.dist_next); });
	}

private:

	static const FitClassifier::FitInfo& nth_fit(const FitClassifierRandomForest::Sample& s, int n)
	{
#if REORDER_FITS
		if (n == 0)
			return s.fits[0];
		if (n == 1)
		{
			if (s.fits[1] < s.fits[2])
				return s.fits[1];
			else
				return s.fits[2];
		}
		if(n == 2)
			if (s.fits[1] < s.fits[2])
				return s.fits[2];
			else
				return s.fits[1];
		assert_break(false);
#else
		return s.fits[n];
#endif
	}

} feature_description;

double calculate_prediction_accuracy(const andres::Marray<double>& probabilities, const std::vector<FitClassifierRandomForest::Sample>& samples, size_t samples_offset)
{
	int correct_predictions = 0;
	double avg_certainty_of_wrong = 0;
	for (int i = 0; i < probabilities.shape(0); ++i)
	{
		double max_p = 0;
		unsigned char max_i = -1;
		for (unsigned char j = 0; j < 2; ++j)
			if (probabilities(i, j) > max_p)
			{
				max_p = probabilities(i, j);
				max_i = j;
			}
		if (max_i == samples[i + samples_offset].label)
			++correct_predictions;
		else
		{
			std::cout << samples[i + samples_offset].model << ": " << "wrong prediction with probabilities " << probabilities(i, 0) << " / " << probabilities(i, 1) << std::endl;
			avg_certainty_of_wrong += (max_p - 0.5) / 0.5;
		}
	}
	avg_certainty_of_wrong /= (probabilities.shape(0) - correct_predictions);
	std::cout << "Average certainty of wrong decisions: " << avg_certainty_of_wrong << std::endl;
	return (double)correct_predictions / probabilities.shape(0);
}

size_t FitClassifierRandomForest::read_samples(
	const std::string& data_dir,
	const std::string& classification_dir,
	SampleSet& sample_set
) {
	const size_t n_samples_old = sample_set.samples.size();

	std::vector<std::string> classification_files;
	polyvec::os::list_txt_files_recursive(classification_files, classification_dir.c_str());	

	for (auto& file : classification_files)
	{
		std::ifstream f(file);
		if (!f.good())
		{
			std::cout << "Cannot open classification file: " << file << std::endl;
			continue;
		}

		std::string model_info;
		int label;
		while (f >> model_info >> label)
		{
			sample_set.samples.emplace_back(Sample::read_from_data(data_dir, model_info, label));
			
#if DUPLICATE_SAMPLES_REORDERED
			auto dup = sample_set.samples.back();
			std::swap(dup.fits[1], dup.fits[2]);
			sample_set.samples.push_back(dup);
#endif
		}

		f.close();
	}	

	prepare_samples(sample_set, sample_set.samples.begin() + n_samples_old, true, data_dir);

	return sample_set.samples.size() - n_samples_old;
}

void FitClassifierRandomForest::prepare_samples(
	SampleSet& sample_set, 
	std::vector<Sample>::iterator begin, 
	bool remove_by_regularities,
	const std::string& data_dir)
{
	//Load the polygon data
	struct PolygonInfo
	{
		size_t index = (size_t)-1; //index in sample set
		std::vector<PolygonCornerInfo> corner_data;
		std::unordered_map<size_t, TangentFitType> corners_affected_by_regularities;
	};
	std::map<std::string, PolygonInfo> model_to_polygon_info;
	size_t removed_due_to_regularities = 0;
	for (auto it = begin; it != sample_set.samples.end();)
	{
		auto& sample = *it;

		auto last_slash = sample.model.find_last_of("/\\");
		auto model_dir = data_dir + sample.model.substr(0, last_slash + 1);

		auto& map_entry = model_to_polygon_info[model_dir];
		if (map_entry.index == (size_t)-1)
		{
			std::cout << "Reconstructing model " << model_dir << std::endl;
			//we have not encountered this model before, load all the data
			auto polygon_file = model_dir + "polygon.txt";
			auto input_ser = model_dir + "input_ser.txt";
			auto polygon_ser = model_dir + "polygon_ser.txt";

			map_entry.corner_data = FitClassifier::PolygonCornerInfo::read(polygon_file);
			auto input_ser_stream = std::ifstream(input_ser);
			if (!input_ser_stream.good())
				throw std::runtime_error("Cannot open input serialization");
			auto input_data = InputData::deserialize(input_ser_stream);
			auto polygon_ser_stream = std::ifstream(polygon_ser);
			if (!polygon_ser_stream.good())
				throw std::runtime_error("Cannot open polygon serialization");
			auto polygon_data = PolygonData::deserialize(polygon_ser_stream, input_data);
			sample_set.polygons.push_back(polygon_data);
			map_entry.index = sample_set.polygons.size() - 1;

#if REMOVE_SAMPLES_AFFECTED_BY_REGULARITIES
			if (remove_by_regularities)
			{
				auto filtered_regs = polygon_data.RE;
				filtered_regs.clear_symmetries();
				//filtered_regs.clear_important_edges();			
				auto regularity_actions = get_regularity_actions<HandledRegularityListOfChanges>(polygon_data.PP, filtered_regs, std::vector<TangentFitType>(), -1);
				for (auto& action : regularity_actions)
					map_entry.corners_affected_by_regularities[action.first] = action.second;
			}
#endif
		}

#if REMOVE_SAMPLES_AFFECTED_BY_REGULARITIES
		if (remove_by_regularities)
		{
			auto affected_entry = map_entry.corners_affected_by_regularities.find(sample.corner_idx);
			if (affected_entry != map_entry.corners_affected_by_regularities.end())
			{
				std::cout << "Removed " << sample.model << std::endl;
				it = sample_set.samples.erase(it);
				++removed_due_to_regularities;
				continue;
			}
		}
#endif

		sample.corner = map_entry.corner_data.at(sample.corner_idx);
		sample.polygon = map_entry.index;
		++it;
	}

	std::cout << "Removed " << removed_due_to_regularities << " samples due to regularities." << std::endl;
}

void FitClassifierRandomForest::train(
	SampleSet& sample_set,
	const std::string& model_uri
) {

	if (sample_set.samples.empty()) {
		std::cerr << "Cannot train Random Forests, no samples loaded." << std::endl;
		throw std::runtime_error("Cannot train Random Forests, no samples loaded.");
	}

	std::random_shuffle(sample_set.samples.begin(), sample_set.samples.end());

	const size_t training_samples = sample_set.samples.size() / 2;
	const size_t shape[] = { training_samples, feature_description.features.size() };
	andres::Marray<double> features(shape, shape + 2);
	andres::Marray<unsigned char> labels(shape, shape + 1);

#if OUTPUT_TRAINING_DATA
	std::ofstream feature_file("features.txt");
	std::ofstream label_file("labels.txt");
	std::ofstream math_file("features_mathematica.txt");

	std::ofstream training_set("training_set.txt");
	std::ofstream test_set("test_set.txt");

	math_file << "trainingSet={";
#endif
	for (int i = 0; i < training_samples; ++i)
	{
		const auto& sample = sample_set.samples[i];
		for (int j = 0; j < feature_description.features.size(); ++j)
			features(i, j) = feature_description.features[j](sample, sample_set.polygons.at(sample.polygon).PP, sample.corner_idx);
		labels(i) = sample.label;

#if OUTPUT_TRAINING_DATA
		if (i != 0)
			math_file << ",";
		math_file << "{";
		for (int j = 0; j < shape[1]; ++j)
		{
			if (j != 0)
			{
				feature_file << " ";
				math_file << ",";
				training_set << " ";
			}
			feature_file << features(i, j);
			math_file << features(i, j);
			training_set << features(i, j);
		}
		math_file << "}";
		feature_file << "\n";
		label_file << int(labels(i)) << "\n";

		training_set << " " << int(labels(i)) << "\n";

#endif
	}
#if OUTPUT_TRAINING_DATA
	math_file << "}";

	for (int i = training_samples; i < sample_set.samples.size(); ++i)
	{
		const auto& sample = sample_set.samples[i];
		for (int j = 0; j < feature_description.features.size(); ++j)
		{
			if (j != 0)
				test_set << " ";
			test_set << feature_description.features[j](sample, sample_set.polygons.at(sample.polygon).PP, sample.corner_idx);
		}
		test_set << " " << int(sample.label) << "\n";		
	}

	feature_file.close();
	label_file.close();
	math_file.close();
	training_set.close();
	test_set.close();
#endif

	std::cout << "Learning forest with " << training_samples << " samples" << std::endl;

	const size_t number_of_trees = 100;
	forest.learn(features, labels, number_of_trees);

	{
		//get feature statistics
		std::vector<size_t> features_used(feature_description.features.size());
		size_t total_features = 0;
		size_t leaves = 0;
		for (size_t iTree = 0; iTree < forest.size(); ++iTree)
		{
			auto& tree = forest.decisionTree(iTree);
			for (size_t iNode = 0; iNode < tree.size(); ++iNode)
			{
				auto& node = tree.decisionNode(iNode);
				if (node.isLeaf())
					++leaves;
				else
				{
					++features_used[node.featureIndex()];
					++total_features;
				}
			}
		}
		std::cout << "Avg leaves per tree: " << (double)leaves / forest.size() << std::endl;
		std::cout << "Feature statistics: " << std::endl;
		for (int i = 0; i < features_used.size(); ++i)
		{
			std::cout << "Feature " << i << ": " << features_used[i] << " / " << total_features << " (" << 100.0 * features_used[i] / total_features << " %)" << std::endl;
		}
	}

	std::cout << "Serializing forest to " << model_uri << std::endl;
	std::ofstream forest_file(model_uri);
	if (!forest_file.good())
		std::cout << "Cannot open forest file. Skipping serialization." << std::endl;
	else
	{
		forest.serialize(forest_file);
		forest_file.close();
	}

	//validate		
	const size_t shape_validation[] = { sample_set.samples.size() - training_samples, feature_description.features.size() };
	const size_t shape_probabilities_valid[] = { shape_validation[0], 2 };
	const size_t shape_probabilities_train[] = { shape[0], 2 };	
	andres::Marray<double> probabilities_train(shape_probabilities_train, shape_probabilities_train + 2);
	forest.predict(features, probabilities_train);
	std::cout << "Ratio of correct predictions on training data set: " << calculate_prediction_accuracy(probabilities_train, sample_set.samples, 0) << std::endl;
	if (shape_validation[0] > 0)
	{
		std::cout << "Validating with " << shape_validation[0] << " samples" << std::endl;
		andres::Marray<double> validation_samples(shape_validation, shape_validation + 2);
		andres::Marray<double> probabilities_valid(shape_probabilities_valid, shape_probabilities_valid + 2);

		for (size_t i = training_samples; i < sample_set.samples.size(); ++i)
		{
			auto& sample = sample_set.samples[i];
			for (int j = 0; j < feature_description.features.size(); ++j)
				validation_samples(i - training_samples, j) = feature_description.features[j](sample, sample_set.polygons.at(sample.polygon).PP, sample.corner_idx);
		}
		forest.predict(validation_samples, probabilities_valid);
		std::cout << "Ratio of correct predictions on validation data set: " << calculate_prediction_accuracy(probabilities_valid, sample_set.samples, training_samples) << std::endl;
	}
}

FitClassifierRandomForest::SampleSet FitClassifierRandomForest::train(
	const std::string& data_dir, 
	const std::string& classification_dir,
	const std::string& model_uri
)
{
	SampleSet samples;
	read_samples(data_dir, classification_dir, samples);
	train(samples, model_uri);
	return samples;
}

void FitClassifierRandomForest::load_from_file(const std::string& file)
{
	std::ifstream f(file);
	if (!f.good())
		throw std::runtime_error("Cannot open forest file \"" + file + "\"");
	forest.deserialize(f);
	f.close();
}

void FitClassifierRandomForest::load_from_default_file()
{
	load_from_file(DEFAULT_FOREST_FILE);	
}

double polyvec::FitClassifierRandomForest::evaluate_fit_probability(const FitClassifier::FitInfo * fits, const FitClassifier::PolygonCornerInfo & corner, const polyfit::mat2x& polygon, size_t polygon_corner) const
{
	const size_t shape_features[] = { 1, feature_description.features.size() };
	const size_t shape_probabilities[] = { 1, 2 };
	andres::Marray<double> features(shape_features, shape_features + 2);
	andres::Marray<double> probabilities(shape_probabilities, shape_probabilities + 2);
	Sample sample(fits, corner);
	for (int j = 0; j < feature_description.features.size(); ++j)
		features(0, j) = feature_description.features[j](sample, polygon, polygon_corner);
	forest.predict(features, probabilities);
	return probabilities(0, 1);
}

double polyvec::FitClassifierGeometric::evaluate_fit_probability(const FitClassifier::FitInfo * fits, const FitClassifier::PolygonCornerInfo & corner, const polyfit::mat2x& polygon, size_t polygon_corner) const
{	
	double prob = 0;
	for (int i = 0; i < 3; ++i)
	{
		double curvature_good_probability = polyvec::Num::smooth_probability_incr(fits[i].error.curvature.r_min, 1.25, 1.75);
		double accuracy_good_probability = polyvec::Num::smooth_probability_decr(
			fits[i].error.accuracy.max_error(), 
			corner.accuracy.max_error(), 
			std::min(1.0, corner.accuracy.max_error() + 0.4));

		prob += curvature_good_probability * accuracy_good_probability;
	}
	return prob / 3.0;
}

FitClassifierConstant::FitClassifierConstant(double constant_probability)
	: constant_probability(constant_probability)
{ }

double FitClassifierConstant::evaluate_fit_probability(const FitClassifier::FitInfo* fits, const FitClassifier::PolygonCornerInfo& corner, const polyfit::mat2x& polygon, size_t polygon_corner) const
{
	return constant_probability;
}