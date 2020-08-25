// Polyvec
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include "drawing.hpp"
#include "training.hpp"

// libc++
#include <fstream>

using namespace polyfit;
using namespace polyvec;

const char* usage =
"Required parameters: <samples file> <image dir> <output dir>\n"
"<samples file>       path to a samples file generated while training.\n"
"<image dir>          must contain the original image files, used for visualization"
"<output dir>         the visualizations will be written to this directory"
;

struct ParsedModel
{
	FitClassifierRandomForest::Sample sample;
	
	std::string model;
	std::string model_dir;
	std::string data_dir;
	int resolution;
	int type;
	int corner_number;

	ParsedModel() : resolution(-1), type(-1), corner_number(-1) { }
	ParsedModel(const polyvec::FitClassifierRandomForest::Sample& sample, const std::string& data_dir)
		: sample(sample), data_dir(data_dir)
	{
		auto last_slash = sample.model.find_last_of("/\\");
		auto model_name = sample.model.substr(0, last_slash);
		
		auto first_dash = model_name.find_first_of('-');
		model = model_name.substr(0, first_dash);

		auto first_slash = sample.model.find_first_of("/\\");
		model_dir = sample.model.substr(0, first_slash);

		resolution = std::stoi(model_name.substr(first_dash + 1));
		type = extract_number(sample.model, "type");
		corner_number = extract_number(sample.model, "corner");
	}

	bool operator<(const ParsedModel& other) const
	{
		if (model != other.model)
			return model < other.model;
		if (resolution != other.resolution)
			return resolution < other.resolution;
		if (type != other.type)
			return type < other.type;
		return corner_number < other.corner_number;
	}

	bool is_same_model(const ParsedModel& other) const
	{
		return model == other.model && resolution == other.resolution && type == other.type;
	}
};

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr, usage);
		return EXIT_FAILURE;
	}

	init_pipeline();

	FitClassifierRandomForest classifier;
	
	const std::string image_dir = std::string(argv[2]) + "/";
	const std::string write_dir = std::string(argv[3]) + "/";

	os::make_dir(write_dir);

	int model_number = -1;

	int found_samples = 0;
	int total_samples = 0;

	std::vector<std::pair<double, double>> acc_errors;
	std::vector<std::pair<double, double>> curv_errors;

	FitClassifierRandomForest::SampleSet samples;
	std::vector<ParsedModel> models;
	load_samples_file(argv[1], samples, [&](const FitClassifierRandomForest::Sample& sample, const std::string& data_dir)
	{
		models.emplace_back(sample, data_dir);
	});

	std::sort(models.begin(), models.end());

	PolygonData polygon_data;
	CurvePrimitiveSequence fits[3];

	std::vector<std::vector<int>> corner_labels;		

	std::vector<CurveFitError> corner_errors; //three entries per corner

	ParsedModel last_model;
	bool has_valid_model = false;

	auto start_new_model = [&](const ParsedModel& model)
	{
		srand(1);
		std::string image_uri = image_dir + model.model + "/" + std::to_string(model.resolution) + ".png";
		auto input_ser = model.data_dir + model.model_dir + "/input_ser.txt";
		auto polygon_ser = model.data_dir + model.model_dir + "/polygon_ser.txt";
		auto input_data = InputData::deserialize(std::ifstream(input_ser));
		polygon_data = PolygonData::deserialize(std::ifstream(polygon_ser), input_data);

		if (input_data.R.size() == 0)
			throw std::runtime_error("Cannot load input from " + input_ser);
		if (polygon_data.P.size() == 0)
			throw std::runtime_error("Cannot load polygon from " + polygon_ser);

		corner_labels.clear();
		corner_labels.resize(polygon_data.PP.cols());

		corner_errors.clear();
		corner_errors.resize(3 * polygon_data.PP.cols());

		CurveSequenceFitter fitter(polygon_data.B, polygon_data.PP, polygon_data.PV, polygon_data.RE);
		fits[0] = fitter.all_with_type(TangentFitType(model.type), true, false, false);
		fits[1] = fitter.all_with_type(TangentFitType(model.type), true, true, false);
		fits[2] = fitter.all_with_type(TangentFitType(model.type), true, false, true);

		for (int i = 0; i < 3; ++i)
			for (auto& p : fits[i].primitives)
				corner_errors[3 * p.corner + i].combine(p.error);
	};

	auto finish_model = [&]()
	{
		if (!has_valid_model)
			return;

		// Count how many pictures we need
		size_t pictures = 0;
		for (auto& e : corner_labels)
			pictures = std::max(pictures, e.size());

		for (int j = 0; j < pictures; ++j)
		{
			++model_number;
			for (int i = 0; i < 3; ++i)
			{
				const std::string fit_uri = misc::sfmt("%s/model-%i_fit-%i.svg", write_dir.c_str(), model_number, i);
				DevicePDF* pdf = new DevicePDF(fit_uri.c_str(), 1, 1);

				draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5));
				draw_raster(polygon_data.B);
				draw_raster(polygon_data.PP);

				for (int ip = 0; ip < fits[i].primitives.size(); ++ip)
				{
					auto& p = fits[i].primitives[ip];
					if (j >= corner_labels[p.corner].size())
						continue;
												
					auto this_color = AlternatingCornerColorFunctor()(ip, p);
					draw::curve(p.curve->get_curve().get(), 7.5, this_color);

					auto pos = p.curve->get_curve()->pos(0.3);
					draw::text(pos, misc::sfmt("C %i, Max Acc Error: %.3f, Min Curv R: %.3f", p.corner, p.error.accuracy.max_error(), p.error.curvature.r_min), draw::font_pdf * 0.5, Style::text());
				}					

				pdf->draw(0, 0);
				delete pdf;
			}

			const std::string uri = misc::sfmt("%s/model-%i_classification.svg", write_dir.c_str(), model_number);
			DevicePDF* pdf = new DevicePDF(uri.c_str(), 1, 1);

			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5));
			draw_raster(polygon_data.B);
			draw_raster(polygon_data.PP);

			for (int iCorner = 0; iCorner < polygon_data.PP.cols(); ++iCorner)
			{
				polyvec::real3 color;
				if (j >= corner_labels[iCorner].size())
					color = colors::dark_gray;
				else
				{
					switch (corner_labels[iCorner][j])
					{
					case 0:
						color = colors::red;
						break;
					case 1:
						color = colors::forest_green;
						break;
					}
				}


				auto iPrev = (iCorner - 1 + polygon_data.PP.cols()) % polygon_data.PP.cols();
				auto iNext = (iCorner + 1 + polygon_data.PP.cols()) % polygon_data.PP.cols();

				auto corner = polygon_data.PP.col(iCorner);
				auto midPrev = 0.5 * (polygon_data.PP.col(iPrev) + polygon_data.PP.col(iCorner));
				auto midNext = 0.5 * (polygon_data.PP.col(iCorner) + polygon_data.PP.col(iNext));

				draw::line(midPrev, corner, Style::outline(color, 30.));
				draw::line(corner, midNext, Style::outline(color, 30.));
			}

			pdf->draw(0, 0);
			delete pdf;
		}
	};				

	for (auto it = models.begin(); it != models.end(); ++it)
	{
		if (!it->is_same_model(last_model))
		{
			finish_model();
			try
			{
				has_valid_model = false;
				start_new_model(*it);
				has_valid_model = true;
				last_model = *it;
			}
			catch (std::exception& x)
			{
				std::cout << "Error: " << x.what() << std::endl;
			}
		}
		if (!has_valid_model)
			continue;
		bool same_measures = true;
		for (int i = 0; i < 3; ++i)
		{
			if (std::abs(corner_errors[3 * it->corner_number + i].accuracy.max_error() - it->sample.fits[i].error.accuracy.max_error()) > 0.001)
			{
				std::cout << "Expected accuracy " << it->sample.fits[i].error.accuracy.max_error() << " for " << it->sample.model << " at corner " << it->corner_number << " in fit " << i << " but got " << corner_errors[3 * it->corner_number + i].accuracy.max_error() << std::endl;
				same_measures = false;
				acc_errors.emplace_back(it->sample.fits[i].error.accuracy.max_error(), corner_errors[3 * it->corner_number + i].accuracy.max_error());
			}
			if (std::abs(corner_errors[3 * it->corner_number + i].curvature.r_min - it->sample.fits[i].error.curvature.r_min) > 0.001)
			{
				std::cout << "Expected curvature " << it->sample.fits[i].error.curvature.r_min << " for " << it->sample.model << " at corner " << it->corner_number << " in fit " << i << " but got " << corner_errors[3 * it->corner_number + i].curvature.r_min << std::endl;
				same_measures = false;
				curv_errors.emplace_back(it->sample.fits[i].error.curvature.r_min, corner_errors[3 * it->corner_number + i].curvature.r_min);
			}
		}

		if (true || same_measures)
		{
			corner_labels[it->corner_number].push_back(it->sample.label);
			++found_samples;
		}
		++total_samples;
	}
	finish_model();
	
	std::cout << "Encountered " << (total_samples - found_samples) << " unmatched samples (" << (100.0 * (total_samples - found_samples) / total_samples) << " %)" << std::endl;

	{
		std::ofstream acc("acc.csv");
		for (auto& p : acc_errors)
			acc << p.first << "," << p.second << std::endl;
	}
	{
		std::ofstream curv("curv.csv");
		for (auto& p : curv_errors)
			curv << p.first << "," << p.second << std::endl;
	}

	return EXIT_SUCCESS;
}