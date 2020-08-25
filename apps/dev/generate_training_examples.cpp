const char* usage =
"./generate_training_examples <image:file> <output:directory>";

// Polyvec
#include <polyvec/api.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/core/macros.hpp>
#include <polyvec/io/image.hpp>
#include <polyvec/io/ascii.hpp>
#include <polyvec/io/draw.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/matrix.hpp> 
#include <polyvec/utils/arguments.hpp>
#include <polyvec/visualize/accuracy.hpp>
#include <polyvec/image-segment/image_segment.hpp>
#include <polyvec/polygon-tracer/iterative-global.hpp>
#include <polyvec/curve-tracer/measure_accuracy_polygon.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/pipeline_helper.hpp>

// todo: remove
#include "drawing.hpp"

// libc++
#include <iostream>
#include <functional>
#include <random>

// if defined, deserializes input and polygon data if available instead of re-calculating them
#define USE_PREVIOUS_CALCULATION_IF_AVAILABLE 1

void draw_polygon_data (const PolygonData& polygon, const std::string& uri);

using namespace polyfit;
using namespace std;
using namespace polyvec;

int main(int argc, char* argv[]) {
	const int variants = 1;

	init_pipeline();
	static std::vector<str> image_addrs;

	if (argc < 3) {
		std::cerr << usage << std::endl;
		return EXIT_FAILURE;
	}

	auto image_uri = argv[1];
	std::cout << "Generating training examples for " << image_uri << std::endl;
	const string write_dir = std::string(argv[2]) + "/" + IO::model_id(argv[1]);

	polyvec::InputData input_data;
	polyvec::PolygonData polygon_data;

	bool deserialized_input = false;
	bool deserialized_polygon = false;

	std::string input_ser_path = write_dir + "/input_ser.txt";
	std::string polygon_ser_path = write_dir + "/polygon_ser.txt";

#if USE_PREVIOUS_CALCULATION_IF_AVAILABLE
	if (polyvec::os::file_exists(input_ser_path))
	{
		auto ifs = std::ifstream(input_ser_path);
		input_data = InputData::deserialize(ifs);
		deserialized_input = true;
	}	
#endif

	if (!deserialized_input)
	{
		input_data = load_input(image_uri);
		auto ofs = std::ofstream(input_ser_path);
		input_data.serialize(ofs);
	}
#if USE_PREVIOUS_CALCULATION_IF_AVAILABLE
	if (polyvec::os::file_exists(polygon_ser_path))
	{
		auto ifs = std::ifstream(polygon_ser_path);
		polygon_data = PolygonData::deserialize(ifs, input_data);
		deserialized_polygon = true;
	}
#endif
	if (!deserialized_polygon)
	{
		polygon_data = extract_polygon(input_data);
		auto ofs = std::ofstream(polygon_ser_path);
		polygon_data.serialize(ofs);
	}

	// preparing output directory
	os::make_dir(write_dir);

	CurveSequenceFitter fitter(polygon_data.B, polygon_data.PP, polygon_data.PV, polygon_data.RE);
	
	std::mt19937 rnd;
	std::uniform_real_distribution<> prob_dist(0.0, 1.0);
	std::uniform_int_distribution<> size_dist(0, 5);

	auto random_type = [&]() 
	{
		auto v = prob_dist(rnd);
		if (v <= 0.6)
			return TANGENT_FIT_LERP;
		else if (v <= 0.9)
			return TANGENT_FIT_LERP_SYM;
		else
			return TANGENT_FIT_CONSTANT;
	};	

	const string polygon_uri = misc::sfmt("%s/polygon.txt", write_dir.c_str());
	FitClassifier::PolygonCornerInfo::write(polygon_data.B, polygon_data.P, polygon_data.PP, polygon_uri);

	const string polygon_data_uri = misc::sfmt("%s/polygon.svg", write_dir.c_str());
	draw_polygon_data(polygon_data, polygon_data_uri.c_str());

	for(int corner = 0; corner < polygon_data.PP.cols(); ++corner)
		for (int type = TANGENT_FIT_LERP; type <= TANGENT_FIT_LERP_SYM; ++type)
			for (int i = 0; i < variants; ++i)
			{
				std::vector<std::pair<int, TangentFitType>> corner_context;

				// generate a random context for the curve
				auto add_before = size_dist(rnd);
				auto add_after = size_dist(rnd);
				if (i == 0)
					add_before = add_after = 0;
				for (int j = 0; j < add_before; ++j)
					corner_context.emplace_back(std::make_pair<int, TangentFitType>((corner - 1 - add_before + j + polygon_data.PP.cols()) % polygon_data.PP.cols(), random_type()));
				corner_context.emplace_back(corner, TangentFitType(type));
				for (int j = 0; j < add_after; ++j)
					corner_context.emplace_back(std::make_pair<int, TangentFitType>((corner + 1 + j) % polygon_data.PP.cols(), random_type()));

				// Edo: The context is not being used anymore?
				const string base_uri = misc::sfmt("%s/corner-%i_type-%i_ex-%i", write_dir.c_str(), corner, type, i);
				const string variant_uri = misc::sfmt("%s_context.txt", base_uri.c_str());
				std::ofstream context_file(variant_uri);
				for (auto& c : corner_context)
					context_file << c.first << " " << c.second << "\n";
				context_file.close();

				for (int fix0 = 0; fix0 < 2; ++fix0)
					for (int fix1 = 0; fix1 < 2; ++fix1)
					{
						const string uri = misc::sfmt("%s_fix-%i%i.svg", base_uri.c_str(), fix0, fix1);
						DevicePDF* pdf = new DevicePDF(uri.c_str(), 1, 1, true);

						draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
						draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5));
						draw_raster(polygon_data.B);
						draw_raster(polygon_data.PP, Style::outline(colors::black, 10.0));

						auto seq = fitter.fit_individual_corners(corner_context, false, fix0 == 1, fix1 == 1);

						CurveFitError error;
						geom::aabb bbox;
						BezierCurve* bezier = nullptr;
						for (auto& p : seq.primitives)
						{
							if (p.corner != corner)
								continue;
							error.combine(p.error);
							draw::curve(p.curve->get_curve().get(), 15, colors::talking_orange);
							auto pbbox = p.curve->get_curve()->get_bounding_box();
							bbox.add(pbbox.min);
							bbox.add(pbbox.max);
							if (!bezier)
								bezier = dynamic_cast<BezierCurve*>(p.curve->get_curve().get());
						}

						auto corner_t = bezier->project(polygon_data.PP.col(corner));

						const string fit_info_uri = misc::sfmt("%s_fix-%i%i_info.txt", base_uri.c_str(), fix0, fix1);
						FitClassifier::FitInfo info(error, (bezier->pos(corner_t) - polygon_data.PP.col(corner)).norm());
						info.write_to_file(fit_info_uri);

						bbox.add(polygon_data.PP.col(corner));
						bbox.add(polygon_data.PP.col((corner + 1) % polygon_data.PP.cols()));
						bbox.add(polygon_data.PP.col((corner - 1 + polygon_data.PP.cols()) % polygon_data.PP.cols()));

						pdf->draw(0, 0, bbox);
						//pdf->draw(0, 0);
						delete pdf;
					}
			}

	return EXIT_SUCCESS;
}

void draw_polygon_data(const PolygonData& polygon, const std::string& uri) {
	DevicePDF* pdf = new DevicePDF(uri.c_str(), 1, 1);

	draw_raster_background(polygon.B, Style::outline(colors::gray * 1.75, 2.5));
	draw_raster(polygon.PP);
	draw_raster_indices(polygon.PP);

	pdf->draw(0, 0);
	delete pdf;
}