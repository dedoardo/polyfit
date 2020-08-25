const char* usage =
	"Supported arguments:\n"
	"	--image      \n"
	"	--output     \n"
	"	--input      \n"
	"	--polygon    \n"
	"	--types      \n"
	"   --classifier   Serialized model used for classifying corner types. (if not set, the default file will be used)\n"
	;

// Polyvec
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/utils/directions.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/curve-tracer/fit_classifier.hpp>

// todo: remove
#include "drawing.hpp"

// libc++
#include <iostream>

using namespace polyfit;
using namespace std;
using namespace polyvec;

int main(int argc, char* argv[]) {

	std::string image_path;
	std::string output_path;
	std::string input_ser_path;
	std::string polygon_ser_path;
	std::string corner_type_path;
	std::string classifier_ser_path;

	for (int i = 1; i < argc; ++i)
	{
		if (strcmp(argv[i], "--image") == 0)
			image_path = argv[++i];
		else if (strcmp(argv[i], "--output") == 0)
			output_path = argv[++i];
		else if (strcmp(argv[i], "--input") == 0)
			input_ser_path = argv[++i];
		else if (strcmp(argv[i], "--polygon") == 0)
			polygon_ser_path = argv[++i];
		else if (strcmp(argv[i], "--types") == 0)
			corner_type_path = argv[++i];
		else if (strcmp(argv[i], "--classifier") == 0)
			classifier_ser_path = argv[++i];
		else
			std::cerr << "Unrecognized argument " << argv[i] << std::endl;
	}

	if (image_path.empty() && input_ser_path.empty())
	{
		std::cout << "--image argument and --input argument not set" << std::endl;
		return -1;
	}
	if (output_path.empty() && polygon_ser_path.empty())
	{
		std::cout << "--output argument and --polygon argument not set" << std::endl;
		return -2;
	}

	init_pipeline();

	InputData input;
	PolygonData polygon;

	if (input_ser_path.empty())
	{
		input = polyvec::load_input(image_path);
		auto ofs = std::ofstream(output_path + "/input_ser.txt");
		input.serialize(ofs);
	}
	else {
		auto ifs = std::ifstream(input_ser_path);
		input = InputData::deserialize(ifs);
	}

	if (polygon_ser_path.empty())
	{
		polygon = polyvec::extract_polygon(input);
		auto ofs = std::ofstream(output_path + "/polygon_ser.txt");
		polygon.serialize(ofs);
	}
	else {
		auto ifs = std::ifstream(polygon_ser_path);
		polygon = PolygonData::deserialize(ifs, input);
	}

	//draw raster
	if(input_ser_path.empty())
	{
		auto uri = output_path + "/raster.svg";
		auto pdf = new DevicePDF(uri.c_str(), 1, 1, true);

		draw_raster_closed(polygon.B);

		pdf->draw(0, 0);
		delete pdf;
	}

	//draw polygon
	if(polygon_ser_path.empty())
	{
		auto uri = output_path + "/polygon.svg";
		auto pdf = new DevicePDF(uri.c_str(), 1, 1, true);

		draw_raster_closed(polygon.PP);

		/*for (int i = 0; i < polygon.PP.cols(); ++i)
		{
			auto dir_prev = (polygon.PP.col(i) - polyfit::CircularAt(polygon.PP, i - 1)).normalized();
			auto dir_next = (polyfit::CircularAt(polygon.PP, i + 1) - polygon.PP.col(i)).normalized();
			auto normal = polyvec::util::normal_dir(dir_prev + dir_next).normalized();
			auto font_size = 0.8;
			auto pos = polygon.PP.col(i) + 0.5 * normal + Eigen::Vector2d::Constant(-0.5 * font_size);
			draw::text(pos, misc::sfmt("%i", i), font_size, Style::text());
		}*/

		pdf->draw(0, 0);
		delete pdf;
	}

	try
	{		
		polyvec::CurveSequenceFitter fitter(polygon.B, polygon.PP, polygon.PV, polygon.RE, -1, std::vector<bool>(), true);
		std::vector<TangentFitType> corner_types;
		CurvePrimitiveSequence fit;
		if (corner_type_path.empty())
		{
			const string polygon_uri = misc::sfmt("%s/polygon.txt", output_path.c_str());
			FitClassifier::PolygonCornerInfo::write(polygon.B, polygon.P, polygon.PP, polygon_uri);

			FitClassifierRandomForest classifier;
			if (classifier_ser_path.empty()) {
				classifier.load_from_default_file();
			}
			else {
				classifier.load_from_file(classifier_ser_path);
			}

			int it = 0;
			fit = fitter.fit_evolutionary_simple(classifier, corner_types, [&](const polyvec::CurveSequenceFitter::EvolutionaryFittingState& state, const polyvec::CurveSequenceFitter::FittingAttempt *fits)
			{
				for (int corner = 0; corner < polygon.PP.cols(); ++corner)
				{
					const string context_uri = misc::sfmt("%s/corner-%i_type-%i_ex-0_context.txt", output_path.c_str(), corner, it);
					std::ofstream f(context_uri);
					f << corner << " " << it << std::endl;
					f.close();
				}

				for (int i = 0; i < 3; ++i)
				{
					bool fix0 = (i & 1) == 1;
					bool fix1 = ((i >> 1) & 1) == 1;

					auto corner_errors = get_corner_errors(polygon.PP.cols(), fits[i].primitive_seq.primitives);

					for (int corner = 0; corner < polygon.PP.cols(); ++corner)
					{
						const string fit_info_uri = misc::sfmt("%s/corner-%i_type-%i_ex-0_fix-%i%i_info.txt", output_path.c_str(), corner, it, fix0, fix1);
						FitClassifier::FitInfo info(corner_errors[corner], fits[i].distances_to_corner[corner]);
						info.write_to_file(fit_info_uri);
					}
				}

				++it;
			});
		}
		else
		{
			auto file = std::ifstream(corner_type_path);
			corner_types.resize(polygon.PP.cols());
			std::vector<std::pair<int, TangentFitType>> corners(corner_types.size());
			for (int i = 0; i < corner_types.size(); ++i)
			{
				int type;
				file >> type;
				corner_types[i] = TangentFitType(type);
				corners[i] = std::make_pair(i, corner_types[i]);
			}
			fit = fitter.fit_individual_corners(corners, true, false, false);
		}


		//draw fit
		{
			auto uri = output_path + "/fit_closed.svg";
			auto pdf = new DevicePDF(uri.c_str(), 1, 1, true);
			draw_curve_primitives_closed(fit.primitives);
			pdf->draw(0, 0);
			delete pdf;

			uri = output_path + "/fit_outline.svg";
			pdf = new DevicePDF(uri.c_str(), 1, 1, true);
			draw_raster_closed(polygon.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon.B, Style::outline(colors::gray * 1.75, 2.5));
			draw_raster(polygon.B);
			//draw_raster(polygon.PP, Style::outline(colors::black, 10.0));
			draw_curve_primitives(fit.primitives, AlternatingCornerColorFunctor());
			pdf->draw(0, 0);

			auto matrix = pdf->get_matrix();		
			delete pdf;

			std::ofstream corners(output_path + "/fit_outline_corners.txt", std::ios::trunc | std::ios::out);
			for (int i = 0; i < polygon.PP.cols(); ++i)
			{
				corners << corner_types[i] << " ";
				for (auto& prim : fit.primitives)
				{
					if (prim.corner != i)
						continue;
					auto curve = prim.curve->get_curve();
					auto samples = std::max(5, (int)std::ceil((curve->pos(0.0) - curve->pos(1.0)).norm() * 5.0));
					for (int s = 0; s <= samples; ++s)
					{
						auto p = curve->pos(1.0 * s / samples);
						double x = matrix.xx * p.x() + matrix.xy * p.y() + matrix.x0;
						double y = matrix.yx * p.x() + matrix.yy * p.y() + matrix.y0;
						corners << x << " " << y << " ";
					}
				}
				corners << "\n";
			}
			corners.close();
		}
	}
	catch (std::exception& x)
	{
		std::cout << "Fitting error: " << x.what() << std::endl;
		return -3;
	}

	return EXIT_SUCCESS;
}