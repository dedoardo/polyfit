// Polyvec
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/eigen_serialization.hpp>
#include "drawing.hpp"

// libc++
#include <fstream>
#include <Windows.h>

using namespace polyfit;
using namespace polyvec;

// Converts training samples to the curves that the current fit produces.

const char* usage =
"Required parameters: <classification dir> <data dir> <output dir data> <output dir labels>\n"
"<data dir>           must contain the training samples that the classifications are based on.\n"
"<output dir data>    the converted samples will be written to this directory\n"
;

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, usage);
		return EXIT_FAILURE;
	}

	const std::string model_uri = argv[1];

	init_pipeline();

	FitClassifierRandomForest classifier;
	
	const std::string data_dir = std::string(argv[1]);
	const std::string write_dir_data = std::string(argv[2]) + "/";

	os::make_dir(write_dir_data);

	std::vector<std::pair<double, double>> acc_errors;
	std::vector<std::pair<double, double>> curv_errors;

	std::vector<std::string> models;


	// Convert fitting data

	WIN32_FIND_DATA find_data;
	HANDLE find_handle;
	std::string search_pattern = data_dir + "/*";
	find_handle = FindFirstFile(search_pattern.c_str(), &find_data);
	if (find_handle == INVALID_HANDLE_VALUE)
		std::cout << "Cannot iterate subfolders in " + data_dir << std::endl;
	else
	{
		do
		{
			std::string filename(find_data.cFileName);
			if (find_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY && filename != "." && filename != "..")
				models.push_back(filename);
		} while (FindNextFile(find_handle, &find_data));
		FindClose(find_handle);
	}

	for (auto& model : models)
	{
		std::cout << "Converting model " << model << std::endl;

		PolygonData polygon_data;
		try
		{
			srand(1);
			auto input_ser = data_dir + "/" + model + "/input_ser.txt";
			auto polygon_ser = data_dir + "/" + model + "/polygon_ser.txt";
			auto input_data = InputData::deserialize(std::ifstream(input_ser));
			polygon_data = PolygonData::deserialize(std::ifstream(polygon_ser), input_data);			
		}
		catch (std::exception& e)
		{
			std::cerr << "Error reading input data for model " << model << ": " << e.what() << std::endl;
			continue;
		}

		auto model_data_dir = data_dir + "/" + model;
		os::make_dir(model_data_dir);

		// Calculate the new fits
		for (int fit_type = TANGENT_FIT_LERP; fit_type <= TANGENT_FIT_LERP_SYM; ++fit_type)
		{
			std::vector<CurveFitError> corner_errors(3 * polygon_data.P.size());
			std::vector<double> dist2corner(3 * polygon_data.P.size(), std::numeric_limits<double>::infinity());
			CurvePrimitiveSequence fits[3];
			CurveSequenceFitter fitter(polygon_data.B, polygon_data.PP, polygon_data.PV, polygon_data.RE);
			fits[0] = fitter.all_with_type(TangentFitType(fit_type), true, false, false);
			fits[1] = fitter.all_with_type(TangentFitType(fit_type), true, true, false);
			fits[2] = fitter.all_with_type(TangentFitType(fit_type), true, false, true);

			for (int i = 0; i < 3; ++i)
				for (auto& p : fits[i].primitives)
				{
					corner_errors[3 * p.corner + i].combine(p.error);

					//find corner dist
					auto corner = polygon_data.PP.col(p.corner);
					double t = p.curve->get_curve()->project(corner);
					double dist = (corner - p.curve->get_curve()->pos(t)).norm();
					if (dist < dist2corner[3 * p.corner + i])
						dist2corner[3 * p.corner + i] = dist;
				}

			std::string write_dir_model_data = write_dir_data + "/" + model;
			os::make_dir(write_dir_model_data);
			os::copy_file(model_data_dir + "/polygon.txt", write_dir_model_data + "/polygon.txt");
			os::copy_file(model_data_dir + "/polygon_ser.txt", write_dir_model_data + "/polygon_ser.txt");
			os::copy_file(model_data_dir + "/input_ser.txt", write_dir_model_data + "/input_ser.txt");
			for (int fix = 0; fix < 3; ++fix)
			{
				int next_primitive = 0;
				for (int i = 0; i < polygon_data.P.size(); ++i)
				{					
					std::stringstream ss;
					ss << write_dir_model_data << "/corner-" << i << "_type-" << fit_type << "_ex-0_fix-";
					ss << (fix & 1 ? "1" : "0");
					ss << (fix & 2 ? "1" : "0");
					ss << "_info.txt";
					std::ofstream corner_txt(ss.str());
					corner_txt << "acc max: " << corner_errors[3 * i + fix].accuracy.max_error() << "\n";
					corner_txt << "acc pos: " << corner_errors[3 * i + fix].accuracy.e_pos << "\n";
					corner_txt << "acc neg: " << corner_errors[3 * i + fix].accuracy.e_neg << "\n";
					corner_txt << "curv: " << corner_errors[3 * i + fix].curvature.r_min << "\n";
					corner_txt << "dist2corner: " << dist2corner[3 * i + fix] << "\n";

					while (next_primitive < fits[fix].primitives.size() && fits[fix].primitives[next_primitive].corner == i)
					{
						auto controls = fits[fix].primitives[next_primitive].curve->get_curve()->get_params();
						corner_txt << Serialized(controls) << "\n";
						++next_primitive;
					}
				}
			}
		}		
	}	

	return EXIT_SUCCESS;
}