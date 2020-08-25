/*
	Generates text file for specific parts of the pipeline.
*/
const char* usage =
	"./polygon_baseline <type> <image> <output>\n"
	"<type> which stage of the pipeline should be exported\n"
	"<image> input image file\n"
	"<output> where the output file should be written\n";

#define STAGE_BOUNDARY "boundary"
#define STAGE_POLYGON "polygon"
#define STAGE_CURVES "curves"
#define STAGE_CORNERS "corners"

// libc++
#include <cstdlib>
#include <cstdio>

// polyvec
#include <polyvec/pipeline_helper.hpp>

// namespaces
using namespace Eigen;
using namespace polyvec;

// Pipeline stage outputs
InputData   input_data;
PolygonData polygon_data;
SplineData  spline_data;

int write_boundary(const char* uri);
int write_polygon(const char* uri);
int write_curves(const char* uri);
int write_corners(const char* uri);

int main(int argc, char* argv[]) {
	if (argc < 4) {
		fprintf(stderr, "%s", usage);
		return EXIT_FAILURE;
	}

	const char* stage = argv[1];
	const char* image_uri = argv[2];
	const char* output_uri = argv[3];

	const char* stages[] = {
		STAGE_BOUNDARY,
		STAGE_POLYGON,
		STAGE_CURVES,
		STAGE_CORNERS
	};

	bool is_stage_valid = false;
	for (int i = 0; i < PF_ARRAY_LEN(stages); ++i) {
		if (strcmp(stage, stages[i]) == 0) {
			is_stage_valid = true;
			break;
		}
	}

	if (!is_stage_valid) {
		fprintf(stderr, "invalid stage %s, supported ones are:\n", stage);
		for (int i = 0; i < PF_ARRAY_LEN(stages); ++i) {
			fprintf(stderr, "%s\n", stages[i]);
		}
		return EXIT_FAILURE;
	}

	init_pipeline();

	input_data = load_input(image_uri);
	polygon_data = extract_polygon(input_data);

	if (strcmp(stage, STAGE_BOUNDARY) == 0) {
		return write_boundary(output_uri);
	}
	else if (strcmp(stage, STAGE_POLYGON) == 0) {
		return write_polygon(output_uri);
	}
	else if (strcmp(stage, STAGE_CURVES) == 0) {
		return write_curves(output_uri);
	}
	else if (strcmp(stage, STAGE_CORNERS) == 0) {
		return write_corners(output_uri);
	}
	else {
		fprintf(stderr, "Stage %s not implemented\n", stage);
		return EXIT_FAILURE;
	}
}

int write_boundary(const char* uri) {
	return EXIT_SUCCESS;
}

int write_polygon(const char* uri) {
	FILE* fp = fopen(uri, "w");
	if (!fp) {
		return EXIT_FAILURE;
	}

	for (Index i = 0; i < polygon_data.PV.size(); ++i) {
		fprintf(fp, "%d\n", (int)polygon_data.PV[i]);
	}

	fclose(fp);
	return EXIT_SUCCESS;
}

int write_curves(const char* uri) {
	return EXIT_SUCCESS;
}

int write_corners(const char* uri) {
	FILE* fp = fopen(uri, "w");
	if (!fp) {
		return EXIT_FAILURE;
	}

	for (size_t i = 0; i < spline_data.tangent_fits.size(); ++i) {
		fprintf(fp, "%d\n", (int)spline_data.tangent_fits[i]);
	}

	fclose(fp);
	return EXIT_SUCCESS;
}