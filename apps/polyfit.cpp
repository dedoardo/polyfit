// Polyvec
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/system.hpp>
#include <polyvec/core/log.hpp>

#include <polyvec/polygon-tracer/multi_polygon.hpp>
#include <polyvec/mc/get_bounding_box.hpp>
#include <polyvec/mc/get_pixel_regions.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
#include <polyvec/mc/junction.hpp>
#include <polyvec/mc/segment_connectivity.hpp>

#include "dev/drawing.hpp"
#include "dev/draw_paper_figures.hpp"
#include <polyvec/io/svg.hpp>

#define EXPORT_DEBUG_SVG 0

// libc++
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;
using namespace polyfit;
using namespace polyvec;
using namespace Eigen;

namespace {
    void draw_all_boundaries(polyfit::mc::RasterImageConnectivity& conn) {
        for (int j = 0; j < conn.n_regions(); ++j) {
            std::vector<int> holes;
            int outer_boundary;
            conn.region_to_polygon(j, outer_boundary, holes);
            draw::polyline(conn.get_polygon_points(outer_boundary).cast<double>(), Style::outline(colors::gray, .75));
        }
    }
}

const char* usage =
"Usage should be one of the following\n"
"<input_image>                           ~ The output will be written to curves_fill.svg"
"<input_image> <output_svg>              ~ The output will be written to <output_svg>"
"<input_image> <classifier> <output_svg> ~ Same as above, but uses the specified classifier";

int main(int argc, char* argv[]) {
    std::string input_image, classifier, output_svg;
    if (argc == 1) {
        fprintf(stderr, usage);
        return EXIT_FAILURE;
    } else if (argc == 2) {
        input_image = argv[1];
        classifier = "trained/random_forest_paper.txt";
        output_svg = "curves_fill.svg";
    } else if (argc == 3) {
        input_image = argv[1];
        classifier = "trained/random_forest_paper.txt";
        output_svg = argv[2];
    } else if (argc >= 4) {
        input_image = argv[1];
        classifier = argv[2];
        output_svg = argv[3];
        if (argc > 4) {
            fprintf(stderr, "Ignoring extra arguments starting from %s\n", argv[4]);
        }
    }

    Log::open(stdout, Log::CHANNEL_DEV);
    
    PF_DEV_F("Input image: %s", input_image.c_str());
    PF_DEV_F("Output svg: %s", output_svg.c_str());
    PF_DEV_F("Classifier: %s", classifier.c_str());

    std::vector<Matrix2Xi> boundaries;
    std::vector<Vector4d>  colors;
    polyfit::mc::get_pixel_regions_from_bmp(input_image.c_str(), boundaries, colors);
    PF_STATUS_F("boundaries %d", boundaries.size());
    PF_STATUS_F("colors %d", colors.size());

    polyfit::mc::RasterImageConnectivity raster_conn = polyfit::mc::RasterImageConnectivity::build(boundaries);
    std::vector<polyfit::mc::JunctionType> junction_types;
    polyfit::mc::find_junctions(raster_conn, junction_types);
    mc::SegmentConnectivity segment_conn = mc::SegmentConnectivity::build(raster_conn, junction_types);

    MultiPolygonTracer tracer(raster_conn, segment_conn, colors, classifier.c_str());
    tracer.trace();
    std::vector<int> boundaries_by_area = tracer.get_boundaries_ordered_by_area();
    PF_DEV_F("Exporting %llu boundaries", boundaries_by_area.size());

    // -----------------------------------------------------------------------------------------------------------------
    {
        DevicePDF* pdf = new DevicePDF(output_svg.c_str(), 1, 1);
        //draw_raster_background(B, Style::outline(colors::gray * 1.75, 2.5));
        for (size_t i = 0; i < boundaries_by_area.size(); ++i) {

            const auto& boundary = tracer.boundaries()[boundaries_by_area[i]];
            draw_curve_primitives(boundary.spline().primitives, PrimitiveTypeColorFunctor());
            //draw_curve_primitives_closed(boundary.spline().primitives, boundary.color());
            PF_VERBOSE_F("Boundary color %f %f %f %f points %d", boundary.color()(0), boundary.color()(1), boundary.color()(2), boundary.color()(3), boundary.boundary_points().size());
        }
        pdf->draw(0, 0);
        delete pdf;
    }

    return EXIT_SUCCESS;
}