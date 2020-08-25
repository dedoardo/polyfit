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

using namespace polyfit;
using namespace polyvec;
using namespace Eigen;

int main(int argc, char* argv[]) {
    if (argc < 3) {
        return EXIT_FAILURE;
    }

    const char* image_uri = argv[1];
    const char* boundary_uri = argv[2];

    std::vector<Matrix2Xi> boundaries;
    std::vector<Vector4d>  colors;
    polyfit::mc::get_pixel_regions_from_bmp(image_uri, boundaries, colors);

    polyfit::mc::RasterImageConnectivity raster_conn = polyfit::mc::RasterImageConnectivity::build(boundaries);

    std::vector<int> regions = raster_conn.regions_sorted_by_area();

    Matrix2Xi P = raster_conn.get_polygon_points(regions[1]);
    FILE* fp = fopen(boundary_uri, "w");
    for (Index i = 0; i < P.cols(); ++i) {
        fprintf(fp, "%.1f %.1f\n", (double)P.col(i)(0), (double)P.col(i)(1));
    }
    fclose(fp);
    return EXIT_SUCCESS;
}