#include <iostream>
//
#include <polyvec/api.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>
//
#include <dev/pdf.hpp>
//
#include <polyvec/mc/get_bounding_box.hpp>
#include <polyvec/mc/get_pixel_regions.hpp>
#include <polyvec/mc/junction.hpp>
#include <polyvec/mc/random_color.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
//
#include <Eigen/Core>
#include <Eigen/Geometry>

NAMESPACE_BEGIN()
void draw_raster_closed(const Eigen::Matrix2Xd& raster,
                        const Eigen::Vector4d& color) {
  ::polyvec::vec<::polyvec::real2> points;
  for (int i = 0; i < raster.cols(); ++i) {
    points.emplace_back(raster.col(i));
  }

  ::polyvec::draw::polygon(points, ::polyvec::Style::fill(color.head<3>()));
}  // end of draw_raster_closed()
NAMESPACE_END()

NAMESPACE_BEGIN(polyvectest)
NAMESPACE_BEGIN(mc)

int junctions(int, char**) {
  std::vector<std::string> image_addrs;
  image_addrs.emplace_back(POLYVEC_TEST_PATH "/mc/abc.png");
  image_addrs.emplace_back(POLYVEC_TEST_PATH "/mc/holes001.png");
  // image_addrs.emplace_back(POLYVEC_TEST_PATH "/mc/chess.png");

  printf("==== Testing raster image connectivity \n");

  polyfit::Log::open("test_dump/log-export-multicolor-traces.txt");

  for (int i = 0; i < (int)image_addrs.size(); ++i) {
    std::vector<Eigen::Matrix2Xi> boundaries;
    std::vector<Eigen::Vector4d> colors;
    polyfit::mc::get_pixel_regions_from_bmp(image_addrs[i], boundaries, colors);

    // Get the bounding box of the whole picture
    Eigen::AlignedBox<int, 2> bbox = polyfit::mc::get_bounding_box(boundaries);

    // Can you build a connectivity
    polyfit::mc::RasterImageConnectivity raster_conn =
        polyfit::mc::RasterImageConnectivity::build(boundaries);
    const bool is_mesh_okay =
        raster_conn.connectivity().check_sanity_slowly(true);
    printf("Mesh okay status: %d \n", (int)is_mesh_okay);

    std::unique_ptr<::polyvec::DevicePDF> pdf;

    //
    // Draw the image by random colors
    //
    pdf.reset(new ::polyvec::DevicePDF(
        ::polyvec::misc::sfmt("test_dump/junctions%02d_%04d.pdf", i, 0).c_str(),
        1, 1));
    //
    std::vector<int> by_area = raster_conn.regions_sorted_by_area();
    //
    for (int j = 0; j < (int)raster_conn.n_regions(); ++j) {
      const int region_id = by_area[j];
      std::vector<int> holes;
      int outer_boundary;
      raster_conn.region_to_polygon(region_id, outer_boundary, holes);
      assert_break(outer_boundary == region_id);

      draw_raster_closed(
          raster_conn.get_polygon_points(outer_boundary).cast<double>(),
          j == 0 ? Eigen::Vector4d(1, 1, 1, 1) : ::polyfit::mc::random_color());

      for (int h : holes) {
        Eigen::Matrix2Xi bdry = raster_conn.get_polygon_points(h);
        bdry.rowwise().reverseInPlace();
        draw_raster_closed(bdry.cast<double>(), Eigen::Vector4d(1, 1, 1, 1));
      }
    }  // end of regions
    // pdf->draw(0, 0);

    //
    // Now get the junction types and draw them
    //

    std::vector<polyfit::mc::JunctionType> junction_types;
    polyfit::mc::find_junctions(raster_conn, junction_types);

    // Draw the junctions
    for (int vid = 0; vid < raster_conn.n_pixel_vertices(); ++vid) {
      if (junction_types[vid] == ::polyfit::mc::JUNCTION_INVALID) {
        continue;
      }

      // Now draw
      Eigen::Vector2i pos = raster_conn.pixel_vertex_index_to_pos(vid);
      auto __draw = [&](const std::string& ss) {
        ::polyvec::draw::text(pos.cast<double>(), ss.c_str(),
                              ::polyvec::draw::font_pdf,
                              ::polyvec::Style::fill(::polyvec::colors::black));
      };

      switch (junction_types[vid]) {
        case polyfit::mc::JUNCTION_3_SHARP:
          __draw("3");
          break;
        case polyfit::mc::JUNCTION_4_SHARP:
          __draw("4");
          break;
        case polyfit::mc::JUNCTION_KOPF:
          __draw("K");
          break;
        default:;
      }
    }  // end of vertices
    pdf->draw(0, 0);
  }  // end of images

  return EXIT_FAILURE;
}

NAMESPACE_END(mc)
NAMESPACE_END(polyvectest)
