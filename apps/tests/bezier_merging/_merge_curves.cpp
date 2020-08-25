#include <iostream>
//
#include <polyvec/api.hpp>
#include <polyvec/curve-tracer/bezier_merging.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/geometry/winding_number.hpp>
//
#include <Eigen/Core>
#include <Eigen/Geometry>

// These are only required by the second test, i.e., _potrace_merge()
#include <dev/pdf.hpp>
#include <polyvec/image-segment/image_segment.hpp>
#include <polyvec/io/image.hpp>
#include <polyvec/utils/potrace.hpp>

NAMESPACE_BEGIN()


::polyvec::BezierCurve build_bezier(const Eigen::Matrix<double, 2, 4> ctrl) {
  ::polyvec::BezierCurve bz;
  bz.set_control_points(ctrl);
  return bz;
};




void draw_control_polygon(const std::vector<Eigen::Matrix2Xd>& curves) {
  using namespace polyvec;
  for (int i = 0; i < (int)curves.size(); ++i) {
    if (curves[i].cols() == 4) {
      draw::line(curves[i].col(0), curves[i].col(1),
                 Style::outline(colors::talking_orange, .1));
      draw::line(curves[i].col(1), curves[i].col(2),
                 Style::outline(colors::talking_orange, .1));
      draw::line(curves[i].col(2), curves[i].col(3),
                 Style::outline(colors::talking_orange, .1));
    } else {
      // draw::line(curves[i].col(0), curves[i].col(1),
      //            Style::outline(colors::forest_green, .1));
    }
  }
}

void draw_curves(const std::vector<Eigen::Matrix2Xd>& curves) {
  using namespace polyvec;
  for (int i = 0; i < (int)curves.size(); ++i) {
    if (curves[i].cols() == 4) {
      Eigen::Vector3d color =
          i % 2 == 0 ? colors::talking_orange : colors::enemys_blood;
      auto bz = build_bezier(curves[i]);
      auto tess = bz.get_tesselation2();
      for (int j = 0; j < (int)tess.cols() - 1; ++j)
        draw::line(tess.col(j), tess.col(j + 1),
                   Style::outline(color, 7));
    } else {
      Eigen::Vector3d color = i % 2 == 0 ? colors::forest_green : colors::green;
      draw::line(curves[i].col(0), curves[i].col(1), Style::outline(color, 7));
    }
  }
}

void draw_curves_fill(const std::vector<Eigen::Matrix2Xd>& curves) {
  using namespace polyvec;
  vec<real2> polygon;
  for (int i = 0; i < (int)curves.size(); ++i) {
    if (curves[i].cols() == 4) {
      auto bz = build_bezier(curves[i]);
      auto tess = bz.get_tesselation2();
      for (int j = 0; j < (int)tess.cols() - 1; ++j)
        polygon.push_back(tess.col(j));
    } else {
      for (int j = 0; j < (int)curves[i].cols() - 1; ++j)
           polygon.push_back(curves[i].col(j));
    } 
  }
  draw::polygon(polygon, Style::fill(colors::dark_gray));
}


void _tiny_merges() {
  using Mat24 = Eigen::Matrix<double, 2, 4>;
  using Mat22 = Eigen::Matrix<double, 2, 2>;

  std::vector<Eigen::Matrix2Xd> curves;
  std::vector<Eigen::Matrix2Xd> curves_merged;
  std::vector<std::vector<int>> out2in;
  bool is_circular;
  double per_curve_cost;

  int dump_counter = 0;
  auto write_results = [&]() {
    // Write the unmerged
    polyvec::VtkCurveWriter writer;

    writer.add_point(curves.back().col(curves.back().cols() - 1));

    for (const Eigen::Matrix2Xd& bzc : curves) {
      if (bzc.cols() == 4)
        writer.add_polyline(build_bezier(bzc).get_tesselation2());
      if (bzc.cols() == 2) writer.add_polyline(bzc);
      for (int i = 0; i < (int)bzc.cols(); ++i) {
        writer.add_point(bzc.col(i));
      }
    }

    writer.dump(polyfit::StringUtils::fmt("test_dump/merge_curves_%04d_%d.vtk",
                                          0, dump_counter));
    writer.clear();

    // Write the merged

    writer.add_point(curves.back().col(curves.back().cols() - 1));

    for (const Eigen::Matrix2Xd& bzc : curves_merged) {
      if (bzc.cols() == 4)
        writer.add_polyline(build_bezier(bzc).get_tesselation2());
      if (bzc.cols() == 2) writer.add_polyline(bzc);
      for (int i = 0; i < (int)bzc.cols(); ++i) {
        writer.add_point(bzc.col(i));
      }
    }

    writer.dump(polyfit::StringUtils::fmt("test_dump/merge_curves_%04d_%d.vtk",
                                          1, dump_counter));
    writer.clear();

    ++dump_counter;
  };

  auto attempt_merge = [&]() {
    polyfit::BezierMerging::MergeRecursivelyOptions options;
    options.distance_tolerance = 1.;
    options.alpha_max = 0.99; // allow a bigger alpha here so that the merge succeed
    options.per_bezier_const =per_curve_cost;
    options.allow_merging_lines = true; // let this work for more testing.

    polyfit::BezierMerging::merge_recursively(
        curves, is_circular, options,  curves_merged, out2in);
    printf("ATTEMP %d \n ", dump_counter);
    for (int i = 0; i < (int)curves_merged.size(); ++i) {
      printf("Curve %d | merged ", i);
      for (int j = 0; j < (int)out2in[i].size(); ++j)
        printf("%d, ", out2in[i][j]);
      printf(" \n ");
    }

    write_results();
  };

  //
  // First try a single bezier
  //
  is_circular = false;
  curves.clear();
  curves.push_back(Mat24());
  curves.back().col(0) << -1, 0;
  curves.back().col(1) << -0.25, 0.75;
  curves.back().col(2) << 0.7, 0.3;
  curves.back().col(3) << 1, 0;
  attempt_merge();

  //
  // Now  single bezier with one line before failes
  //
  is_circular = false;
  per_curve_cost = 0.2;
  curves.clear();
  //
  curves.push_back(Mat22());
  curves.front().col(0) << -2, -1;
  curves.front().col(1) << -1, 0;
  //
  curves.push_back(Mat24());
  curves.back().col(0) << -1, 0;
  curves.back().col(1) << -0.25, 0.75;
  curves.back().col(2) << 0.7, 0.3;
  curves.back().col(3) << 1, 0;
  attempt_merge();

  //
  // Reduce the cost and try again.
  //
  is_circular = false;
  per_curve_cost = 0.10;
  attempt_merge();

  //
  //  Now  single bezier with one line before and one after
  //
  is_circular = false;
  per_curve_cost = 0.2;
  curves.clear();
  //
  curves.push_back(Mat22());
  curves.back().col(0) << -1.5, -0.5;
  curves.back().col(1) << -1, 0;
  //
  curves.push_back(Mat24());
  curves.back().col(0) << -1, 0;
  curves.back().col(1) << -0.25, 0.75;
  curves.back().col(2) << 0.7, 0.3;
  curves.back().col(3) << 1, 0;
  //
  curves.push_back(Mat22());
  curves.back().col(0) << 1, 0;
  curves.back().col(1) << 1.5, -0.5;
  //
  attempt_merge();

  //
  // Reduce the cost and try again.
  //
  is_circular = false;
  per_curve_cost = 0.025;
  attempt_merge();

  //
  // Same as before just transform
  //
  is_circular = false;
  per_curve_cost = 0.2;
  curves.clear();
  curves.push_back(Mat22());
  curves.back().col(0) << -1.25, -0.25;
  curves.back().col(1) << -1, 0;
  //
  curves.push_back(Mat24());
  curves.back().col(0) << -1, 0;
  curves.back().col(1) << -0.25, 0.75;
  curves.back().col(2) << 0.7, 0.3;
  curves.back().col(3) << 1, 0;
  //
  curves.push_back(Mat22());
  curves.back().col(0) << 1, 0;
  curves.back().col(1) << 1.1, -0.1;
  {
    Eigen::Matrix2d A;
    A << -3, 1, 2, -4;
    Eigen::Vector2d b;
    b << -1, 4;
    for (Eigen::Matrix2Xd& bz : curves) {
      bz = ((A * bz).colwise() + b).eval();
    }
  }
  attempt_merge();

  //
  // Reduce the cost and try again.
  //
  is_circular = false;
  per_curve_cost = 0.025;
  attempt_merge();

  //
  // Now add two beziers
  //
  is_circular = false;
  per_curve_cost = 0.2;
  curves.clear();
  //
  curves.push_back(Mat22());
  curves.back().col(0) << -1.25, -0.25;
  curves.back().col(1) << -1, 0;
  //
  curves.push_back(Mat24());
  curves.back().col(0) << -1, 0;  //
  curves.back().col(1) << -0.6, 0.4;
  curves.back().col(2) << -0.2, 0.6;  //
  curves.back().col(3) << 0, 0.6;     //
  curves.push_back(Mat24());
  curves.back().col(0) << 0, 0.6;    //
  curves.back().col(1) << 0.3, 0.6;  //
  curves.back().col(2) << 0.65, 0.35;
  curves.back().col(3) << 1, 0;  //
  //
  curves.push_back(Mat22());
  curves.back().col(0) << 1, 0;
  curves.back().col(1) << 1.1, -0.1;

  {
    Eigen::Matrix2d A;
    A << -3, 1, 2, -4;
    Eigen::Vector2d b;
    b << -1, 4;

    for (Eigen::Matrix2Xd& bz : curves) {
      bz = ((A * bz).colwise() + b).eval();
    }
  }
  attempt_merge();

  //
  // Reduce cost and run again
  //
  is_circular = false;
  per_curve_cost = 0.025;
  attempt_merge();
}

void _potrace_output_merges(const std::string input_address,
                            const std::vector<double>& percurve_penalty) {
  using namespace polyfit;
  static int dump_id = 0;

  // Read the image
  IO::Image I;
  if (!IO::read_image(input_address.c_str(), I)) {
    printf("image not found %s", input_address.c_str());
    return;
  }

  // Get the boundary
  std::vector<mat2x> boundaries;
  std::vector<vec4> colors;
  ImageSegment::expand_and_cleanup(I);
  ImageSegment::extract_closed_regions(I, boundaries, colors);
  const mat2x& R =
      boundaries[ImageSegment::find_binary_color_region(boundaries)];

  // Run the potrace pipeline to get some curves
  std::vector<Eigen::Matrix2Xd> curves = PotraceUtil::run_pipeline(
      R, false /* does not really matter */, false /*absolutely don't merge */);

  std::unique_ptr<::polyvec::DevicePDF> pdf;

  //
  // Draw the image by random colors
  //
  pdf.reset(new ::polyvec::DevicePDF(
      ::polyvec::misc::sfmt("test_dump/merge_curves_potrace_%04d.pdf", dump_id).c_str(),
      (int)(2 + percurve_penalty.size()), 2));

  draw_curves_fill({R});
  pdf->draw(0, 0);

  draw_curves(curves);
  pdf->draw(1, 0);

  draw_curves_fill(curves);
  pdf->draw(1, 1);

  for (int j = 0; j < (int)percurve_penalty.size(); ++j) {
    std::vector<Eigen::Matrix2Xd> curves_merged;
    std::vector<std::vector<int>> out2in;

    polyfit::BezierMerging::MergeRecursivelyOptions options;
    options.distance_tolerance = 0.3;
    options.alpha_max = 0.8;
    options.per_bezier_const =percurve_penalty[j];
    options.allow_merging_lines = true; // let this work for more testing.

    polyfit::BezierMerging::merge_recursively(
        curves, true, options, curves_merged, out2in);

    draw_curves(curves_merged);
    pdf->draw(j+2, 0);

    draw_curves_fill(curves_merged);
    pdf->draw(j+2, 1);
  }

  ++dump_id;
}

void _polyfit_output_merges(const std::string input_address,
                            const std::vector<double>& percurve_penalty) {
  using namespace polyfit;
  static int dump_id = 0;

  std::unique_ptr<::polyvec::DevicePDF> pdf;

  //
  // Draw the image by random colors
  //
  pdf.reset(new ::polyvec::DevicePDF(
      ::polyvec::misc::sfmt("test_dump/merge_curves_polyfit_%04d.pdf", dump_id).c_str(),
      (int)(1 + percurve_penalty.size()), 2));

  FILE *fl = fopen(input_address.c_str(),"r");
  std::vector<Eigen::Matrix2Xd> curves = BezierMerging::read_curve_sequence(fl);
  fclose(fl);

  draw_curves(curves);
  pdf->draw(0, 0);

  draw_curves_fill(curves);
  pdf->draw(0, 1);

  for (int j = 0; j < (int)percurve_penalty.size(); ++j) {
    std::vector<Eigen::Matrix2Xd> curves_merged;
    std::vector<std::vector<int>> out2in;

    polyfit::BezierMerging::MergeRecursivelyOptions options;
    options.distance_tolerance = 0.3;
    options.alpha_max = 0.8;
    options.per_bezier_const =percurve_penalty[j];
    options.allow_merging_lines = true; // let this work for more testing.

    polyfit::BezierMerging::merge_recursively(
        curves, true, options, curves_merged, out2in);

    draw_curves(curves_merged);
    pdf->draw(j+1, 0);

    draw_curves_fill(curves_merged);
    pdf->draw(j+1, 1);
  }

  ++dump_id;
}

NAMESPACE_END()

NAMESPACE_BEGIN(polyvectest)
NAMESPACE_BEGIN(BezierMerging)

int merge_curves(int, char**) {
  // ::_polyfit_output_merges(POLYVEC_TEST_PATH "/bezier_merging/merge_failure_01.txt", {3});
  // return EXIT_FAILURE;

  // MINY TESTS
  ::_tiny_merges();

  // POLYFIT FAILURES
  ::_polyfit_output_merges(POLYVEC_TEST_PATH "/bezier_merging/merge_failure_01.txt", {3});

  // POTRACE
  #define PREFIX "/home/hooshi/code/raster2vector-2019/polyvec_data/binary-perfect"
  std::vector<double> params_to_test = {1, 2, 3 , 4};
  // ::_potrace_output_merges(PREFIX "/castle/32.png",params_to_test);
  ::_potrace_output_merges(PREFIX "/castle/64.png", params_to_test);
  // ::_potrace_output_merges(PREFIX "/apple/32.png", params_to_test);
  ::_potrace_output_merges(PREFIX "/apple/64.png", params_to_test);
  // ::_potrace_output_merges(PREFIX "/plane/32.png", params_to_test);
  ::_potrace_output_merges(PREFIX "/plane/64.png", params_to_test);
  ::_potrace_output_merges(PREFIX "/tabletennis/32.png", params_to_test);
  ::_potrace_output_merges(PREFIX "/tabletennis/64.png", params_to_test);
  #undef PREFIX

  return EXIT_FAILURE;
}

NAMESPACE_END(BezierMerging)
NAMESPACE_END(polyvectest)
