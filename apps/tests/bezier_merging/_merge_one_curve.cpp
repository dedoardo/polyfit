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

NAMESPACE_BEGIN()
::polyvec::BezierCurve build_bezier(const Eigen::Matrix<double, 2, 4> ctrl) {
  ::polyvec::BezierCurve bz;
  bz.set_control_points(ctrl);
  return bz;
};
NAMESPACE_END()

NAMESPACE_BEGIN(polyvectest)
NAMESPACE_BEGIN(BezierMerging)

int merge_one_curve(int, char**) {
  using Mat24 = Eigen::Matrix<double, 2, 4>;
  using Mat22 = Eigen::Matrix<double, 2, 2>;

  std::vector<Mat22> line_begin;
  std::vector<Mat22> line_end;
  std::vector<Mat24> beziers;
  Mat24 merged_bezier;
  polyfit::BezierMerging::PotraceParameters potparams;
  bool doable;

  int dump_counter = 0;
  auto write_results = [&](const bool doable) {
    // Write the unmerged
    polyvec::VtkCurveWriter writer;

    if (line_begin.size()) {
      writer.add_polyline(line_begin[0]);
      writer.add_point(line_begin[0].col(0));
    }

    if (line_end.size()) {
      writer.add_polyline(line_end[0]);
      writer.add_point(line_end[0].col(0));
      writer.add_point(line_end[0].col(1));
    } else {
      writer.add_point(beziers.back().col(3));
    }

    for (const Mat24& bzc : beziers) {
      writer.add_polyline(build_bezier(bzc).get_tesselation2());
      writer.add_point(bzc.col(0));
      writer.add_point(bzc.col(1));
      writer.add_point(bzc.col(2));
    }

    writer.dump(polyfit::StringUtils::fmt("test_dump/merge_one_curve_%04d_%d.vtk",
                                          0, dump_counter));
    writer.clear();

    // Write the merged

    Eigen::MatrixXd todraw(2, 3);

    if (doable) {
      todraw.col(0) = potparams.get_yy().control_points.col(0);
      todraw.col(1) = potparams.get_oy();
      todraw.col(2) = potparams.get_yy().control_points.col(3);
      writer.add_polyline(todraw);
      //
      writer.add_point(potparams.get_yy().control_points.col(0));
      writer.add_point(potparams.get_yy().control_points.col(1));
      writer.add_point(potparams.get_yy().control_points.col(2));
      writer.add_point(potparams.get_yy().control_points.col(3));
      //
      writer.add_polyline(
          build_bezier(potparams.get_yy().control_points).get_tesselation2());
    }

    //
    writer.dump(polyfit::StringUtils::fmt("test_dump/merge_one_curve_%04d_%d.vtk",
                                          1, dump_counter));
    writer.clear();

    ++dump_counter;
  };

  auto attempt_merge = [&]() {
    polyfit::BezierMerging::merged_curve(line_begin, beziers, line_end, stdout,
                                         doable, merged_bezier);

    if (doable) {
      potparams = ::polyfit::BezierMerging::PointParameters({merged_bezier})
                      .as_potrace_params();
      const double tolerance = 1.0;
      const double alphamax = 0.9;
      double penalty = 0, potrace_penalty=0;
      bool is_acceptable = false;
      ::polyfit::BezierMerging::is_merge_acceptable(beziers, merged_bezier, tolerance,alphamax,
                                              potrace_penalty, is_acceptable);
      ::polyfit::BezierMerging::merge_penalty(beziers, merged_bezier, penalty);
      printf("Case %d | acceptable %d penalty_potrace %.8g penalty %.8g \n", dump_counter,
             is_acceptable, potrace_penalty, penalty);
    } else {
      printf("Case %d | cannot fit merged bezier \n", dump_counter);
    }

    write_results(doable);
  };

  //
  // First try a single bezier
  //
  line_begin.clear();
  line_end.clear();
  beziers.push_back(Mat24());
  beziers.back().col(0) << -1, 0;
  beziers.back().col(1) << -0.25, 0.75;
  beziers.back().col(2) << 0.7, 0.3;
  beziers.back().col(3) << 1, 0;
  attempt_merge();

  //
  // Now  single bezier with one line before failes
  //
  line_end.clear();
  //
  line_begin.clear();
  line_begin.push_back(Mat22());
  line_begin.front().col(0) << -2, -1;
  line_begin.front().col(1) << -1, 0;
  //
  beziers.clear();
  beziers.push_back(Mat24());
  beziers.back().col(0) << -1, 0;
  beziers.back().col(1) << -0.25, 0.75;
  beziers.back().col(2) << 0.7, 0.3;
  beziers.back().col(3) << 1, 0;
  attempt_merge();

  //
  //  Now  single bezier with one line before and one after
  //
  line_end.clear();
  line_end.push_back(Mat22());
  line_end.front().col(0) << 1, 0;
  line_end.front().col(1) << 1.5, -0.5;
  //
  line_begin.clear();
  line_begin.push_back(Mat22());
  line_begin.front().col(0) << -1.5, -0.5;
  line_begin.front().col(1) << -1, 0;
  //
  beziers.clear();
  beziers.push_back(Mat24());
  beziers.back().col(0) << -1, 0;
  beziers.back().col(1) << -0.25, 0.75;
  beziers.back().col(2) << 0.7, 0.3;
  beziers.back().col(3) << 1, 0;
  attempt_merge();

  //
  // Same as before just transform
  //
  line_end.clear();
  line_end.push_back(Mat22());
  line_end.front().col(0) << 1, 0;
  line_end.front().col(1) << 1.1, -0.1;
  //
  line_begin.clear();
  line_begin.push_back(Mat22());
  line_begin.front().col(0) << -1.25, -0.25;
  line_begin.front().col(1) << -1, 0;
  //
  beziers.clear();
  beziers.push_back(Mat24());
  beziers.back().col(0) << -1, 0;
  beziers.back().col(1) << -0.25, 0.75;
  beziers.back().col(2) << 0.7, 0.3;
  beziers.back().col(3) << 1, 0;
  {
    Eigen::Matrix2d A;
    A << -3, 1, 2, -4;
    Eigen::Vector2d b;
    b << -1, 4;
    line_end.front() = ((A * line_end.front()).colwise() + b).eval();
    line_begin.front() = ((A * line_begin.front()).colwise() + b).eval();

    for (Mat24& bz : beziers) {
      bz = ((A * bz).colwise() + b).eval();
    }
  }
  attempt_merge();

  //
  // Now add two beziers
  //
  line_end.clear();
  line_end.push_back(Mat22());
  line_end.front().col(0) << 1, 0;
  line_end.front().col(1) << 1.1, -0.1;
  //
  line_begin.clear();
  line_begin.push_back(Mat22());
  line_begin.front().col(0) << -1.25, -0.25;
  line_begin.front().col(1) << -1, 0;
  //
  beziers.clear();
  beziers.push_back(Mat24());
  beziers.back().col(0) << -1, 0;  //
  beziers.back().col(1) << -0.6, 0.4;
  beziers.back().col(2) << -0.2, 0.6;  //
  beziers.back().col(3) << 0, 0.6;     //
  beziers.push_back(Mat24());
  beziers.back().col(0) << 0, 0.6;    //
  beziers.back().col(1) << 0.3, 0.6;  //
  beziers.back().col(2) << 0.65, 0.35;
  beziers.back().col(3) << 1, 0;  //
  {
    Eigen::Matrix2d A;
    A << -3, 1, 2, -4;
    Eigen::Vector2d b;
    b << -1, 4;
    line_end.front() = ((A * line_end.front()).colwise() + b).eval();
    line_begin.front() = ((A * line_begin.front()).colwise() + b).eval();

    for (Mat24& bz : beziers) {
      bz = ((A * bz).colwise() + b).eval();
    }
  }
  attempt_merge();

  //
  // This was a problem with potrace incoming curves
  //
  line_end.clear();
  //
  line_begin.clear();
  //
  beziers.clear();
  beziers.push_back(Mat24());
  beziers.push_back(Mat24());
  beziers.push_back(Mat24());
  beziers.push_back(Mat24());
  beziers.push_back(Mat24());
 beziers[0] <<   21, 22.5238,  23.119,    23.5,      30,      30 , 30.8333 ,   33.5;
  beziers[1] << 23.5, 23.9952, 24.0476,      29,    33.5 ,36.9667  ,    37,      37;
 beziers[2] <<   29,  31.75, 34.225,   34.5,     37,     37, 36.775  , 36.5;
 beziers[3] << 34.5, 34.775, 35.225,   35.5 ,  36.5, 36.225 ,30.825 ,  24.5;
 beziers[4] <<  35.5,  35.881, 35.6429,    34.5,    24.5, 15.7381, 12.0476 ,      9;
  attempt_merge();

  //
  // This was a problem with potrace incoming curves
  //
  line_end.clear();
  //
  line_begin.clear();
  //
  beziers.clear();
  beziers.push_back(Mat24());
  beziers.push_back(Mat24());
  beziers[0] <<   29,  31.75, 34.225,   34.5,     37,     37, 36.775,   36.5;
  beziers[1] << 34.5, 34.775, 35.225,   35.5,   36.5, 36.225, 30.825,   24.5;
  attempt_merge();

   

  return EXIT_FAILURE;
}

NAMESPACE_END(BezierMerging)
NAMESPACE_END(polyvectest)
