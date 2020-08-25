#include <iostream>
//
#include <polyvec/api.hpp>
#include <polyvec/curve-tracer/bezier_merging.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/utils/num.hpp>
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

int potrace_params(int, char**) {
  // Creat the points and then test the mapping
  Eigen::Matrix<double, 2, 3> points0;
  Eigen::Matrix<double, 2, 4> yy;
  Eigen::Matrix2d A;
  Eigen::Vector2d b;
  Eigen::Matrix2d Amanu;
  Eigen::Vector2d bmanu;
  polyvec::VtkCurveWriter writer;

  points0.col(0) << 2, 10;
  points0.col(1) << 5, 7;
  points0.col(2) << 4, 12;

  //
  // First Just Create a bezier and draw it
  //
  const double alpha_manu = 0.8;
  const double beta_manu = 0.3;

  yy.col(0) = points0.col(0);
  yy.col(1) = ::polyvec::Num::lerp<Eigen::Vector2d>(points0.col(0),
                                                    points0.col(1), alpha_manu);
  yy.col(2) = ::polyvec::Num::lerp<Eigen::Vector2d>(points0.col(2),
                                                    points0.col(1), beta_manu);
  yy.col(3) = points0.col(2);

  writer.add_polyline(points0);
  writer.add_point(yy.col(0));
  writer.add_point(yy.col(1));
  writer.add_point(yy.col(2));
  writer.add_point(yy.col(3));
  writer.add_polyline(build_bezier(yy).get_tesselation2());
  writer.dump("test_dump/potrace_params_00.vtk");
  writer.clear();

  //
  // Now create the potrace parameterization
  //
  ::polyfit::BezierMerging::PointParameters point_params{yy};
  assert_break(point_params.is_potracable());
  {
    double alpha, beta;
    Eigen::Vector2d oo;
    bool is;
    point_params.is_potracable(is, oo, alpha, beta);
    assert_break(is == true);
    assert_break(std::abs(alpha - alpha_manu) < 1e-8);
    assert_break(std::abs(beta - beta_manu) < 1e-8);
    assert_break((oo - points0.col(1)).norm() < 1e-8);
  }
  ::polyfit::BezierMerging::PotraceParameters potrace_params =
      point_params.as_potrace_params();

  writer.add_polyline(
      build_bezier(potrace_params.get_yy().control_points).get_tesselation2());
  writer.dump("test_dump/potrace_params_01.vtk");
  writer.clear();

  //
  // Now draw the closest equiparam bezier
  //

  writer.add_polyline(
      build_bezier(potrace_params.equiparamed().get_yy().control_points)
          .get_tesselation2());
  writer.dump("test_dump/potrace_params_02.vtk");
  writer.clear();

  //
  // Now check that the area works
  //
  Eigen::Matrix2Xd tess;
  tess = build_bezier(yy).get_tesselation2();
  double area_tess;
  {
    bool is_ccw;
    polyvec::WindingNumber::compute_orientation(tess, is_ccw, area_tess);
  }
  double area_analytic = potrace_params.get_areay();
  double area_eqp_analytic = potrace_params.equiparamed().get_areay();

  printf("Testing area \n");
  printf("areas: tess: %.6f, ana: %.6f, anaeq: %.6f \n", area_tess,
         area_analytic, area_eqp_analytic);

  //
  // Now check that the tangent things works
  //
  printf("Testing inverse tangent solver \n");
  {
    ::polyvec::BezierCurve bzcurve = build_bezier(yy);
    for (int i = 1; i < 99 ; ++i) {
      double t_solved;
      bool success;

      const double t = 1. / 99 * i;
      const Eigen::Vector2d tan = bzcurve.dposdt(t).normalized();
      potrace_params.t_for_tangent(tan, t_solved, success);

      assert_break( success );
      assert_break( std::abs(t-t_solved) < 1e-10 );
    }
  }

  return EXIT_FAILURE;
}

NAMESPACE_END(BezierMerging)
NAMESPACE_END(polyvectest)
