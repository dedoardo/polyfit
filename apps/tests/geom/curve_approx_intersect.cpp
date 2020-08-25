#include <polyvec/geometry/line.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>

using namespace polyvec;

namespace {
void run_test(GlobFitCurve& curve, const Eigen::Vector2d& shooting_plane_src,
              const Eigen::Vector2d& shooting_dir) {
  static int n_called = 0;

  // Find the begin and end of shooting
  double shooting_t_min = +std::numeric_limits<double>::max();
  double shooting_t_max = -std::numeric_limits<double>::max();
  Eigen::Matrix2Xd tess2 = curve.get_tesselation2();

  for (int i = 0; i < (int)tess2.cols(); ++i) {
    const double t_proj = LineUtils::project_t(
        tess2.col(i), shooting_plane_src, shooting_plane_src + shooting_dir);
    shooting_t_min = std::min(shooting_t_min, t_proj);
    shooting_t_max = std::max(shooting_t_max, t_proj);
  }

  // Extend slightly
  shooting_t_min -= 0.2;
  shooting_t_max += 0.2;

  // Now start shooting
  const int n_rays_shot = 200;
  VtkCurveWriter writer;
  const Eigen::Vector2d ray_dir(shooting_dir[1], -shooting_dir[0]);
  writer.add_polyline(tess2);
  std::vector<double> line_at;
  std::vector<double> curve_at;

  for (int i = 0; i < n_rays_shot; ++i) {
    const double t_shooting_line =
        shooting_t_min +
        i * (shooting_t_max - shooting_t_min) / (n_rays_shot - 1);
    const Eigen::VectorXd ray_base =
        t_shooting_line * shooting_dir + shooting_plane_src;

    curve.approximate_intersect(ray_base, ray_base + ray_dir, line_at,
                                 curve_at);

    if (!line_at.size()) {
      writer.add_point(ray_base);
    } else {
      for (int j = 0; j < (int)line_at.size(); ++j) {
        writer.add_line(ray_base, line_at[j] * ray_dir + ray_base);
        writer.add_point(curve.pos(curve_at[j]));
      }
    }
  }

  writer.dump(polyvec_str("test_dump/approx_intersect" << n_called << ".vtk"));

  ++n_called;
}

}  // namespace

namespace polyvectest {
namespace geom {
int curve_approx_intersect(int /*argc*/, char** /*argv*/) {
  Eigen::Matrix<double, 2, 4> pts;
  pts.col(0) << 0, 0;
  pts.col(1) << 2, 2;
  pts.col(2) << -1, 3;
  pts.col(3) << 1.5, -0.5;

  BezierCurve bz;
  bz.set_control_points(pts);

  run_test( bz, Eigen::Vector2d(5, 1),  Eigen::Vector2d(0, 1) );
  run_test( bz, Eigen::Vector2d(1, -5),  Eigen::Vector2d(1, 0) );
  run_test( bz, Eigen::Vector2d(5, 1),  Eigen::Vector2d(1, 1).normalized() );
  run_test( bz, Eigen::Vector2d(5, 1),  Eigen::Vector2d(0.25, 0.75).normalized() );
  run_test( bz, Eigen::Vector2d(5, 1),  Eigen::Vector2d(-0.25, 0.75).normalized() );

  return 0;
}  // end of function
}  // namespace geom
}  // namespace polyvectest
