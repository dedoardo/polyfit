#include <polyvec/geometry/line.hpp>

// Header
#include <polyvec/geom.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/num.hpp>
#include <polyvec/utils/directions.hpp>

using namespace polyfit;

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(LineUtils)

double project_t(const vec2& point, const vec2& line_src,
                 const vec2& line_dest) {
  const vec2 point_v = point - line_src;
  const vec2 line_dir = (line_dest - line_src);
  const double vd_dd = point_v.dot(line_dir) / (line_dir.squaredNorm());
  return vd_dd;
}

Eigen::Vector2d project_point(
    const Eigen::Vector2d& point,
    const Eigen::Vector2d& line_src,
    const Eigen::Vector2d& line_dst
) {
    return Num::lerp(line_src, line_dst, project_t(point, line_src, line_dst));
}

vec2 line_at(const vec2& src, const vec2& dst, const double t) {
  return src + (dst - src) * t;
}

// 
bool intersect(const vec2& lhs_src, const vec2& lhs_dst, const vec2& rhs_src,
               const vec2& rhs_dst, double& lhs_at, double& rhs_at) {
  const double tol = 1e-15;

  const Eigen::Vector2d tang_line0 = (lhs_dst - lhs_src);
  const Eigen::Vector2d tang_norm_line0 = (tang_line0).normalized();
  const Eigen::Vector2d n_line0(-tang_norm_line0.y(), tang_norm_line0.x());

  const Eigen::Vector2d tang_line1 = (rhs_dst - rhs_src);
  const Eigen::Vector2d tang_norm_line1 = (tang_line1).normalized();
  const Eigen::Vector2d n_line1(-tang_norm_line1.y(), tang_norm_line1.x());

  const double denom_line0 = tang_line0.dot(n_line1);
  const double denom_line1 = tang_line1.dot(n_line0);

  bool ret;

  if (std::abs(denom_line0) < tol) {
    ret = false;
  } else {
    ret = true;
    const double nom_line0 = n_line1.dot(rhs_src - lhs_src);
    const double nom_line1 = n_line0.dot(lhs_src - rhs_src);
    lhs_at = nom_line0 / denom_line0;
    rhs_at = nom_line1 / denom_line1;
  }

  return ret;
}

void intersect_polyline(const vec2& lhs_src, const vec2& lhs_dst,
                        const Eigen::Ref<const Eigen::Matrix2Xd> rhs,
                        std::vector<double>& line_at, 
                        std::vector<int>& rhs_seg_no,
                        std::vector<double>& rhs_at) {
                          
  const double tol = 1e-15;
  double rhs_at_now;
  double line_at_now;

  line_at.clear();
  rhs_seg_no.clear();
  rhs_at.clear();
  

  for (int i = 0; i < (int)(rhs.cols() - 1); ++i) {
    const bool does_inersect =
        intersect(lhs_src, lhs_dst, rhs.col(i), rhs.col(i + 1), line_at_now, rhs_at_now);

    if (does_inersect && (rhs_at_now > -tol) && (rhs_at_now < 1. + tol)) {
      rhs_seg_no.push_back(i);
	  if (rhs_at_now < 0)
		  rhs_at_now = 0;
	  if (rhs_at_now > 1)
		  rhs_at_now = 1;
      rhs_at.push_back(rhs_at_now);
      line_at.push_back( line_at_now ) ;
    }
  }

}

// Returns the normal from the tangent averaging the neighboring edges 
Eigen::Vector2d path_normal_avg_corner(const Eigen::Matrix2Xd& path, const Eigen::Index corner) {
	assert_break(corner < path.cols());
	const Eigen::Vector2d t0 = (path.col(corner) - CircularAt(path, corner - 1)).normalized();
	const Eigen::Vector2d t1 = (CircularAt(path, corner + 1) - path.col(corner));
	return polyvec::util::normal_dir(misc::lerp(t0, t1, .5)).normalized();
}

// Returns the outward normal orthogonal to the path segment starting at the input corner
Eigen::Vector2d path_normal_avg_edge_from(const Eigen::Matrix2Xd& path, const Eigen::Index segment_src) {
	assert_break(segment_src < path.cols());
	const Eigen::Vector2d np = (CircularAt(path, segment_src) - CircularAt(path, segment_src - 1)).normalized();
	const Eigen::Vector2d nn = (CircularAt(path, segment_src) - CircularAt(path, segment_src + 1)).normalized();
	return misc::lerp(np, nn, .5).normalized();
}

// there is an extra norm here..
double distance_from_point(const vec2& p0, const vec2& p1, const vec2& p) {
	return abs(signed_distance_from_point(p0,p1,p));
}

double signed_distance_from_point(const vec2& p0, const vec2& p1, const vec2& p) {
	const vec2 p01 = (p1 - p0).normalized(); // 
	const vec2 n(-p01.y(), p01.x());
	return (p - p0).dot(n);
}

// This predicate is subject to precision errors
int which_side_is_the_point_on(
	const Eigen::Vector2d& src,
	const Eigen::Vector2d& dst,
	const Eigen::Vector2d& p
) {
	return Num::sign(
		(dst(1) - src(1)) * (p(0) - src(0)) - (dst(0) - src(0)) * (p(1) - src(1))
	);
}

NAMESPACE_END(LineUtils)
NAMESPACE_END(polyvec)