#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

// All functions assume that the path points are oriented clockwise
NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(LineUtils)

// does not check bounds for t on either line
double project_t(
    const polyfit::vec2& point, 
    const polyfit::vec2& line_src,
    const polyfit::vec2& line_dest
);

// Identical behavior to project_t, but returns the location of the point rather than the distance along the line
Eigen::Vector2d project_point(
    const Eigen::Vector2d& point,
    const Eigen::Vector2d& line_src,
    const Eigen::Vector2d& line_dst
);

polyfit::vec2 line_at(const polyfit::vec2& src, const polyfit::vec2& dst, const double t);

// does not check bounds on either line
bool intersect ( const polyfit::vec2& lhs_src, const polyfit::vec2& lhs_dst, const polyfit::vec2& rhs_src,
                 const polyfit::vec2& rhs_dst, double& lhs_at, double& rhs_at );

// Does not check bounds for t on the line
void intersect_polyline(const polyfit::vec2& lhs_src, const polyfit::vec2& lhs_dst,
                        const Eigen::Ref<const Eigen::Matrix2Xd> rhs,
                        std::vector<double>& line_at, 
                        std::vector<int>& rhs_seg_no,
                        std::vector<double>& rhs_at);

double distance_from_point(const polyfit::vec2& line_src, const polyfit::vec2& line_dst, const polyfit::vec2& pt);

double signed_distance_from_point(const polyfit::vec2& line_src, const polyfit::vec2& line_dst, const polyfit::vec2& pt);

int which_side_is_the_point_on (
	const Eigen::Vector2d& src,
	const Eigen::Vector2d& dst,
	const Eigen::Vector2d& p
);

NAMESPACE_END(LineUtils)
NAMESPACE_END(polyvec)

