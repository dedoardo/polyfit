#pragma once

#include <polyvec/core/types.hpp>
#include <polyvec/core/macros.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(ErrorMetrics)

enum EdgeProperties {
	EDGE_FIRST_IN_PATH = 1 << 0,
	EDGE_LAST_IN_PATH = 1 << 1,
	EDGE_HAS_INFLECTION = 1 << 2
};
		
enum Type {
	SMOOTHNESS = 0,
	CONTINUITY,
	ACCURACY
};

struct Options {
	double smoothness_weight;
	double smoothness_limit;
	double inflection_weight;
	double inflection_limit;
	double continuity_weight;
	double continuity_limit;
	double accuracy_weight;
};

double smoothness(const double rad0, const double rad1, const int flags);
double smoothness(const double rad0);
double accuracy(const double d_min, const double d_max, const int flags);
double accuracy_implicit(const double d_min, const double d_max, const int flags);
double continuity(const double rad0, const double rad1, const int flags);
bool   accuracy_within_bounds(const double d_min, const double d_max, const vec2& edge_dir);
bool   accuracy_within_bounds_relaxed(const double d_min, const double d_max);

double calculate_inner(double d_min, double d_max, const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3);
double calculate_first(double d_min, double d_max, const vec2 p0, const vec2 p1, const vec2 p2);
double calculate_last(double d_min, double d_max, const vec2 p1, const vec2 p2, const vec2 p3);

vec3 calculate_inner_sep(double d_min, double d_max, const vec2 p0, const vec2 p1, const vec2 p2, const vec2 p3);
vec3 calculate_first_sep(double d_min, double d_max, const vec2 p0, const vec2 p1, const vec2 p2);
vec3 calculate_last_sep(double d_min, double d_max, const vec2 p1, const vec2 p2, const vec2 p3);

NAMESPACE_END(ErrorMetrics)
NAMESPACE_END(polyfit)