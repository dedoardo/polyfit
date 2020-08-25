/*
	This is meant to replace core/options.hpp 
*/
#pragma once

// polyvec0
#include <polyvec/core/macros.hpp>
#include <polyvec/geometry/angle.hpp>

NAMESPACE_BEGIN(polyfit)

struct Options {
	double short_edge_penalty = .5;
    double edge_penalty = .0;

    // Raster::distance_bounds_from_points_with_slack
    double midpoint_accuracy_trigger_projection = 0.5;

    // PolygonTracer::regularized
    double regularity_raster_parallel_min_feature_ratio = 0.;
    double regularity_raster_parallel_proximity_scale = 1.;
    double regularity_raster_parallel_min_length_ratio = 0.;

    // Regularity::construct_regularity_graph
    double regularity_continuation_distance_max_32x32 = 1.;
    
    double regularity_continuation_angle_max = PF_RAD(135);
    double regularity_continuation_angle_separation_max = PF_RAD(45);

    double regularity_continuation_angle_polygon_max = PF_RAD(30);
    double regularity_continuation_angle_intersection_max = PF_RAD(110);

    // PostProcessor
    double post_processor_max_short_edge_length = 2.2360679775; // sqrt(5). diagonal of a 2by1 rectangle

	static Options* get();
};
#define SQUARED_WITH_EPS(v) ((v * v) + PF_EPS)

NAMESPACE_END(polyfit)