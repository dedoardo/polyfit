#pragma once

#include <vector>

#include <polyvec/api.hpp>

#include <polyvec/mc/mesh_connectivity.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
#include <polyvec/mc/junction.hpp>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )

// Stores the connectivity of the segments in a multicolored images
//
// USED IN: SH_export_traces.cpp, SH_mc_fit_polygons.cpp
//
// Add extra info if needed
class SegmentConnectivity {
    // =====
    // Functions
    // =====
public:
    // A dummy index for the outer polygon.
    constexpr static int infinite_region = RasterImageConnectivity::infinite_region;

    // Number of segments in the image
    int n_segments();

    // The underlying raster image connectivity
    RasterImageConnectivity & raster_image_connectivity();

    // ==== Map segment and pixel half-edges
    // The half-edge of raster_image_connectivity() that lies at the beginning of
    // the segment
    int segment_to_pixelhe_begin(const int segment_id);
    // The half-edge of raster_image_connectivity() that lies at the end of
    // the segment (non-inclusive)
    int segment_to_pixelhe_end(const int segment_id);
    // For a half-edge in raster_image_connectivity(), returns the segment to which 
    // the half-edge belongs.
    // TODO: also return if the half-edge is on the front or back side of the segment.
    int pixelhe_to_segment(const int pixelhe_id);

    // Returns all the points on a polygon
    // Treats kopf junctions as two points, so no need to worry about them
    void get_segment_points(const int segment_id, Eigen::Matrix2Xi &points, bool &is_closed);

    // Returns the regions which are sharing the specific segment
    void get_incident_polygons_for_segment(const int segment_id, int& region0, int& region1);

private:
    SegmentConnectivity() = default;

    //
    // Static functions
    //
public:
    // constructor
    // Input:
    //     RasterImageConnectivity: the underlying raster image connectivity   
    //     junction_types: junction type for each points (size n+1*m+1for an nxm image)
    static SegmentConnectivity build ( RasterImageConnectivity &, const std::vector<JunctionType> & junction_types);

private:
    // The segment connectivities
    RasterImageConnectivity *_raster_image_connectivity;

    // Segment he to pixel he mapping
    int _n_segments;
    std::vector<int> _pixelhe_to_segment;
    std::vector<int> _segment_to_pixelhe_begin;
    std::vector<int> _segment_to_pixelhe_end;
};


NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
