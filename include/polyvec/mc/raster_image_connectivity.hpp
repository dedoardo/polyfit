#pragma once

#include <vector>

#include <polyvec/api.hpp>

#include <polyvec/mc/mesh_connectivity.hpp>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )

// Stores the connectivity of the regions in a multicolored images
//
// TEST: apps_test/mc/_raster_image_connectivity.cpp
// USED IN: SH_export_traces.cpp, SH_mc_fit_polygons.cpp
//
// Add extra info if needed
class RasterImageConnectivity {
    // =====
    // Functions
    // =====
public:
    enum PolygonType {
      // A polygon representing the outer boundary of a region
      POLYGON_OUTER_BOUNDARY = 0,
      // A polygon representing a hole in a region
      POLGYON_HOLE,
      // nothing
      POLYGON_INVALID
    };
 
    // The index for the polygon on the outerside of the bounding box
    constexpr static int infinite_region = Mesh_connectivity::invalid_index;

    // Return the position of each pixel point from its index
    // There are a total of ( (n+1)*(m+1) vertices for an nxm image)
    Eigen::Vector2i pixel_vertex_index_to_pos ( const int index ) ;
    int pixel_vertex_pos_to_index ( const Eigen::Vector2i& pos ) ;

    // Number of pixel vertices including the deactivated ones
    // There are a total of ( (n+1)*(m+1) vertices for an nxm image)
    int n_pixel_vertices() ;

    // Dimensions of the vertices in the images 
    // i.e., n+1 and m+1 for an nxm image.
    const Eigen::Vector2i& pixel_vertex_dims() const;

    // connectivity of pixels as a half-edge data structure
    Mesh_connectivity& connectivity();

    // ====== Polygons and regions and mapping between them
    // number of polygons
    int n_polygons();
    // number of regions
    int n_regions();
    // how many polygons are associated with a region 
    // (one for the outer boundary, and possibly many for holes)
    int n_region_polygons(const int region_id);
    // Get the region index of a polygon
    // Note: the outer boundary has no region and will just return
    // the infinite region index.
    int polygon_to_region(const int polygon_id);
    // Get the polygons of a region
    // outerboundary is the outer_boundary polygon
    // holes are the hole polygons
    void region_to_polygon(const int region_id, int &outer_boundary, std::vector<int>& holes);
    // Get the type of a polygon (hole or outer boundary)
    PolygonType get_polygon_type(const int polygon_id);

    // Returning all the points on a polygon
    Eigen::Matrix2Xi get_polygon_points(const int pid);

    // Return the region indices sorted by area (descending order)
    // The area of the holes are not subtracted 
    // useful for drawing the image
    std::vector<int> regions_sorted_by_area();

    // for debugging
    void dump_as_vtk ( const std::string filename ) ;

    // Kill constructor (I mean make it private not delete)
private:
    RasterImageConnectivity() = default;

    //
    // Static functions
    //
public:

    // constructor
    // input:
    // region_boundaries: the regions returned from get_region_boudaries_from_bmp() 
    static RasterImageConnectivity build ( const std::vector<Eigen::Matrix2Xi>& region_boundaries );

private: 
    // Give numbers to the holes
    static void enumerate_holes ( Mesh_connectivity&, std::vector<int> &half_edge_on_hole );

    // =====
    // Fields
    // =====
private:

    // The boundig box of pixel_vertices -- min must be 0,0
    Eigen::Vector2i _pixel_vertex_dims;

    // The connectivity of the edges of pixels
    Mesh_connectivity _connectivity;

    // Region to polygon mapping 
    // a polygon is either the outer boundary or a hole in a region
    int _n_regions;
    std::vector<int> _polygon_to_region;
    std::vector<int> _region_to_polygon_val;
    std::vector<int> _region_to_polygon_xval;
};


NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )

