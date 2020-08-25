// polyvec
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
#include <polyvec/mc/get_bounding_box.hpp>
#include <polyvec/mc/bit_vector.hpp>
#include <polyvec/mc/associate_holes.hpp>

// libc++
#include <cstdio>
#include <memory>

// eigen
#include <Eigen/Geometry>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )


constexpr int RasterImageConnectivity::infinite_region;

Eigen::Vector2i
RasterImageConnectivity::pixel_vertex_index_to_pos ( const int index ) {
    return Eigen::Vector2i ( index % pixel_vertex_dims().x(), index / pixel_vertex_dims().x() );
}

int
RasterImageConnectivity::pixel_vertex_pos_to_index ( const Eigen::Vector2i& pos ) {
    return pos.x() + pos.y() * pixel_vertex_dims().x();
}


int
RasterImageConnectivity::n_pixel_vertices() {
    return pixel_vertex_dims().x() * pixel_vertex_dims().y();
}

const Eigen::Vector2i&
RasterImageConnectivity::pixel_vertex_dims() const {
    return _pixel_vertex_dims;
}

Mesh_connectivity&
RasterImageConnectivity::connectivity() {
    return _connectivity;
}

int
RasterImageConnectivity::n_polygons() {
    return connectivity().n_active_faces();
}

int
RasterImageConnectivity::n_regions() {
    return _n_regions;
}

int
RasterImageConnectivity::n_region_polygons ( const int region_id ) {
    assert_break ( region_id < n_regions() );
    const int ans= _region_to_polygon_xval[region_id+1]-_region_to_polygon_xval[region_id];
    assert_break ( ans >= 1 );
    return ans;
}

int
RasterImageConnectivity::polygon_to_region ( const int polygon_id ) {
    assert_break ( polygon_id < n_polygons() );
    return _polygon_to_region[ polygon_id ];
}

void
RasterImageConnectivity::region_to_polygon ( const int region_id, int& outer_boundary, std::vector<int>& holes ) {
    outer_boundary = _region_to_polygon_val[ _region_to_polygon_xval[ region_id ] ];
    assert_break ( outer_boundary == region_id );

    holes.resize ( 0 );
    holes.reserve ( std::max ( 0, n_region_polygons ( region_id )-1 ) );

    for ( int j= _region_to_polygon_xval[ region_id ]+1 ; j < _region_to_polygon_xval[ region_id+1 ] ; ++j ) {
        holes.push_back (  _region_to_polygon_val [j] );
    }
}

RasterImageConnectivity::PolygonType
RasterImageConnectivity::get_polygon_type ( const int polygon_id ) {
    if ( polygon_id < n_regions() ) {
        return POLYGON_OUTER_BOUNDARY;
    } else if ( polygon_id < n_polygons() ) {
        return POLGYON_HOLE ;
    } else {
        assert_break ( 0 );
        return POLYGON_INVALID;
    }
}

Eigen::Matrix2Xi
RasterImageConnectivity::get_polygon_points ( const int pid ) {
    // Lambda to extract the points on polygon
    std::vector<int> points;

    Mesh_connectivity::Face_iterator face = connectivity().face_at ( pid );
    Mesh_connectivity::Half_edge_iterator he_end = face.half_edge();
    Mesh_connectivity::Half_edge_iterator he = he_end;

    int safe_guard = 0;
    const int safe_gaurd_cap = 10000;

    do {
        Eigen::Vector2i pos =  pixel_vertex_index_to_pos ( he.origin().index() );
        points.push_back ( pos.x() );
        points.push_back ( pos.y() );
        assert_break ( safe_guard < safe_gaurd_cap );
        ++safe_guard;
        he = he.next();
    } while ( !he.is_equal ( he_end ) );

    return Eigen::Matrix2Xi::ConstMapType ( points.data(), 2, points.size() /2 );
}

std::vector<int>
RasterImageConnectivity::regions_sorted_by_area() {
    std::vector<int> by_area(n_regions());
    std::vector<double> areas(n_regions());
    for (int i = 0; i < (int)by_area.size(); ++i) {
      by_area[i] = i;
      bool _1;
      Eigen::Matrix2Xi bdry = get_polygon_points(i);
      polyvec::WindingNumber::compute_orientation(bdry.cast<double>(), _1,
                                                  areas[i]);
    }
    std::sort(by_area.begin(), by_area.end(),
              [areas](int i, int j) { return areas[i] > areas[j]; });

    return by_area;
}



// dump as vtk for debugging
void
RasterImageConnectivity::dump_as_vtk ( const std::string filename ) {

    // Defrag the mesh
    Mesh_connectivity::Defragmentation_maps defrag;
    connectivity().compute_defragmention_maps ( defrag );


    FILE* fl = fopen ( filename.c_str(), "w" );
    assert_break ( fl && "FILE should be open" );

    /*
     * Write the vtk file.
     */

    // write the header
    fprintf ( fl, "# vtk DataFile Version 2.0\n" );
    fprintf ( fl, "Shayan's output mesh\n" );
    fprintf ( fl, "ASCII\n" );
    fprintf ( fl, "DATASET UNSTRUCTURED_GRID\n" );
    fprintf ( fl, "\n" );

    // write the vertices
    fprintf ( fl, "POINTS %d float\n", connectivity().n_active_vertices() );

    for ( int vnidx = 0; vnidx < connectivity().n_active_vertices(); vnidx++ ) {
        int voidx = defrag.new2old_vertices[vnidx];
        Eigen::Vector2i pos = pixel_vertex_index_to_pos ( voidx );
        fprintf ( fl, "%d %d 0 \n",  pos.x(), pos.y() );
    }

    fprintf ( fl, "\n" );

    //
    // write the faces
    //

    // count their total number of vertices.
    int total_vert_duplicated_per_face = 0;

    for ( int fn = 0; fn < connectivity().n_active_faces(); ++fn ) {
        int fo = defrag.new2old_faces[fn];
        Mesh_connectivity::Face_iterator face = connectivity().face_at ( fo );
        Mesh_connectivity::Half_edge_iterator he_end = face.half_edge();
        Mesh_connectivity::Half_edge_iterator he = face.half_edge();

        do {
            ++total_vert_duplicated_per_face;
            he = he.next();
        } while ( !he.is_equal ( he_end ) );
    }

    int face_verts_cache[4096]; // this must be more than maximum number of edges per face.
    fprintf ( fl, "CELLS %d %d \n", connectivity().n_active_faces(), connectivity().n_active_faces() + total_vert_duplicated_per_face );

    for ( int fn = 0; fn < connectivity().n_active_faces(); ++fn ) {
        int fo = defrag.new2old_faces[fn];
        Mesh_connectivity::Face_iterator face = connectivity().face_at ( fo );
        int n_verts = 0;
        Mesh_connectivity::Half_edge_iterator he_end = face.half_edge();
        Mesh_connectivity::Half_edge_iterator he = face.half_edge();

        do {
            face_verts_cache[n_verts] = defrag.old2new_vertices[he.origin().index()];
            ++n_verts;
            he = he.next();
        } while ( !he.is_equal ( he_end ) );

        fprintf ( fl, "%d ", n_verts );

        for ( int voffset = 0; voffset < n_verts; ++voffset ) {
            //
            fprintf ( fl, "%d ", face_verts_cache[voffset] );
        }

        fprintf ( fl, "\n" );
    }

    fprintf ( fl, "\n" );

    // write the face types
    fprintf ( fl, "CELL_TYPES %d \n", connectivity().n_active_faces() );

    for ( int f = 0; f < connectivity().n_active_faces(); f++ ) {
        fprintf ( fl, "7 \n" ); // VTK POLYGON
    }

    fprintf ( fl, "\n" );


    // Now write the pixel indices
    fprintf ( fl, "POINT_DATA %d \n", connectivity().n_active_vertices() );

//
    fprintf ( fl, "SCALARS %s int 1 \n", "pixel_point_index" );
    fprintf ( fl, "LOOKUP_TABLE default \n" );

    for ( int i = 0; i < connectivity().n_active_vertices(); i++ ) {
        int voidx = defrag.new2old_vertices[i];
        fprintf ( fl, "%d \n", voidx );
    }

// Write region colors?
    fprintf ( fl, "CELL_DATA %d \n", connectivity().n_active_faces() );

//
    fprintf ( fl, "SCALARS %s int 1 \n", "region_index" );
    fprintf ( fl, "LOOKUP_TABLE default \n" );

    for ( int i = 0; i < connectivity().n_active_faces(); i++ ) {
        fprintf ( fl, "%d \n", i );
    }

    fclose ( fl );
}

void
RasterImageConnectivity::enumerate_holes ( Mesh_connectivity& pixconn, std::vector<int>& half_edge_on_hole ) {
    assert_break ( pixconn.n_active_half_edges() == pixconn.n_total_half_edges() );
    assert_break ( pixconn.n_active_faces() == pixconn.n_total_faces() );

    // Start marking all the half-edges
    Bit_vector is_half_edge_marked ( pixconn.n_active_half_edges() );

    // A lambda to traverse a hole and mark all the half edges on it
    auto traverse_bounday = [&] ( const int seed_he ) {
        Mesh_connectivity::Half_edge_iterator he_end = pixconn.half_edge_at ( seed_he );
        Mesh_connectivity::Half_edge_iterator he = he_end;
        assert_break ( he.face().is_equal ( pixconn.hole() ) );

        int safe_guard = 0;
        const int safe_gaurd_cap = 10000;

        do {
            is_half_edge_marked.set ( he.index() );
            assert_break ( safe_guard < safe_gaurd_cap );
            ++safe_guard;
            he = he.next();
        } while ( !he.is_equal ( he_end ) );
    };

    // Loop over all half-edges, if not marked and boundary, it is a hole
    half_edge_on_hole.resize ( 0 );

    for ( int heid = 0 ; heid < pixconn.n_total_half_edges() ; ++heid ) {
        Mesh_connectivity::Half_edge_iterator he = pixconn.half_edge_at ( heid );

        if ( he.face().is_equal ( pixconn.hole() ) && ( !is_half_edge_marked[heid] ) ) {
            half_edge_on_hole.push_back ( heid );
            traverse_bounday ( heid );
        }
    }
}

// constructor
RasterImageConnectivity
RasterImageConnectivity::build ( const std::vector<Eigen::Matrix2Xi>& region_boundaries ) {

    RasterImageConnectivity ans;

    //
    // First find the bounding box of the image
    //
    {
        Eigen::AlignedBox<int, 2> bbox = get_bounding_box ( region_boundaries );
        // Make sure min is 0,0. This must have been done in get_pixel_regions()
        assert_break ( bbox.min() == Eigen::Vector2i ( 0, 0 ) );
        // Set the dims to the max
        ans._pixel_vertex_dims = bbox.max()+Eigen::Vector2i ( 1, 1 );
    }

    //
    // Now bould the connectivity by treating each region as a BIG polygon
    //
    {
        std::vector<int> polygon_verts;
        std::vector<int> polygon_xadj = {0};

        for ( int i = 0 ; i < ( int ) region_boundaries.size() ; ++i ) {
            polygon_xadj.push_back ( polygon_xadj.back() );
            const Eigen::Matrix2Xi& region_boundary = region_boundaries[i];

            for ( int j = 0 ; j < ( int ) region_boundary.cols() ; ++j ) {
                polygon_verts.push_back ( ans.pixel_vertex_pos_to_index ( region_boundary.col ( j ) ) );
                ++polygon_xadj.back();
            }
        }

        ans._connectivity.build_from_polygons ( ans.n_pixel_vertices(), polygon_verts, polygon_xadj );
    }

    //
    // Now get all the holes and add them
    // Even the outer boundary -- it should be fine
    //
    {
        Mesh_connectivity& pixconn = ans.connectivity();

        // Set number of regions
        ans._n_regions = ( int ) region_boundaries.size();

        // A lambda to modify the connectivity
        auto associate_face_and_he = [&] ( const int seed_he, const int parent_face ) {
            Mesh_connectivity::Half_edge_iterator he_end = pixconn.half_edge_at ( seed_he );
            Mesh_connectivity::Half_edge_iterator he = he_end;
            Mesh_connectivity::Face_iterator face = pixconn.face_at ( parent_face );
            assert_break ( he.face().is_equal ( pixconn.hole() ) );

            int safe_guard = 0;
            const int safe_gaurd_cap = 10000;

            // Make the association
            face.data().half_edge = he.index();

            do {
                he.data().face = face.index();
                assert_break ( safe_guard < safe_gaurd_cap );
                ++safe_guard;
                he = he.next();
            } while ( !he.is_equal ( he_end ) );
        };

        // Get all the holes
        std::vector<int> he_on_holes;
        enumerate_holes ( pixconn,  he_on_holes );

        for ( int heid: he_on_holes ) {
            Mesh_connectivity::Face_iterator face =  pixconn.add_face();
            associate_face_and_he ( heid, face.index() );
        }

    } // End of creating the hole polygons


    // Now associate holes and boundary (polygons to regions)
    {
        Mesh_connectivity& pixconn = ans.connectivity();

        // No extract the points
        // Get the holes and polygons as CCW paths
        // And also find the areas
        std::vector<Eigen::Matrix2Xd> outer_boundaries;
        std::vector<double> outer_boundary_areas;
        std::vector<Eigen::Matrix2Xd> holes;
        std::vector<double> hole_areas;

        for ( int i = 0 ; i < ans.n_regions() ; ++i ) {
            outer_boundaries.push_back ( ans.get_polygon_points ( i ).cast<double>() );
            bool is_ccw;
            double area;
            polyvec::WindingNumber::compute_orientation ( outer_boundaries.back(), is_ccw, area );
            assert_break ( is_ccw );
            outer_boundary_areas.push_back ( area );
        }

        for ( int i = ans.n_regions()  ; i < ans.n_polygons() ; ++i ) {
            holes.push_back ( ans.get_polygon_points ( i ).cast<double>() );
            holes.back().rowwise().reverseInPlace(); // make CCW
            bool is_ccw;
            double area;
            polyvec::WindingNumber::compute_orientation ( holes.back(), is_ccw, area );
            assert_break ( is_ccw );
            hole_areas.push_back ( area );
        }

        //
        // Create lambdas that say if a polygon is inside another
        //  assumes that out and in are both pixel boundaries
        // and are ccw
        // Only works if the area of out is bigger than in
        //
        auto is_inside = [&] ( const Eigen::Matrix2Xd &out, const Eigen::Matrix2Xd &in ) {
            const Eigen::Vector2d mid = (in.col ( 0 ) + in.col ( 1 )) /2.;
            const Eigen::Vector2d tang = in.col ( 1 ) -  in.col ( 0 );
            const Eigen::Vector2d normal ( -tang.y(), tang.x() );
            const Eigen::Vector2d definitely_inside_in = normal/2. + mid;

            bool is_trustable;
            double winding;
            ::polyvec::WindingNumber::compute_winding ( out, definitely_inside_in, winding, is_trustable );
            assert_break ( is_trustable );
            return std::abs ( winding ) > 1e-5;
        };

        // Now call the association function
        std::vector<std::vector<int>> boundary_holes;
        associate_holes (
            ( int ) outer_boundaries.size(),
            ( int ) holes.size(),
        [&] ( const int bout, const int hin ) {
            if ( ( int ) outer_boundary_areas[bout] <= ( int ) hole_areas[hin] ) {
                return false;
            } else {
                return  is_inside ( outer_boundaries[bout], holes[hin] );
            }
        },
        [&] ( const int hout, const int hin ) {
            if ( ( int ) hole_areas[hout] <= ( int ) hole_areas[hin] ) {
                return false;
            } else {
                return  is_inside ( holes[hout], holes[hin] );
            }
        },
        boundary_holes );

        // Now build the adj and xadj
        ans._polygon_to_region.resize ( ans.n_polygons(), infinite_region );
        ans._region_to_polygon_xval.resize ( 1, 0 );
        ans._region_to_polygon_val.resize ( 0 );

        for ( int i = 0 ; i < ans.n_regions() ; ++i ) {
            const int region_id = i;
            ans._region_to_polygon_xval.push_back ( ans._region_to_polygon_xval.back() );

            // Don't forget yourself :)
            ans._polygon_to_region[region_id] = region_id;
            ++ans._region_to_polygon_xval.back();
            ans._region_to_polygon_val.push_back ( region_id );

            for ( int  j = 0 ; j < ( int ) boundary_holes[region_id].size() ; ++j ) {
                const int polygon_id = ans.n_regions() + boundary_holes[region_id][j];
                ans._polygon_to_region[polygon_id] = region_id;
                ++ans._region_to_polygon_xval.back();
                ans._region_to_polygon_val.push_back ( polygon_id );
            }
        }
    } // ALL DONE WITH REGION HOLE ASSOCIATION

    return ans;
}


NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
