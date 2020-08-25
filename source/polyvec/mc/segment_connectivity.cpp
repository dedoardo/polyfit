#include <cstdio>
#include <memory>
#include <queue>

#include <polyvec/mc/segment_connectivity.hpp>
#include <polyvec/mc/bit_vector.hpp>

#include <Eigen/Geometry>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )


constexpr int SegmentConnectivity::infinite_region;


int
SegmentConnectivity::n_segments() {
    return _n_segments;
}

RasterImageConnectivity&
SegmentConnectivity::raster_image_connectivity() {
    return *_raster_image_connectivity;
}

int
SegmentConnectivity::segment_to_pixelhe_begin ( const int segment_id ) {
    return _segment_to_pixelhe_begin[ segment_id ];
}

int
SegmentConnectivity::segment_to_pixelhe_end ( const int segment_id ) {
    return _segment_to_pixelhe_end[ segment_id ];
}

int
SegmentConnectivity::pixelhe_to_segment ( const int pixelhe_id ) {
    return _pixelhe_to_segment[ pixelhe_id ];
}

void
SegmentConnectivity::get_segment_points ( const int segment_id, Eigen::Matrix2Xi& points, bool& is_closed ) {
    Mesh_connectivity& pixconn = raster_image_connectivity().connectivity();
    //
    const int  heid_begin =  segment_to_pixelhe_begin ( segment_id );
    const int  heid_end =  segment_to_pixelhe_end ( segment_id );
    //
    Mesh_connectivity::Half_edge_iterator he_begin = pixconn.half_edge_at ( heid_begin );
    Mesh_connectivity::Half_edge_iterator he_end = pixconn.half_edge_at ( heid_end );
    Mesh_connectivity::Half_edge_iterator he;

    // Set is closed
    is_closed = he_begin.is_equal ( he_end );

    // Count the number of vertices
    int n_verts = 0;
    he = he_begin;

    do {
        ++n_verts;
        he = he.next();
    } while ( !he.is_equal ( he_end ) );

    if ( !is_closed ) {
        ++n_verts;
    }

    // Now store the vert positions
    points.resize ( 2, n_verts );
    int offset = 0;
    he = he_begin;

    do {
        Eigen::Vector2i pos = raster_image_connectivity().pixel_vertex_index_to_pos ( he.origin().index() );
        points.col ( offset ) = pos;
        ++offset;
        he = he.next();
    } while ( !he.is_equal ( he_end ) );

    if ( !is_closed ) {
        Eigen::Vector2i pos = raster_image_connectivity().pixel_vertex_index_to_pos ( he.origin().index() );
        points.col ( offset ) = pos;
    }

}

void SegmentConnectivity::get_incident_polygons_for_segment(const int segment_id, int& region0, int& region1) {
    int he_begin = segment_to_pixelhe_begin(segment_id);
    Mesh_connectivity::Half_edge_iterator he_it = _raster_image_connectivity->connectivity().half_edge_at(he_begin);
    region0 = he_it.face().index();
    region1 = he_it.twin().face().index();
}

SegmentConnectivity
SegmentConnectivity::build ( RasterImageConnectivity& raster_image_connectivity, const std::vector<JunctionType>& junction_types ) {
    assert_break ( ( int ) junction_types.size() == raster_image_connectivity.n_pixel_vertices() );
    Mesh_connectivity& pixconn = raster_image_connectivity.connectivity();

    // Mappings
    std::vector<int> pixelhe_to_segment ( pixconn.n_active_half_edges(), infinite_region );
    std::vector<int> segment_to_pixelhe_begin ( 0 );
    std::vector<int> segment_to_pixelhe_end ( 0 );
    int n_segments = 0;

    // Detect 3 or 4 junction
    auto is_34junction = [&junction_types] ( const int jt ) {
        return ( junction_types[jt] == JUNCTION_3_SHARP ) || ( junction_types[jt] == JUNCTION_4_SHARP ) ;
    };

    //
    // Push all he that begin at a junction into a stack
    //
    std::queue<int> to_inspect;

    for ( int heid = 0 ; heid < pixconn.n_total_half_edges() ; ++heid ) {
        Mesh_connectivity::Half_edge_iterator he = pixconn.half_edge_at ( heid );
        Mesh_connectivity::Vertex_iterator vert = he.origin();
        Mesh_connectivity::Face_iterator face = he.face();
        assert_break ( vert.is_active() );


        bool is_okay = true;
        is_okay = is_okay && ( !face.is_equal ( pixconn.hole() ) ); // must check before calling last one
        is_okay = is_okay && is_34junction ( vert.index() ) ;
        is_okay = is_okay && ( raster_image_connectivity.get_polygon_type ( face.index() ) == RasterImageConnectivity::POLYGON_OUTER_BOUNDARY );

        if ( is_okay ) {
            to_inspect.push ( he.index() );
        }
    }

    //
    // Traverse half edges and create a segment for them
    //
    while ( !to_inspect.empty() ) {
        const int candidate = to_inspect.front();
        to_inspect.pop();

        //
        // Non-hole half-edges can be visited twice. Prevent it
        if ( pixelhe_to_segment[ candidate ] != infinite_region ) {
            continue;
        }

        //
        // Now create the segment
        const int segment_id = n_segments;
        ++n_segments;
        //
        Mesh_connectivity::Half_edge_iterator he = pixconn.half_edge_at ( candidate );
        Mesh_connectivity::Half_edge_iterator he_begin = he;
        bool should_continue = false;

        int safe_gaurd_counter = 0;
        const int safe_gaurd_max = 10000;

        do {
            // Set maps
            assert_break ( raster_image_connectivity.get_polygon_type ( he.face().index() ) == RasterImageConnectivity::POLYGON_OUTER_BOUNDARY );
            pixelhe_to_segment[ he.index() ] = segment_id;
            pixelhe_to_segment[ he.twin().index() ] = segment_id;
            // Advance
            he = he.next();
            // Should continue ?
            should_continue = true;
            should_continue = should_continue && ( !he.is_equal ( he_begin ) );
            should_continue = should_continue &&  ( !is_34junction ( he.origin().index() ) );
            // Safe gaurd
            ++safe_gaurd_counter;
            assert_break ( safe_gaurd_counter < safe_gaurd_max );
        } while ( should_continue );

        //
        segment_to_pixelhe_begin.push_back ( he_begin.index() );
        segment_to_pixelhe_end.push_back ( he.index() );
    } // done creating segments

    //
    // Now push all non-visited half-edges to a queue
    //
    to_inspect = std::queue<int>();

    for ( int heid = 0 ; heid < pixconn.n_total_half_edges() ; ++heid ) {
        Mesh_connectivity::Half_edge_iterator he = pixconn.half_edge_at ( heid );
        Mesh_connectivity::Vertex_iterator vert = he.origin();
        Mesh_connectivity::Face_iterator face = he.face();

        assert_break ( vert.is_active() );

        bool is_okay = true;
        is_okay = is_okay && ( !face.is_equal ( pixconn.hole() ) ); // must check before next one
        is_okay = is_okay && ( raster_image_connectivity.get_polygon_type ( face.index() ) == RasterImageConnectivity::POLYGON_OUTER_BOUNDARY );
        is_okay = is_okay && ( pixelhe_to_segment[heid] == infinite_region );

        if ( is_okay ) {
            assert_break ( !is_34junction ( vert.index() ) );
            to_inspect.push ( he.index() );
        }
    }

    //
    // Now extract closed loops
    //
    while ( !to_inspect.empty() ) {
        // Get a non-marked candidate
        const int candidate = to_inspect.front();
        to_inspect.pop();

        if ( pixelhe_to_segment[ candidate ] != infinite_region ) {
            continue;
        }

        // Create a segment
        const int segment_id = n_segments;
        ++n_segments;
        //
        // Now create the segment
        Mesh_connectivity::Half_edge_iterator he = pixconn.half_edge_at ( candidate );
        Mesh_connectivity::Half_edge_iterator he_begin = he;
        bool should_continue = false;

        int safe_gaurd_counter = 0;
        const int safe_gaurd_max = 10000;

        do {
            // Set maps
            assert_break ( raster_image_connectivity.get_polygon_type ( he.face().index() ) == RasterImageConnectivity::POLYGON_OUTER_BOUNDARY );
            assert_break ( !is_34junction ( he.origin().index() ) );
            pixelhe_to_segment[ he.index() ] = segment_id;
            pixelhe_to_segment[ he.twin().index() ] = segment_id;
            // Advance
            he = he.next();
            // Should continue ?
            should_continue = true;
            should_continue = should_continue && ( !he.is_equal ( he_begin ) );
            // Safe gaurd
            ++safe_gaurd_counter;
            assert_break ( safe_gaurd_counter < safe_gaurd_max );
        } while ( should_continue );

        //
        segment_to_pixelhe_begin.push_back ( he_begin.index() );
        segment_to_pixelhe_end.push_back ( he.index() );
    } // done creating segments


    //
    // Make sure all half-edges were visited
    //
    for ( int heid = 0 ; heid < pixconn.n_total_half_edges() ; ++heid ) {
        assert_break ( pixelhe_to_segment[heid] != infinite_region );
    }

    //
    // Create output
    //
    SegmentConnectivity buildee;
    buildee._n_segments = n_segments;
    buildee._pixelhe_to_segment = pixelhe_to_segment;
    buildee._segment_to_pixelhe_end = segment_to_pixelhe_end;
    buildee._segment_to_pixelhe_begin = segment_to_pixelhe_begin;
    buildee._raster_image_connectivity = &raster_image_connectivity;
    return buildee;
}


NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )