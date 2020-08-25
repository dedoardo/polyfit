#include <set>

#include <Eigen/Geometry>

#include <polyvec/mc/junction.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>

NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )

std::string
junction_type_as_str ( const JunctionType t ) {
    switch ( t ) {
#define TAKECAREOF(XX) case XX: return #XX
        TAKECAREOF ( JUNCTION_KOPF );
        TAKECAREOF ( JUNCTION_3_SHARP );
        TAKECAREOF ( JUNCTION_4_SHARP );
        TAKECAREOF ( JUNCTION_INVALID );

    default:
        assert_break ( 0 );
#undef TAKECAREOF
    }
}

void
find_junctions ( /*const*/ RasterImageConnectivity& image, std::vector<JunctionType>& junction_types ) {
    junction_types.resize ( image.connectivity().n_total_vertices() );
    std::fill ( junction_types.begin(), junction_types.end(), JUNCTION_INVALID );

    // Counter number of outgoing half-edges
    std::vector<int> n_outgoing_he ( image.connectivity().n_total_vertices() );

    for ( int i = 0 ; i < image.connectivity().n_total_half_edges() ; ++ i ) {
        Mesh_connectivity::Half_edge_iterator hiter = image.connectivity().half_edge_at ( i );
        ++n_outgoing_he[ hiter.origin().index() ];
    }

    for ( int vid = 0 ; vid < image.n_pixel_vertices() ; ++vid ) {
        Mesh_connectivity::Vertex_iterator viter = image.connectivity().vertex_at ( vid );

        if ( viter.is_active() ) {
            JunctionType jtype;
            is_junction ( image, vid, n_outgoing_he[vid], jtype );
            junction_types[viter.index()] = jtype ;
        } // end of active
    } // end of vertcies
}

void
is_junction ( /*const*/ RasterImageConnectivity& image, const int vertex_id, const int n_outgoing_he, JunctionType& junction_type ) {
    Mesh_connectivity::Vertex_iterator viter = image.connectivity().vertex_at ( vertex_id );

    if ( !viter.is_active() ) {
        junction_type = JUNCTION_INVALID;
    } else {
        int n_visited_edges = 0;

        Mesh_connectivity::Vertex_ring_iterator ring_iter = image.connectivity().vertex_ring_at ( vertex_id );
        Mesh_connectivity::Half_edge_iterator hiter = viter.half_edge();
        ( void ) hiter;
        Mesh_connectivity::Half_edge_iterator hiter_twin = hiter.twin();
        ( void ) hiter_twin;

        do {
            const int polygon_id = ring_iter.half_edge().face().index();
            // note that the outer boundary will return the infinite_region index
            const int region_id = image.polygon_to_region ( polygon_id );
            const int other_vertex_id = ring_iter.half_edge().origin().index();
            ( void ) other_vertex_id; // (unused);
            ++n_visited_edges;
        } while ( ring_iter.advance() );

        // Now find the type
        if ( n_visited_edges == 2 ) {
            if ( n_outgoing_he == 2 ) {
                junction_type = JUNCTION_INVALID;
            } else {
                junction_type = JUNCTION_KOPF;
            }
        } else if ( n_visited_edges == 3 ) {
            junction_type = JUNCTION_3_SHARP;
        }  else if ( n_visited_edges == 4 ) {
            junction_type = JUNCTION_4_SHARP;
        } else {
            assert_break ( 0 );
        }
    }

} // end of active


NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
