#include <cassert>

#include <polyvec/mc/associate_holes.hpp>
#include <polyvec/mc/bit_vector.hpp>


NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( mc )


void
associate_holes (
    const int n_outer_boundaries,
    const int n_holes,
    const std::function<bool ( const int out, const int in ) > does_boundary_contain_hole,
    const std::function<bool ( const int out, const int in ) > does_hole_contain_hole,
    std::vector<std::vector<int>>& per_boundary_holes
) {

    //
    // BUILD HOLE PARENTS
    //
    // hole_parent[i] is the index of holes that contain hole i
    std::vector<std::vector<int>> hole_parents ( n_holes );

    for ( int i = 0 ; i < n_holes ; ++i ) {
        for ( int j = 0 ; j < n_holes ; ++j ) {
            if ( does_hole_contain_hole ( j, i ) ) {
                hole_parents[i].push_back ( j );
            }
        }
    }

    //
    // A lambda that returns the biggest parent of a hole
    // that is contained inside a polygon
    //
    auto find_contained_hole_ancestor = [hole_parents, does_boundary_contain_hole] ( const int boundy_id , const int seed_hole_id ) -> int {
        assert ( does_boundary_contain_hole ( boundy_id, seed_hole_id ) );
        int candiate_hole_id = seed_hole_id ;
        bool should_continue;

        do {
            should_continue = false;

            for ( const int parent_id : hole_parents[candiate_hole_id] ) {
                if ( does_boundary_contain_hole ( boundy_id, parent_id ) ) {
                    candiate_hole_id = parent_id;
                    should_continue = true;
                    break;
                }
            }
        } while ( should_continue );

        return candiate_hole_id;
    };


    //
    // Now for each boundary find its holes
    //
    Bit_vector is_hole_visited ( n_holes );
    per_boundary_holes.resize ( n_outer_boundaries );

    for ( int boundary_id = 0 ; boundary_id < n_outer_boundaries ; ++boundary_id )  {
        is_hole_visited.reset_all();
        per_boundary_holes[boundary_id].resize ( 0 );

        for ( int hole_id = 0 ; hole_id < n_holes ; ++hole_id ) {
            if ( does_boundary_contain_hole ( boundary_id, hole_id ) ) {
                const int ancestor_hole = find_contained_hole_ancestor ( boundary_id,  hole_id);

                if ( !is_hole_visited[ancestor_hole] ) {
                    per_boundary_holes[boundary_id].push_back ( ancestor_hole );
                    is_hole_visited.set( ancestor_hole );
                }
            }
        }

        // Just sort it for unification
        std::sort(per_boundary_holes[boundary_id].begin() , per_boundary_holes[boundary_id].end());
    }

}

NAMESPACE_END ( mc )
NAMESPACE_END ( polyfit )
