#include <polyvec/mc/associate_holes.hpp>
#include <iostream>


NAMESPACE_BEGIN()
NAMESPACE_END()

NAMESPACE_BEGIN ( polyvectest )
NAMESPACE_BEGIN ( mc )

// void associate_holes(
//   const int n_outer_boundaries,
//   const int n_holes,
//   const std::function<bool (const int out, const int in)> does_boundary_contain_hole,
//   const std::function<bool (const int out, const int in)> does_hole_contain_hole,
//   std::vector<std::vector<int>>& per_boundary_holes
//   );

void associate_holes(int , char**) {
    // se associate_holes.svg
    const int n_outer_boundaries = 5;
    const int n_holes = 4;

    std::vector< std::vector<int> > in_bdry_holes  = {
        {},
        {1},
        {},
        {0},
        {0,1,3}
    };
    auto does_boundary_contain_hole = [in_bdry_holes] ( const int out, const int in ) {
        return ( std::find ( in_bdry_holes[out].begin(), in_bdry_holes[out].end(), in ) != in_bdry_holes[out].end() );
    };

    std::vector< std::vector<int> > in_hole_holes = {
        {},
        {},
        {0, 1, 3},
        {0, 1}
    };
    auto does_hole_contain_hole = [in_hole_holes] ( const int out, const int in ) {
        return ( std::find ( in_hole_holes[out].begin(), in_hole_holes[out].end(), in ) != in_hole_holes[out].end() );
    };

    std::vector<std::vector<int>> per_boundary_holes;
    ::polyfit::mc::associate_holes ( n_outer_boundaries, n_holes, does_boundary_contain_hole, does_hole_contain_hole, per_boundary_holes );

    std::cout << "p0: " << Eigen::RowVectorXi::ConstMapType( per_boundary_holes[0].data(), per_boundary_holes[0].size() ) << std::endl;
    std::cout << "p1: " << Eigen::RowVectorXi::ConstMapType( per_boundary_holes[1].data(), per_boundary_holes[1].size() )  << std::endl;
    std::cout << "p2: " << Eigen::RowVectorXi::ConstMapType( per_boundary_holes[2].data(), per_boundary_holes[2].size() )  << std::endl;
    std::cout << "p3: " << Eigen::RowVectorXi::ConstMapType( per_boundary_holes[3].data(), per_boundary_holes[3].size() )  << std::endl;
    std::cout << "p4: " << Eigen::RowVectorXi::ConstMapType( per_boundary_holes[4].data(), per_boundary_holes[4].size() )  << std::endl;
}

NAMESPACE_END ( mc )
NAMESPACE_END ( polyvectest )