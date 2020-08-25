#pragma once

#include <polyvec/api.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// Associates the holes to the outer boundaries in a multicolored image
// Input:
//   n_outer_boundaries, number of outer boundaries
//   n_holes, number of holes
//   does_boundary_contain_hole, indicates if a outer boundary contains a hole
//   does_hole_contain_hole, indicates if a hole contains a hole
// Output:
//  per_boundary_holes, the biggest holes for each outer boundary
// The test is here: 
//   apps/apps_test/mc/_associate_holes.cpp
//   apps/apps_test/mc/associate-holes.svg
void associate_holes(
  const int n_outer_boundaries,
  const int n_holes,
  const std::function<bool (const int out, const int in)> does_boundary_contain_hole,
  const std::function<bool (const int out, const int in)> does_hole_contain_hole,
  std::vector<std::vector<int>>& per_boundary_holes
  );

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

