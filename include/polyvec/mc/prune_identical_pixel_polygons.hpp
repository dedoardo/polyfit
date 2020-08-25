#pragma once

#include <vector>

#include <polyvec/api.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// NOTE: This is not being used. Can be thrown away.
//
// Given a set of regions, prunes the identical ones.
// out2in is the mapping between output/input regions.
void prune_identical_pixel_polygons(
  const std::vector<Eigen::Matrix2Xi>& ,
  std::vector<Eigen::Matrix2Xi>& out , 
  std::vector<int>& out2in);
  

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

