#pragma once

#include <polyvec/api.hpp>

#include <Eigen/Core>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

// ======== PUBLIC API 

// Read the pixel regions along with their color from a bmp or png file
// also works with png, the name is misleading.
void get_pixel_regions_from_bmp(
  const std::string fname,
  std::vector<Eigen::Matrix2Xi>& regions,
  std::vector<Eigen::Vector4d>& colors);

// ======== Private API -- exposed for testing

// processes the regions
// since it discards the CW regions, out2in is the mapping between
// input/output regions.
void post_process_pixel_regions(
  const std::vector<Eigen::Matrix2Xd> &in,
  std::vector<Eigen::Matrix2Xd> &out,
  std::vector<int> & out2in ); 

NAMESPACE_END(mc)
NAMESPACE_END(polyfit)
