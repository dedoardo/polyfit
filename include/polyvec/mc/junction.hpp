#pragma once

#include <polyvec/api.hpp>

#include <Eigen/Core>

#include <string>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

class RasterImageConnectivity;

// Different types of junctions in an image
enum JunctionType {
  // This is not really important as it would be treated just like a 
  // normal point adjacent to two colored regions.
  JUNCTION_KOPF = 0,
  // A point where three colors meet.
  JUNCTION_3_SHARP ,
  // A point where four colors meet.
  JUNCTION_4_SHARP,
  // nothing
  JUNCTION_INVALID
};

// convert junction type to a string
std::string junction_type_as_str(const JunctionType);

// Get an image and then find the junctions and their type
// the indices are the same pixel_vertex_index of the RasterImageConnectivity
// TEST: apps_test/mc/_junctions.cpp
void find_junctions(/*const*/ RasterImageConnectivity&, std::vector<JunctionType>& junction_types);

// ======== PRIVATE API -- only used internally

void is_junction(/*const*/ RasterImageConnectivity&, const int vertex_id, const int n_outgoing_he, JunctionType& junction_type);


NAMESPACE_END(mc)
NAMESPACE_END(polyfit)
