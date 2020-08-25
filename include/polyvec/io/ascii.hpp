#pragma once

#include <polyvec/core/types.hpp>
#include <polyvec/core/macros.hpp>

#include <string>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(IO)

void read_raster(const std::string& addr, Eigen::Matrix2Xd& points);
bool read_polygon(const std::string& addr, std::vector<Eigen::Index>& indices);
void read_polygon_g0(const std::string& addr, std::vector<Eigen::Index>& indices);
void read_manual_polygon(const std::string& raster_addr, const std::string& polygon_addr, Eigen::Matrix2Xd& raster, Eigen::Matrix2Xd& polygon, std::vector<Eigen::Index>& polygon_vertices);

void read_baseline_boundary(const std::string& uri, mat2x& B);
void read_baseline_polygon(const std::string& uri, vecXi& P);

bool write_baseline_boundary (const std::string& uri, const mat2x& B);
bool write_baseline_polygon(const std::string& uri, const vecXi& P);

// creates a unique string from a path, considering the filename and last-directory
std::string model_id(const std::string& image_uri);

// writes the samples transposed into the file
bool write_points(const std::string& uri, const Eigen::MatrixXd& P);

// reads space separated error parameters from the first line in the file
bool read_parameter_limits(const std::string& uri,
	double& lim_accuracy,
	double& lim_smoothness,
	double& lim_continuity,
	double& lim_inflection);

NAMESPACE_END(IO)
NAMESPACE_END(polyfit)