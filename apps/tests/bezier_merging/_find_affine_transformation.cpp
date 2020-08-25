#include <iostream>
//
#include <polyvec/api.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/curve-tracer/bezier_merging.hpp>
//
#include <Eigen/Core>
#include <Eigen/Geometry>

NAMESPACE_BEGIN()
NAMESPACE_END()

NAMESPACE_BEGIN(polyvectest)
NAMESPACE_BEGIN(BezierMerging)

int find_affine_transformation(int, char**) {

  // Creat the points and then test the mapping
  Eigen::Matrix<double, 2,3> points0;
  Eigen::Matrix<double, 2,3> points1;
  Eigen::Matrix2d A;
  Eigen::Vector2d b;
  Eigen::Matrix2d Amanu;
  Eigen::Vector2d bmanu;

  points0.col(0) << 2,10;
  points0.col(1) << 5,7;
  points0.col(2) << 4,12;

  //
  // Identity test
  // 
  Amanu = Eigen::Matrix2d::Identity();
  bmanu = Eigen::Vector2d::Zero();
  points1 = (Amanu * points0).colwise() + bmanu;

  ::polyfit::BezierMerging::find_affine_transformation(points0, points1, A, b);
  assert_break( (A - Amanu).norm() < 1e-12 );
  assert_break( (b - bmanu).norm() < 1e-12 );

  //
  // Translation test
  // 
  Amanu = Eigen::Matrix2d::Identity();
  bmanu = Eigen::Vector2d(2, -1.48);
  points1 = (Amanu * points0).colwise() + bmanu;

  ::polyfit::BezierMerging::find_affine_transformation(points0, points1, A, b);
  assert_break( (A - Amanu).norm() < 1e-12 );
  assert_break( (b - bmanu).norm() < 1e-12 );

  //
  // A test
  // 
  Amanu << 5, 12, -1, 30;
  bmanu = Eigen::Vector2d::Zero();
  points1 = (Amanu * points0).colwise() + bmanu;

  ::polyfit::BezierMerging::find_affine_transformation(points0, points1, A, b);
  assert_break( (A - Amanu).norm() < 1e-12 );
  assert_break( (b - bmanu).norm() < 1e-12 );

  //
  // together
  // 
  Amanu << 5, 12, -1, 30;
  bmanu << 9, -3;
  points1 = (Amanu * points0).colwise() + bmanu;

  ::polyfit::BezierMerging::find_affine_transformation(points0, points1, A, b);
  assert_break( (A - Amanu).norm() < 1e-12 );
  assert_break( (b - bmanu).norm() < 1e-12 );

  return EXIT_FAILURE;
}

NAMESPACE_END(BezierMerging)
NAMESPACE_END(polyvectest)
