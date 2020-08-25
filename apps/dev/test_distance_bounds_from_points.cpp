// Polyvec
#include <polyvec/core/types.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/path.hpp>

// libc
#include <stdio.h>

using namespace Eigen;
using namespace polyfit;
// using namespace polyvec;

bool test_distance_bounds_from_points(
	// Connecting first from last point
	const mat2x& P,

	// Expected value
	const double d_min, const double d_max
) {
	printf("----  TEST DISTANCE BOUNDS -----\n");
	for (Index i = 0; i < P.cols(); ++i) {
		printf("%f %f\n", P(0, i), P(1, i));
	}

	mat2 E;
	E.col(0) = P.col(0);
	E.col(1) = P.col(P.cols() - 1);

	vec2 error = PathUtils::distance_bounds_from_points(P, E, vec2i(0, P.cols() - 1), true);
	printf("error %f %f\n", error(0), error(1));

	if (abs(d_min - error.minCoeff()) > PF_EPS || abs(d_max - error.maxCoeff()) > PF_EPS) {
		printf("FAIL! expected d_min %f\n", d_min);
		printf("FAIL! expected d_max %f\n", d_max);
		// return 0;
	}

	return 1;
}

/*
	TODO: transform the polylines
*/
#define TEST(expected_min, expected_max) if (!test_distance_bounds_from_points(B, expected_min, expected_max)) { return EXIT_FAILURE; }
int main(int argc, char* argv[]) {
	mat2x B;

	// 
	B.resize(2, 4);
	B.col(0) = vec2(0., 0.);
	B.col(1) = vec2(1., 0.);
	B.col(2) = vec2(1., 1.);
	B.col(3) = vec2(2., 1.);
	TEST(.25, .25); // 1 - 3/4

	//
	B.resize(2, 6); 
	B.col(0) = vec2(0., 0.);
	B.col(1) = vec2(1., 0.);
	B.col(2) = vec2(2., 0.);
	B.col(3) = vec2(2., 1.);
	B.col(4) = vec2(3., 1.);
	B.col(5) = vec2(4., 1.);
	TEST(0.375, 0.375); // 1 - 5/8

	//
	B.resize(2, 8);
	B.col(0) = vec2(0., 0.);
	B.col(1) = vec2(1., 0.);
	B.col(2) = vec2(2., 0.);
	B.col(3) = vec2(3., 0.);
	B.col(4) = vec2(3., 1.);
	B.col(5) = vec2(4., 1.);
	B.col(6) = vec2(5., 1.);
	B.col(7) = vec2(6., 1.);
	TEST(0.416667, 0.416667); // 1 - 7/12

	// 
	B.resize(2, 6);
	B.col(0) = vec2(0., 0.);
	B.col(1) = vec2(1., 0.);
	B.col(2) = vec2(1., 1.);
	B.col(3) = vec2(1., 2.);
	B.col(4) = vec2(1., 3.);
	B.col(5) = vec2(2., 3.);
	TEST(0.615385, 0.615385); // 1 - 7/12

	// 
	B.resize(2, 5);
	B.col(0) = vec2(0., 0.);
	B.col(1) = vec2(1., 0.);
	B.col(2) = vec2(2., 0.);
	B.col(3) = vec2(2., 1.);
	B.col(4) = vec2(2., 2.);
	TEST(0., 1.); // .5 + .5

	// 
	B.resize(2, 7);
	B.col(0) = vec2(0., 0.);
	B.col(1) = vec2(1., 0.);
	B.col(2) = vec2(2., 0.);
	B.col(3) = vec2(3., 0.);
	B.col(4) = vec2(3., 1.);
	B.col(5) = vec2(3., 2.);
	B.col(6) = vec2(3., 3.);
	TEST(0., 1.5); // .5 + 1

	return EXIT_SUCCESS;
}