// polyvec
#include <polyvec/core/types.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/geometry/path.hpp>

// libc++
#include <cstdlib>

#define TEST(expr) if (!(expr)) { fprintf(stderr, "failed (%s:%d) " #expr, __FILE__, __LINE__); return EXIT_FAILURE; }

using namespace polyvec;
using namespace polyfit;

int main(int argc, char* argv[]) {
	{
		mat2x P(2, 4);
		P.col(0) = vec2(0., 0.);
		P.col(1) = vec2(1., 0.);
		P.col(2) = vec2(1., 1.);
		P.col(3) = vec2(0., 1.);

		bool is_ccw;
		WindingNumber::compute_orientation(P, is_ccw);
		printf("---------\n");
		printf("(winding) orientation %d\n", is_ccw ? -1 : +1);
		//printf("(path)    orientation %d\n", PathUtils::compute_orientation(P));

		std::vector<int> C;
		PathUtils::compute_convexities(P, C);
		for (int i = 0; i < (int)C.size(); ++i) {
			printf("cvx(%d) = %d\n", i, C[i]);
		}

		TEST(C[0] == 1);
		TEST(C[1] == 1);
		TEST(C[2] == 1);
		TEST(C[3] == 1);
	}

	{
		mat2x P(2, 3);
		P.col(0) = vec2(0., 0.);
		P.col(1) = vec2(1., 1.);
		P.col(2) = vec2(2., 0.);

		bool is_ccw;
		WindingNumber::compute_orientation(P, is_ccw);
		printf("---------\n");
		printf("(winding) orientation %d\n", is_ccw ? -1 : +1);
		//printf("(path)    orientation %d\n", PathUtils::compute_orientation(P));

		std::vector<int> C;
		PathUtils::compute_convexities(P, C);
		for (int i = 0; i < (int)C.size(); ++i) {
			printf("cvx(%d) = %d\n", i, C[i]);
		}

		TEST(C[0] == 1);
		TEST(C[1] == 1);
		TEST(C[2] == 1);
	}

	{
		mat2x P(2, 3);
		P.col(0) = vec2(2., 0.);
		P.col(1) = vec2(1., 1.);
		P.col(2) = vec2(0., 0.);

		bool is_ccw;
		WindingNumber::compute_orientation(P, is_ccw);
		printf("--------\n");
		printf("(winding) orientation %d\n", is_ccw ? -1 : +1);

		std::vector<int> C;
		PathUtils::compute_convexities(P, C);
		for (int i = 0; i < (int)C.size(); ++i) {
			printf("cvx(%d) = %d\n", i, C[i]);
		}

		TEST(C[0] == 1);
		TEST(C[1] == 1);
		TEST(C[2] == 1);
	}

	{
		mat2x P(2, 5);
		P.col(0) = vec2(0., 0.);
		P.col(1) = vec2(2., 0.);
		P.col(2) = vec2(2., 1.);
		P.col(3) = vec2(1., .5);
		P.col(4) = vec2(0., 1.);

		bool is_ccw;
		WindingNumber::compute_orientation(P, is_ccw);
		printf("--------\n");
		printf("(winding) orientation %d\n", is_ccw ? -1 : +1);

		std::vector<int> C;
		PathUtils::compute_convexities(P, C);
		for (int i = 0; i < C.size(); ++i) {
			printf("cvx(%d) = %d\n", i, C[i]);
		}

		TEST(C[0] == 1);
		TEST(C[1] == 1);
		TEST(C[2] == -1);
		TEST(C[3] == 1);
		TEST(C[4] == 1);
	}

	{
		mat2x P(2, 5);
		P.col(0) = vec2(0., 1.);
		P.col(1) = vec2(1., .5);
		P.col(2) = vec2(2., 1.);
		P.col(3) = vec2(2., 0.);
		P.col(4) = vec2(0., 0.);

		bool is_ccw;
		WindingNumber::compute_orientation(P, is_ccw);
		printf("--------\n");
		printf("(winding) orientation %d\n", is_ccw ? -1 : +1);

		std::vector<int> C;
		PathUtils::compute_convexities(P, C);
		for (int i = 0; i < C.size(); ++i) {
			printf("cvx(%d) = %d\n", i, C[i]);
		}

		TEST(C[0] == 1);
		TEST(C[1] == 1);
		TEST(C[2] == -1);
		TEST(C[3] == 1);
		TEST(C[4] == 1);
	}

	return EXIT_FAILURE;
}