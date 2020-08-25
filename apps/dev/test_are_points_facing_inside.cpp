#include <polyvec/core/types.hpp>
#include <polyvec/regularity/are_points_facing_inside.hpp>
#include <polyvec/utils/matrix.hpp>

using namespace polyvec;
using namespace polyfit;

#define TEST(v) if (!(v)) { return EXIT_FAILURE; }

int main() {
    {
        mat2x P(2, 6);
        P.col(0) = vec2(0., 0.);
        P.col(1) = vec2(1., 1.);
        P.col(2) = vec2(0., 2.);
        P.col(3) = vec2(3., 2.);
        P.col(4) = vec2(2., 1.);
        P.col(5) = vec2(3., 0.);

        TEST(are_points_facing_inside(P, 1, 4));
    }
    
    {
        mat2x P(2, 8);
        P.col(0) = vec2(0., 0.);
        P.col(1) = vec2(0., 2.);
        P.col(2) = vec2(1., 2.);
        P.col(3) = vec2(1., 1.);
        P.col(4) = vec2(2., 1.);
        P.col(5) = vec2(2., 2.);
        P.col(6) = vec2(3., 2.);
        P.col(7) = vec2(3., 0.);
        
        TEST(are_points_facing_inside(P, 0, 2));
        TEST(are_points_facing_inside(P, 0, 3));
        TEST(are_points_facing_inside(P, 0, 4));
        TEST(!are_points_facing_inside(P, 0, 5));
        TEST(!are_points_facing_inside(P, 0, 6));
    }

    return EXIT_SUCCESS;
}
