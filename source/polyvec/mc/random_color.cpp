#include <cstdlib>

#include <polyvec/mc/random_color.hpp>

// Not that random -- but also not too bad
// https://stackoverflow.com/questions/43044/algorithm-to-randomly-generate-an-aesthetically-pleasing-color-palette

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(mc)

Eigen::Vector4d
random_color()
{    
    Eigen::Vector3i mix(255,255,255);
    Eigen::Vector3i rgb(rand() % 256,rand() % 256,rand() % 256);
    Eigen::Vector4d  ans;
    ans <<  (mix+rgb).cast<double>()/2./255. , 1. ;

    return ans;
}


NAMESPACE_END(mc)
NAMESPACE_END(polyfit)

