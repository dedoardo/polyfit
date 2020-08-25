#include <polyvec/curve-tracer/curve_fit_error.hpp>

using namespace polyvec;

void CurveFitError::combine(const CurveFitError& other)
{
	accuracy.combine(other.accuracy);
	curvature.combine(other.curvature);
}

void CurveFitError::reset()
{
	accuracy.reset();
	curvature.r_max = 0;
	curvature.r_min = std::numeric_limits<double>::infinity();
}