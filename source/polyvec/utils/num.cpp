// Polyvec
#include <polyvec/utils/num.hpp>
#include <polyvec/api.hpp>

// Eigen
#include <Eigen/LU>

// std
#include <cstdlib>

using namespace std;

NAMESPACE_BEGIN ( polyvec )
NAMESPACE_BEGIN ( Num )

double
determinant(Eigen::Ref<const Eigen::MatrixXd> mat) {
  return mat.determinant();
}
	
bool test_and_calculate_interval_overlap(
    const double i0, const double i1,
    const double j0, const double j1,
    double& overlap
) {
    if (j0 > i1 || i0 > j1) {
        return false;
    }

    overlap = min(i1, j1) - max(i0, j0);
    return overlap > -PF_EPS;
}

Eigen::MatrixXd 
solve_linear_system (
    const Eigen::Ref<const Eigen::MatrixXd> LHS,
    const Eigen::Ref<const Eigen::MatrixXd> RHS )  {
      assert_break(LHS.rows() == LHS.cols());
      assert_break(RHS.rows() == LHS.cols());
      assert_break(RHS.cols() == 1);
      assert_break(std::abs(LHS.determinant()) > 1e-12 );
      return LHS.lu().solve( RHS );
}

double smooth_probability_incr(double x, double zero_up_to, double one_beyond)
{
	if (x <= zero_up_to)
		return 0;
	if (x >= one_beyond)
		return 1;
	return 0.5 + 0.5 * std::cos(M_PI * (x - one_beyond) / (one_beyond - zero_up_to));
}

double smooth_probability_decr(double x, double one_up_to, double zero_beyond)
{
	return 1 - smooth_probability_incr(x, one_up_to, zero_beyond);
}

NAMESPACE_END ( Num )
NAMESPACE_END ( polyvec )

