#ifndef POLYVEC_DERIV_TEST_IS_INCLUDED
#define POLYVEC_DERIV_TEST_IS_INCLUDED

#include <cstring>
#include <iomanip>
#include <iostream>

#include <Eigen/Core>

namespace polyvec{
  
namespace derivtest {
//
// Finite difference test
//

inline void
run(const Eigen::VectorXd &state, const Eigen::VectorXd &delta,
    std::function<Eigen::VectorXd(const Eigen::VectorXd &)> eval,
    std::function<Eigen::VectorXd(const Eigen::VectorXd &)> eval_jacobian,
    std::ostream &cout = std::cout, int n_halving = 5, double h = 0.5) {
  char dump[4096];
  double diff0(0), diff1(0);

  cout << std::setw(8) << "n" << std::setw(16) << "DIFF0" << std::setw(16)
       << "DIFF1" << std::endl;

  Eigen::VectorXd value0 = eval(state);
  for (int ih = 0; ih < n_halving; ++ih) {
    h /= 2.;

    Eigen::VectorXd hdelta = h * delta;
    Eigen::VectorXd value1exact = eval(state + hdelta);

    eval(state);
    Eigen::VectorXd linear_correction = eval_jacobian(hdelta);
    Eigen::VectorXd value1lin = value0 + linear_correction;

    double diff0new = (value0 - value1exact).norm();
    double diff1new = (value1lin - value1exact).norm();
    sprintf(dump, "%8d %8.4e(%4.1f) %8.4e(%4.1f)  \n", ih, diff0new,
            diff0 / diff0new, diff1new, diff1 / diff1new);
    cout << dump;
    cout.flush();
    diff0 = diff0new;
    diff1 = diff1new;
  }
} // All done

inline void
run(const Eigen::VectorXd &state, const Eigen::VectorXd &delta,
    std::function<Eigen::VectorXd(const Eigen::VectorXd &)> eval,
    std::function<Eigen::VectorXd(const Eigen::VectorXd &)> eval_jacobian,
    std::function<Eigen::VectorXd(const Eigen::VectorXd &)> eval_hessian,
    std::ostream &cout = std::cout, int n_halving = 5, double h = 0.5) {
  char dump[4096];
  double diff0(0), diff1(0), diff2(0);

  cout << std::setw(8) << "n" << std::setw(16) << "DIFF0" << std::setw(16)
       << "DIFF1" << std::setw(16) << "DIFF2\n";

  Eigen::VectorXd value0 = eval(state);
  for (int ih = 0; ih < n_halving; ++ih) {
    h /= 2.;

    Eigen::VectorXd hdelta = h * delta;
    Eigen::VectorXd value1exact = eval(state + hdelta);

    eval(state);
    Eigen::VectorXd value1lin = value0 + eval_jacobian(hdelta);
    Eigen::VectorXd value1quad = value1lin + 0.5 * eval_hessian(hdelta);

    double diff0new = (value0 - value1exact).norm();
    double diff1new = (value1lin - value1exact).norm();
    double diff2new = (value1quad - value1exact).norm();
    sprintf(dump, "%8d %8.4e(%4.0f) %8.4e(%4.0f)  %8.4e(%4.0f)  \n", ih,
            diff0new, diff0 / diff0new, diff1new, diff1 / diff1new, diff2new,
            diff2 / diff2new);
    cout << dump;
    cout.flush();
    diff0 = diff0new;
    diff1 = diff1new;
    diff2 = diff2new;
  }
} // All done

} // namespace derivtest
} // namespace polyvec

#endif /* CONTRIB_DERIV_TEST_IS_INCLUDED */
