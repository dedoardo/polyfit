#include <sstream>

#include <polyvec/api.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/curve_bezier.hpp>
#include <polyvec/utils/deriv_test.hpp>
#include <polyvec/io/vtk_curve_writer.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/geometry/smooth_curve.hpp>

using namespace polyvec;

namespace {

Eigen::MatrixXd
gaussian(const Eigen::MatrixXd & in)
{
  Eigen::MatrixXd out(in.rows(), in.cols());
  const double sqrtpiinv = 1. / std::sqrt(2 * 3.14);

  for(unsigned i = 0; i < in.rows(); ++i)
    {
      for(unsigned j = 0; j < in.cols(); ++j)
        {
          out(i, j) = sqrtpiinv * std::exp(-in(i, j) * in(i, j) / 2.);
        }
    }
  return out;
}


Eigen::MatrixXd
randg(unsigned int i, unsigned int j)
{
  return gaussian(Eigen::MatrixXd::Random(i, j));
}

void
test_tangent() {
  
    double h0 = 0.2;
    int n_halving = 5;
    Eigen::Vector2d drdt0  = randg(2,1) ;
    Eigen::Vector2d delta  = randg(2,1) ;
    Eigen::Matrix2d dtang_ddrdt;
    {
      Eigen::Vector2d dummy_tangent;
      SmoothCurveUtil::tangent_derivatives(drdt0, dummy_tangent, dtang_ddrdt);
    }

    derivtest::run (
        drdt0,
        delta,
    [&] ( const Eigen::VectorXd & drdt_perturbed_in ) -> Eigen::Vector2d {
      Eigen::Vector2d ans;
      SmoothCurveUtil::tangent_eval(drdt_perturbed_in, ans);
      return ans; 
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Vector2d {
        return dtang_ddrdt * hdelta;
    },
    std::cout,
    n_halving,
    h0 );
} // All done with test_tangent

void
test_curvature() {
  
    double h0 = 0.2;
    int n_halving = 5;
    Eigen::Matrix<double, 4, 1> drdt0_and_drdtdt0  = randg(4,1) ;
    Eigen::Matrix<double, 4, 1> delta  = randg(4,1) ;
    Eigen::RowVector2d dcurvature_ddrdt, dcurvature_ddrdtdt;
    {
      double dummy_curvature;
      SmoothCurveUtil::curvature_derivatives(drdt0_and_drdtdt0.head<2>(), drdt0_and_drdtdt0.tail<2>(), dummy_curvature, dcurvature_ddrdt, dcurvature_ddrdtdt);
    }

    derivtest::run (
        drdt0_and_drdtdt0,
        delta,
    [&] ( const Eigen::VectorXd & drdt0_and_drdtdt0_pert ) -> Eigen::Matrix<double,1,1> {
      Eigen::Matrix<double,1,1> ans;
      SmoothCurveUtil::curvature_eval(drdt0_and_drdtdt0_pert.head<2>(), drdt0_and_drdtdt0_pert.tail<2>(), ans(0));
      return ans; 
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Matrix<double,1,1> {
        return dcurvature_ddrdt* hdelta.head<2>() +  dcurvature_ddrdtdt* hdelta.tail<2>();
    },
    std::cout,
    n_halving,
    h0 );
} // All done with test_tangent

void
test_smooth_abs(const double x_in) {
  
    double h0 = 0.1;
    int n_halving = 7;
    Eigen::Matrix<double, 1, 1> x0;
    x0 << x_in ;
    Eigen::Matrix<double, 1, 1> delta  = randg(1,1) ;
    Eigen::Matrix<double, 1, 1> dsmoothabsdx;
    {
      double dummy_abs;
      SmoothCurveUtil::smoothabs_derivatives(x0(0), dummy_abs, dsmoothabsdx(0));
    }

    derivtest::run (
        x0,
        delta,
    [&] ( const Eigen::VectorXd & x_pert ) -> Eigen::Matrix<double,1,1> {
      Eigen::Matrix<double,1,1> ans;
      SmoothCurveUtil::smoothabs_eval(x_pert(0), ans(0));
      return ans; 
    },
    [&] ( const Eigen::VectorXd & hdelta ) -> Eigen::Matrix<double,1,1> {
        return dsmoothabsdx* hdelta;
    },
    std::cout,
    n_halving,
    h0 );
} // All done with test_tangent

} // end of anonymus

namespace polyvectest {
namespace geom {
int
test_curvature_and_tangent ( int /*argc*/, char** /*argv*/ ) {

printf("Testing tangent derivative \n");
test_tangent();
printf("Testing curvature derivative \n");
test_curvature();
printf("Testing smooth abs derivative at 0\n");
test_smooth_abs(0);
printf("Testing smooth abs derivative at 0.015\n");
test_smooth_abs(0.015);
printf("Testing smooth abs derivative at -0.015\n");
test_smooth_abs(-0.015);

return 0;

}
} // geom
} // polyvectest
