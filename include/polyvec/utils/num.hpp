#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

#include <algorithm>

#ifndef logic_xor // https://benpfaff.org/writings/clc/logical-xor.html
#define logic_xor(a, b) ((!(a)) != (!(b)))
#endif

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(Num)
	template <typename T>
	T lerp(const T& v0, const T& v1, const double t) {
		return v0 + (v1 - v0) * t;
	}

    // hlsl ftw
	template <typename T>
	T saturate(const T& v) {
		return std::max(std::min(v, (T)1.), (T)0.);
	}

	inline double sign(const double v) {
		return (0 < v) - (v < 0);
	}

	inline double sign_with_tol(const double v) {
		const double tol = 1e-10;
		if(v > tol ) return 1;
		if ( v < -tol) return -1;
		return 0;
	}

    double determinant(Eigen::Ref<const Eigen::MatrixXd> mat);
		
    // Returns true if the two intervals are overlapping and the value is written to overlap
    bool test_and_calculate_interval_overlap(
        const double i0, const double i1,
        const double j0, const double j1,
        double& overlap
    );

    // put this here prevent huge compilation times
	Eigen::MatrixXd solve_linear_system(
			const Eigen::Ref<const Eigen::MatrixXd> LHS, 
			const Eigen::Ref<const Eigen::MatrixXd> RHS) ;

	//Implements a smooth transition function that returns
	//  0  iff x <= zero_up_to
	//  1  iff x >= one_beyond
	//  a smooth transition otherwise
	double smooth_probability_incr(double x, double zero_up_to, double one_beyond);

	//Implements a smooth transition function that returns
	//  1  iff x <= one_up_to
	//  0  iff x >= zero_beyond
	//  a smooth transition otherwise
	double smooth_probability_decr(double x, double one_up_to, double zero_beyond);
NAMESPACE_END(Num)
NAMESPACE_END(polyvec)
