// polyvec
#include <polyvec/regularity/check_opposite_points_winding_number.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/utils/num.hpp>

using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

bool check_opposite_points_winding_number (
	const mat2x& P,
	const Vertex i,
	const Vertex j
) {
	const vec2 pi = P.col(i);
	const vec2 pj = P.col(j);

	const int n_samples = (int)(std::ceil((pj - pi).norm()) + 1);
	const double t_min = .025;
	const double t_step = 1.f / n_samples;
	double t = t_step;

	// todo the winding order should be assumed as input
	bool is_ccw;
	WindingNumber::compute_orientation(P, is_ccw);
	const double winding_number_sign = is_ccw ? 1. : -1.;

	while (t < 1. + PF_EPS) {
		const vec2 point = Num::lerp(pi, pj, std::max(std::min(t, 1 - t_min), t_min));
		double winding_number;

		// this is true if the sample is withing epsilon from and edge. Since the samples
		// are offset by a constant amount it is true when the two edges have dot product +-1
		bool winding_number_trustable;

		WindingNumber::compute_winding(P, point, winding_number, winding_number_trustable);
		if (isnan(winding_number)) {
			winding_number_trustable = false;
		}

		if (winding_number_trustable && abs(winding_number) < PF_EPS) {
			return false;
		}

		t += t_step;
	}

	return true;
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)