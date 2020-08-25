// polyvec
#include <polyvec/regularity/check_same_line_continuation.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/raster.hpp>

#include <algorithm>

#define MAX_DISTANCE_TO_LINE (1. + PF_EPS)

using namespace std;
using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(Regularity)

bool check_same_line_continuation(
	const mat2x& B,
	const vecXi& P,
	const Vertex i,
	const Vertex j,
	double& weight
) {
	weight = 0.;

	const vec2 l0 = B.col(CircularAt(P, i - 1));
	const vec2 l1 = B.col(CircularAt(P, j + 1));
	const vec2 pi = B.col(P(i));
	const vec2 pj = B.col(P(j));

	if (!AngleUtils::have_opposite_convexity(l0, pi, pj, l1)) {
		return false;
	}

	const double di = LineUtils::distance_from_point(l0, l1, pi);
	const double dj = LineUtils::distance_from_point(l0, l1, pj);

	weight = fmax(di, dj) / MAX_DISTANCE_TO_LINE;
	return di < MAX_DISTANCE_TO_LINE &&
		dj < MAX_DISTANCE_TO_LINE;
}

bool check_line_continuity(
	const vec2 pi0,
	const vec2 pi1,
	const vec2 pj0,
	const vec2 pj1,
	double& weight
) {
	const vec2 d0 = (pi1 - pi0).normalized();
	const vec2 d1 = (pj0 - pi1).normalized();
	const vec2 d2 = (pj1 - pj0).normalized();

	// are all the vectors aligned?
	if (abs(d0.dot(d1) - 1.) > PF_EPS ||
		abs(d1.dot(d2) - 1.) > PF_EPS) {
		return false;
	}

	// are they axis aligned?
	// todo: 1 check should be enough
	const double devh0 = AngleUtils::deviation_from_horizontal(d0);
	const double devh1 = AngleUtils::deviation_from_horizontal(d1);
	const double devh2 = AngleUtils::deviation_from_horizontal(d2);

	if (max(devh0, max(devh1, devh2)) < 1e-1) {
		return true;
	}

	const double devv0 = AngleUtils::deviation_from_vertical(d0);
	const double devv1 = AngleUtils::deviation_from_vertical(d1);
	const double devv2 = AngleUtils::deviation_from_vertical(d2);

	if (max(devv0, max(devv1, devv2)) < 1e-1) {
		return true;
	}

	return false;
}

NAMESPACE_END(Regularity)
NAMESPACE_END(polyfit)