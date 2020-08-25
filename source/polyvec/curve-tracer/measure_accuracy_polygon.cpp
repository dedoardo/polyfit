// polyvec
#include <polyvec/curve-tracer/measure_accuracy_polygon.hpp>
#include <polyvec/curve-tracer/find-fitting-points.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/curve-tracer/spline.hpp>
#include <polyvec/utils/directions.hpp>

// remove todo
#include <polyvec/geom.hpp> // nop
#include <polyvec/debug.hpp>


// libc++
#include <cstdlib> // min, max

using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

void measure_accuracy_signed_extended_polygon_fit_symmetric(
	const mat2x& B,       // raster boundary
	const vecXi& P,        // polygon
	const Vertex corner,      // corner
	AccuracyMeasurement& m,// result
	const bool circular
) {
	double len_prev = 0.;
	if (circular || corner > 0) {
		len_prev = (B.col(CircularAt(P, corner - 1)) - B.col(CircularAt(P, corner))).norm();
	}

	double len_next = 0.;
	if (circular || corner < P.size() - 1) {
		len_next = (B.col(CircularAt(P, corner + 1)) - B.col(CircularAt(P, corner))).norm();
	}

	int e_start = corner;
	double t_start = 0.;
	if (circular || corner > 0) {
		e_start = Circular(P, corner - 1);
		t_start = (.5 * min(len_prev, len_next)) / len_prev;
	}

	// i don't think this works for non-circular segments (todo)
	int e_end = corner;
	double t_end = 0.;
	if (circular || corner < P.size() - 1) {
		t_end = (.5 * min(len_prev, len_next)) / len_next;
	}

	// Fitting data
	const vec2 p0 = B.col(P(e_start));
	const vec2 p1 = B.col(P(e_end));
	const vec2 p2 = B.col(CircularAt(P, e_end + 1)); // todo circularity

	mat2xi V_fit;
	mat2x P_fit;

	// (1) find fitting points
	find_pixel_centers_for_subpath(B, P, e_start, t_start, e_start, 1., V_fit, &P_fit, circular);

	// (1) compute fitting normals
	mat2x N_fit(2, P_fit.cols());
	for (int i = 0; i < P_fit.cols(); ++i) {
		N_fit.col(i) = polyvec::util::normal_dir(B.col(V_fit(0, i)) - B.col(V_fit(1, i)));
	}

	// (1) construct edge segment
	GlobFitCurve_Line edge0;
	edge0.set_points(
		p0 + (p1 - p0) * t_start,
		p1
	);

	// (1) test accuracy
	AccuracyMeasurement m0 = measure_accuracy(edge0, P_fit, N_fit);

	// (2) find fitting points
	V_fit.resize(2, 0);
	P_fit.resize(2, 0);
	find_pixel_centers_for_subpath(B, P, e_end, 0., e_end, t_end, V_fit, &P_fit, circular);

	// (2) compute fitting normals 
	N_fit.resize(2, P_fit.cols());
	for (int i = 0; i < P_fit.cols(); ++i) {
		N_fit.col(i) = polyvec::util::normal_dir(B.col(V_fit(0, i)) - B.col(V_fit(1, i)));
	}

	// (2) construct edge segment
	GlobFitCurve_Line edge1;
	edge1.set_points(
		p1,
		p1 + (p2 - p1) * t_end
	);

	// (2) test accuracy
	AccuracyMeasurement m1 = measure_accuracy(edge1, P_fit, N_fit);
	
	new (&m) AccuracyMeasurement;
	m.combine(m0);
	m.combine(m1);
}

void measure_accuracy_signed_extended_polygon_fit_asymmetric(
	const mat2x& B,       // raster boundary
	const vecXi& P,        // polygon
	const Vertex corner,      // corner
	AccuracyMeasurement& m,// result
	const bool circular
) {
	double t_start = .5;
	int e_start = corner;
	if (circular || corner > 0) {
		e_start = Circular(P, corner - 1);
	}

	double t_end = .5;
	int e_end = corner;

	// Fitting data
	const vec2 p0 = B.col(P(e_start));
	const vec2 p1 = B.col(P(e_end));
	const vec2 p2 = B.col(CircularAt(P, e_end + 1)); // todo circularity

	mat2xi V_fit;
	mat2x P_fit;

	// (1) find fitting points
	find_pixel_centers_for_subpath(B, P, e_start, t_start, e_start, 1., V_fit, &P_fit, circular);

	// (1) compute fitting normals
	mat2x N_fit(2, P_fit.cols());
	for (int i = 0; i < P_fit.cols(); ++i) {
		N_fit.col(i) = polyvec::util::normal_dir(B.col(V_fit(1, i)) - B.col(V_fit(0, i)));
	}

	// (1) construct edge segment
	GlobFitCurve_Line edge0;
	edge0.set_points(
		p0 + (p1 - p0) * t_start,
		p1
	);

	// (1) test accuracy
	AccuracyMeasurement m0 = measure_accuracy(edge0, P_fit, N_fit);

	// (2) find fitting points
	V_fit.resize(2, 0);
	P_fit.resize(2, 0);
	find_pixel_centers_for_subpath(B, P, e_end, 0., e_end, t_end, V_fit, &P_fit, circular);

	// (2) compute fitting normals 
	N_fit.resize(2, P_fit.cols());
	for (int i = 0; i < P_fit.cols(); ++i) {
		N_fit.col(i) = polyvec::util::normal_dir(B.col(V_fit(1, i)) - B.col(V_fit(0, i)));
	}

	// (2) construct edge segment
	GlobFitCurve_Line edge1;
	edge1.set_points(
		p1,
		p1 + (p2 - p1) * t_end
	);

	// (2) test accuracy
	AccuracyMeasurement m1 = measure_accuracy(edge1, P_fit, N_fit);

	new (&m) AccuracyMeasurement;
	m.combine(m0);
	m.combine(m1);
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)