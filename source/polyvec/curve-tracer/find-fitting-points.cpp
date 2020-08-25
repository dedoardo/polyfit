#include <polyvec/curve-tracer/find-fitting-points.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/utils/num.hpp>

// todo: remove
#include <polyvec/debug.hpp>

using namespace polyvec;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(CurveTracer)

int find_pixel_centers_for_subpath(
	const mat2x& B,     // raster boundary
	const vecXi& V,     // vertices of the polygonal approximation
	const size_t e_beg, // index of the first edge (source endpoint)
	const double t_beg, // t along the first edge where the subpath begins,
	const size_t e_end, // index of the last edge (source endpoint)
	const double t_end, // t along the first edge where the subpath begin
	mat2xi& V_fit,   	// output vertices
	mat2x* P_fit,		// output midpoints
	bool circular
) {
	if (V_fit.size() == 0) {
		V_fit.resize(2, 0);
	} 

	if (P_fit && P_fit->size() == 0) {
		P_fit->resize(2, 0);
	}

	// the midpoints are a subset of those contains between these two boundary vertices
	const Vertex v_beg = V(e_beg);
	const Vertex v_end = V((e_end + 1) % V.size());

	Vertex v = GeomRaster::find_next_grid_aligned(B, v_beg, circular);

	if (v == -1) {
		return 0;
	}

	// current edge being tested
	Vertex e = e_beg;
	Vertex en = (e_beg + 1) % V.size();

	while (v != v_end) {
		// edge endpoints
		const vec2 e_src = B.col(V(e));
		const vec2 e_dst = B.col(V(en));

		// finding the first pair of (v, vn) which are consecutive and grid aligned
		Vertex vn = GeomRaster::find_next_grid_aligned(B, (v + 1) % B.cols(), circular);

		if (vn == -1) {
			break;
		}

		// consecutive grid aligned points should be on consecutive pixels
		if ((B.col(v) - B.col(vn)).norm() > 1. + PF_EPS_LARGE) {
			PF_LOGF("error! this distance should be shorter than 1 pixel B(%d) = %f %f B(%d) = %f %f",
			v, B.col(v)(0), B.col(v)(1), vn, B.col(vn)(0), B.col(vn)(1));
			PF_ABORT;
		}

		// candidate midpoint
		const vec2 p_mid = (B.col(v) + B.col(vn)) * .5;

		// parametric value of the closest point on the curve
		const double t_mid = LineUtils::project_t(p_mid, e_src, e_dst);
		PF_ASSERT(!isnan(t_mid));

		// if we are on the first or last segment we should check the start/end t
		bool valid = true;
		if(e == e_beg)
			valid &= t_mid >= t_beg - PF_EPS;
		if(e == e_end)
			valid &= t_mid <= t_end + PF_EPS;

		if (valid) {
			MatrixUtils::append(V_fit, vec2i(v, vn));

			if (P_fit) {
				MatrixUtils::append(*P_fit, p_mid);
			}
		} 

		// the first point which exceeds the bound on the last edge is the termination criteria
		if (t_mid > t_end && e == e_end) {
			break;
		}

		// moving to the next edge once the next grid point corresponds to one of the polygon's vertices
		if (vn == V(en))
		{
			e = en;
			en = (en + 1) % V.size();
		}

		v = vn;
	}

	return (int)V_fit.size();
}

NAMESPACE_END(CurveTracer)
NAMESPACE_END(polyfit)