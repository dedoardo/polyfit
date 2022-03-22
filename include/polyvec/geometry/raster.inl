NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(GeomRaster)

// this is wrong
PV_INLINE bool is_grid_aligned(const vec2& p, const double tol) {
	return fmod(p(0), .5) < tol &&
		fmod(p(1), .5) < tol;
}

PV_INLINE Vertex find_next_grid_aligned(const mat2x& P, const Vertex v, const bool circular) {
	Vertex d = 0;
	while (d < P.cols()) {
		if (!circular && v + d >= P.cols()) {
			return -1;
		}

		Vertex vn = Circular(P, v + d);

		if (is_grid_aligned(P.col(vn))) {
			return vn;
		}

		++d;
	}

	return -1;
}

PV_INLINE bool are_overlapping(const vec2& p0, const vec2& p1, const double tol) {
	return (p1 - p0).cwiseAbs().maxCoeff() < tol;
}

PV_INLINE bool is_edge_degenerate(const vec2& src, const vec2& dst) {
    return (src - dst).cwiseAbs().maxCoeff() < 1. + PF_EPS &&
           (src - dst).cwiseAbs().minCoeff() < 0. + PF_EPS;
}

#define PIXEL_DIAGONAL_LEN (1.41421356237 * .5)
#define PIXEL_DIAGONAL_LEN_SQ (PIXEL_DIAGONAL_LEN * PIXEL_DIAGONAL_LEN)
PV_INLINE vec2 distance_bounds_from_points_with_slack(
    const mat2x& B,	  // raster boundary
    const vec2 p_src, // segment source point
    const vec2 p_dst, // segment destination point
    const int  v_src, // segment source vertex
    const int  v_dst  // segment destination vertex
) {
    PF_ASSERT(v_src >= 0);
    PF_ASSERT(v_dst >= 0);

    // This happens when the raster region has 1-pixel wide diagonals
    if (v_src != v_dst && GeomRaster::are_overlapping(p_src, p_dst)) {
        return vec2(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    }

    // If the intersection between the normal and the line farther than this 
    // distance from the midpoints, we trigger the projection test.
    // This should really be 0.5
    // todo: adding an epsilon would avoid extra computations in case the distance
    // is exactly 0.5, but it's still correct.
    const double max_distance_intersect = Options::get()->midpoint_accuracy_trigger_projection;

    // Finding how voxel midpoints to check for accuracy
    const int n_midpoints_to_test = CircularDist(B, v_src, v_dst);

    // Keeping track of the maximum violation
    vec2 error_max(0., 0.);

    for (int i = 0; i < n_midpoints_to_test; ++i) {
        // Calculating the position of the midpoint to test. 
        const int v = Circular(B, v_src + i);
        int vn = Circular(B, v + 1);

        // The boundary is allowed to contain midpoints, in which case we safely skip them
        // and move to the next point
        const double v_vn_dist_sq = (B.col(v) - B.col(vn)).squaredNorm();
        // fprintf(stderr, "v %d vn %d v_vn_dist_sq %f\n", v, vn, v_vn_dist_sq);
        if (abs(v_vn_dist_sq) < 1. - PF_EPS) {
            vn = Circular(B, v + 2);

            // Not expecting more than one point which is not pixel aligned
            //fprintf(stderr, "v %d vn %d v_vnn_dist_sq %f\n", v, vn, (B.col(v) - B.col(vn)).squaredNorm());
            //fprintf(stderr, "v %f %f\n", B(0, v), B(1, v));
            //fprintf(stderr, "vn %f %f\n", B(0, vn), B(1, vn));
            PF_ASSERT((B.col(v) - B.col(vn)).squaredNorm() >= PIXEL_DIAGONAL_LEN_SQ - PF_EPS_MEDIUM);
            ++i;
        }

        const vec2 p = B.col(v);
        const vec2 pn = B.col(vn);
        const vec2 pmid = .5 * (p + pn);

        // Let's avoid a normalization, but we expect the points to be axis-aligned.
        // and the polyline to be oriented clockwise
        vec2 p_normal = pn - p;
        p_normal = vec2(p_normal(1), -p_normal(0)).normalized();

        // Testing the distance along the normal at which the edge is intersected
        double d_intersect, d_intersect_1;
        if (polyvec::LineUtils::intersect(pmid, pmid + p_normal, p_src, p_dst, d_intersect, d_intersect_1)) {
            if (d_intersect_1 <= -PF_EPS_MEDIUM || d_intersect_1 >= 1 + PF_EPS_MEDIUM) {
                error_max(0) = error_max(1) = INFINITY;
                break;
            }

            int violation_sign = d_intersect >= 0 ? +1 : -1;
            d_intersect = abs(d_intersect);

            if (d_intersect < max_distance_intersect) {
                if (violation_sign > 0) {
                    error_max(0) = ::std::max(error_max(0), d_intersect);
                }
                else {
                    error_max(1) = ::std::max(error_max(1), d_intersect);
                }

                continue;
            }
        }

        // Otherwise we test the distance from the point offset along the normal to its projection on 
        // the line.
        const vec2 pmid_test_out = pmid + max_distance_intersect * p_normal;
        const vec2 pmid_test_in = pmid - max_distance_intersect * p_normal;
        const vec2 pmid_test_out_proj = polyvec::LineUtils::project_point(pmid_test_out, p_src, p_dst);
        const vec2 pmid_test_in_proj = polyvec::LineUtils::project_point(pmid_test_in, p_src, p_dst);

        // todo: should this be an assert?
        // todo: is this really happening
        //if (!geom::point_on_line(edge_segment, pmid_test_out_proj) ||
        //	!geom::point_on_line(edge_segment, pmid_test_in_proj)) {
        //	return vec2(FLT_MAX, FLT_MAX);
        //}

        // Taking the L1 of the distance
        const double d_project_out = (pmid_test_out_proj - pmid_test_out).cwiseAbs().maxCoeff();
        const double d_project_in = (pmid_test_in_proj - pmid_test_in).cwiseAbs().maxCoeff();

        // Testing the sign of the greater of the two errors and updating the respecting value tracker
        int violation_sign = d_project_out < d_project_in ? +1 : -1;
        if (violation_sign > 0) {
            error_max(0) = ::std::max(error_max(0), abs(.5 + d_project_out));
        }
        else {
            error_max(1) = ::std::max(error_max(1), abs(.5 + d_project_in));
        }
    }

    return vec2(error_max.minCoeff(), error_max.maxCoeff());
}


PV_INLINE double get_resolution_scaling_factor_32x32(
    const mat2x& B
) {
    vec2 B_min = vec2(INFINITY, INFINITY);
    vec2 B_max = vec2(-INFINITY, -INFINITY);

    for (Eigen::Index i = 0; i < B.cols(); ++i) {
        B_min = B_min.cwiseMin(B.col(i));
        B_max = B_max.cwiseMax(B.col(i));
    }
   
    const double R = (B_max - B_min).maxCoeff();
    return R / 32;
}

PV_INLINE double raster_corner_size(const mat2x& raster_boundary, Vertex v) {	
	auto corner = raster_boundary.col(v);
	
	auto get_equal_dim = [&](int offset) -> int
	{
		polyfit::vec2 d = CircularAt(raster_boundary, v + offset) - corner;
		for (int i = 0; i < 2; ++i)
			if (std::abs(d(i)) < PF_EPS)
				return i;
		return -1;
	};

	auto pos_dir = get_equal_dim(1);
	auto neg_dir = get_equal_dim(-1);
	if (pos_dir == neg_dir)
		return 0; //this is not a corner
	int length = 1;
	int offset_pos = 1;
	int offset_neg = 1;
	while (true)
	{
		while ((CircularAt(raster_boundary, v + offset_pos) - corner).squaredNorm() < length * length - PF_EPS)
			++offset_pos;
		while ((CircularAt(raster_boundary, v - offset_neg) - corner).squaredNorm() < length * length - PF_EPS)
			++offset_neg;
		if (get_equal_dim(offset_pos + 1) != pos_dir)
			return length;
		if (get_equal_dim(-offset_neg - 1) != neg_dir)
			return length;
		++length;
	}
}

PV_INLINE bool are_points_axis_aligned(const vec2& p0, const vec2& p1, const double tol) {
	return (p0 - p1).cwiseAbs().minCoeff() < tol;
}

NAMESPACE_END(GeomRaster)
NAMESPACE_END(polyfit)