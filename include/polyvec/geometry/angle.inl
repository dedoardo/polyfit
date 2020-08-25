PV_INLINE bool lie_in_the_same_quadrant(
	const vec2& p0,
	const vec2& p1
) {
	// I have the impression this could be done with half the number of tests
	const vec2i x0_sign(p0(0) > -PF_EPS, p0(0) < PF_EPS);
	const vec2i x1_sign(p1(0) > -PF_EPS, p1(0) < PF_EPS);
	const vec2i y0_sign(p0(1) > -PF_EPS, p0(1) < PF_EPS);
	const vec2i y1_sign(p1(1) > -PF_EPS, p1(1) < PF_EPS);

	return x0_sign == x1_sign && y0_sign == y1_sign;
}

PV_INLINE double deviation_from_horizontal(const vec2& d) {
	return deviation_from_horizontal(atan2(d.y(), d.x()));
}

PV_INLINE double deviation_from_vertical(const vec2& d) {
	return deviation_from_vertical(atan2(d.y(), d.x()));
}

PV_INLINE double deviation_from_horizontal(double angle) {
	const double angle1 = fmod(angle + 8. * M_PI, M_PI);
	return std::min(std::abs(angle1), abs(M_PI - angle1));
}

PV_INLINE double deviation_from_vertical(double angle) {
	const double angle1 = fmod(angle + 8. * M_PI, M_PI);
	return  std::fabs(angle1 - M_PI_2);
}

PV_INLINE double spanned_shortest(const vec2 & p0, const vec2 & p1, const vec2 & p2) {
    const vec2 d0 = (p0 - p1).normalized();
    const vec2 d1 = (p2 - p1).normalized();
    return ::std::max(0., ::std::min(acos(d0.dot(d1)), M_PI));
}

PV_INLINE double spanned_shortest_between(const vec2 & d0, const vec2 & d1) {
	return acos(d0.normalized().dot(d1.normalized()));
}

// https://stackoverflow.com/questions/14066933/direct-way-of-computing-clockwise-angle-between-2-vectors
PV_INLINE double spanned_clockwise_between(const vec2& d0, const vec2& d1) {
	const double dot = d0.dot(d1);
	const double det = d0(0) * d1(1) - d0(1) * d1(0);
	const double x = atan2(det, dot);
	return x > 0 ? x : (M_PI * 2 + x);
}

PV_INLINE bool have_opposite_convexity(const vec2 & p0, const vec2 & p1, const vec2 & p2, const vec2 & p3) {
	vec3 v01 = vec3::Zero(), v23 = vec3::Zero(), v12 = vec3::Zero();
	v01.segment(0, 2) = p1 - p0;
	v23.segment(0, 2) = p3 - p2;
	v12.segment(0, 2) = p2 - p1;

	const int sign_xprev = (int)polyvec::Num::sign(v12.cross(v01).z());
	const int sign_xnext = (int)polyvec::Num::sign(v12.cross(v23).z());
	return sign_xprev * sign_xnext == 0 ? false : sign_xprev == sign_xnext;
}

PV_INLINE bool have_opposite_convexity_with_tol(const vec2 & p0, const vec2 & p1, const vec2 & p2, const vec2 & p3) {
	vec3 v01 = vec3::Zero(), v23 = vec3::Zero(), v12 = vec3::Zero();
	v01.segment(0, 2) = p1 - p0;
	v23.segment(0, 2) = p3 - p2;
	v12.segment(0, 2) = p2 - p1;

	const int sign_xprev = (int)polyvec::Num::sign_with_tol(v12.cross(v01).z());
	const int sign_xnext = (int)polyvec::Num::sign_with_tol(v12.cross(v23).z());
	return sign_xprev * sign_xnext == 0 ? false : sign_xprev == sign_xnext;
}


PV_INLINE int count_inflections(const mat2x & points, bool circular) {
    if (points.cols() < 4) {
        return 0;
    }

	int count = 0;
	const Eigen::Index offset = circular ? 0 : 1;
	for (Eigen::Index i = offset; i < points.cols() - offset; ++i) {
		if (AngleUtils::have_opposite_convexity(CircularAt(points, i - 1), CircularAt(points, i), CircularAt(points, i + 1), CircularAt(points, i + 2))) {
			++count;
		}
	}

	return count;
}

PV_INLINE bool is_flat(const vec2 pp, const vec2 p, const vec2 pn) {
	return (1. - abs((pn - p).normalized().dot((p - pp).normalized()))) < 1e-4;
}

PV_INLINE int convexity(const vec2 pp, const vec2 p, const vec2 pn, int orientation) {
	const vec2 tn = (pn - p).normalized();
	const vec2 tp = (p - pp).normalized();
	const vec2 n = orientation * vec2(-tp.y(), tp.x());
	double dcvx = -tn.dot(n);

	if (abs(dcvx) < 1e-5) {
		return 0;
	}
	else if (dcvx > 0) {
		return +1;
	}
	else {
		return -1;
	}
}

PV_INLINE int number_of_neighbors_with_different_convexity(const Eigen::Matrix2Xd& polygon, int corner)
{
	// Edo: Shouldn't the neighbors be -1 0 +1 ?
	int conv[3];
	for (int i = 0; i < 3; ++i)
		conv[i] = convexity(CircularAt(polygon, corner + i - 2), CircularAt(polygon, corner + i - 1), CircularAt(polygon, corner + i), 1);

	int result = 0;
	if (conv[0] != conv[1])
		++result;
	if (conv[2] != conv[1])
		++result;

	return result;
}

PV_INLINE void neighborhood_convexity(const Eigen::Matrix2Xd& P, const int corner, const int neighborhood_size, int* convexity_nhb) {
	for (int nhb = -neighborhood_size; nhb <= neighborhood_size; ++nhb) {
		const int nhb_idx = nhb + neighborhood_size;

		convexity_nhb[nhb_idx] = convexity(
			CircularAt(P, corner + nhb - 1), CircularAt(P, corner + nhb), CircularAt(P, corner + nhb + 1), 1
		);
	}
}

PV_INLINE int count_visually_degenerate_edges(const mat2x& B, const vecXi& P, const bool circular) {
	int count = 0;

	for (int i = 0; i < P.size(); ++i) {
		if (!circular && (i >= P.cols() - 1)) {
			continue;
		}

		const Eigen::Index v0 = P(i);
		const Eigen::Index v1 = CircularAt(P, i + 1);
		const vec2 p0 = B.col(v0);
		const vec2 p1 = B.col(v1);
		const bool is_short = (p1 - p0).squaredNorm() < 1. + PF_EPS * 2;
		const bool is_aa = PathUtils::are_axis_aligned(p0, p1);

		if (is_aa && is_short) {
			++count;
		}
	}

	return count;
}

PV_INLINE int count_visually_inconsistent_edges(const mat2x & P, bool circular) {
	Vertex off = circular ? 0 : 1; (void)off;
	int count = 0;

	for (int i = 0; i < P.cols(); ++i) {
		if (!circular && (i == 0 || i >= P.cols() - 1)) {
			continue;
		}

		const vec2 p0p = CircularAt(P, i - 1);
		const vec2 p0 = P.col(i);
		const vec2 p1 = CircularAt(P, i + 1);
		const vec2 p1n = CircularAt(P, i + 2);

		const bool is_short_aa = PathUtils::are_axis_aligned(p0, p1) && (p1 - p0).norm() < 1. + 1e-4;
		const bool is_inflection = AngleUtils::have_opposite_convexity(p0p, p0, p1, p1n);

		if (is_short_aa && is_inflection) {
			++count;
		}
	}

	return count;
}
