// Header
#include <polyvec/geometry/path.hpp>

// polyvec
#include <polyvec/api.hpp>
#include <polyvec/core/constants.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/debug.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/directions.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/core/log.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/angle.hpp>
#include <polyvec/geometry/line.hpp>
#include <polyvec/geometry/winding_number.hpp>
#include <polyvec/polygon-tracer/error-metrics.hpp>

// c++ stl
#include <algorithm>
#include <limits>

// 
#include <Eigen/Geometry>

using namespace Eigen;
using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(PathUtils)

void order_clockwise(mat2x& P) {
	bool ccw;
	WindingNumber::compute_orientation(P, ccw);
	PF_LOGF("orientation %s", ccw ? "counter-clockwise" : "clockwise");

	if (ccw) {
		Reverse(P);
	}
}

bool DegenerateClassifierAxisAligned1Px::should_skip_edge (const vec2& p0, const vec2& p1) {
    const vec2 dabs = (p1 - p0).cwiseAbs();
    return dabs.maxCoeff() < 1. + PF_EPS && dabs.minCoeff() < PF_EPS;
    //return dabs.minCoeff() < PF_EPS;
}

bool DegenerateClassifier2Px::should_skip_edge(const vec2& p0, const vec2& p1) {
	const vec2 dabs = (p1 - p0).cwiseAbs();
	return dabs.maxCoeff() < 2. + PF_EPS;
}

template <typename TDegenerateClassifier>
int find_next_non_degenerate_corner(
    const Eigen::Matrix2Xd& P,
    const int vertex,          // Current vertex
    const int vertex_next,     // Next expected vertex
    const bool circular
) {
    if (!circular && vertex_next == P.cols() - 1) {
        return vertex;
    }

    const int direction = CircularDist(P, vertex, vertex_next) < CircularDist(P, vertex_next, vertex) ? +1 : -1;
    
    const vec2 p0 = P.col(vertex);
    const vec2 p1 = P.col(vertex_next);
    const vec2 p2 = CircularAt(P, vertex_next + direction);

    // If this actually happens, the conditionals below should be made into a for loop statement maybe?
    const bool is_next_degenerate = TDegenerateClassifier::should_skip_edge(p0, p1);
    
    const bool is_next_flat = AngleUtils::spanned_shortest(p0, p1, p2) > PV_FLAT_ANGLE_MIN;

    if (is_next_degenerate || is_next_flat) {
        return Circular(P, vertex_next + direction);
    } else {
        return vertex_next;
    }
}

//explicit instantiations
template int find_next_non_degenerate_corner<DegenerateClassifierAxisAligned1Px>(const Eigen::Matrix2Xd& P, const int vertex, const int vertex_next, const bool circular);
template int find_next_non_degenerate_corner<DegenerateClassifier2Px>(const Eigen::Matrix2Xd& P, const int vertex, const int vertex_next, const bool circular);

void compute_hull_with_overlaps(const Eigen::Matrix2Xd& points, Eigen::Matrix2Xd& boundary, const int direction) {
	std::vector<int> convexities;
	PathUtils::compute_convexities(points, convexities);

	boundary.resize(Eigen::NoChange, 0);

	for (Index i = 0; i < (Index)convexities.size(); ++i) {
		if (convexities[i] == 0) {
			continue;
		}

		const Vector2d p = points.col(i);
		const Vector2d p_next = CircularAt(points, i + 1);
		const Vector2d p_prev = CircularAt(points, i - 1);

		const Vector2d cp = (p - p_next) * .5 + (p - p_prev) * .5;
		MatrixUtils::append(boundary, p + cp * convexities[i] * direction);
	}
}

double bounding_box_diagonal(const Eigen::Matrix2Xd& points) {
    Vector2d min(INFINITY, INFINITY), max(-INFINITY, -INFINITY);
    for (Index i = 0; i < points.cols(); ++i) {
        min = min.cwiseMin(points.col(i));
        max = max.cwiseMax(points.col(i));
    }
    return (max - min).norm();
}

void compute_convexities(const Eigen::Matrix2Xd& points, std::vector<int>& convexities) {
	convexities.resize(points.cols(), 0);

	bool is_ccw;
	WindingNumber::compute_orientation(points, is_ccw);
	const int orientation = is_ccw ? -1 : +1;

	for (Index i = 0; i < points.cols(); ++i) {
		int j = (i + 1) % (int)points.cols();
		int k = (i - 1 + (int)points.cols()) % (int)points.cols();

		Vector2d tangent_i = (points.col(j) - points.col(i)).normalized();
		Vector2d tangent_k = (points.col(i) - points.col(k)).normalized();
		Vector2d normal_k = orientation * Vector2d(-tangent_k[1], tangent_k[0]);

		double dconvexity = -tangent_i.dot(normal_k);
		int iconvexity;

		if (abs(dconvexity) < 1e-5) {
			iconvexity = 0;
		}
		else if (dconvexity > 0) {
			iconvexity = 1;
		}
		else {
			iconvexity = -1;
		}

		convexities[i] = iconvexity;
	}
}
		
// Starting from a corner, it counts the number of flat vertices before reaching another corner
int count_flat_vertices_before_transition(const std::vector<int>& convexities, const Eigen::Index from, const Eigen::Index direction) {
	assert_break(convexities[from] != 0);
	assert_break(direction == +1 || direction == -1);

	int distance = 0;
	while (++distance < convexities.size()) {
		const int next_index = (from + (distance * direction) + (int)convexities.size() * 2) % (int)convexities.size();
				
		if (convexities[next_index] != 0) {
			break;
		}
	}

	return distance;
}

int count_transitions_between_vertices(const std::vector<int>& convexities, const Eigen::Index src, const Eigen::Index dst) {
	int vertex = Circular(convexities, src + 1);
	int transitions = 0;

	while (vertex != dst) {
		transitions += convexities[vertex] != 0 ? 1 : 0;
		vertex = Circular(convexities, vertex + 1);
	}

	return transitions;
}

bool is_edge_inflected(const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, const Eigen::Vector2d& v2, const Eigen::Vector2d& v3) {
	const Vector2d vfrom2d = v1 - v0;
	const Vector2d vto2d = v3 - v2;
	const Vector2d sep2d = v2 - v1;

	const Vector3d vfrom(vfrom2d.x(), vfrom2d.y(), 0.);
	const Vector3d vto(vto2d.x(), vto2d.y(), 0.);
	const Vector3d sep(sep2d.x(), sep2d.y(), 0.);

	const Vector3d cross_from = sep.cross(vfrom);
	const Vector3d cross_to = sep.cross(vto);

	int sfrom = misc::sign(cross_from.z());
	int sto = misc::sign(cross_to.z());
	return sfrom * sto == 0 ? false : sfrom == sto;
}
		
Vertex end_of_convex_region(const std::vector<int>& C, const Vertex v0, const int d, bool circular) {
	PF_ASSERT(v0 >= 0 && v0 < (Vertex)C.size());

	if (!circular && v0 == C.size() - 1) {
		return v0;
	}
			
	int c_last = C[v0];
	Vertex v = Circular(C, v0 + d);

	while (true) {
		if (!circular && v == C.size() - 1) {
			return v;
		}

		if (C[v] != 0) {
			if (c_last == C[v]) {
				break;
			}
			c_last = C[v];
		}

		v = Circular(C, v + d);
	}



	return v;
}

bool contains_closed(const int size, Vertex v0, Vertex v1, Vertex v, bool circular) {
	Vertex d01 = circular ? CircularDist(size, v0, v1) : (v1 - v0);
	Vertex d1 = circular ? CircularDist(size, v, v1) : (v1 - v);
	return d1 <= d01;
}

bool contains_closed(const mat2x& P, Vertex v0, Vertex v1, Vertex v, bool circular) {
	return contains_closed(P.cols(), v0, v1, v, circular);	
}

bool contains_open(const mat2x& P, Vertex v0, Vertex v1, Vertex v, bool circular) {
	Vertex d01 = circular ? CircularDist(P, v0, v1) : (v1 - v0);
	Vertex d1 = circular ? CircularDist(P, v, v1) : (v1 - v);
	return d1 > 0 && d1 < d01;
}

bool has_equal_distance(const mat2x& P, Vertex v0, Vertex v, Vertex v1, bool circular) {
	if (!circular && (v < v0 || v > v1)) {
		return false;
	}

	return CircularDist(P, v0, v) == CircularDist(P, v, v1);
}

Eigen::Index next_transition_vertex(const std::vector<int>& convexities, const Eigen::Index from, const int direction, const bool circular) {
	Index vertex = Circular(convexities, from + direction);

	while (convexities[vertex] == 0) {
		vertex = Circular(convexities, vertex + direction);
	}

	return vertex;
}

Vertex distance_to_next_corner (const std::vector<int>& convexities, const Eigen::Index from, const int direction) {
    return direction > 0 ?
        CircularDist(convexities, from, next_transition_vertex(convexities, from, direction)) :
        CircularDist(convexities, next_transition_vertex(convexities, from, direction), from);
}

bool are_axis_aligned(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1) {
	return abs((p0.x() - p1.x())) < constants::Eps ||
		abs((p0.y() - p1.y())) < constants::Eps;
}

double shortest_angle_spanned(const vec2& p0, const vec2& p1, const vec2& p2) {
	return acos((p0 - p1).normalized().dot((p2 - p1).normalized()));
}

// todo: is a cross product really needed? 
bool are_consecutive_angles_opposite(const vec2& p0, const vec2& p1, const vec2& p2, const vec2& p3) {
	vec3 v01 = vec3::Zero(), v23 = vec3::Zero(), v12 = vec3::Zero();
	v01.segment(0, 2) = p1 - p0;
	v23.segment(0, 2) = p3 - p2;
	v12.segment(0, 2) = p2 - p1;

	const int sign_xprev = misc::sign(v12.cross(v01).z());
	const int sign_xnext = misc::sign(v12.cross(v23).z());
	return sign_xprev * sign_xnext == 0 ? false : sign_xprev == sign_xnext;
}

/* OLD VERSION USING WYKOBY primitives
vec2 distance_bounds_from_points(const mat2x& points, const mat2& p, vec2i vp) {
	Vertex voxel_src = vp(0);
	Vertex voxel_dst = vp(1);
	vec2 from = p.col(0);
	vec2 to = p.col(1);

	// This happens when the raster region has 1-pixel wide diagonals
	if (voxel_src != voxel_dst && GeomRaster::are_overlapping(from, to)) {
		return vec2(FLT_MAX, FLT_MAX);
	}

	const double WEAK_THR = .5;

	// This is perfectly valid. The reason is that generated nodes store the respective voxel
	// with the `-` sign in front.
	voxel_src = abs(voxel_src);
	voxel_dst = abs(voxel_dst);

	//
	// Finding how voxel midpoints to check for accuracy
	geom::segment edge_segment(from, to);
	const double edge_segment_len = edge_segment.length() + constants::Eps;
	::polyvec::index n_voxels = points.cols();
	::polyvec::index num_nodes = -1;

	if (voxel_dst < voxel_src) {
		num_nodes = voxel_dst + n_voxels - voxel_src;
	} else {
		num_nodes = voxel_dst - voxel_src;
	}

	//
	// Keeping track of the maximum violation
	double max_in_accuracy_violation = -1.;
	double max_out_accuracy_violation = -1.;

	for (int i = 0; i < num_nodes; ++i) {
		::polyvec::index v = (voxel_src + i) % n_voxels;
		::polyvec::index next_v = (v + 1) % n_voxels;

		// only considering pixel points, not raster points,
		// todo: this shouldn't really be here
		if ((points.col(v) - points.col(next_v)).norm() < 1. - PF_EPS) {
			next_v = (v + 2) % n_voxels;
			PF_ASSERT((points.col(v) - points.col(next_v)).norm() >= 1. - PF_EPS);
			++i;
		}

		Vector2d cur = points.col(v);
		Vector2d next = points.col(next_v);
		Vector2d midpoint = 0.5 * (cur + next);
		Vector2d normal = next - cur;

		normal = polyvec::util::normal_dir(normal).normalized();		

		// WEAK ACCURACY (<0.5) - FIRST TEST
		// We first check along the midpoint normal
		double d = std::numeric_limits<double>::quiet_NaN();

		const vec2 n_abs = normal.cwiseAbs();

		if (geom::ray_intersect(midpoint, normal, edge_segment, d)) {
			//PF_ASSERT(abs(d - d2) < 1e-4 ||abs(d + d2) < 1e-4);
					
			if (std::abs(d) >= WEAK_THR) {
				d = std::numeric_limits<double>::quiet_NaN();				
			}
		}

		if (std::isnan(d))
		{
			// STRONG ACCURACY (>=0.5 to 0.5+slack) - SECOND TEST
			// Accuracy is checked as distance from the offset midpoint and its projection on the segment
			real2 off_midpoint_out = midpoint + 0.5 * normal;
			real2 off_midpoint_in = midpoint - 0.5 * normal;

			double dist_out = -1;
			double dist_in = -1;

			real2 out_proj = geom::project(off_midpoint_out, edge_segment);
			real2 in_proj = geom::project(off_midpoint_in, edge_segment);

			// Is the point on the line ?
			if (!geom::point_on_line(edge_segment, out_proj) && !geom::point_on_line(edge_segment, in_proj)) {
				//printf("edge    (%f %f) (%f %f)\n", edge_segment.src(0), edge_segment.src(1), edge_segment.dst(0), edge_segment.dst(1));
				//printf("out proj %f %f\n", out_proj(0), out_proj(1));
				//printf("in  proj %f %f\n", in_proj(0), in_proj(1));
				return vec2(FLT_MAX, FLT_MAX);
			}

			real2 d_out = (out_proj - off_midpoint_out).cwiseAbs();
			real2 d_in = (in_proj - off_midpoint_in).cwiseAbs();

			dist_out = ::std::max(d_out.x(), d_out.y());
			dist_in = ::std::max(d_in.x(), d_in.y());

			assert_break(dist_out >= 0.);
			assert_break(dist_in >= 0.);

			if (dist_out < dist_in) {
				d = abs(.5 + dist_out);
			}
			else {
				d = -abs(.5 + dist_in);
			}
		}

		if (d > 0) {
			max_out_accuracy_violation = max(max_out_accuracy_violation, abs(d));
		}
		else {
			max_in_accuracy_violation = max(max_in_accuracy_violation, abs(d));
		}
	}

	double accuracy_in = 0.;
	double accuracy_out = 0.;

	if (max_out_accuracy_violation > 0) {
		accuracy_out = max_out_accuracy_violation;
	}

	if (max_in_accuracy_violation > 0) {
		accuracy_in = max_in_accuracy_violation;
	}

	return vec2(std::min(accuracy_in, accuracy_out), std::max(accuracy_in, accuracy_out));
}

bool edge_is_within_relaxed_distance_bounds(const mat2x& P, const int vsrc, const int vdst) {
	mat2 e;
	e.col(0) = P.col(vsrc);
	e.col(1) = P.col(vdst);
	const vec2 d_error = PathUtils::distance_bounds_from_points(P, e, {vsrc, vdst});
	return ErrorMetrics::accuracy_within_bounds_relaxed(d_error.minCoeff(), d_error.maxCoeff());


}
*/

// returns true if exp is the first vertex encountered traversing the polygon P from i in direction d
bool points_in_same_partition(const mat2x& P, Vertex i, Vertex exp, Vertex fail, int d) {
	Vertex v = i + d;
	while (v != fail && v != i) {
		if (v == exp) {
			return true;
		}

		v = Circular(P, v + d);
	}

	return false;
}

Vertex opposite(const Eigen::Index& size, const Eigen::Index& v, const Eigen::Index& vsym) {
	if (CircularDist(size, vsym, v) < CircularDist(size, v, vsym )) {
		return Circular(size, vsym - 1);
	}
	else {
		return Circular(size, vsym + 1);
	}
}

void order(const mat2x& P, std::vector<Vertex>& V, Vertex vstart, bool circular) {
	// todo: test 
	PF_ASSERT(circular);

	std::vector<Vertex> Vord;

	Vertex prev = vstart;
	while (Vord.size() < V.size()) {
		Vertex dmin = numeric_limits<Vertex>::max();
		Vord.emplace_back(-1);

		// finding the next closest vertex in the path
		for (int i = 0; i < V.size(); ++i) {
			const Vertex d = CircularDist(P, vstart, V[i]);
					
			const bool exists = find(Vord.begin(), Vord.end(), V[i]) != Vord.end();

			if (!exists && d < dmin) {
				dmin = d;
				Vord.back() = V[i];
			}
		}

		PF_ASSERT(Vord.back() != -1);
	}

	V = std::move(Vord);
}

double angle_at_corner(
	const mat2x& B,
	const vecXi& P,
	const int v,
	const bool circular
) {
	if (!circular && (v == 0 || v == P.size() - 1)) {
		return 0.;
	}

	return AngleUtils::spanned_shortest(
		B.col(CircularAt(P, v - 1)),
		B.col(CircularAt(P, v)),
		B.col(CircularAt(P, v + 1))
	);
}

int count_degenerate_edges_exclusive(
	const mat2x& P, // raster boundary
	const mat2xi& E, // exceptions
	const bool circular
) {
	Vertex off = circular ? 0 : 1; (void)off;
	int count = 0;

	for (int i = 0; i < P.cols(); ++i) {
		if (!circular && (i == 0 || i >= P.cols() - 1)) {
			continue;
		}

		// testing if any of the two edges is contained in the intervals to exclude
		bool exclude = false;
		for (int j = 0; j < E.cols(); ++j) {
			if (PathUtils::contains_open(P, E(0, j), E(1, j), i, circular)) {
				exclude = true;
				break;
			}
		}

		if (exclude) {
			continue;
		}

		// testing inflection + length 1
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

Eigen::Index find_closest_corner_before_strict(
	const mat2x& B,      // list of points
	const vecXi& P,      // polygon vertices
	const Eigen::Index q // corner limit (index to B)
) {
	Index min_d = numeric_limits<Index>::max();
	Index min_v = -1;

	for (Index i = 0; i < P.size(); ++i) {
		const Index d = CircularDist(B, P(i), q);

		if (d > 0 && d < min_d) {
			min_d = d;
			min_v = i;
		}
	}

	return min_v;
}

// Given a list of points and a list of vertices with the same orientation, it finds the first 
// vertex of P which lies after the queried point, excluding the point itself.
Eigen::Index find_closest_corner_after_strict(
	const mat2x& B,      // list of points
	const vecXi& P,      // polygon vertices
	const Eigen::Index q // corner limit (index to B)
) {
	Index min_d = numeric_limits<Index>::max();
	Index min_v = -1;

	for (Index i = 0; i < P.size(); ++i) {
		const Index d = CircularDist(B, q, P(i));
		if (d > 0 && d < min_d) {
			min_d = d;
			min_v = i;
		}
	}

	return min_v;
}

NAMESPACE_END(PathUtils)
NAMESPACE_END(polyvec)