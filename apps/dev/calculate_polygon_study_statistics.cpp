#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/system.hpp>
#include <polyvec/utils/matrix.hpp>
#include <polyvec/geometry/path.hpp>
#include <polyvec/geometry/raster.hpp>
#include <polyvec/geometry/line.hpp>

#include "../dev/drawing.hpp"

using namespace std;
using namespace Eigen;
using namespace polyvec;
using namespace polyfit;

const char* results_small [] = {
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-all-points/tmp/benjamin",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-all-points/tmp/gabriele",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-all-points/tmp/polina",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-all-points/tmp/daniele",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-all-points/tmp/jerry"
};

const char* results_big[] = {
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-only-corners/tmp/helgen",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-only-corners/tmp/mikhail",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-only-corners/tmp/matteo",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-only-corners/tmp/fan",
    "D:/data/polyvec/svn/binary-tracing-study/results/2019-12-2-only-corners/tmp/paolo"
};

enum PointType {
    POINT_TYPE_FLAT,
    POINT_TYPE_CORNER,
    POINT_TYPE_CORNER_NHB,
    POINT_TYPE_MIDPOINT
};

void read_points(
    const string& file,
    Matrix2Xd& points
) {
    ifstream fs(file);
    assert(fs.good());

    string line;
    while (getline(fs, line)) {
        double x, y;
        sscanf(line.c_str(), "%lf %lf", &x, &y);
        MatrixUtils::append(points, Vector2d(x, y));
    }
}

void read_and_classify_vertices(
    const string& file,
    const Matrix2Xd& points,
    VectorXi& vertices,
    vector<PointType>& types
) {
    ifstream fs(file);
    assert(fs.good());

    string line;
    while (getline(fs, line)) {
        int v;
        sscanf(line.c_str(), "%d", &v);
        MatrixUtils::append(vertices, v);
    }

    vector<int> turns;
    PathUtils::compute_convexities(points, turns);
    
    types.clear();
    types.resize(vertices.size());
    for (Index i = 0; i < vertices.size(); ++i) {
        if (vertices(i) >= points.cols()) {
            const int v_prev = vertices(i) - points.cols();
            const int v_next = (vertices(i) + 1) % points.cols();

            const int turn_prev = turns[v_prev];
            const int turn_next = turns[v_next];

            if (turn_prev != 0 && turn_next ) {
                const int len_prev = CircularDist(points, PathUtils::next_transition_vertex(turns, v_prev, -1), v_prev);
                const int len_next = CircularDist(points, v_next, PathUtils::next_transition_vertex(turns, v_next, +1));

                if (len_prev == len_next) {
                    types[i] = POINT_TYPE_MIDPOINT;
                } else {
                    types[i] = POINT_TYPE_CORNER_NHB;
                }
            } else if (turn_prev != 0 || turn_next != 0) {
                types[i] = POINT_TYPE_CORNER_NHB;
            } else {
                types[i] = POINT_TYPE_FLAT;
            }

        } else {
            const int turn = turns[vertices(i)];
            const int turn_prev = CircularAt(turns, vertices(i) - 1);
            const int turn_next = CircularAt(turns, vertices(i) + 1);
            
            if (turn != 0) {
                types[i] = POINT_TYPE_CORNER;
            }else if (turn_prev != 0 && turn_next != 0) {
                types[i] = POINT_TYPE_MIDPOINT;
            } else if (turn_prev != 0 || turn_next != 0) {
                types[i] = POINT_TYPE_CORNER_NHB;
            } else {
                types[i] = POINT_TYPE_FLAT;
            }
        }
    }
}

// Murdered version of the original function which exports more information and doesn't
// terminate at the first pixel failing
#define PIXEL_DIAGONAL_LEN (1.41421356237 * .5)
#define PIXEL_DIAGONAL_LEN_SQ (PIXEL_DIAGONAL_LEN * PIXEL_DIAGONAL_LEN)
PV_INLINE vec4 distance_bounds_from_points(
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
        return vec4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX);
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
    vec2 error_max_manhattan(0., 0.);

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
                error_max_manhattan(0) = error_max_manhattan(1) = INFINITY;
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
        const double pmid_test_out_proj_t = polyvec::LineUtils::project_t(pmid_test_out, p_src, p_dst);
        const double pmid_test_in_proj_t = polyvec::LineUtils::project_t(pmid_test_in, p_src, p_dst);

        const double d_project_out = (pmid_test_out_proj - pmid_test_out).cwiseAbs().maxCoeff();
        const double d_project_in = (pmid_test_in_proj - pmid_test_in).cwiseAbs().maxCoeff();

        // Testing the sign of the greater of the two errors and updating the respecting value tracker
        int violation_sign;
        if (d_project_out < d_project_in + PF_EPS && pmid_test_out_proj_t > 0. && pmid_test_out_proj_t < 1.) {
            violation_sign = +1;
        } else {
            violation_sign = -1;
        }

        //int violation_sign = d_project_out < d_project_in ? +1 : -1;
        if (violation_sign > 0) {
            error_max(0) = ::std::max(error_max(0), abs(.5 + d_project_out));
            error_max_manhattan(0) = std::max(error_max(0), (pmid_test_out_proj - pmid_test_out).cwiseAbs().sum());
        }
        else {
            error_max(1) = ::std::max(error_max(1), abs(.5 + d_project_in));
            error_max_manhattan(1) = ::std::max(error_max(1), (pmid_test_in_proj - pmid_test_in).cwiseAbs().sum());
        }
    }

    return vec4(error_max.minCoeff(), error_max.maxCoeff(), error_max_manhattan.minCoeff(), error_max_manhattan.maxCoeff());
}

int main() {
    const int n_results_small = sizeof(results_small) / sizeof(results_small[0]);
    const int n_results_big = sizeof(results_big) / sizeof(results_big[0]);

    const int n_shapes_small = 7;
    const int n_shapes_big = 14;

    int n_corners = 0;
    int n_corner_nhbs = 0;
    int n_midpoints = 0;
    int n_total_vertices = 0;

    for (int i = 0; i < n_results_small; ++i) {
        for (int shape = 0; shape < n_shapes_small; ++shape) {
            char raster_txt[1024], polygon_txt[1024];
            sprintf(raster_txt, "%s/shape_%03d/image.png.raster.txt", results_small[i], shape);
            sprintf(polygon_txt, "%s/shape_%03d/image.png.polygon.1", results_small[i], shape);

            Matrix2Xd points;
            read_points(raster_txt, points);

            VectorXi vertices;
            vector<PointType> point_types;
            read_and_classify_vertices(polygon_txt, points, vertices, point_types);

            for (Index vidx = 0; vidx < vertices.size(); ++vidx) {
                n_corners += point_types[vidx] == POINT_TYPE_CORNER ? 1 : 0;
                n_corner_nhbs += point_types[vidx] == POINT_TYPE_CORNER_NHB ? 1 : 0;
                n_midpoints += point_types[vidx] == POINT_TYPE_MIDPOINT ? 1 : 0;
            }
            n_total_vertices += vertices.size();
        }
    }

    const float corner_pct = (float)n_corners / n_total_vertices;
    const float corner_nhb_pct = (float)n_corner_nhbs / n_total_vertices;
    const float midpoint_pct = (float)n_midpoints / n_total_vertices;
    const int n_remaining = n_total_vertices - n_corners - n_corner_nhbs - n_midpoints;
    const float n_remaining_pct = (float)n_remaining / n_total_vertices;
    printf("small study\n");
    printf("corner percentage     %f (%d / %d)\n", corner_pct, n_corners, n_total_vertices);
    printf("corner nhb percentage %f (%d / %d)\n", corner_nhb_pct, n_corner_nhbs, n_total_vertices);
    printf("midpoint percentage   %f (%d / %d)\n", midpoint_pct, n_midpoints, n_total_vertices);
    printf("remaining percentage  %f (%d / %d)\n", n_remaining_pct, n_remaining, n_total_vertices);
    printf("\n\n\n");

    int n_edges_with_wrong_classification = 0;
    int n_edges_with_opposite_sides = 0;
    int n_edges_within_distance_1 = 0;
    int n_total_edges = 0;

    for (int i = 0; i < n_results_big; ++i) {
        for (int shape = 0; shape < n_shapes_big; ++shape) {
            char raster_txt[1024], polygon_txt[1024];
            sprintf(raster_txt, "%s/shape_%03d/image.png.raster.txt", results_big[i], shape);
            sprintf(polygon_txt, "%s/shape_%03d/image.png.polygon.1", results_big[i], shape);

            Matrix2Xd points;
            read_points(raster_txt, points);
            
            VectorXi vertices;
            vector<PointType> point_types;
            read_and_classify_vertices(polygon_txt, points, vertices, point_types);
            
            // eigen < 3.4
            vector<int> vertices_buf(vertices.size());
            for (Index j = 0; j < vertices.size(); ++j) {
                vertices_buf[j] = vertices(j);
            }

            sort(vertices_buf.begin(), vertices_buf.end(), 
                [&points](const int& v0_weird, const int& v1_weird) {
                    const int v0 = v0_weird >= points.cols() ? (v0_weird - points.cols()) : v0_weird;
                    const int v1 = v1_weird >= points.cols() ? (v1_weird - points.cols()) : v1_weird;
                    return v0 < v1;
                });

            for (size_t j = 0; j < vertices_buf.size(); ++j) {
                vertices(j) = vertices_buf[j];
            }

            vector<int> turns;
            PathUtils::compute_convexities(points, turns);

            for (Index k = 0; k < vertices.size(); ++k) {
                const int v_src_weird = vertices(k);
                const int v_dst_weird = CircularAt(vertices, k + 1);

                int v_src = v_src_weird >= points.cols() ? ((v_src_weird + 1) % points.cols()) : v_src_weird;
                int v_dst = v_dst_weird >= points.cols() ? (v_dst_weird - points.cols()) : v_dst_weird;
                
                // shit can be ordered randomly 
                //bool flip = CircularDist(points, v_src, v_dst) > CircularDist(points, v_dst, v_src);
                bool flip = false;
                if (flip) {
                    swap(v_src, v_dst);
                }

                vec2 p_src;
                if (v_src_weird >= points.cols()) {
                    p_src = .5 * (points.col(v_src_weird - points.cols()) + points.col((v_src_weird + 1) % points.cols()));
                } else {
                    p_src = points.col(v_src_weird);
                }

                vec2 p_dst;
                if (v_dst_weird >= points.cols()) {
                    p_dst = .5 * (points.col(v_dst_weird - points.cols()) + points.col((v_dst_weird + 1) % points.cols()));
                } else {
                    p_dst = points.col(v_dst_weird);
                }

                if (flip) {
                    swap(p_src, p_dst);
                }

                vec2 err = distance_bounds_from_points(points, p_src, p_dst, v_src, v_dst).segment(0, 2);
                if (err.maxCoeff() > .5 + PF_EPS) {
                    ++n_edges_with_wrong_classification;

                    if (err(0) > .5 + PF_EPS && err(1) > .5 + PF_EPS) {
                        ++n_edges_with_opposite_sides;
                    }
                }

                vec2 dirs[4] = {
                    Vector2d(1., 0.),
                    Vector2d(-1., 0.),
                    Vector2d(0., 1.),
                    Vector2d(0., -1.)
                };

                static int draw = 0;
                if (!draw && shape == 6 && i == 3) {
                    DevicePDF* pdf = new DevicePDF("D:/data/polyvec/out/castle-debug.svg", 1, 1);
                    draw_raster(points);
                    draw_raster_indices(points);
                    for (Index j = 0; j < vertices.size(); ++j) {
                        const int v_src_weird = vertices(j);
                        const int v_dst_weird = CircularAt(vertices, j + 1);

                        vec2 p_src;
                        if (v_src_weird >= points.cols()) {
                            p_src = .5 * (points.col(v_src_weird - points.cols()) + points.col((v_src_weird + 1) % points.cols()));
                        }
                        else {
                            p_src = points.col(v_src_weird);
                        }

                        vec2 p_dst;
                        if (v_dst_weird >= points.cols()) {
                            p_dst = .5 * (points.col(v_dst_weird - points.cols()) + points.col((v_dst_weird + 1) % points.cols()));
                        }
                        else {
                            p_dst = points.col(v_dst_weird);
                        }

                        draw::line(p_src, p_dst, Style::outline(colors::red, .75));
                    }

                    pdf->draw(0, 0);
                    delete pdf;

                    draw = 1;
                }

                if (shape == 6 && i == 3 && v_src == 8 && v_dst == 22) {
                    printf("moo");
                }

                double max_manhattan_distance = 0.;
                int v = v_src;
                while (v != v_dst) {
                    double max_distance_pixel = 0.;

                    const vec2 p_corner = points.col(v);
                    for (int dir = 0; dir < 4; ++dir) {
                        const vec2 p_dir = dirs[dir];
                        double d_ray, d_line;
                        if (LineUtils::intersect(p_corner, p_corner + p_dir * 1.5, p_src, p_dst, d_ray, d_line)) {
                            if (d_ray < -PF_EPS || d_ray > 1. + PF_EPS) {
                                continue;
                            }

                            max_distance_pixel = max(max_distance_pixel, d_ray * 1.5);
                        }
                    }

                    max_manhattan_distance = max(max_manhattan_distance, max_distance_pixel);
                    v = Circular(points, v + 1);
                }

                if (max_manhattan_distance < 1. + PF_EPS_MEDIUM) {
                    ++n_edges_within_distance_1;
                }
            }

            n_total_edges += vertices.size();
        }
    }

    const float edges_within_1_pct = (float)n_edges_within_distance_1 / n_total_edges;
    const float edges_wrong_classification_pct = (float)n_edges_with_wrong_classification / n_total_edges;
    const float edges_opposite_sides_pct = (float)n_edges_with_opposite_sides / n_total_edges;

    printf("study big\n");
    printf("edges within manhattan 1 percentage %f (%d / %d)\n", edges_within_1_pct, n_edges_within_distance_1, n_total_edges);
    printf("edges with wrong classification percentage %f (%d / %d)\n", edges_wrong_classification_pct, n_edges_with_wrong_classification, n_total_edges);
    printf("edges with opposite sides pct %f (%d / %d)\n", edges_opposite_sides_pct, n_edges_with_opposite_sides, n_total_edges);
    return EXIT_SUCCESS;
}