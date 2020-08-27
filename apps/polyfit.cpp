const char* usage =
    "<image_uri>      Path to the image file to be vectorized\n"
    "<classifier_uri> Trained corner classification model\n"
    "<csv_file>       Path to the file where the CSV data will be written.";

// Polyvec
#include <polyvec/pipeline_helper.hpp>
#include <polyvec/utils/string.hpp>
#include <polyvec/utils/system.hpp>
#include <polyvec/core/log.hpp>

#include <polyvec/polygon-tracer/multi_polygon.hpp>
#include <polyvec/mc/get_bounding_box.hpp>
#include <polyvec/mc/get_pixel_regions.hpp>
#include <polyvec/mc/raster_image_connectivity.hpp>
#include <polyvec/mc/junction.hpp>
#include <polyvec/mc/segment_connectivity.hpp>

#include "dev/drawing.hpp"
#include "dev/draw_paper_figures.hpp"
#include <polyvec/io/svg.hpp>

#define EXPORT_DEBUG_SVG 0

// libc++
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;
using namespace polyfit;
using namespace polyvec;
using namespace Eigen;

namespace {
    void draw_all_boundaries(polyfit::mc::RasterImageConnectivity& conn) {
        for (int j = 0; j < conn.n_regions(); ++j) {
            std::vector<int> holes;
            int outer_boundary;
            conn.region_to_polygon(j, outer_boundary, holes);
            draw::polyline(conn.get_polygon_points(outer_boundary).cast<double>(), Style::outline(colors::gray, .75));
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Expected <input image> <serialized classifier> <output directory>");
        return EXIT_FAILURE;
    }

    const char* image_uri = argv[1];
    const char* classifier_uri = argv[2];
    const char* write_dir = argv[3];

    os::make_dir(write_dir);
    Log::open(stdout, Log::CHANNEL_DEV);

    std::vector<Matrix2Xi> boundaries;
    std::vector<Vector4d>  colors;
    polyfit::mc::get_pixel_regions_from_bmp(image_uri, boundaries, colors);
    PF_STATUS_F("boundaries %d", boundaries.size());
    PF_STATUS_F("colors %d", colors.size());

    polyfit::mc::RasterImageConnectivity raster_conn = polyfit::mc::RasterImageConnectivity::build(boundaries);
    std::vector<polyfit::mc::JunctionType> junction_types;
    polyfit::mc::find_junctions(raster_conn, junction_types);
    mc::SegmentConnectivity segment_conn = mc::SegmentConnectivity::build(raster_conn, junction_types);

    MultiPolygonTracer tracer(raster_conn, segment_conn, colors, classifier_uri);
    tracer.trace();
    std::vector<int> boundaries_by_area = tracer.get_boundaries_ordered_by_area();
    PF_DEV_F("Exporting %llu boundaries", boundaries_by_area.size());

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "input.svg").c_str(), 1, 1);

        for (const int boundary_id : boundaries_by_area) {
            draw_raster_closed(tracer.boundaries()[boundary_id].raster_points_as_double(), tracer.boundaries()[boundary_id].color());
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        const string polygons_uri = StringUtils::join_path(write_dir, "polygons.svg");
        DevicePDF* pdf = new DevicePDF(polygons_uri.c_str(), 1, 1);
        draw_all_boundaries(raster_conn);

        for (size_t i = 0; i < tracer.boundaries().size(); ++i) {
            const ImageBoundary& boundary = tracer.boundaries()[i];
            draw_raster(boundary.polygon_points(), colors::talking_orange);
            draw_raster_indices(boundary.polygon_points());

            const auto& P = boundary.polygon_points();
            vec2 com = vec2::Zero();
            for (Index j = 0; j < boundary.polygon_points().cols(); ++j) {
                com += boundary.polygon_points().col(j);
            }

            com /= boundary.polygon_points().cols();
            //draw::text(com, to_string(i), draw::font_pdf, Style::text());
        }

        pdf->draw(0, 0);
        delete pdf;
    }


    // -----------------------------------------------------------------------------------------------------------------
    {
        const string curves_color_uri = StringUtils::join_path(write_dir, "curves_closed.svg");
        DevicePDF* pdf = new DevicePDF(curves_color_uri.c_str(), 1, 1);
        //draw_raster_background(B, Style::outline(colors::gray * 1.75, 2.5));
        for (size_t i = 0; i < boundaries_by_area.size(); ++i) {

            const auto& boundary = tracer.boundaries()[boundaries_by_area[i]];
            //draw_curve_primitives(boundary.spline().primitives, boundary.color());
            draw_curve_primitives_closed(boundary.spline().primitives, boundary.color());
            PF_DEV_F("Boundary color %f %f %f %f points %d", boundary.color()(0), boundary.color()(1), boundary.color()(2), boundary.color()(3), boundary.boundary_points().size());
        }
        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "boundary-graph.svg").c_str(), 1, 1);
        draw_all_boundaries(raster_conn);
        for (const auto& boundary : tracer.boundaries()) {
            for (const auto& e : boundary.edges()) {
                draw::line(boundary.boundary_points().col(e.v0), boundary.boundary_points().col(e.v1), Style::outline(colors::forest_green, 5.));
            }
        }
        pdf->draw(0, 0);
        delete pdf;
    }
    if (0)
    {
        for (int j = 0; j < raster_conn.n_polygons(); ++j) {
            const string svg_uri = StringUtils::join_path(write_dir, StringUtils::fmt("region-%d.svg", j).c_str());
            DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);
            draw_all_boundaries(raster_conn);
            std::vector<int> holes;
            draw::polyline(raster_conn.get_polygon_points(j).cast<double>(), Style::outline(colors::red, 2.5));

            pdf->draw(0, 0);
            delete pdf;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        for (auto& boundary : tracer.boundaries()) {
            const string svg_uri = StringUtils::join_path(write_dir, StringUtils::fmt("regularity-graph-%d.svg", boundary.id()).c_str());
            DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);

            const auto& PP = boundary.polygon_points();
            const auto& P = boundary.polygon_vertices();
            const auto& B = boundary.boundary_points();

            draw_raster_closed(B, real3(colors::gray * 1.9));
            draw_raster_background(B, Style::outline(colors::gray * 1.75, 2.5), false);
            draw_raster(PP);

            draw_raster_polygon_indices(B, P);

            const vec3 color_parallel = colors::forest_green;
            for (auto& r : boundary.regularity_graph().parallels()) {
                const vec2 p0 = .5 * (PP.col(r.v00) + PP.col(r.v01));
                const vec2 p1 = .5 * (PP.col(r.v10) + PP.col(r.v11));
                draw::line(p0, p1, Style::outline(color_parallel, 7.5, LineType::Dash));

                if (r.aligned_00_11) {
                    draw::line(PP.col(r.v00), PP.col(r.v11), Style::outline(color_parallel, 3.5, LineType::Dash));
                }

                if (r.aligned_01_10) {
                    draw::line(PP.col(r.v01), PP.col(r.v10), Style::outline(color_parallel, 3.5, LineType::Dash));
                }
            }
			for (auto& r : boundary.regularity_graph().continuations()) {
				const vec2 p0 = PP.col(r.v0);
				const vec2 p1 = PP.col(r.v1);
				draw::point(p0, .1, Style::fill(colors::red));
				draw::point(p1, .1, Style::fill(colors::red));
				draw::line(p0, p1, Style::outline(colors::red, 2.5));
			}
            pdf->draw(0, 0);
            delete pdf;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            const string svg_uri = StringUtils::join_path(write_dir, StringUtils::fmt("tangent-fits-%d.svg", boundary.id()).c_str());
            DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);

            for (Index i = 0; i < boundary.polygon_vertices().size(); ++i) {
                const vec2 pp = CircularAt(boundary.polygon_points(), i - 1);
                const vec2 p = boundary.polygon_points().col(i);
                const vec2 pn = CircularAt(boundary.polygon_points(), i + 1);

                vec3 color;
                if (boundary.tangents_fits()[i] == TANGENT_FIT_CONSTANT) {
                    color = colors::red;
                }
                else if (boundary.tangents_fits()[i] == TANGENT_FIT_LERP_SYM) {
                    color = colors::talking_orange;
                }
                else if (boundary.tangents_fits()[i] == TANGENT_FIT_LERP) {
                    color = colors::forest_green;
                }

                draw::line((pp + p) * .5, p, Style::outline(color, .5));
                draw::line((pn + p) * .5, p, Style::outline(color, .5));
            }

            pdf->draw(0, 0);
            delete pdf;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            const string svg_uri = StringUtils::join_path(write_dir, StringUtils::fmt("curves-%d.svg", boundary.id()).c_str());
            DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);

            draw_all_boundaries(raster_conn);
            //draw_curve_primitives(boundary.spline().primitives, TangentFitColorFunctor(boundary.tangents_fits()));
            draw_raster(boundary.polygon_points());
            draw_curve_primitives(boundary.spline().primitives, AlternatingColorFunctor());
            for (size_t i = 0; i < boundary.spline().primitives.size(); ++i) {
                const auto curve = boundary.spline().primitives[i].curve->get_curve();
                const vec2 l = curve->pos(.5);
                const string label = StringUtils::fmt("%d(%d)", i, boundary.spline().primitives[i].corner);
                draw::text(l, label, draw::font_pdf / 3, Style::text());
                draw::point(curve->pos(0.), .075, Style::fill(colors::black));
            }

            pdf->draw(0, 0);
            delete pdf;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "boundaries.svg").c_str(), 1, 1);

        for (size_t i = 0; i < tracer.boundaries().size(); ++i) {
            const ImageBoundary& boundary = tracer.boundaries()[i];
            draw_raster(boundary.boundary_points());
            draw_raster_indices(boundary.boundary_points());

            vec2 com = vec2::Zero();
            for (Index j = 0; j < boundary.polygon_points().cols(); ++j) {
                com += boundary.polygon_points().col(j);
            }

            com /= boundary.polygon_points().cols();
            draw::text(com, to_string(i), draw::font_pdf * 4, Style::text());
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "convexities.svg").c_str(), 1, 1);

        for (size_t i = 0; i < tracer.boundaries().size(); ++i) {
            const ImageBoundary& boundary = tracer.boundaries()[i];
            draw_raster(boundary.boundary_points());

            for (size_t j = 0; j < boundary.boundary_points().cols(); ++j) {
                draw::text(boundary.boundary_points().col(j), to_string(boundary.convexities()[j]), draw::font_pdf, Style::text());
            }
        }

        pdf->draw(0, 0);
        delete pdf;
    }


    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        for (size_t i = 0; i < tracer.boundaries().size(); ++i) {
            const ImageBoundary& boundary = tracer.boundaries()[i];
            DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, StringUtils::fmt("boundary-graph-original-%llu.svg", i).c_str()).c_str(), 1, 1);

            const auto color = colors::random();
            draw_all_boundaries(raster_conn);
            for (size_t j = 0; j < boundary.edges_original().size(); ++j) {
                const vec2 p0 = boundary.boundary_points().col(boundary.edges_original()[j].v0);
                const vec2 p1 = boundary.boundary_points().col(boundary.edges_original()[j].v1);
                draw::line(p0, p1, Style::outline(color, 1., LineType::Solid));
            }

            pdf->draw(0, 0);
            delete pdf;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            const string svg_uri = StringUtils::join_path(write_dir, StringUtils::fmt("polygon-%d.svg", boundary.id()).c_str());
            DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);

            draw_all_boundaries(raster_conn);
            draw_raster(boundary.polygon_points(), colors::talking_orange);
            draw_raster_polygon_indices(boundary.boundary_points(), boundary.polygon_vertices());

            pdf->draw(0, 0);
            delete pdf;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        for (size_t i = 0; i < tracer.boundaries().size(); ++i) {
            const ImageBoundary& boundary = tracer.boundaries()[i];
            DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, StringUtils::fmt("boundary-graph-%llu.svg", i).c_str()).c_str(), 1, 1);

            const auto color = colors::random();
            draw_all_boundaries(raster_conn);
            draw_raster_indices(boundary.boundary_points());
            for (size_t j = 0; j < boundary.edges().size(); ++j) {
                    const vec2 p0 = boundary.boundary_points().col(boundary.edges()[j].v0);
                    const vec2 p1 = boundary.boundary_points().col(boundary.edges()[j].v1);
                    draw::line(p0, p1, Style::outline(color, 1., LineType::Solid));
                }

            pdf->draw(0, 0);
            delete pdf;
            }
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        const string curves_uri = StringUtils::join_path(write_dir, "curves_outline.svg");
        DevicePDF* pdf = new DevicePDF(curves_uri.c_str(), 1, 1);
        draw_all_boundaries(raster_conn);
        for (size_t i = 0; i < boundaries_by_area.size(); ++i) {
            const auto& boundary = tracer.boundaries()[boundaries_by_area[i]];
            //draw_raster(boundary.polygon_points(), colors::talking_orange);
            //draw_curve_primitives(boundary.spline().primitives, ConstantColorFunctor(boundary.color().segment(0, 3)));
            draw_curve_primitives(boundary.spline().primitives, AlternatingColorFunctor());
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        const string svg_uri = StringUtils::join_path(write_dir, "segments.svg");
        DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);
        for (size_t i = 0; i < tracer.segments().size(); ++i) {
            const Eigen::Matrix2Xd pts = tracer.segments()[i].cast<double>();
            const auto color = colors::random4();
            //draw_raster(pts, color);
            for (Index j = 0; j < pts.cols() - 1; ++j) {
                draw::line(pts.col(j), pts.col(j + 1), Style::outline(color.segment(0, 3), 6.));
            }

            //draw::text(pts.col(pts.cols() / 2), to_string(i), draw::font_pdf, Style::text());
            draw::text(pts.col(1), to_string(i), draw::font_pdf / 4, Style::text());

            //draw::point(pts.col(0), .25, Style::fill(color.segment(0, 3)));
            //draw::point(pts.col(1), .25, Style::fill(color.segment(0, 3)));
            //draw::point(pts.col(pts.cols() - 1), .25, Style::fill(color.segment(0, 3)));
            //draw::point(pts.col(pts.cols() - 2), .25, Style::fill(color.segment(0, 3)));
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        const string svg_uri = StringUtils::join_path(write_dir, "flat_corners.svg");
        DevicePDF* pdf = new DevicePDF(svg_uri.c_str(), 1, 1);

        for (size_t i = 0; i < tracer.boundaries().size(); ++i) {
            const ImageBoundary& boundary = tracer.boundaries()[i];

            const auto& P = boundary.polygon_points();
            draw_raster(boundary.polygon_points(), colors::talking_orange);

            for (Index i = 0; i < P.cols(); ++i) {
                const vec2 pp = CircularAt(P, i - 1);
                const vec2 p = P.col(i);
                const vec2 pn = CircularAt(P, i + 1);
                const double angle = AngleUtils::spanned_shortest(pp, p, pn);

                const bool is_prev_aa = (pp - p).cwiseAbs().minCoeff() < PF_EPS;
                const bool is_next_aa = (pn - p).cwiseAbs().minCoeff() < PF_EPS;

                if (angle > PF_RAD(160) && !is_prev_aa && !is_next_aa) {
                    draw::point(P.col(i), .25, Style::fill(colors::red));
                }
            }
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "input.svg").c_str(), 1, 1);

        for (const int boundary_id : boundaries_by_area) {
            draw_raster_closed(tracer.boundaries()[boundary_id].raster_points_as_double(), tracer.boundaries()[boundary_id].color());
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "junctions.svg").c_str(), 1, 1);
        draw_all_boundaries(raster_conn);
        
        for (const JunctionInfo& info : tracer.junctions()) {
            const vec2 point = tracer.boundary_vertex_pos_from_id(info.vertex);
            draw::point(point, .025, Style::fill(colors::red));
            draw::text(point, to_string(info.vertex), draw::font_pdf, Style::text());
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    // -----------------------------------------------------------------------------------------------------------------
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "half_edges.svg").c_str(), 1, 1);
        draw_all_boundaries(raster_conn);

        for (const HalfEdge& he : tracer.half_edges()) {
            if (he.id_twin == -1) {
                continue;
            }
            
            const auto edge = BoundaryGraph::unpack_edge_id(he.edge);
            const vec2 p = tracer.boundaries()[he.boundary].polygon_points().col(edge(0));
            const vec2 pn = tracer.boundaries()[he.boundary].polygon_points().col(edge(1));

            draw::line(p, pn, Style::outline(colors::red, 1.5));
        }

        pdf->draw(0, 0);
        delete pdf;
    }
    
    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "highlight_flat_angles.svg").c_str(), 1, 1);
        const double highlight_limit = PF_RAD(165) - PF_EPS;
        draw_all_boundaries(raster_conn);

        for (const ImageBoundary& boundary : tracer.boundaries()) {
            const auto& P = boundary.polygon_points();
            for (Index i = 0; i < P.cols(); ++i) {
                const vec2 pp = CircularAt(P, i - 1);
                const vec2 p = P.col(i);
                const vec2 pn = CircularAt(P, i + 1);
                
                const double angle = AngleUtils::spanned_shortest(pp, p, pn);
                vec3 color = colors::black;
                if (angle > highlight_limit) {
                    color = colors::red;
                }

                draw::line(.5 * (pp + p), p, Style::outline(color, 7.5));
                draw::line(.5 * (p + pn), p, Style::outline(color, 7.5));
            }
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    if (0)
    {
        DevicePDF* pdf = new DevicePDF(StringUtils::join_path(write_dir, "highlight_short_edges.svg").c_str(), 1, 1);
        const double highlight_limit = std::max(raster_conn.pixel_vertex_dims().cast<double>().norm() / 20, 3.);
        draw_all_boundaries(raster_conn);

        for (const ImageBoundary& boundary : tracer.boundaries()) {
            const auto& P = boundary.polygon_points();
            for (Index i = 0; i < P.cols(); ++i) {
                const vec2 pp = CircularAt(P, i - 1);
                const vec2 p = P.col(i);
                const vec2 pn = CircularAt(P, i + 1);
                const vec2 pnn = CircularAt(P, i + 2);

                vec3 color = colors::black;
                const double ratio_len_prev = (p - pp).squaredNorm() / (pn - p).squaredNorm();
                const double ratio_len_next = (pnn - pn).squaredNorm() / (pn - p).squaredNorm();
                if ((ratio_len_prev > 4. + PF_EPS || ratio_len_next > 4. + PF_EPS) && (pn - p).norm() < highlight_limit) {
                    color = colors::red;
                }

                draw::line(p, pn, Style::outline(color, 7.5));
            }
        }

        pdf->draw(0, 0);
        delete pdf;
    }

    if (0)
    {
        FILE* fp = fopen(StringUtils::join_path(write_dir, "image.png.raster.txt").c_str(), "w");
        
        const auto& boundary = tracer.boundaries()[boundaries_by_area[1]]; 
        for (Index i = 0; i < boundary.boundary_points().cols(); ++i) {
            fprintf(fp, "%.1f %.1f\n", boundary.boundary_points().col(i)(0), boundary.boundary_points().col(i)(1));
        }

        fclose(fp);
    }

    if (0)
    {
        FILE* fp = fopen(StringUtils::join_path(write_dir, "image.png.polygon.1").c_str(), "w");
        const auto& boundary = tracer.boundaries()[boundaries_by_area[1]];
        for (Index i = 0; i < boundary.polygon_vertices().size(); ++i) {
            fprintf(fp, "%d\n", boundary.polygon_vertices()(i));
        }
        fclose(fp);
    }
    
#if 0
    int canvas_width = 0;
    int canvas_height = 0;
    for (const auto& boundary : tracer.boundaries()) {
        canvas_width = max(canvas_width, (int)boundary.boundary_points().row(0).maxCoeff());
        canvas_height = max(canvas_height, (int)boundary.boundary_points().row(1).maxCoeff());
    }

    draw_paper_figures::Data paper_figure_data;
    
    auto init_paper_figure_data_from_boundary = [&paper_figure_data, &tracer](const int boundary_id) {
        auto& boundary = tracer.boundaries()[boundary_id];
        paper_figure_data.B = &boundary.boundary_points();
        paper_figure_data.P = &boundary.polygon_vertices();
        paper_figure_data.PP = &boundary.polygon_points();
        paper_figure_data.E = &boundary.edges();
        paper_figure_data.curves = &boundary.spline().primitives;
        paper_figure_data.tangent_fits = &boundary.tangents_fits();
        paper_figure_data.curves_attempts = &boundary.fitting_attempts();
		paper_figure_data.color = boundary.color();
    };

    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_grid.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::grid(paper_figure_data);
        }
    }
    
    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_raster.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::raster_boundary(paper_figure_data);
        }
    }

    {    
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_graph.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::graph(paper_figure_data);
        }
    }

    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_polygon.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::polygon(paper_figure_data);
        }
    }
    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_polygon_curves.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::polygon_curves(paper_figure_data);
        }
    }
    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_curves_primitives_classification.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::curve_primitives_classification(paper_figure_data);
        }
    }
    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_curves_primitives_optimized.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::curve_primitives(paper_figure_data);
        }
    }
    //{
    //    io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_curves_fill_classification.svg").c_str(), canvas_width, canvas_height);
    //    for (const ImageBoundary& boundary : tracer.boundaries()) {
    //        init_paper_figure_data_from_boundary(boundary.id());
    //        draw_paper_figures::curve_fill_classification(paper_figure_data);
    //    }
    //}

    {
        io::SvgCanvas svg(StringUtils::join_path(write_dir, "paper_figure_curves_fill_optimized.svg").c_str(), canvas_width, canvas_height);
        for (const ImageBoundary& boundary : tracer.boundaries()) {
            init_paper_figure_data_from_boundary(boundary.id());
            draw_paper_figures::curve_fill(paper_figure_data);
        }
    }
#endif

    return EXIT_SUCCESS;
}