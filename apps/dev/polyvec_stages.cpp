/*
	Writes a vector file (.svg) for the type of trace specified
	<type>   type of trace (see the implementation)
	<image>  address of the input image
	<output> address of the directory where the files should be written 
*/
const char* usage =
	"./polyvec-stages <image> <classifier> [<stage>=<svg>, ...]\n"
	"<image>      - Image file.\n"
	"<classifier> - Serialized classifier used for inferring the corner type.\n"
	"<stage>      - One of (input, graph, polygon, polygon_closed, regularity_graph, curves, curves_closed, ...).\n"
	"<vector>     - Output svg address.\n"
	"";	

// polyvec
#include <polyvec/pipeline_helper.hpp>
#include "../dev/drawing.hpp"
#include "../dev/draw_paper_figures.hpp"
#include <polyvec/io/svg.hpp>

// todo: remove
#include <polyvec/regularity/check_opposite_points_winding_number.hpp>

// libc++
#include <cstdlib>

using namespace std;
using namespace polyvec;
using namespace polyfit;

struct ArgExport {
	string stage;
	string svg_uri;
};

struct CandidateMatch {
	// Opposite corners
	int  i;
	int  j;

	// Squared distance between the two corners
	double dist_sq;
};

auto state_color_prob = [](const CurveSequenceFitter::EvolutionaryFittingState& state, int iCorner) -> real3
{
	if (state.corner_good_probability.empty())
		return colors::dark_gray;
	auto prob = state.corner_good_probability[iCorner];
	if (prob > 0.75)
		return colors::forest_green;
	if (prob > 0.5)
		return colors::forest_green * 1.8;
	return colors::red;
};

int main(int argc, char** argv) {
	if (argc < 2) {
		fprintf(stderr, usage);
		return EXIT_FAILURE;
	}

	std::vector<ArgExport> args;

	for (int i = 3; i < argc; ++i) {
		const string text = string(argv[i]);
		const size_t token_split = text.find_first_of("=");
		if (token_split == std::string::npos) {
			fprintf(stderr, "skipping unrecognized argument %s", argv[i]);
			continue;
		}

		ArgExport arg;
		arg.stage = text.substr(0, token_split);
		arg.svg_uri = text.substr(token_split + 1);
		args.emplace_back(move(arg));
	}

	Log::open((char*)stdout, Log::CHANNEL_DEV);
	//Log::open((char*)stdout, Log::CHANNEL_DEV | Log::CHANNEL_STATUS | Log::CHANNEL_VERBOSE);

	const char* image_uri = argv[1];
	const char* classifier_uri = argv[2];

	init_pipeline();
	InputData input_data = load_input(image_uri);
	PolygonData polygon_data = extract_polygon(input_data);
	std::vector<AttemptInfo> attempts;
	SplineData spline_data = extract_curves_with_classifier(polygon_data, classifier_uri, &attempts);

    PostProcessor post_processor = post_process(polygon_data, spline_data, classifier_uri);    

    const auto& B = polygon_data.B;
    //const int canvas_width = B.row(0).maxCoeff() - B.row(0).minCoeff();
    //const int canvas_height = B.row(1).maxCoeff() - B.row(1).minCoeff();
    const int canvas_width = B.row(0).maxCoeff();
    const int canvas_height = B.row(1).maxCoeff();
    draw_paper_figures::Data paper_figure_data;
    paper_figure_data.B = &polygon_data.B;
    paper_figure_data.P = &polygon_data.P;
    paper_figure_data.PP = &polygon_data.PP;
    paper_figure_data.E = &polygon_data.E;
    paper_figure_data.curves = &spline_data.curves.primitives;
    paper_figure_data.tangent_fits = &spline_data.tangent_fits;
    paper_figure_data.curves_attempts = &attempts;

	for (ArgExport& arg : args) {
		printf("export stage %s svg %s\n", arg.stage.c_str(), arg.svg_uri.c_str());

		if (arg.stage.compare("input") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(input_data.R);
			pdf->draw(0, 0);
			delete pdf;
		}
        else if (arg.stage.compare("boundary_graph") == 0) {
            DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);

			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
            draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);
			//draw_raster(polygon_data.B, Style::outline(colors::gray * 1.5, 5));
			draw_raster(polygon_data.PP, Style::outline(colors::gray, 15.0));
            draw_raster_indices(polygon_data.B);

            for (int i = 0; i < polygon_data.E.size(); ++i) {
                const auto& e = polygon_data.E[i];
                const vec2 psrc = polygon_data.B.col(e.v0);
                const vec2 pdst = polygon_data.B.col(e.v1);
                draw::line(psrc, pdst, Style::outline(colors::forest_green, 1.5));
            }

#if KEEP_ERASED_ENTRIES
			for (auto& e : ErasedEntries<BoundaryGraph::Edge>::entries) {
				const vec2 psrc = polygon_data.B.col(e.v0);
				const vec2 pdst = polygon_data.B.col(e.v1);
				draw::line(psrc, pdst, Style::outline(vec3(0.85, 0.37, 0.01), 3.0));
			}
#endif

            pdf->draw(0, 0);
            delete pdf;
        }
        else if (arg.stage.compare("boundary") == 0) {
            DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);			
            draw_raster_indices(polygon_data.B);
            pdf->draw(0, 0);
            delete pdf;
        }
		else if (arg.stage.compare("polygon") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5));			
			draw_raster(polygon_data.PP);
            draw_raster_indices(polygon_data.PP);
			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("polygon_closed") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.PP);
			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("opposite_corners") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			
			draw_raster(polygon_data.PP);

			std::vector<int> C;
			PathUtils::compute_convexities(polygon_data.PP, C);

			std::vector<CandidateMatch> candidates;

			for (int i = 0; i < (int)polygon_data.PP.cols(); ++i) {
				double dist_sq_min = INFINITY;
				int i_opp_min = -1;
				for (int i_opp = i + 1; i_opp < polygon_data.PP.cols(); ++i_opp) {
					if (C[i] != -1 || C[i_opp] != -1) {
						continue;
					}

					if (ShortestDist(polygon_data.PP, i, i_opp) < 2) {
						continue;
					}

					if (!Regularity::check_opposite_points_winding_number(polygon_data.PP, i, i_opp)) {
						continue;
					}

					const double dist_sq = (polygon_data.PP.col(i) - polygon_data.PP.col(i_opp)).squaredNorm();

					if (dist_sq < dist_sq_min) {
						dist_sq_min = dist_sq;
						i_opp_min = i_opp;
					}
				}

				if (i_opp_min != -1) {
					CandidateMatch match;
					match.i = i;
					match.j = i_opp_min;
					match.dist_sq = dist_sq_min;
					candidates.emplace_back(move(match));
				}
			}

			std::vector<size_t> candidates_delete;
			for (int i = 0; i < candidates.size(); ++i) {
				for (int j = 0; j < candidates.size(); ++j) {
					const auto& ci = candidates[i];
					const auto& cj = candidates[j];

					if (ci.i == cj.i || ci.i == cj.j || ci.j == cj.i || ci.j == cj.j) {
						if (ci.dist_sq > cj.dist_sq) {
							candidates_delete.emplace_back(i);
							break;
						}
					}
				}
			}

			EraseOrdered(candidates, candidates_delete);

			for (auto& c : candidates) {
				const vec2 pi = polygon_data.PP.col(c.i);
				const vec2 pj = polygon_data.PP.col(c.j);
				draw::line(pi, pj, Style::outline(colors::red, 3.5));
			}

			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("regularity_graph") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);
			draw_raster(polygon_data.PP);            

			draw_raster_polygon_indices(polygon_data.B, polygon_data.P);

			const vec3 color_parallel = colors::forest_green;
			const vec3 color_closure_aligned = colors::red;
            const vec3 color_symmetry = colors::calm_blue;

			for (auto& r : polygon_data.RE.parallels()) {
				const vec2 p0 = .5 * (polygon_data.PP.col(r.v00) + polygon_data.PP.col(r.v01));
				const vec2 p1 = .5 * (polygon_data.PP.col(r.v10) + polygon_data.PP.col(r.v11));
				draw::line(p0, p1, Style::outline(color_parallel, 7.5, LineType::Dash));

				if (r.aligned_00_11) {
					draw::line(polygon_data.PP.col(r.v00), polygon_data.PP.col(r.v11), Style::outline(color_parallel, 3.5, LineType::Dash));
				}

				if (r.aligned_01_10) {
					draw::line(polygon_data.PP.col(r.v01), polygon_data.PP.col(r.v10), Style::outline(color_parallel, 3.5, LineType::Dash));
				}
			}
			for (auto& r : polygon_data.RE.continuations()) {
				const vec2 p0 = polygon_data.PP.col(r.v0);
				const vec2 p1 = polygon_data.PP.col(r.v1);
				draw::point(p0, .1, Style::fill(color_closure_aligned));
				draw::point(p1, .1, Style::fill(color_closure_aligned));
				draw::line(p0, p1, Style::outline(color_closure_aligned, 2.5));
#if PV_REGULARITY_CONTINUATION_STORE_CURVE
				draw::curve(reg.data.continuation.curve, 2.5, colors::forest_green);
#endif
			}			

#if 0
            for (auto& r : polygon_data.RE.vertex_symmetries()) {
                draw::line(polygon_data.PP.col(r.v0), polygon_data.PP.col(r.v1), Style::outline(color_symmetry, 4., LineType::Dash));
            }

            for (auto& r : polygon_data.RE.edge_symmetries()) {
                const vec2 p00 = polygon_data.PP.col(r.v0);
                const vec2 p01 = CircularAt(polygon_data.PP, r.v0 + 1);
                const vec2 p10 = polygon_data.PP.col(r.v1);
                const vec2 p11 = CircularAt(polygon_data.PP, r.v1 + 1);
                draw::line(.5 * (p00 + p01), .5 * (p10 + p11), Style::outline(color_symmetry, 4., LineType::Dash));
            }
#endif
			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("curves_raster") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);
			draw_curve_primitives(spline_data.curves.primitives, TangentFitColorFunctor(spline_data.tangent_fits));
			//draw_curve_primitives(spline_data.curves.primitives, AlternatingColorFunctor());

			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("curves_polygon") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);
			draw_raster(polygon_data.PP);			
			//draw_raster_polygon_indices(polygon_data.B, polygon_data.P);
			//draw_curve_primitives(spline_data.curves.primitives, TangentFitColorFunctor(spline_data.tangent_fits));
			draw_curve_primitives(spline_data.curves.primitives, AlternatingColorFunctor());
			
			draw_curves_curvature(spline_data.curves.primitives);

			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("curves_closed") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_curve_primitives_closed(spline_data.curves.primitives);
			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("corners_errors") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);

			CurveSequenceFitter fitter(polygon_data.B, polygon_data.PP, polygon_data.PV, polygon_data.RE);
			CurvePrimitiveSequence spline = fitter.all_with_type(TANGENT_FIT_LERP, true, false, false);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);
			draw_raster(polygon_data.PP);
			draw_curve_primitives(spline.primitives, AlternatingColorFunctor());
			draw_raster_indices(polygon_data.PP);
	
			pdf->draw(0, 0);
			delete pdf;
		}
		else if (arg.stage.compare("corner_convexity") == 0) {
			DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);

			std::vector<int> C;
			PathUtils::compute_convexities(polygon_data.PP, C);
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5));
			draw_raster(polygon_data.PP);

			for (size_t i = 0; i < C.size(); ++i) {
				draw::text(polygon_data.PP.col(i), to_string(C[i]), draw::font_pdf, Style::text());
			}

			pdf->draw(0, 0);
			delete pdf;
		}
        else if (arg.stage.compare("post_processor_immovable_corners") == 0) {
            DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
			draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
			draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);
            draw_raster(polygon_data.PP);

            const auto& immovable = post_processor.get_immovable_corners();
            for (size_t i = 0; i < immovable.size(); ++i) {
                if (immovable[i]) {
                    draw::point(polygon_data.B.col(i), .25, Style::fill(colors::red));
                }
            }

            draw_raster_indices(polygon_data.B);

            pdf->draw(0, 0);
            delete pdf;
        }
        else if (arg.stage.compare("post_processor_deleted_edges") == 0) {
            DevicePDF* pdf = new DevicePDF(arg.svg_uri.c_str(), 1, 1);
            const auto& E = post_processor.get_original_edges();
            const auto& E_delete = post_processor.get_permanently_deleted_edges();

            draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5));

            for (size_t e_delete : E_delete) {
                draw::line(
                    polygon_data.B.col(E[e_delete].v0), polygon_data.B.col(E[e_delete].v1),
                    Style::outline(colors::red, 2.5)
                );
            }

            pdf->draw(0, 0);
            delete pdf;
        }
		else if (arg.stage.compare("fit_attempts") == 0) {
			for (int i_attempt = 0; i_attempt < attempts.size(); ++i_attempt)
			{
				auto& attempt = attempts[i_attempt];
				for (int j = 0; j < 3; ++j)
				{
					const string uri = misc::sfmt("%s-%i-%i.svg", arg.svg_uri, i_attempt, j);
					DevicePDF* pdf = new DevicePDF(uri.c_str(), 1, 1);
					
					draw_raster_closed(polygon_data.B, real3(colors::gray * 1.9));
					draw_raster_background(polygon_data.B, Style::outline(colors::gray * 1.75, 2.5), false);					
					draw_raster(polygon_data.PP, Style::outline(colors::gray, 6.0));

					draw_curve_primitives(attempts[i_attempt].attempt[j].primitive_seq.primitives, AlternatingCornerColorFunctor(), true);
					draw_curves_curvature(attempts[i_attempt].attempt[j].primitive_seq.primitives);

					pdf->draw(0, 0);
					delete pdf;
				}

				{
					const string uri = misc::sfmt("%s-%i-classification.svg", arg.svg_uri, i_attempt);
					DevicePDF* pdf = new DevicePDF(uri.c_str(), 1, 1);

					for (int iCorner = 0; iCorner < polygon_data.PP.cols(); ++iCorner)
					{
						polyvec::real3 color = state_color_prob(attempt.state, iCorner);


						auto iPrev = (iCorner - 1 + polygon_data.PP.cols()) % polygon_data.PP.cols();
						auto iNext = (iCorner + 1 + polygon_data.PP.cols()) % polygon_data.PP.cols();

						auto corner = polygon_data.PP.col(iCorner);
						auto midPrev = 0.5 * (polygon_data.PP.col(iPrev) + polygon_data.PP.col(iCorner));
						auto midNext = 0.5 * (polygon_data.PP.col(iCorner) + polygon_data.PP.col(iNext));

						draw::line(midPrev, corner, Style::outline(color, 30.));
						draw::line(corner, midNext, Style::outline(color, 30.));

						int line = 0;
						draw::text(corner + Eigen::Vector2d(0.0, 0 * 0.5 * draw::font_pdf), misc::sfmt("%i:", iCorner), 0.5 * draw::font_pdf, Style::text());
						if (!attempt.state.corner_good_probability.empty())
							draw::text(corner + Eigen::Vector2d(0.4, (line++) * 0.5 * draw::font_pdf), misc::sfmt("Prob: %.3f", attempt.state.corner_good_probability[iCorner]), 0.5 * draw::font_pdf, Style::text());
						draw::text(corner + Eigen::Vector2d(0.4, (line++) * 0.5 * draw::font_pdf), misc::sfmt("Angle: %.3f", attempt.state.polygon_corner_angles[iCorner]), 0.5 * draw::font_pdf, Style::text());
					}

					for (auto& s : attempt.state.get_connected_sequences())
					{
						if (s.last_incl() == polygon_data.PP.cols())
							continue; //circular polygon
						auto c1 = polygon_data.PP.col(s.last_incl() % polygon_data.PP.cols());
						auto c2 = polygon_data.PP.col((s.last_incl() + 1) % polygon_data.PP.cols());
						auto mid = 0.5 * (c1 + c2);
						Eigen::VectorXd normal = (c2 - c1).normalized();
						std::swap(normal.x(), normal.y());
						normal.y() *= -1;
						draw::line(mid + normal, mid - normal, Style::outline(colors::black, 3.));
					}

					pdf->draw(0, 0);
					delete pdf;
				}
			}
		}
        else if (arg.stage.compare("paper_figure_points") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::grid(paper_figure_data);
            draw_paper_figures::raster_boundary(paper_figure_data);

            //std::cout << "B dims " << paper_figure_data.B->cols() << "x" << paper_figure_data.B->cols() << std::endl;
            for (Eigen::Index i = 0; i < paper_figure_data.B->cols(); ++i) {
                const auto pt = paper_figure_data.B->col(i);
                const auto pt_next = .5 * (pt + paper_figure_data.B->col((i + 1) % paper_figure_data.B->cols()));
                draw::point(pt, .18, Style::fill(real3(1., 0., 0.)));
                draw::point(pt_next, .18, Style::fill(real3(1., 0., 0.)));
            }
        }
        else if (arg.stage.compare("paper_figure_grid") == 0) { 
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::grid(paper_figure_data);
            draw_paper_figures::raster_boundary(paper_figure_data);
            draw_paper_figures::polygon_curves(paper_figure_data);
            draw_paper_figures::curve_primitives(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_raster") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::grid(paper_figure_data);
            draw_paper_figures::raster_boundary(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_graph") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::graph(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_polygon") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::polygon(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_polygon_curves") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::polygon_curves(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_curves_primitives_classification") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::curve_primitives_classification(paper_figure_data);
        }else if (arg.stage.compare("paper_figure_curves_primitives_optimized") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::curve_primitives(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_curves_fill_classification") == 0) {
            //io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            //draw_paper_figures::curve_fill_classification(paper_figure_data);
        }
        else if (arg.stage.compare("paper_figure_curves_fill_optimized") == 0) {
            io::SvgCanvas svg(arg.svg_uri.c_str(), canvas_width, canvas_height);
            draw_paper_figures::curve_fill(paper_figure_data);
        }
        else {
			fprintf(stderr, "Unrecognized export stage %s", arg.stage.c_str());
			continue;
		}
	}

	return EXIT_SUCCESS;
}
