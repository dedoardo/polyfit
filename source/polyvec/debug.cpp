// Header
#include <polyvec/debug.hpp>

// C++ STL
#include <cstdarg>
#include <random>

// Polyvec
#include <polyvec/utils/matrix.hpp>
#include <polyvec/utils/system.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/curve-tracer/spline.hpp>

#define _STYLE const ::Style& style
#define _LABEL const str& _name
#define _STR   const str&
#define _F2    const real2&
#define _F3    const real3&

using namespace std;
using namespace polyfit;

namespace polyvec {
    namespace {
        FILE* info_fp = stdout;
        FILE* warning_fp = stderr;
        FILE* error_fp = stderr;
        FILE* verbose_fp = stdout;
        FILE* dbg_fp = nullptr;

        void print ( FILE* fp, const char* tag, const char* fun, int segment, const char* fmt, va_list args ) {
            str mod = misc::sfmt ( "%-14s", "" );
            str fun_s = str ( fun );
            auto mod_it = fun_s.find_first_of ( "::" );

            if ( mod_it != fun_s.npos ) {
                mod = misc::sfmt ( " :%-12s", str ( fun, mod_it ).c_str() );
            }

            str routine = fun;
            auto routine_it = fun_s.find_last_of ( "::" );

            if ( routine_it != fun_s.npos ) {
                routine = str ( fun + routine_it );
            }

#if 0

            if ( strcmp ( tag, "Info" ) == 0 ) {
                os::text_bright_green();
            }

            if ( strcmp ( tag, "Verb" ) == 0 ) {
                os::text_green();
            }

            if ( strcmp ( tag, "Warn" ) == 0 ) {
                os::text_orange();
            }

            if ( strcmp ( tag, "Err " ) == 0 ) {
                os::text_red();
            }

#endif

            if ( fp != nullptr ) {
                str new_fmt = misc::sfmt ( "%s%s> %s @ %s:%d\n", tag, mod.c_str(), fmt, routine.c_str(), segment );
                vfprintf ( fp, new_fmt.c_str(), args );
            }

            // os::text_reset();
        }
    }

    namespace dbg {
        void set_targets ( FILE* _info_fp, FILE* _warning_fp, FILE* _error_fp, FILE* _verbose_fp ) {
            info_fp = _info_fp;
            warning_fp = _warning_fp;
            error_fp = _error_fp;
            verbose_fp = _verbose_fp;
        }

#define GEN_LOG(name, tag)\
    void name (const char* fun, int segment, const char* fmt, ...) {\
        va_list args; va_start(args, fmt); \
        print(name##_fp, tag, fun, segment, fmt, args);\
        va_end(args); \
    }

        GEN_LOG ( verbose, "Verb" );
        GEN_LOG ( info, "Info" );
        GEN_LOG ( warning, "Warn" );
        GEN_LOG ( error, "Err " );
    }

    Style Style::text ( const real3& color ) {
        Style s;
        s.line_type = LineType::Solid;
        s.line_thickness = 1.;
        s.line_color = color;
        s.fill_color = real4(color(0), color(1), color(2), 1.0);
        return s;
    }

    Style Style::fill(const real3& color) {
        return fill(real4(color(0), color(1), color(2), 1.0));
    }

    Style Style::fill ( const real4& color ) {
        Style s;
        memset ( &s, 0, sizeof ( s ) );
        s.line_type = LineType::Solid;
        s.line_thickness = 0.;
        s.line_color = color.segment(0, 3);
        s.fill_color = color;
        return s;
    }

    Style Style::outline ( const real3& color, double thickness, LineType line_type ) {
        Style s;
        memset ( &s, 0, sizeof ( s ) );
        s.line_type = line_type;
        s.line_thickness = thickness;
        s.line_color = color;
        s.fill_color = real4(color(0), color(1), color(2), 1.0);
        return s;
    }

	void Scene::begin_closed_shape(const Style & style)
	{
		draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_set_source_rgba(ctx, style.fill_color.x(), style.fill_color.y(), style.fill_color.z(), style.fill_color.w()); });
	}

	void Scene::close_shape()
	{
		draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_close_path(ctx); });
		draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_fill(ctx); });
	}

	void Scene::begin_outline(const Style & style)
	{
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_set_source_rgba(_context, style.line_color.x(), style.line_color.y(), style.line_color.z(), 1.0); });
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_set_line_width(_context, style.line_thickness * thickness_scale); });

        if (style.line_type == LineType::Dash) {
            draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { double dash = thickness_scale * 25;  cairo_set_dash(_context, &dash, 1, 0.); });
        }
	}

	void Scene::finish_outline()
	{
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_stroke(_context); });
        draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_set_dash(_context, 0, 0, 0.);  });
	}

	void Scene::polygon ( const std::vector<real2>& pts, const Style& style ) {
		begin_closed_shape(style);

		draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_move_to(ctx, pts[0].x() * scale.x(), pts[0].y() * scale.y()); });
		bounding_box.add(pts[0]);

		for (index i = 0; i < (int)pts.size(); ++i) {
			draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_line_to(ctx, pts[(i + 1) % pts.size()].x() * scale.x(), pts[(i + 1) % pts.size()].y() * scale.y()); });
			bounding_box.add(pts[i]);
		}

		close_shape();
    }

    void Scene::polyline(const Eigen::Matrix2Xd& pts, const Style& style, bool is_closed) {
		begin_outline(style);
        
        draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_move_to(ctx, pts.col(0).x() * scale.x(), pts.col(0).y() * scale.y()); });        
        bounding_box.add(pts.col(0));

        for (index i = 0; i < (int)pts.cols(); ++i) {
            if (!is_closed && i == pts.cols() - 1) {
                continue;
            }
            draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_line_to(ctx, CircularAt(pts, i + 1)(0) * scale.x(), CircularAt(pts, i + 1)(1) * scale.y()); });
            bounding_box.add(pts.col(i));
        }

        if (is_closed) {
            draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_close_path(ctx); });
        }
		finish_outline();
    }

	void Scene::text ( const real2& l, const str& text, double size, const Style& style ) {
        texts.push_back ( { l, text, size, style } );
    }

    void Scene::line ( const real2& src, const real2& dst, const Style& style ) {
		begin_outline(style);

		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_move_to ( _context, src.x(), src.y() ); });
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_line_to ( _context, dst.x(), dst.y() ); });
		
		finish_outline();

		bounding_box.add(src);
		bounding_box.add(dst);
    }

    void Scene::point ( const real2& l, double radius, const Style& style ) {
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_set_source_rgba(_context, style.fill_color.x(), style.fill_color.y(), style.fill_color.z(), 1.0); });
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_arc ( _context, l.x(), l.y(), radius, 0, constants::PI2 ); });
		draw_calls.push_back([=](cairo_t* _context, double thickness_scale, const real2& scale) { cairo_fill ( _context ); });
		bounding_box.add(l);
    }

    namespace {
        unordered_map<str, ptr<Scene>> scenes;
        Scene*                         scene_active = nullptr;
    }

    namespace draw {
        Scene* scene ( const str& name ) {
            return scenes[name].get();
        }

        Scene& scene_cur() {
            assert_break ( scene_active );
            return *scene_active;
        }

        void clear() {
            assert_break ( scene_active );
			new (&scene_active->bounding_box) geom::aabb;
            scene_active->draw_calls.clear();
            scene_active->texts.clear();
        }

        void scene_open ( const str& name ) {
            if ( !misc::container_in ( scenes, name ) ) {
                scenes.insert ( make_pair ( name, ptr<Scene> ( new Scene ) ) );
            }

            scene_active = scenes[name].get();
        }

        void scene_close ( const str& name ) {
            assert_break ( misc::container_in ( scenes, name ) );
            scenes.erase ( scenes.find ( name ) );

            if ( scenes.empty() ) {
                scene_active = nullptr;
            } else {
                scene_active = scenes[""].get();
            }
        }

        // Primitives
        void polygon ( const std::vector<real2>& pts, const Style& style ) {
            scene_cur().polygon ( pts, style );
        }

        void polyline(const Eigen::Matrix2Xd& pts, const Style& style, bool is_closed) {
            scene_cur().polyline(pts, style, is_closed);
        }

        void text ( const real2& l, const str& text, double size, const Style& style ) {
            scene_cur().text ( l, text, size, style );
        }

        void line ( const real2& src, const real2& dst, const Style& style ) {
            scene_cur().line ( src, dst, style );
        }

        void point ( const real2& l, double radius, const Style& style ) {
            scene_cur().point ( l, radius, style );
        }

		void circle(const real2& c, double r, const Style& style) {
			// ehm, really a cosine?
			double a = 0;
			const double a_step = .05;
			const real2 p0 = c + real2(cos(a), sin(a)) * r;
			real2 p = p0;
			a += a_step;

			while (a < 2 * M_PI + PF_EPS) {
				const real2 pn = c + real2(cos(a), sin(a)) * r;
				draw::line(p, pn, style);
				a += a_step;
				p = pn;
			}

			draw::line(p, p0, style);
		}

		void add_draw_commands(GlobFitCurve* curve, polyvec::Scene* scene, bool initial_move)
		{
			auto bbox = curve->get_bounding_box();
			scene->bounding_box.add(bbox.min);
			scene->bounding_box.add(bbox.max);
			auto line = dynamic_cast<GlobFitCurve_Line*>(curve);
			if(line)
			{
				auto start = line->pos(Curve::START);
				auto end = line->pos(Curve::END);
				if(initial_move)
					scene->draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_move_to(ctx, start.x(), start.y()); });
				scene->draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_line_to(ctx, end.x(), end.y()); });				
				return;
			}

			auto bezier = dynamic_cast<BezierCurve*>(curve);
			if (bezier)
			{
				auto controls = bezier->get_control_points();
				if (initial_move)
					scene->draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_move_to(ctx, controls(0, 0), controls(1, 0)); });
				scene->draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_curve_to(ctx, controls(0, 1), controls(1, 1), controls(0, 2), controls(1, 2), controls(0, 3), controls(1, 3)); });
				return;
			}

			//generic implementation
			const double t_step = .025;
			double t = t_step;

			auto start = curve->pos(Curve::START);
			if (initial_move)
				scene->draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_move_to(ctx, start.x(), start.y()); });

			while (t < Curve::END + constants::Eps) {
				auto p = curve->pos(misc::saturate(t));
				scene->draw_calls.push_back([=](cairo_t* ctx, double thickness_scale, const real2& scale) { cairo_line_to(ctx, p.x(), p.y()); });
				t += t_step;
			}
		}

		void curve(GlobFitCurve* curve, double thickness, const real3& bez_color) {
			auto& scene = draw::scene_cur();
			scene.begin_outline(Style::outline(bez_color, thickness));
			add_draw_commands(curve, &scene, true);
			scene.finish_outline();
		}

		void curves_fill(const std::vector<GlobFitCurve*>& curves, const real4& _color) {
            auto color = _color;
            if (color(3) < PF_EPS) {
                color = Eigen::Vector4d::Ones();
            }

            draw::scene_cur().begin_closed_shape(Style::fill(color));

			for (size_t icurve = 0; icurve < curves.size(); ++icurve) {
				GlobFitCurve* curve = curves[icurve];

				add_draw_commands(curve, &draw::scene_cur(), icurve == 0);
			}

			draw::scene_cur().close_shape();
		}
    }

#ifdef _MSC_VER
#include <Windows.h>
#endif

    namespace dbg {
        void vs_debug ( const char* fun, int segment, const char* fmt, ... ) {
#ifdef _MSC_VER
            va_list args, args_cp;
            va_start ( args, fmt );
            va_copy ( args_cp, args );
            size_t len = vsnprintf ( nullptr, 0, fmt, args );
            str msg;
            msg.resize ( len );
            va_end ( args );
            vsnprintf ( &msg[0], len, fmt, args_cp );
            msg = misc::sfmt ( "%s %s:%d\n", msg.c_str(), fun, segment );
            OutputDebugStringA ( msg.c_str() );
            va_end ( args );
#endif
        }
    }
}