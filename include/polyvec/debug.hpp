#ifndef polyvec_debug_h_
#define polyvec_debug_h_

// Polyvec
#include <polyvec/api.hpp>
#include <polyvec/geom.hpp>

#include <cairo.h>

namespace polyvec {
    struct Node;
    struct Edge;
    class EdgeError;
    class Polygon;
    struct Segment;
    struct SubSegment;
    struct PolyPath;
    struct Debugging;
    struct Plottable2D;
    struct CornerGenerator;
    class PolyGraph;
    struct CurvePrimitive;
    class GlobFitCurve;
    class Spline;

    class Polygon;

//
// Logging / Debug
    namespace dbg {
        void set_targets ( FILE* info_fp = stdout, FILE* warning_fp = stderr, FILE* error_fp = stderr,
                           FILE* verbose_fp = nullptr );

        void verbose ( const char* fun, int segment, const char* fmt, ... );
        void info ( const char* fun, int segment, const char* fmt, ... );
        void warning ( const char* fun, int segment, const char* fmt, ... );
        void error ( const char* fun, int segment, const char* fmt, ... );
        void vs_debug ( const char* fun, int segment, const char* fmt, ... );

#define FMT(fmt, ...) __FUNCTION__, __LINE__, fmt, __VA_ARGS__
#define STR(str) __FUNCTION__, __LINE__, str
    }

	//
	// Drawing
    namespace colors {
        inline const real3 from_packed ( const std::uint32_t packed ) {
            return real3{ ( double ) ( ( packed >> 24 ) & 0xFF ) / 255,
                          ( double ) ( ( packed >> 16 ) & 0xFF ) / 255,
                          ( double ) ( ( packed >> 8 ) & 0xFF ) / 255 };
        }

        static const real3 white = from_packed ( 0xFFFFFFFF );
        static const real3 black = from_packed ( 0x000000FF );
        const real3 ivory = from_packed ( 0xFFFFF0FF );
        const real3 beige = from_packed ( 0xF5F5DCFF );
        const real3 wheat = from_packed ( 0xF5DeB3FF );
        const real3 tan = from_packed ( 0xD2B48CFF );
        const real3 khaki = from_packed ( 0xC3B091FF );
        const real3 silver = from_packed ( 0xC0C0C0FF );
        const real3 gray = from_packed ( 0x808080FF );
        const real3 dark_gray = from_packed ( 0x3F3F3FFF );
        const real3 charcoal = from_packed ( 0x464646FF );
        const real3 navy_blue = from_packed ( 0x000080FF );
        const real3 royal_blue = from_packed ( 0x084C9EFF );
        const real3 medium_blue = from_packed ( 0x0000CDFF );
        const real3 azure = from_packed ( 0x007FFFFF );
        const real3 cyan = from_packed ( 0x00FFFFFF );
        const real3 aquamarine = from_packed ( 0x7FFFD4FF );
        const real3 teal = from_packed ( 0x008080FF );
        const real3 forest_green = from_packed ( 0x228B22FF );
        const real3 olive = from_packed ( 0x808000FF );
        const real3 chartreuse = from_packed ( 0x7FFF00FF );
        const real3 lime = from_packed ( 0xBFFF00FF );
        const real3 golden = from_packed ( 0xFFD700FF );
        const real3 golden_rod = from_packed ( 0xDAA520FF );
        const real3 coral = from_packed ( 0xFF7F50FF );
        const real3 salmon = from_packed ( 0xFA8072FF );
        const real3 hot_pink = from_packed ( 0xFC9FC9FF );
        const real3 fuchsia = from_packed ( 0xFF77FFFF );
        const real3 puce = from_packed ( 0xCC8899FF );
        const real3 mauve = from_packed ( 0xE0B0FFFF );
        const real3 lavendere = from_packed ( 0xB57EDCFF );
        const real3 plum = from_packed ( 0x843179FF );
        const real3 indigo = from_packed ( 0x4B0082FF );
        const real3 maroon = from_packed ( 0x800000FF );
        const real3 crimson = from_packed ( 0xDC143CFF );
        const real3 talking_orange = from_packed ( 0xEF9928FF );
        const real3 shy_purple = from_packed ( 0xd042edFF );
		const real3 calm_blue = from_packed(0x5cade0FF);
		const real3 blue = from_packed ( 0x5cade0FF );
        const real3 red = from_packed ( 0xf70045FF );
        const real3 green = from_packed ( 0x37e83aFF );
        const real3 white_whale = from_packed ( 0xdee9fcff );
        const real3 turtle_purple = from_packed ( 0xef94e3ff );
        const real3 enemys_blood = from_packed ( 0xf70045FF );

        inline const real3 random() {
            return real3 ( .2 + .6 * ( double ) rand() / RAND_MAX, .2 + .6 * ( double ) rand() / RAND_MAX,
                           .2 + .6 * ( float ) rand() / RAND_MAX );
        }

        inline const real4 random4() {
            return real4(.2 + .6 * (double)rand() / RAND_MAX, .2 + .6 * (double)rand() / RAND_MAX,
                .2 + .6 * (float)rand() / RAND_MAX, 1.0);
        }
    }

    enum class LineType {
        Solid,
        Dash,
    };

    struct Style {
        double   line_thickness;
        real3    line_color;
        LineType line_type;
        real4    fill_color;

        static Style text ( const real3& color = colors::black );
        static Style fill ( const real4& color );
        static Style fill (const real3& color);
        static Style outline ( const real3& color, double thickness, LineType line_type = LineType::Solid );
    };

    struct Scene {
        struct Text {
            real2 loc;
            str    text;
            double size;
            Style  style;
        };        

        str          name;
        std::vector<Text>    texts;
		std::vector<std::function<void(cairo_t* context, double thickness_scale, const real2& scale)>> draw_calls;
		geom::aabb bounding_box;

		void begin_closed_shape(const Style& style);
		void close_shape();
		void begin_outline(const Style& style);
		void finish_outline();
        void polygon ( const std::vector<real2>& pts, const Style& style );
        void polyline(const Eigen::Matrix2Xd& pts, const Style& style, bool is_closed);
        void text ( const real2& l, const str& text, double size, const Style& style );
        void line ( const real2& src, const real2& dst, const Style& style );
        void point ( const real2& l, double radius, const Style& style );
    };

    namespace draw {
        static const real2 O = real2 ( 0., 0. );
        static const double font_pdf    = .45;
        static const double font_tiny   = 3.5;
        static const double font_small  = 1.;
        static const double font_normal = 12.;
        static const double font_big    = 16.;
        static const double font_title  = 18;

        Scene* scene ( const str& name = "" );
        Scene& scene_cur ();
        void   clear();
        void   scene_open ( const str& name );
        void   scene_close ( const str& name );

        // Primitives
        void polygon ( const std::vector<real2>& pts, const Style& style );
        void polyline(const Eigen::Matrix2Xd& pts, const Style& style, bool is_closed = true);
        void text ( const real2& l, const str& text, double size, const Style& style );
        void line ( const real2& src, const real2& dst, const Style& style );
        void point ( const real2& l, double radius, const Style& style );
		void circle(const real2& c, double radius, const Style& style);
		void curve(GlobFitCurve* curve, double thickness, const real3& bez_color = colors::turtle_purple);
		void curves_fill(const std::vector<GlobFitCurve*>& curves, const real4& color);
    }

}

#endif // _polyvec_debug_h_