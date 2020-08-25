// Header
#include <polyvec/io/pdf.hpp>
#include <polyvec/api.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/utils/system.hpp>

// cairo
#include <cairo-svg.h>

namespace polyvec {
    DevicePDF::DevicePDF( const char* addr, unsigned int rows, unsigned int cols, bool svg_only)
		: svg_only(svg_only)
	{
        assert_break ( addr );

        // Open pdf
        const double aspect_ratio = ( double ) rows / cols;
        const double pdf_width_pt = 3072;
        const double pdf_height_pt = pdf_width_pt * aspect_ratio;

        _addr_pdf = addr;
        _addr_svg = addr;

        _surface = cairo_svg_surface_create ( _addr_svg.c_str(), pdf_width_pt, pdf_height_pt );
        assert_break ( _surface );

        _context = cairo_create ( _surface );
        assert_break ( _context );

        // Document bounds
        double x, y;
        cairo_clip_extents ( _context, &x, &y, &_surface_w, &_surface_h );
        assert_break ( x < constants::Eps && y < constants::Eps );

        // Redirect draw calls to _scene
        draw::scene_open ( addr );
        _scene = &draw::scene_cur();

        // Grid
        _rows = rows;
        _cols = cols;

        _addr_pdf = addr;
    }

    DevicePDF::~DevicePDF() {
        assert_break ( _surface );
        assert_break ( _context );
        assert_break ( _scene );

        cairo_show_page ( _context );
        cairo_surface_flush ( _surface );
        cairo_surface_destroy ( _surface );
        cairo_destroy ( _context );
        draw::scene_close ( _addr_pdf.c_str() );

        _context = nullptr;
        _surface = nullptr;
        _scene = nullptr;
    }

    void DevicePDF::draw ( unsigned int row, unsigned int  col, const geom::aabb& override_bounding_box) {
        assert_break ( _surface );
        assert_break ( _context );
        assert_break ( _scene );

        assert_break ( row < _rows && col < _cols );
        
		auto bbox = override_bounding_box;
		if (!bbox.valid())
			bbox = _calculate_scene_box();

        double scale = _fit_scene_to_cell ( row, col, bbox );        

        _draw ( scale, scale );
    }

    void DevicePDF::draw_in_box ( const double centerx, const double centery, const double w, const double h ) {
        assert_break ( _surface );
        assert_break ( _context );
        assert_break ( _scene );

        assert_break ( centerx >= 0 && centerx <= 1. );
        assert_break ( centery >= 0 && centery <= 1. );
        assert_break ( w > 0 && w <= 1. );
        assert_break ( h > 0 && h <= 1. );

		auto bbox = _calculate_scene_box();

        // Nothing drawn
        if ( bbox.width() < 0 || bbox.height() < 0 ) {
            return;
        }

        const double cell_w = _surface_w * w;
        const double cell_h = _surface_h * h;

        // Fitting scene box to cell without stretching
        const double scale_factor_x = cell_w / bbox.width();
        const double scale_factor_y = cell_h / bbox.height();

        const double scalex = std::min ( scale_factor_x, scale_factor_y );
        const double scaley = std::min ( scale_factor_x, scale_factor_y );

        // top-left corner of the grid
        const double translatex = _surface_w * ( centerx - .5 * h )- bbox.min.x() * scalex;
        const double translatey = _surface_h * ( centery - .5 * h )- bbox.min.y() * scaley;

        cairo_identity_matrix ( _context );
        cairo_translate ( _context, translatex, translatey );
        cairo_scale ( _context, scalex, scaley );

        _draw ( scalex, scaley );
    }

    void DevicePDF::draw_in_box_contour ( const double centerx, const double centery, const double w, const double h ) {
        const double pad = .1;
        draw_in_box ( centerx, centery, w * ( 1. -pad ), h * ( 1. -pad ) );
        //draw::box ( { 0., 0., 10., 10 }, colors::black );
        draw_in_box ( centerx, centery, w, h );
    }

	cairo_matrix_t DevicePDF::get_matrix() const
	{
		cairo_matrix_t   matrix;
		cairo_get_matrix(_context, &matrix);
		return  matrix;
	}

    void DevicePDF::_draw ( const double scalex, const double scaley ) {
        const double no_scale = 1. / std::min ( scalex, scaley );

		for (auto& call : _scene->draw_calls)
			call(_context, no_scale, real2::Constant(1.0));

        for ( const auto& t : _scene->texts ) {
            cairo_set_source_rgba ( _context, t.style.fill_color.x(), t.style.fill_color.y(), t.style.fill_color.z(), 1.0 );
            cairo_set_font_size ( _context, t.size );

            // bottom left corner -> top left corner
            cairo_text_extents_t text_box;
            cairo_text_extents ( _context, t.text.c_str(), &text_box );
            cairo_move_to ( _context, t.loc.x(), t.loc.y() + t.size );

            cairo_show_text ( _context, t.text.c_str() );
        }

        draw::clear();
    }

    geom::aabb DevicePDF::_calculate_scene_box ( ) {
		auto bbox = _scene->bounding_box;
		
		for (const auto& t : _scene->texts) {			
			cairo_set_font_size(_context, t.size);
			cairo_text_extents_t text_box;
			cairo_text_extents(_context, t.text.c_str(), &text_box);
			
			bbox.add(t.loc);
			bbox.add(t.loc + Eigen::Vector2d(text_box.width, text_box.height));
		}

		return bbox;
    }

    double DevicePDF::_fit_scene_to_cell ( unsigned int row, unsigned int col, geom::aabb bounding_box) {


		bounding_box.min -= Eigen::Vector2d::Constant(_border);
		bounding_box.max += Eigen::Vector2d::Constant(_border);

        // Nothing drawn
        if ( bounding_box.width() < 0 || bounding_box.height() < 0 ) {
            return 0;
        }

        const double cell_w = ( _surface_w ) / _cols;
        const double cell_h = ( _surface_h ) / _rows;

        // Fitting scene box to cell without stretching
        const double scale_factor_x = cell_w / bounding_box.width();
        const double scale_factor_y = cell_h / bounding_box.height();

        const double scale = std::min ( scale_factor_x, scale_factor_y );

		cairo_identity_matrix(_context);
		
		//move local origin to cell center
		cairo_translate(_context, (col + 0.5) * cell_w, (row + 0.5) * cell_h);

		//apply scale
		cairo_scale ( _context, scale, scale);

		//move center of bounding box to local origin
		auto bbox_center = 0.5 * (bounding_box.min + bounding_box.max);
		cairo_translate ( _context, -bbox_center.x(), -bbox_center.y() );		

		return scale;
    }
}
