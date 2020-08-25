#ifndef _polyvec_pdf_h_
#define _polyvec_pdf_h_

// polyvec
#include <polyvec/debug.hpp>

// cairo
#include <cairo.h>

namespace polyvec {
    // Layouts draw::* calls into a grid and writes them to pdf
    class DevicePDF {
    public:
        // Opens a pdf for writing, successive draw::* calls will be redirected to the file
        DevicePDF ( const char* pdf_addr, unsigned int rows, unsigned int cols, bool svg_only = false );

        // Writes and closes the document
        ~DevicePDF();

        // Draws adjusting elements in the canvas to fit the specific grid cell
        void draw ( unsigned int row, unsigned int col, const geom::aabb& override_bounding_box = geom::aabb());

        // Draws the element in the specified box in global normalized coordinates
        void draw_in_box ( const double centerx, const double centery, const double w, const double h );
        void draw_in_box_contour ( const double centerx, const double centery, const double w, const double h );

		//Returns the model transform matrix used in the last call to draw*().
		cairo_matrix_t get_matrix() const;

    private:
        void _draw ( const double scalex, const double scaley );
        geom::aabb _calculate_scene_box ();
        double _fit_scene_to_cell ( unsigned int row, unsigned int col, geom::aabb bounding_box);

    private:
        cairo_t*         _context;
        cairo_surface_t* _surface;		

        Scene*  _scene;
        str     _addr_svg;
        str     _addr_pdf;

        const double _border = 0.5;
        double       _surface_w;
        double       _surface_h;
        unsigned int _rows;
        unsigned int _cols;

		bool svg_only;
    };
}

#endif // _polyvec_pdf_h_