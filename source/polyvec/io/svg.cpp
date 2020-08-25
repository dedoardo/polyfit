// polyvec
#include <polyvec/io/svg.hpp>
#include <polyvec/api.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/core/log.hpp>

#include <cairo-svg.h>

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(io)

SvgCanvas::SvgCanvas(const char* fname, const unsigned int width, const unsigned int height) :
    _file(fname) {
    const int canvas_scaling = 72;
    PF_DEV_F("Creating svg canvas %s %d %d", fname, width, height);
    _surface = cairo_svg_surface_create(fname, (width +1)* canvas_scaling, (height +1)* canvas_scaling);
    PF_ASSERT(_surface);
    _context = cairo_create(_surface);
    PF_ASSERT(_context);

    draw::scene_open(fname);
    _scene = &draw::scene_cur();
}

SvgCanvas::~SvgCanvas() {
    cairo_identity_matrix(_context);
    cairo_scale(_context, 72., 72.);

    for (auto& call : _scene->draw_calls) {
        call(_context, 1., Eigen::Vector2d::Constant(1.));
    }
    
    cairo_show_page(_context);
    cairo_surface_flush(_surface);
    cairo_surface_destroy(_surface);
    cairo_destroy(_context);
    draw::scene_close(_file.c_str());
}

NAMESPACE_END(io)
NAMESPACE_END(polyvec)