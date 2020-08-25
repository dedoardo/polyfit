#pragma once

// polyvec
#include <polyvec/core/macros.hpp>
#include <polyvec/debug.hpp>

// cairo
#include <cairo.h>

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(io)

class SvgCanvas {
public:
    SvgCanvas(const char* fname, const unsigned int width, const unsigned int height);
    ~SvgCanvas();

private:
    const std::string _file;
    Scene* _scene;
    cairo_t*         _context;
    cairo_surface_t* _surface;
};

NAMESPACE_END(io)
NAMESPACE_END(polyvec)