#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/io/image.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(ImageSegment)

void extract_closed_regions (const IO::Image& I, std::vector<mat2x>& boundaries, std::vector<vec4>& colors);

// miscellaneous cleanup function that expands the image by 1 pixel in every direction guessing a
// reasonable color for the background. It also remove small individual pixels with non-matching alpha.
// this should really not be called on pixel art inputs
void expand_and_cleanup(IO::Image& I);

uint32_t find_background_color(IO::Image& I);

int find_binary_color_region(const std::vector<mat2x>& boundaries);

NAMESPACE_END(ImageSegment)
NAMESPACE_END(polyfit)