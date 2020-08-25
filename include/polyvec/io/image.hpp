#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(IO)

struct Image {
	uint8_t* pixels = nullptr;
	int width = -1;
	int height = -1;
	int channels = -1;

	void reset(uint8_t* pixels, int width, int height);
	~Image();
};

// reads an image with any of the formats supported by stb and returns the rgba array of pixels
bool read_image(const char* uri, Image& I);

NAMESPACE_END(IO)
NAMESPACE_END(polyfit)