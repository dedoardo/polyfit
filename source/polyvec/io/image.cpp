#include <polyvec/io/image.hpp>
#include <polyvec/core/log.hpp>

// stb
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

// libc++
#include <unordered_map>
#include <algorithm>

using namespace std;

namespace polyfit {
    namespace IO {
        uint32_t pack ( const uint8_t* xyzw ) {
            uint32_t ret = 0x0;
            ret |= ( uint32_t ) xyzw[0] << 0;
            ret |= ( uint32_t ) xyzw[1] << 8;
            ret |= ( uint32_t ) xyzw[2] << 16;
            ret |= ( uint32_t ) xyzw[3] << 24;
            return ret;
        }

        void unpack ( uint32_t xyzw, uint8_t& x, uint8_t& y, uint8_t& z, uint8_t& w ) {
            x = ( uint8_t ) ( ( xyzw >> 0 ) & 0xff );
            y = ( uint8_t ) ( ( xyzw >> 8 ) & 0xff );
            z = ( uint8_t ) ( ( xyzw >> 16 ) & 0xff );
            w = ( uint8_t ) ( ( xyzw >> 24 ) & 0xff );
        }

        vec4i unpack ( uint32_t xyzw ) {
            return {
                ( int ) ( ( xyzw >> 0 ) & 0xff ),
                ( int ) ( ( xyzw >> 8 ) & 0xff ),
                ( int ) ( ( xyzw >> 16 ) & 0xff ),
                ( int ) ( ( xyzw >> 24 ) & 0xff )
            };
        }

        void Image::reset ( uint8_t* pixels, int width, int height ) {
            if ( this->pixels ) {
                free ( this->pixels );
                this->pixels = nullptr;
            }

            this->pixels = pixels;
            this->width = width;
            this->height = height;
        }

        Image::~Image() {
            if ( pixels ) {
                free ( pixels );
                pixels = nullptr;
            }
        }

        bool read_image ( const char* uri, Image& I ) {
            new ( &I ) Image;
            int channels;
            I.pixels = stbi_load ( uri, &I.width, &I.height, &channels, 4 );
            I.channels = 4;

            if ( !I.pixels ) {
                PF_LOGF ( "Failed to open %s as %s", uri, stbi_failure_reason() );
                return false;
            }

            return true;
        }
    }
}