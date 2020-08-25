// Header
#include <polyvec/misc.hpp>

// libc
#include <cstdarg>

namespace polyvec {

namespace misc {
    bool is_pointer_aligned ( intptr_t v ) {
        return ( v & ( sizeof ( intptr_t ) - 1 ) ) == 0;
    }

    hash64 hash ( const real2& v ) {
        return hash_combine ( ( uint64_t& ) v.x(), ( uint64_t& ) v.y() );
    }

    // https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
    hash64 hash_combine ( hash64 lhs, hash64 rhs ) {
        lhs ^= rhs + 0x9e3779b9 + ( lhs << 6 ) + ( lhs >> 2 );
        return lhs;
    }

    // MongoDB's FNV1a implementation
    hash64 hash_bytes ( const void* buf, size_t len, uint64_t hval ) {
        const unsigned char* bp = ( const unsigned char* ) buf; /* start of buffer */
        const unsigned char* be = bp + len; /* beyond end of buffer */

        while ( bp < be ) {
            hval ^= ( uint64_t ) * bp++;
            hval += ( hval << 1 ) + ( hval << 4 ) + ( hval << 5 ) +
                    ( hval << 7 ) + ( hval << 8 ) + ( hval << 40 );
        }

        return ( hval );
    }

    bool same ( const real2& lhs, const real2& rhs ) {
        return ( lhs - rhs ).norm() < constants::SmallEps;
    }

    bool same ( double lhs, double rhs ) {
        return ::abs ( lhs - rhs ) < constants::SmallEps;
    }

    bool same_flipped ( double lhs, double rhs ) { // ehm.. less abs?
        return ::abs ( ::abs ( lhs ) - ::abs ( rhs ) ) < constants::Eps && ( sign ( lhs ) == -sign ( rhs ) );
    }

    double max0 ( double v ) {
        return ::std::max ( 0., v );
    }

    double min1 ( double v ) {
        return ::std::min ( 1., v );
    }

    double saturate ( double v ) {
        return max0 ( min1 ( v ) );
    }

    int sign ( double v ) {
        return ( 0 < v ) - ( v < 0 );
    }

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

    int4 unpack ( uint32_t xyzw ) {
        return {
            ( int ) ( ( xyzw >> 0 ) & 0xff ),
            ( int ) ( ( xyzw >> 8 ) & 0xff ),
            ( int ) ( ( xyzw >> 16 ) & 0xff ),
            ( int ) ( ( xyzw >> 24 ) & 0xff )
        };
    }

    int min_index ( const std::initializer_list<double>& vals ) {
        assert_break ( vals.size() );

        int min_idx = -1;
        double min_val = FLT_MAX;
        int i = 0;

        for ( const double val : vals ) {
            if ( val < min_val ) {
                min_val = val;
                min_idx = i;
            }

            ++i;
        }

        return min_idx;
    }

    int      max_index ( const std::initializer_list<double>& vals ) {
        assert_break ( vals.size() );

        int max_idx = -1;
        double max_val = -FLT_MAX;
        int i = 0;

        for ( const double val : vals ) {
            if ( val > max_val ) {
                max_val = val;
                max_idx = i;
            }

            ++i;
        }

        return max_idx;

    }

    Eigen::MatrixXd
normal_dist ( const Eigen::MatrixXd& in ) {
    Eigen::MatrixXd out ( in.rows(), in.cols() );
    const double sqrtpiinv = 1. / std::sqrt ( 2 * 3.14 ); // don't use for places where need accurate answer

    for ( unsigned i = 0; i < in.rows(); ++i ) {
        for ( unsigned j = 0; j < in.cols(); ++j ) {
            out ( i, j ) = sqrtpiinv * std::exp ( -in ( i, j ) * in ( i, j ) / 2. );
        }
    }

    return out;
}

Eigen::MatrixXd
randn ( unsigned int i, unsigned int j ) {
    return normal_dist ( Eigen::MatrixXd::Random ( i, j ) );
}

Eigen::MatrixXd
reshaped ( const Eigen::MatrixXd& in, int m, int n ) {
    assert ( in.size() == m*n );
    return Eigen::Map< const Eigen::MatrixXd > ( in.data(), m, n );
}


    str path_join ( const str& lhs, const str& rhs ) {
        return lhs + rhs;

        // TODO: handle `..`
#if 0
        int lhs_sep = lhs.back() == '/' || lhs.back() == '\\';
        int rhs_sep = rhs.front() == '/' || rhs.front() == '\\';

        if ( ( lhs_sep + rhs_sep ) > 0 && ( lhs_sep * rhs_sep ) == 0 ) {
            return lhs + rhs;
        }

        if ( lhs_sep + rhs_sep == 0 ) {
            return lhs + "/" + rhs;
        }

        return lhs + str ( rhs.data() + 1 );
#endif

    }

    str path_new_ext ( const str& path, const str& new_ext ) {
        str base = path.substr ( 0, path.find_last_of ( '.' ) + 1 );
        return base + new_ext;
    }

    str path_replace_ext ( const str& path, const str& new_ending ) {
        str base = path.substr ( 0, path.find_last_of ( '.' ) );
        return base + new_ending;
    }

    str sfmt ( const char* fmt, ... ) {
        va_list args, args_cp;
        va_start ( args, fmt );
        va_copy ( args_cp, args );
        int len = vsnprintf ( nullptr, 0, fmt, args );
        va_end ( args );
        str ret;
        ret.resize ( len + 1 );
        assert_break ( vsnprintf ( &ret[0], len + 1, fmt, args_cp ) );
        va_end ( args_cp );
        return ret;
    }

    str sfmt ( const real2& v ) {
        return sfmt ( "(%f %f)", v.x(), v.y() );
    }

    str sfmt ( const real3& v ) {
        return sfmt ( "(%f %f %f)", v.x(), v.y(), v.z() );
    }
}
}