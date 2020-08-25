#ifndef polyvec_misc_h_
#define polyvec_misc_h_

// C++ STL
#include <utility>
#include <type_traits>
#include <map>
#include <iomanip>
#include <unordered_map>
#include <unordered_set>
#include <fstream>

// Polyvec
#include <polyvec/api.hpp>

namespace polyvec {

/*
    Behold adventurers, you are about to enter the realm of unnecessary string copies.

    Only misplaced trust in the compiler's string pool will leave your eyes and soul unharmed.

    Try to look for comfort in its expressiveness..
*/

namespace misc {
    //
    // Misc
    bool     is_pointer_aligned ( intptr_t v );
    hash64   hash               ( const real2& v );
    hash64   hash_combine       ( hash64 lhs, hash64 rhs ); // Let's not xor and spend days debugging a random collision..
    hash64   hash_bytes         ( const void* buf, size_t len, uint64_t init = ( ( uint64_t ) 0xcbf29ce484222325ULL ) );
    bool     same               ( const real2& lhs, const real2& rhs );
    bool     same               ( double lhs, double rhs );
    bool     same_flipped       ( double lhs, double rhs );
    double   max0               ( double v );
    double   min1               ( double v );
    double   saturate           ( double v );
    int      sign               ( double v );
    uint32_t pack               ( const uint8_t* xyzw );
    void     unpack             ( uint32_t xyzw, uint8_t& x, uint8_t& y, uint8_t& z, uint8_t& w );
    int4     unpack             ( uint32_t xyzw );
    int      min_index          ( const std::initializer_list<double>& vals );
    int      max_index          ( const std::initializer_list<double>& vals );


//
// numeric functions
//
Eigen::MatrixXd
normal_dist ( const Eigen::MatrixXd& in );

Eigen::MatrixXd
randn ( unsigned int i, unsigned int j );

Eigen::MatrixXd
reshaped ( const Eigen::MatrixXd& in, int m, int n );



    template <typename T>
    constexpr T lerp ( const T lo, const T hi, const double t ) {
        return lo + ( hi - lo ) * t;
    }

    // path
    str path_join ( const str& lhs, const str& rhs );
    str path_new_ext ( const str& path, const str& new_ext );
    str path_replace_ext ( const str& path, const str& new_ending );

    // high overhead formatting functions
    str sfmt ( const char* fmt, ... );
    str sfmt ( const real2& v );
    str sfmt ( const real3& v );

    template <typename T>
    str sfmt ( const std::vector<T>& vec ) {
        if ( vec.empty() ) {
            return "";
        }

        std::stringstream ss;

        for ( size_t i = 0; i < vec.size() - 1; ++i ) {
            ss << std::setprecision ( 6 ) << vec[i] << " ";
        }

        ss << vec.back();
        return ss.str();
    }

    inline str sfmt ( const double* vals, index n ) {
        std::stringstream ss;

        for ( index i = 0; i < n - 1; ++i ) {
            ss << std::setprecision ( 6 ) << vals[i] << " ";
        }

        ss << vals[n - 1];
        return ss.str();
    }


    // Validate typed enumerator ( from ::None to ::Count)
    template <typename EnumT>
    bool enum_valid ( const EnumT& val ) {
        using EnumBaseT = typename std::underlying_type<EnumT>::type;
        static_assert ( std::is_integral<EnumBaseT>::value, "" );
        constexpr std::int64_t MinEnumV = ( std::int64_t ) EnumT::None;
        constexpr std::int64_t MaxEnumV = ( std::int64_t ) EnumT::Count;
        return ( int64_t ) val > MinEnumV && ( int64_t ) val < MaxEnumV;
    }

    template <typename EnumT>
    EnumT enum_match ( const std::string& name, std::initializer_list<std::pair<EnumT, const char*>> matches ) {
        for ( const auto& match : matches ) {
            if ( strcmp ( name.c_str(), match.second ) == 0 ) {
                return match.first;
            }
        }

        return EnumT::None;
    }


    // load/store vectors to file
    template <typename T>
    T read_element ( const std::string& str ) {
        return ( T ) std::stoll ( str );
    }

    template <>
    inline std::string read_element ( const std::string& str ) {
        return str;
    }

    template <typename T>
    void write_vec ( const std::string& path, const std::vector<T>& vec ) {
        std::ofstream fs ( path );

        if ( !fs.good() ) {
            return;
        }

        for ( const T& el : vec ) {
            fs << el;
        }
    }

    template <typename V, typename ...Vs>
    typename std::map<Vs...>::const_iterator find ( const std::map<Vs...>& c, const V& val ) {
        return c.find ( val );
    }

    template <typename V, typename ...Vs>
    typename std::unordered_map<Vs...>::const_iterator find ( const std::unordered_map<Vs...>& c, const V& val ) {
        return c.find ( val );
    }

    template <typename V, typename ...Vs>
    typename std::unordered_set<Vs...>::const_iterator find ( const std::unordered_set<Vs...>& c, const V& val ) {
        return c.find ( val );
    }

    template <typename V, typename ...Vs>
    typename std::vector<Vs...>::const_iterator find ( const std::vector<Vs...>& c, const V& val ) {
        return std::find ( c.begin(), c.end(), val );
    }

    // todo: remove double lookup
    // find/in/diff stl containers
    template <typename C, typename V>
    bool container_in ( const C& c, const V& val ) {
        return misc::find ( c, val ) != c.end();
    }
}

struct hashpair {
    template <typename U, typename V>
    inline size_t operator() ( const std::pair<U, V>& x ) const {
        return misc::hash_combine ( misc::hash_bytes ( &x.first, sizeof ( U ) ), misc::hash_bytes ( &x.second, sizeof ( V ) ) );
    }
};

}

#endif // polyvec_misc_h_
