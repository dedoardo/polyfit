/*
	This file contains the primitive type definition, macros and includes shared by other project files
*/
#ifndef polyvec_api_h_
#define polyvec_api_h_

// libc
#define _USE_MATH_DEFINES
#include <math.h>

// C++ STL
#include <stdint.h>
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <limits>

#ifndef FLT_MAX
#define FLT_MAX std::numeric_limits<float>::max()
#endif

// Eigen
#include <Eigen/Core>

// Platform specific macros
#ifdef _MSC_VER

#define break_here throw std::exception()
#define polyvec_inline __forceinline

#define polyvec_win32 1
#define polyvec_posix 0

#else
#define break_here do{assert(0); throw std::runtime_error("error");}while(0)
#define polyvec_inline __attribute__((always_inline)) inline

#define polyvec_win32 0
#define polyvec_posix 1

#endif

#ifndef break_here
#define throw std::exception()
#endif

// Generic macros
#ifndef param_unused
#define param_unused(x) (void)(x)
#endif

#define result_unused(x) (void)(x)

#ifndef assert_break
#define assert_break(x) { if (!(x)) { fprintf(stderr, "\nFAILED ASSERT " #x  "\n      LOCATION %s:%d", __FILE__, __LINE__); break_here; } }
#define assert_break_msg(x, msg) { if (!(x)) { fprintf(stderr, "\nFAILED ASSERT %s"  "\n      LOCATION %s:%d", msg,__FILE__, __LINE__); break_here; } }
#endif

#define array_len(x) (sizeof(x)/sizeof((x)[0]))
#ifndef logic_xor // https://benpfaff.org/writings/clc/logical-xor.html
#define logic_xor(a, b) ((!(a)) != (!(b)))
#endif

#define polyvec_str(__X)   static_cast<std::ostringstream&>(std::ostringstream().flush() << __X).str()
#define polyvec_c_str(__X) polyvec_str(__X).c_str()

//
// Type aliases
namespace polyvec {

    using index = std::int64_t;
    using indexu = std::size_t;


    using str = std::string;
    using hash64 = uint64_t; // let's not confuse this with the other bunch of integers

    template <typename T>
    using ptr = std::unique_ptr<T>;

    template <typename K, typename V>
    using umap = std::unordered_map<K, V>;

    template <typename V>
    using uset = std::unordered_set<V>;

// Matrix aliases
    using byte4 = Eigen::Matrix<uint8_t, 4, 1>;
    using real = double;
    using int2 = Eigen::Vector2i;
    using int3 = Eigen::Vector3i;
    using int4 = Eigen::Vector4i;
    using real2 = Eigen::Vector2d;
    using real4 = Eigen::Vector4d;
    using mat2 = Eigen::Matrix2d;
    using mat3 = Eigen::Matrix3d;
    using mat4 = Eigen::Matrix4d;
    using real2xN = Eigen::Matrix2Xd;
    using real1xN = Eigen::VectorXd;

// Matrix Constants
    extern const index index_null;
    extern const hash64 hash64_null;
    extern const index index_max;
    extern const int   int_max;
    extern const uint32_t uint32_max;
    extern const int32_t int32_max;
    extern const int64_t int64_max;
    extern const double real_max;
    extern const real2 real2_0;
    extern const real2 real2_1;
    extern const real2 real2_lo; // lowest
    extern const real2 real2_hi; // max
#define            real2_xy(xy) (real2((xy), (xy))) // avoid evaluating twice?

    using real3 = Eigen::Vector3d;
    extern const real3 real3_0;
    extern const real3 real3_1;
    extern const real3 real3_lo; // lowest
    extern const real3 real3_hi; // max

//
// Constants
    namespace constants {
        static constexpr double PI          = M_PI;
        static constexpr double PI_inv      = M_1_PI;
        static constexpr double PI_half     = M_PI_2;
        static constexpr double PI_half_inv = 1. / PI_half;
        static constexpr double PI2         = PI * 2;
        static constexpr double PI3_4       = 3 * M_PI_4;
        static constexpr double PI_4        = M_PI_4;
        static constexpr double PIDiv180    = 0.017453292519943295769236907684886;

        static constexpr double SmallEps = 1e-4;
        static constexpr double BigEps   = 1e-07;
        static constexpr double Eps      = 1e-10;
    }



//
// Indexing
// The code works with stricly typed indices that do not allow implicit conversions to other types.
// Warnings should be set to level4/Wall and stop compilation.
//
// Add any overload when required, *keep* in mind that to mantain
// strict typing:
// 1. IF conversion operator to `index`        NO non-explicit constructor from `index`
// 2. IF non-explicit constructor from `index` NO  conversion operator to `index`
// (1) is better as allows having to overload all integer operations
//
// DO NOT COMPARE VALIDITY with -1 as it's a perfectly valid index
//
#define define_typed_index(t, i)\
    struct t{\
        const static t null;\
        i _v;\
        explicit t(i v = index_null) : _v(v) { }\
        operator std::int64_t()const { return _v; }\
        explicit operator bool()const { return _v != null; }\
        t& operator=(i v) { _v = v; return *this; }\
        t& operator=(int v) { _v = (i)v; return *this; }\
        bool operator <(i rhs) const { return (i)_v < rhs;}\
        t& operator++() { ++_v; return *this; }\
        t& operator--() { --_v; return *this; }\
    };\
    inline t operator-(t lhs, t rhs) { return (t)(lhs._v - rhs._v); }\
    inline t operator+(t lhs, t rhs) { return (t)(lhs._v + rhs._v); }\
    inline t operator-(t lhs, int rhs) { return (t)(lhs._v - rhs); }\
    inline t operator+(t lhs, int rhs) { return (t)(lhs._v + rhs); }\
    inline t operator%(t lhs, t rhs) { return (t)(lhs._v % rhs._v); }\
    static_assert(sizeof(i) == sizeof(t), "Replace compiler with less edgy cousin");

#define define_typed_index_hash(t)\
    struct _std_hash_index_##t {\
        hash64 operator() (t##_idx i) const {\
            return (hash64)i._v;\
        }\
    };\
    struct _std_hash_index_pair_##t {\
        hash64 operator() (const std::pair<t##_idx, t##_idx>& p) const {\
            return hash_index_pair(p.first._v, p.second._v);\
        }\
    };

    hash64 hash_index_pair ( index lhs, index rhs );
    hash64 hash_index_tuple ( index e0, index e1, index e2 );

// Since there are many different buffers around and they also may be circular or not
// each index has different types to avoid mistakes (naming helps too). It can be used
// as a normal index, but won't be implicitly converted to other types
    define_typed_index ( vertex_idx, index );
    define_typed_index ( corner_idx, index );
    define_typed_index ( node_idx, index );
    define_typed_index ( edge_idx, index );

    define_typed_index_hash ( vertex );
    define_typed_index_hash ( corner );
    define_typed_index_hash ( node );
    define_typed_index_hash ( edge );

    struct _std_hash_node {
        hash64 operator() ( const std::tuple<corner_idx, corner_idx, corner_idx>& t ) const {
            return hash_index_tuple ( std::get<0> ( t )._v, std::get<1> ( t )._v, std::get<2> ( t )._v );
        }
    };

}

// Prevent astyle from indenting namespace; while not touching already written files
#define NAMESPACE_BEGIN(NAME__) namespace NAME__ {
#define NAMESPACE_END(NAME__)                    }

#endif // polyvec_api_h_