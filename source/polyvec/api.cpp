// Header
#include <polyvec/api.hpp>

// C++ STL
#include <limits>

// Polyvec
#include <polyvec/misc.hpp>
#include <polyvec/geom.hpp>
#include <polyvec/debug.hpp>

using namespace std;

namespace polyvec {

// Constants
const hash64 hash64_null = 0xffffffffffffffff;
const index index_max = numeric_limits<index>::max();
const int int_max = numeric_limits<int>::max();
const uint32_t uint32_max = numeric_limits<uint32_t>::max();
const int32_t int32_max = numeric_limits<int32_t>::max();
const int64_t int64_max = numeric_limits<int64_t>::max();
const double real_max = numeric_limits<double>::max();

const index  index_null = 0xefffffffffffffff;
const real2 real2_0 = real2::Zero();
const real2 real2_1 = real2::Ones();
const real2 real2_lo = real2 ( numeric_limits<double>::lowest(), numeric_limits<double>::lowest() );
const real2 real2_hi = real2 ( numeric_limits<double>::max(), numeric_limits<double>::max() );

const real3 real3_0 = real3::Zero();
const real3 real3_1 = real3::Ones();
const real3 real3_lo = real3 ( numeric_limits<double>::lowest(), numeric_limits<double>::lowest(), numeric_limits<double>::lowest() );
const real3 real3_hi = real3 ( numeric_limits<double>::max(), numeric_limits<double>::max(), numeric_limits<double>::max() );


}