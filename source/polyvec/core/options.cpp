// polyvec
#include <polyvec/core/options.hpp>
#include <polyvec/core/macros.hpp>

#include <memory>

using namespace std;

namespace polyvec {

// Global options are set at startups, local-threads can be se and subsequent calls
// to VectorOptions::get() will return the local one.
namespace {
    unique_ptr<VectorOptions>              options_global;
    thread_local unique_ptr<VectorOptions> options_local;
}

VectorOptions* VectorOptions::get() {
    if ( options_local ) {
        return options_local.get();
    }

    PF_ASSERT( options_global );
    return options_global.get();
}

void VectorOptions::make_global() const {
    options_global = unique_ptr<VectorOptions> ( new VectorOptions ( *this ) );
}

void VectorOptions::make_thread_local() const {
    options_local = unique_ptr<VectorOptions> ( new VectorOptions ( *this ) );
}
}