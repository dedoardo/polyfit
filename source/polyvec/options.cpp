// polyvec
#include <polyvec/options.hpp>

NAMESPACE_BEGIN(polyfit)

Options* Options::get() {
    static Options options_default;
    return &options_default;
}

NAMESPACE_END(polyfit)