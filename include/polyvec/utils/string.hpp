#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>

#include <string>

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(StringUtils)

// sprintfs into a c++ string
std::string fmt(const char* fmt, ...);

// replaces the occurrences of seq in s with the variable number of arguments
// without exceeding the length of the token. pads the tail with spaces and clamps the formatted output
int replace_bounded(char* s, const char* seq, const char* fmt, ...);

// concatenates the two addresses making sure that only one path separator is present
// for the directory + path an extra optional extension argument can be checked which
// prevents duplicate '.'
std::string join_path(const char* p0, const char* p1, const char* ext = nullptr);

NAMESPACE_END(StringUtils)
NAMESPACE_END(polyfit)

