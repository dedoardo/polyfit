#pragma once

#define PV_INLINE inline

#define PV_LEN(x) (sizeof(x)/sizeof((x)[0]))

// Prevent astyle from indenting namespace; while not touching already written files
#define NAMESPACE_BEGIN(NAME__) namespace NAME__ {
#define NAMESPACE_END(NAME__)                    }

#include <exception>

#ifndef _DEBUG
#define PF_ABORT throw std::exception()
#else
#include <intrin.h>
#define PF_ABORT __debugbreak()
#endif
#define PF_ASSERT(expr) { if (!(expr)) { fprintf(stderr, "%s: %d - Failed expression: " #expr, __FILE__, __LINE__); PF_ABORT; }}

#define PF_EPS 1e-4
#define PF_EPS_MEDIUM 1e-3
#define PF_EPS_LARGE 1e-1

#define PF_ARRAY_LEN(x) (sizeof(x)/sizeof((x)[0]))

#define PV_UNUSED(x) (void)(x)
