#pragma once

// std
#include <vector>
#include <memory>

// Potrace
#if defined(POLYVEC_USE_POTRACE)
#include <potrace/potracelib.h>
#else
struct potrace_path_t {};
struct potrace_curve_t {};
struct potrace_param_t {};
inline void potrace_path_free ( potrace_path_t* ) {}
inline void potrace_param_free ( potrace_param_t* ) {}
// Deprecated
struct potrace_state_t {};
inline void potrace_state_free ( potrace_state_t* ) {}
#endif

// Polyvec
#include <polyvec/api.hpp>


NAMESPACE_BEGIN ( polyfit )
NAMESPACE_BEGIN ( PotraceUtil )

// Unique pointers for potrace types
#define DEFINE_POTRACE_UNIQUE_PTR(CPP_TYPE_NAME, POTRACE_TYPE_NAME, DELETER_FUNC) \
    class CPP_TYPE_NAME: \
    public  std::unique_ptr<POTRACE_TYPE_NAME, void ( * ) ( POTRACE_TYPE_NAME* )> { \
    public: \
        CPP_TYPE_NAME(POTRACE_TYPE_NAME *in=nullptr): \
         std::unique_ptr<POTRACE_TYPE_NAME, void ( * ) ( POTRACE_TYPE_NAME* )> (in, DELETER_FUNC){} \
    }

// Define the types
// param_unique_ptr, path_unique_ptr, and state_unique_ptr
// to hold the types potrace_param_t, potrace_path_t, and potrace_state_t,
// respectively.
DEFINE_POTRACE_UNIQUE_PTR ( param_unique_ptr, potrace_param_t, potrace_param_free );
DEFINE_POTRACE_UNIQUE_PTR ( path_unique_ptr, potrace_path_t, potrace_path_free );
DEFINE_POTRACE_UNIQUE_PTR ( state_unique_ptr, potrace_state_t, potrace_state_free );

#undef DEFINE_POTRACE_UNIQUE_PTR

// Convert a potrace path to a set of control points.
// Input: potrace_curve, curves in potrace formate
// Returns: control points (2 for line, 4 for bezier)
std::vector<Eigen::Matrix2Xd> decode_potrace_curve ( potrace_curve_t* potrace_curve );

// My non-finished attempt in extracting only the fitting part of the potrace pipeline
// std::vector<Eigen::Matrix2Xd> fit_polygon(const Eigen::Matrix2Xd &polygon, const double alphamax);

// Run the full potrace pipeline on a pixel boundary path.
// Returns: the traces curves (4 control points for each bezier, 2 for line)
// Input: path, the pixel boundary path
//        should_adjust_path, a potrace option. Indicates whether the polygon points are allowed to
//          not be on pixel path points
//        optimize_curve, whether to merge beziers 
std::vector<Eigen::Matrix2Xd> run_pipeline ( const Eigen::Matrix2Xd& path, bool should_adjust_path, bool optimize_curve );

NAMESPACE_END ( PotraceUtil )
NAMESPACE_END ( polyfit )
