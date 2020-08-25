/*
    This file contains the list of configurable options. Note that there are no
    default values and all value should be set to meaningful values, tipically through
    command line options. An example argv is stored in the root directory of the repository
*/
#ifndef polyvec_options_h_
#define polyvec_options_h_

#include <polyvec/geometry/angle.hpp>

namespace polyvec {
	// experiments
	enum class AccuracyMetric {
		OneSided,
		Implicit
	};


	// Vectorization options
    struct VectorOptions {
        static VectorOptions* get();

        void make_thread_local() const;
        void make_global() const;

		// Options for polygon tracing
		const double accuracy_thr = 1.0;
		const AccuracyMetric accuracy_metric = AccuracyMetric::Implicit;
		const double error_accuracy_weight = 1.;

		const double error_smoothness_weight = 1.;
		const double error_smoothness_limit = PF_RAD(135);

		const double error_continuity_weight = .25;
		const double error_continuity_limit = PF_RAD(60);

		const double error_inflection_limit = PF_RAD(90);
		const double error_inflection_penalty = 0.1;		        	
			
		// Deprecated (still referenced in the code but not actively used)
		const double spline_fit_distance_corner_g0_thr = std::numeric_limits<double>::signaling_NaN();
		const double spline_fit_distance_point_g0_thr = 0.75;
		const double spline_fit_curvature_g0_thr = 1.0;		
    };
}

#endif // polyvec_options_h_