#pragma once

#define POLYVEC_SEQUENCE_REMOVAL_NO_PRIORITY 0
#define POLYVEC_FITTER_EXPORT_TANGENTS_T 1

#define USE_CURVATURE_VARIATION_MINIMIZING_INITIALIZATION 1

#define REGULARIZE_BEZIER 1

#define ALLOW_CORNERS_TO_MOVE 1

#define G2_FOR_LINES 0

#define POLYVEC_LOG_CURVE_FITTER 0

#define ADAPTIVE_POINT_WEIGHT 0
const double ADAPTIVE_ACCURACY_CHECK_THRESHOLD = 1.0;

#define POST_ACCURACY_CHECK 1

namespace polyvec {
	// Spline fitting objective weights
	struct FitOptions {
		double tangents_weight;
		double points_weight;
		double curvature_variation_weight;
		double g2_weight;

		bool fix_endtangents;
		bool keep_endpoints_on_bisectors;
		bool add_soft_g2_constraint;

		bool consider_parallel;

		enum
		{
			DENSE_TANGENT_SAMPLES,
			CURVATURE_VARIATION,
		} fairness_formulation;

		FitOptions(bool final_fit)
		{
			tangents_weight = final_fit ? 1000 : 1000.0;
			points_weight = final_fit ? 10.0 : 10.0;
			curvature_variation_weight = 100.0;
			g2_weight = 100.0;

			fix_endtangents = !final_fit;

			// The endpoints must stay on the bisector if
			//   - the primitives fit the same corner (in order to preserve the appearance of the fit)
			//   - both primitives are Bezier curves
			keep_endpoints_on_bisectors = !final_fit;


			add_soft_g2_constraint = final_fit;

			fairness_formulation = (final_fit ? CURVATURE_VARIATION : DENSE_TANGENT_SAMPLES);

			consider_parallel = final_fit;
		}
	};

	const FitOptions final_fit = FitOptions(true);
	const FitOptions classification_fit = FitOptions(false);
}