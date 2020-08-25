#ifndef DEADLINE_CODE_GLOBFIT_OBJECTIVE_IS_INCLUDED
#define DEADLINE_CODE_GLOBFIT_OBJECTIVE_IS_INCLUDED

// polyvec

#include <polyvec/curve-tracer/curve_parametrizations.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/curve-tracer/curve_objective.hpp>

namespace polyvec {  


// =============================================================
//                        Fit Points
// =============================================================
//    class GlobFitObjective_FitPointsProject: public GlobFitObjective {
//    public:
//        using PointWeightFeeder = std::function<void ( int, Eigen::Vector2d&, double& ) >;
//
//        void set_points_to_fit ( PointWeightFeeder, const int n_points );
//        void get_points_to_fit ( PointWeightFeeder& feeder, int& n_points ) const;
//
//        int n_equations();
//        GlobFitObjectiveType get_type();
//        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;
//
//        GlobFitObjective_FitPointsProject() = default;
//    private:
//        PointWeightFeeder _point_feeder = PointWeightFeeder ( nullptr );
//        int _n_points_to_fit = 0;
//
//        bool _is_all_data_provided();
//    };
//
//// =============================================================
////                        FIT TANGENTS
//// =============================================================
//
//    class GlobFitObjective_FitTangentsProject: public GlobFitObjective {
//    public:
//        using PointTangentWeightFeeder = std::function<void ( int, Eigen::Vector2d&, Eigen::Vector2d&, double& ) >;
//
//        void set_points_to_fit ( PointTangentWeightFeeder, const int n_points );
//        void get_points_to_fit ( PointTangentWeightFeeder& feeder, int& n_points ) const;
//
//        int n_equations();
//        GlobFitObjectiveType get_type();
//        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;
//
//        GlobFitObjective_FitTangentsProject() = default;
//    private:
//        PointTangentWeightFeeder _tangent_feeder = PointTangentWeightFeeder ( nullptr );
//        int _n_points_to_fit = 0;
//
//        bool _is_all_data_provided();
//    };

// =============================================================
//                          FIX POSITION
// =============================================================

    class GlobFitObjective_FitPointLength : public GlobFitObjective {
    public:
		void set_params(GlobFitCurveParametrization* curve, const double t_curve, const Eigen::Vector2d& pos);        

        int n_equations() const;
        GlobFitObjectiveType get_type() const;
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
        double _t_curve = -1;
        Eigen::Vector2d _fixed_pos = Eigen::Vector2d::Constant ( std::numeric_limits<double>::max() );

        bool _is_all_data_provided() const;
    };

// =============================================================
//                          FIX TANGENT
// =============================================================
#if 1
	class GlobFitObjective_FitTangentLength : public GlobFitObjective {
	public:
		void set_params(GlobFitCurveParametrization* curve);
		void add_tangent(const double t_curve, const Eigen::Vector2d& tangent);

		// Default is using the cross product formulation (1 equation)
		// Use this to change to the entry-wise difference formultation (2 equations)
		bool get_using_cross_product_formulation() const;
		void set_using_cross_product_formulation(const bool);

		int n_equations() const;
		GlobFitObjectiveType get_type() const;
		void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams) const;

	private:
		bool  _use_cross_product_formulation = true;

		std::vector<double> ts;
		std::vector<Eigen::Vector2d> tangents;		

		bool _is_all_data_provided() const;
	};
#else

class GlobFitObjective_FitTangentLength : public GlobFitObjective {
    public:
        void set_t_for_fixed_tangent ( const double t_curve );
        void set_fixed_tangent ( const Eigen::Vector2d& );

        // Default is using the cross product formulation (1 equation)
        // Use this to change to the entry-wise difference formultation (2 equations)
        bool get_using_cross_product_formulation();
        void set_using_cross_product_formulation ( const bool );

        int n_equations();
        GlobFitObjectiveType get_type();
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
        bool  _use_cross_product_formulation = true;
        double _t_curve=-1;
        Eigen::Vector2d _fixed_tangent = Eigen::Vector2d::Constant ( std::numeric_limits<double>::max() );

        bool _is_all_data_provided();
    };

#endif

// =============================================================
//                          SAME POSITION
// =============================================================

    class GlobFitObjective_SamePosition : public GlobFitObjective {
    public:
		void set_params(GlobFitCurveParametrization* curve1, GlobFitCurveParametrization* curve2);

        int n_equations() const;
        GlobFitObjectiveType get_type() const;
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
        bool _is_all_data_provided() const;
    };

// =============================================================
//                          SAME TANGENT
// =============================================================

    class GlobFitObjective_SameTangent : public GlobFitObjective {
    public:
        enum FormulationType { FORMULATION_VECTOR_DIFF=0, FORMULATION_VECTOR_CROSS };

		void set_params(GlobFitCurveParametrization* curve1, GlobFitCurveParametrization* curve2);

        int n_equations() const;
        GlobFitObjectiveType get_type() const;
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

        FormulationType get_formulation() const;
        void set_formulation(FormulationType);

    private:
        FormulationType _formulation_type = FORMULATION_VECTOR_CROSS;
        bool _is_all_data_provided() const;
    };

// =============================================================
//                          SAME CURVATURE
// =============================================================

    class GlobFitObjective_SameCurvature : public GlobFitObjective {
    public:
        void set_params ( GlobFitCurveParametrization* curve1, const double t1, GlobFitCurveParametrization* curve2, const double t2 );

        int n_equations() const;
        GlobFitObjectiveType get_type() const;
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
		double t1 = -1, t2 = -1;
        bool _is_all_data_provided() const;
    };


// =============================================================
//                BEZIER SMALL DEVIATION FROM INIT
// =============================================================

    class GlobFitObjective_BezierFairness : public GlobFitObjective {
    public:
        // This will make an internal copy, must be called
        void set_initial_bezier ( const Eigen::Matrix2Xd& control_points );

        void set_target_length ( const double length );

		void set_curve(GlobFitCurveParametrization*);

        //
        // Description of each term.
        //
        // weight_end_tangent_length_dev:
        //   prevents the distance from the first control point and the second (also third and fourth)
        //   one to deviate much from the initial guess.
        // weight_end_tangent_angle_dev
        //   prevents the the angle of the line between the first and second (also fourth and third)
        //   control points to deviate from the initial value.
        // weight_end_points_dev
        //   prevents the position of the first and last control points to deviate much from their initial values.
        // weight_length
        //   Regularizer on the length
        // weight_prevent_control_point_overlap
        // thr_prevent_control_point_overlap
        //   prevents the first and second (also third and fourth) control points to overlap.
        // weight_not_too_far_control_points
        //   If the line passing from the first-second and fourth-third control points intersect on the
        //   (second-third) control points side, encourages the control points not to get further than
        //   than that point. This objective can prevent inflection and loops.
        void set_sub_weights (
            const double  weight_end_tangent_length_dev,
            const double  weight_end_tangent_angle_dev,
            const double  weight_end_points_dev,
            const double  weight_length,
            const double weight_prevent_control_point_overlap,
            const double thr_prevent_control_point_overlap,
            const double weight_not_too_far_control_points );

        int n_equations() const;
        GlobFitObjectiveType get_type() const;
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
        bool _is_initial_bezier_set = false;

        double _weight_end_tangent_length_dev = -1;
        double _weight_end_tangent_angle_dev = -1;

        double _weight_end_points_dev = -1;

        double _length_target = -1.;
        double _weight_length = -1;

        double _weight_prevent_control_point_overlap = -1.;
        double _thr_prevent_control_point_overlap = 1.;

        double _weight_not_too_far_control_points_beg = -1;
        double _weight_not_too_far_control_points_end = -1;

		BezierCurve _initial_bezier;
        bool _is_all_data_provided() const;
    };


// =============================================================
//                Point should lie on line
// =============================================================

    class GlobFitObjective_PointLieOnLine : public GlobFitObjective {
    public:
		void set_params(GlobFitCurveParametrization* curve, const double t_curve, const Eigen::Vector2d& point, const Eigen::Vector2d& tangent);

        int n_equations() const;
        GlobFitObjectiveType get_type() const;
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
        double _t_curve = -1;
        Eigen::Vector2d _line_point = Eigen::Vector2d::Constant ( std::numeric_limits<double>::max() );
        Eigen::Vector2d _line_tangent = Eigen::Vector2d::Constant ( std::numeric_limits<double>::max() );
        Eigen::Matrix2d _one_minus_ttT =  Eigen::Matrix2d::Constant ( std::numeric_limits<double>::max() );

        bool _is_all_data_provided() const;
    };

// =============================================================
//                Curvature Variation Minimization
// =============================================================

	class GlobFitObjective_CurvatureVariation : public GlobFitObjective {
	public:
		void set_params(GlobFitCurveParametrization* curve);

		GlobFitObjectiveType get_type()  const { return GLOBFIT_OBJECTIVE_CURVATURE_VARIATION; }

		int n_equations() const;
		void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams) const;
	private:
		int _M; //half number of samples
		double _h; //sample spacing

		//Calculates the curvature variation at a given position of the curve and fills the gradient
		//with respect to the parametrization.
		double calculate_curvature_variation(double t, Eigen::Block<Eigen::MatrixXd, 1> gradient) const;
	};

// =============================================================
//                Regularization
// =============================================================

	class GlobFitObjective_Barrier : public GlobFitObjective
	{
	public:
		void set_params(double soft_limit, double hard_limit);

	//protected:
		template <typename T>
		double compute_barrier(double value, const T& value_jacobian, Eigen::MatrixXd& out_barrier_jacobian) const
		{
			out_barrier_jacobian.resize(1, value_jacobian.cols());
			out_barrier_jacobian.setZero();

			double residual = value - soft_limit;

			if ((residual > 0) != ((hard_limit - soft_limit) > 0))
				return 0; //not exceeding the soft limit
			else
			{				
				out_barrier_jacobian = 2 * residual * value_jacobian;
				return residual * residual;
			}
		}

	//private:
		double soft_limit, hard_limit;
	};

	//Make sure that tangent lengths do not degenerate and no inflection is introduced
	class GlobFitObjective_ParameterBound : public  GlobFitObjective_Barrier {
	public:
		// Adds a penalty for parameters exceeding the soft_limit. The penalty is scaled such that it becomes 10,000 when
		// the parameter approaches hard_limit
		void set_params(GlobFitCurveParametrization* curve, int parameter, double soft_limit, double hard_limit);

		GlobFitObjectiveType get_type() const { return GLOBFIT_OBJECTIVE_PARAMETER_BOUND; }

		int n_equations() const { return 1; }
		void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams) const;

	//private:
		int parameter;		
	};

	class GlobFitObjective_LineRegularization : public GlobFitObjective_Barrier {
	public:
		void set_params(GlobFitCurveParametrization* curve);

		GlobFitObjectiveType get_type() const { return GLOBFIT_OBJECTIVE_LINE_REGULARIZATION; }

		int n_equations() const { return 1; }
		void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams) const;
	};



	class GlobFitObjective_LineAngle : public  GlobFitObjective {
	public:

		void set_params(GlobFitLineParametrization* line, double angle);

		GlobFitObjectiveType get_type() const { return GLOBFIT_OBJECTIVE_LINE_ANGLE; }

		int n_equations() const { return 1; }
		void compute_objective_and_jacobian(Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams) const;

		double angle;
	};

// =============================================================
//        Regularize the scale and transition only Beziers
// =============================================================

   /* class GlobFitObjective_ScaleAndTranslationOnlyRegularizer : 
    public GlobFitObjective {
    public:

        void set_sub_weights (
            const double  weight_for_scale,
            const double weight_for_similarity,
            const double  weight_for_translation);

        int n_equations();
        GlobFitObjectiveType get_type();
        void compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const;

    private:
        double _sqrt_weight_for_scale=-1;
        double _sqrt_weight_for_similarity = -1;
        double _sqrt_weight_for_translation =-1;
        bool _is_all_data_provided();
    };*/

} // end of polyvec


#endif // DEADLINE_CODE_GLOBFIT_OBJECTIVE_IS_INCLUDED
