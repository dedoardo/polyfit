// Author: A very tired Shayan Hoshyari (Mr. Chain rule).
// Header
#include <polyvec/curve-tracer/curve.hpp>

// polyvec
#include <polyvec/geometry/line.hpp>
#include <polyvec/curve-tracer/curve_objectives.hpp>
#include <polyvec/curve-tracer/curve_bezier.hpp>
#include <polyvec/geometry/smooth_curve.hpp>

#define MODIFIED_G2 0

using namespace polyfit;

namespace polyvec {

// =============================================================
//                           FIT POINTS
// =============================================================

//    void GlobFitObjective_FitPointsProject::set_points_to_fit ( PointWeightFeeder feeder, const int n_points ) {
//        _point_feeder = feeder;
//        _n_points_to_fit = n_points;
//    }
//
//    void GlobFitObjective_FitPointsProject::get_points_to_fit ( PointWeightFeeder& feeder, int& n_points ) const {
//        feeder = _point_feeder;
//        n_points = _n_points_to_fit;
//    }
//
//    int GlobFitObjective_FitPointsProject::n_equations() {
//        return int ( _n_points_to_fit ) * 2;
//    }
//
//    GlobFitObjectiveType GlobFitObjective_FitPointsProject::get_type() {
//        return GLOBFIT_OBJECTIVE_FIT_POINTS;
//    }
//
//    void GlobFitObjective_FitPointsProject::compute_objective_and_jacobian (
//        Eigen::VectorXd& obj,
//        Eigen::MatrixXd& dobj_dcurveparams ) const {
//
//        assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );
//
//        dobj_dcurveparams.setZero ( n_equations(), _n_sum_curve_params() );
//        obj.setZero ( n_equations() );
//
//        GlobFitCurve& curve = *_cached_curves.front();
//
//        //
//        // Points to fit
//        //
//        for ( int ptid = 0; ptid < ( int ) _n_points_to_fit; ++ptid ) {
//
//            Eigen::Vector2d point_to_fit;
//            double weight_to_fit, sqrt_weight_to_fit;
//            _point_feeder ( ptid, point_to_fit, weight_to_fit );
//            sqrt_weight_to_fit = sqrt ( weight_to_fit );
//
//            const double t_proj = curve.project ( point_to_fit );
//            Eigen::Vector2d projection = curve.pos ( t_proj );
//            const Eigen::VectorXd dtprojdparams = curve.dtprojectdparams ( t_proj, point_to_fit );
//            const Eigen::Matrix2Xd dpointprojdparams = curve.dposprojectdparams ( t_proj, dtprojdparams );
//
//            const Eigen::Vector2d distance_vector = projection - point_to_fit;
//
//            obj.segment<2> ( 2 * ptid ) = sqrt_weight_to_fit * distance_vector;
//            dobj_dcurveparams.middleRows<2> ( 2 * ptid ) = sqrt_weight_to_fit * dpointprojdparams;
//        }
//
//        // force user to call set_curves() each time;
//        _forget_curves();
//    }
//
//    bool GlobFitObjective_FitPointsProject::_is_all_data_provided() {
//        bool answer = true;
//        answer = answer && GlobFitObjective::_is_all_data_provided();
//        answer = answer && ( _cached_curves.size() == 1 );
//        answer = answer && ( _n_points_to_fit > 0 );
//        answer = answer && ( _point_feeder );
//        return answer;
//    }
//
//// =============================================================
////                           FIT TANGENTS PROJECT
//// =============================================================
//
//    void GlobFitObjective_FitTangentsProject::set_points_to_fit (
//        PointTangentWeightFeeder feeder,
//        const int n_points ) {
//        _tangent_feeder = feeder;
//        _n_points_to_fit = n_points;
//    }
//
//    void GlobFitObjective_FitTangentsProject::get_points_to_fit (
//        PointTangentWeightFeeder& feeder,
//        int& n_points ) const {
//        feeder = _tangent_feeder;
//        n_points = _n_points_to_fit;
//    }
//
//    int GlobFitObjective_FitTangentsProject::n_equations() {//
//        return _n_points_to_fit * 2;
//    }
//
//    GlobFitObjectiveType GlobFitObjective_FitTangentsProject::get_type() {//
//        return GLOBFIT_OBJECTIVE_FIT_TANGENTS;
//    }
//
//    void GlobFitObjective_FitTangentsProject::compute_objective_and_jacobian (
//        Eigen::VectorXd& obj,
//        Eigen::MatrixXd& jac ) const {
//
//        assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );
//
//        jac.setZero ( n_equations(), _n_sum_curve_params() );
//        obj.setZero ( n_equations() );
//
//        GlobFitCurve& curve = *_cached_curves.front();
//
//        //
//        // Points to fit
//        //
//        for ( int ptid = 0; ptid < ( int ) _n_points_to_fit; ++ptid ) {
//
//            Eigen::Vector2d point_to_fit, tangent_to_fit;
//            double weight_to_fit, sqrt_weight_to_fit;
//            _tangent_feeder ( ptid, point_to_fit, tangent_to_fit, weight_to_fit );
//            sqrt_weight_to_fit = sqrt ( weight_to_fit );
//            tangent_to_fit.normalize();
//
//            const double t_proj = curve.project ( point_to_fit );
//            Eigen::Vector2d dposdt = curve.dposdt ( t_proj );
//            const Eigen::VectorXd dtdparams = curve.dtprojectdparams ( t_proj, point_to_fit );
//            const Eigen::Matrix2Xd dposdtdbezierparams = curve.dposdtprojectdparams ( t_proj, dtdparams );
//
//            const double tol = 1e-15;
//
//#if 1
//            // Actual tangent and its derivative -- must normalizer derivative
//            Eigen::Vector2d tang;
//            Eigen::Matrix2d dtang_ddposdt;
//            SmoothCurveUtil::tangent_derivatives ( dposdt, tang, dtang_ddposdt );
//# else
//            // just for testing -- don't normalize
//            param_unused ( tol );
//            param_unused ( dposdt_norm2 );
//            param_unused ( dposdt_norm );
//            const Eigen::Vector2d tang = dposdt;
//            const Eigen::Matrix2d dtang_ddposdt = Eigen::Matrix2d::Identity();
//#endif
//            const Eigen::Vector2d diff = tang - tangent_to_fit;
//
//            // Chain rule is your friend.
//            obj.segment<2> ( 2 * ptid ) = sqrt_weight_to_fit * diff;
//            jac.middleRows<2> ( 2 * ptid ) = sqrt_weight_to_fit * dtang_ddposdt * dposdtdbezierparams;
//        }
//
//        // force user to call set_curves() each time;
//        _forget_curves();
//    }
//
//
//    bool GlobFitObjective_FitTangentsProject::_is_all_data_provided() {
//        bool answer = true;
//        answer = answer && GlobFitObjective::_is_all_data_provided();
//        answer = answer && ( _cached_curves.size() == 1 );
//        answer = answer && ( _n_points_to_fit > 0 );
//        answer = answer && ( _tangent_feeder );
//        return answer;
//    }

// =============================================================
//                          FIX POSITION
// =============================================================

	void GlobFitObjective_FitPointLength::set_params(GlobFitCurveParametrization * curve, const double t_curve, const Eigen::Vector2d & pos)
	{
		_t_curve = t_curve;
		_fixed_pos = pos;
		_curves.clear();
		_curves.push_back(curve);
	}

	int GlobFitObjective_FitPointLength::n_equations() const {
        return 2;
    }

    GlobFitObjectiveType GlobFitObjective_FitPointLength::get_type() const {
        return GLOBFIT_OBJECTIVE_FIX_POSITION;
    }

	void GlobFitObjective_FitPointLength::compute_objective_and_jacobian (
		Eigen::VectorXd& obj,
		Eigen::MatrixXd& dobj_dcurveparams) const {

		assert_break_msg(_is_all_data_provided(), "Some data is not provided yet or is wrong");

		dobj_dcurveparams.resize(n_equations(), _n_sum_curve_params());
		obj.resize(n_equations());

		GlobFitCurveParametrization& curve0 = *_curves.front();
		Eigen::Vector2d pos_curve0 = curve0.get_curve()->pos(_t_curve);

		obj = (pos_curve0 - _fixed_pos);
		dobj_dcurveparams = curve0.dposdparams(_t_curve);
    }

    bool GlobFitObjective_FitPointLength::_is_all_data_provided() const {
        bool answer = true;
        const double tol = 1e-12;
        answer = answer && GlobFitObjective::_is_all_data_provided();
        answer = answer && ( _curves.size() == 1 );
        answer = answer && ( _t_curve >= 0 - tol );
        answer = answer && ( _t_curve <= GlobFitCurve::t_end + tol );
        answer = answer && ( _fixed_pos.norm() < std::numeric_limits<double>::max() );
        return answer;
    }

// =============================================================
//                          FIX TANGENT
// =============================================================

	void GlobFitObjective_FitTangentLength::set_params(GlobFitCurveParametrization* curve)
	{
		_curves.clear();
		_curves.push_back(curve);
	}

	void GlobFitObjective_FitTangentLength::add_tangent(const double t_curve, const Eigen::Vector2d& tangent)
	{
		ts.push_back(t_curve);
		tangents.push_back(tangent);
	}

	bool GlobFitObjective_FitTangentLength::get_using_cross_product_formulation() const {
		return _use_cross_product_formulation;
	}

	void GlobFitObjective_FitTangentLength::set_using_cross_product_formulation(const bool in) {
		_use_cross_product_formulation = in;
	}


	int GlobFitObjective_FitTangentLength::n_equations() const {
		// Using the cross product formulation
		if (_use_cross_product_formulation) {
			return 1 * ts.size();
			// Using the coordinatewise difference formultation
		}
		else {
			return 2 * ts.size();
		}
	}

	GlobFitObjectiveType GlobFitObjective_FitTangentLength::get_type() const {
		return GLOBFIT_OBJECTIVE_FIX_TANGENT;
	}

	void GlobFitObjective_FitTangentLength::compute_objective_and_jacobian (
		Eigen::VectorXd& obj,
		Eigen::MatrixXd& dobj_dcurveparams) const {

		assert_break_msg(_is_all_data_provided(), "Some data is not provided yet or is wrong");

		dobj_dcurveparams.resize(n_equations(), _n_sum_curve_params());
		obj.resize(n_equations());

		GlobFitCurveParametrization& curve0 = *_curves.front();

		Eigen::RowVectorXd dt_dparams(1, curve0.n_params());

		for (int i = 0; i < ts.size(); ++i)
		{
			double t = ts[i];			
			dt_dparams.setZero();
			
			auto& tangent = tangents[i];

			Eigen::Vector2d dposdt = curve0.get_curve()->dposdt(t);
			auto dposdt_dparams = curve0.dposdtdparams(t) + curve0.get_curve()->dposdtdt(t) * dt_dparams;

			// Using the cross product formulation
			if (_use_cross_product_formulation) {
				obj(i) = (dposdt.x() * tangent.y() - dposdt.y() * tangent.x());

				Eigen::RowVector2d dobj_ddposdt;
				dobj_ddposdt << tangent.y(), -tangent.x();

				dobj_dcurveparams.row(i) = dobj_ddposdt * dposdt_dparams;
			}
			// Using the difference formulation
			else {
				const double tol = 1e-4;
				double dposdt_norm2 = dposdt.squaredNorm();
				if (dposdt_norm2 < tol)
				{
					obj.row(2 * i + 0).setZero();
					obj.row(2 * i + 1).setZero();
					dobj_dcurveparams.row(2 * i + 0).setZero();
					dobj_dcurveparams.row(2 * i + 1).setZero();
				}
				else
				{
					double dposdt_norm = sqrt(dposdt_norm2);

					const Eigen::Vector2d tang = dposdt / dposdt_norm;
					const Eigen::Matrix2d dtang_ddposdt =
						Eigen::Matrix2d::Identity() / dposdt_norm
						- 1. / (dposdt_norm * dposdt_norm2) * dposdt * dposdt.transpose();
					auto dtang_dparams = dtang_ddposdt * dposdt_dparams;

					const Eigen::Vector2d diff = tang - tangent;

					// Chain rule is your friend.
					obj.block<2, 1>(2 * i, 0) = diff;
					dobj_dcurveparams.block<2, Eigen::Dynamic>(2 * i, 0, 2, dtang_dparams.cols()) = dtang_dparams;
				}
			}
		}	
	}

	bool GlobFitObjective_FitTangentLength::_is_all_data_provided() const {
		bool answer = true;
		answer = answer && GlobFitObjective::_is_all_data_provided();
		answer = answer && (_curves.size() == 1);		
		return answer;
	}
    
// =============================================================
//                          SAME POSITION
// =============================================================

	void GlobFitObjective_SamePosition::set_params(GlobFitCurveParametrization * curve1, GlobFitCurveParametrization * curve2)
	{
		_curves.clear();
		_curves.push_back(curve1);
		_curves.push_back(curve2);
	}

	int GlobFitObjective_SamePosition::n_equations() const {
        return 2;
    }

    GlobFitObjectiveType GlobFitObjective_SamePosition::get_type() const {
        return GLOBFIT_OBJECTIVE_SAME_POSITION;
    }

    void GlobFitObjective_SamePosition::compute_objective_and_jacobian (
        Eigen::VectorXd& obj,
        Eigen::MatrixXd& dobj_dcurveparams ) const {

        assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );

        dobj_dcurveparams.resize ( n_equations(), _n_sum_curve_params() );
        obj.resize ( n_equations() );

        GlobFitCurveParametrization& curve0 = *_curves.front();
        GlobFitCurveParametrization& curve1 = *_curves.back();

        Eigen::Vector2d pos_curve0 = curve0.get_curve()->pos ( 1 );
        Eigen::Vector2d pos_curve1 = curve1.get_curve()->pos ( 0 );

        obj = ( pos_curve0 - pos_curve1 );
        dobj_dcurveparams.leftCols ( curve0.n_params() ) = curve0.dposdparams ( 1 );
        dobj_dcurveparams.rightCols ( curve1.n_params() ) = - curve1.dposdparams ( 0 );
    }

    bool GlobFitObjective_SamePosition::_is_all_data_provided() const {
        bool answer = true;
        const double tol = 1e-12;
        answer = answer && GlobFitObjective::_is_all_data_provided();
        answer = answer && ( _curves.size() == 2 );
        return answer;
    }

// =============================================================
//                          SAME TANGENT
// =============================================================

	void GlobFitObjective_SameTangent::set_params(GlobFitCurveParametrization * curve1, GlobFitCurveParametrization * curve2)
	{
		_curves.clear();
		_curves.push_back(curve1);
		_curves.push_back(curve2);
	}

	int GlobFitObjective_SameTangent::n_equations() const {
        switch ( _formulation_type ) {
        case FORMULATION_VECTOR_CROSS:
            return 1;

        case FORMULATION_VECTOR_DIFF:
            return 2;
        }
    }

    GlobFitObjective_SameTangent::FormulationType GlobFitObjective_SameTangent::get_formulation() const {
        return _formulation_type;
    }

    void GlobFitObjective_SameTangent::set_formulation ( FormulationType in ) {
        _formulation_type=in;
    }


    GlobFitObjectiveType GlobFitObjective_SameTangent::get_type() const {
        return GLOBFIT_OBJECTIVE_SAME_TANGENT;
    }

    void GlobFitObjective_SameTangent::compute_objective_and_jacobian (
        Eigen::VectorXd& obj,
        Eigen::MatrixXd& dobj_dcurveparams ) const {

        assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );

        dobj_dcurveparams.resize ( n_equations(), _n_sum_curve_params() );
        obj.resize ( n_equations() );

        GlobFitCurveParametrization& curve0 = *_curves.front();
        GlobFitCurveParametrization& curve1 = *_curves.back();

        Eigen::Vector2d dposdt_curve0 = curve0.get_curve()->dposdt ( 1 );
        Eigen::Vector2d dposdt_curve1 = curve1.get_curve()->dposdt ( 0 );

        switch ( _formulation_type ) {
        case FORMULATION_VECTOR_CROSS: {
            obj ( 0 ) = ( dposdt_curve0.x() * dposdt_curve1.y() - dposdt_curve0.y() * dposdt_curve1.x() );

            Eigen::RowVector2d dobj_ddposdt_curve0;
            Eigen::RowVector2d dobj_ddposdt_curve1;
            dobj_ddposdt_curve0 << dposdt_curve1.y(), -dposdt_curve1.x();
            dobj_ddposdt_curve1 << -dposdt_curve0.y(), dposdt_curve0.x();

            dobj_dcurveparams.leftCols ( curve0.n_params() ) =
                dobj_ddposdt_curve0 * curve0.dposdtdparams ( 1 );
            dobj_dcurveparams.rightCols ( curve1.n_params() ) =
                dobj_ddposdt_curve1 * curve1.dposdtdparams ( 0 );

            break;
        }

        case FORMULATION_VECTOR_DIFF: {
            const double tol = 1e-15;
            //
            Eigen::Vector2d tang_curve0 ;
            Eigen::Matrix2d dtang_curve0_ddposdt;
            SmoothCurveUtil::tangent_derivatives ( dposdt_curve0, tang_curve0, dtang_curve0_ddposdt );
            //
            Eigen::Vector2d tang_curve1 ;
            Eigen::Matrix2d dtang_curve1_ddposdt;
            SmoothCurveUtil::tangent_derivatives ( dposdt_curve1, tang_curve1, dtang_curve1_ddposdt );

            obj = tang_curve0 - tang_curve1;

            dobj_dcurveparams.leftCols ( curve0.n_params() ) =
                dtang_curve0_ddposdt * curve0.dposdtdparams ( 1 );
            dobj_dcurveparams.rightCols ( curve1.n_params() ) =
                -dtang_curve1_ddposdt * curve1.dposdtdparams ( 0 );
            break;
        }

        }
    }

    bool GlobFitObjective_SameTangent::_is_all_data_provided() const {
        bool answer = true;
        const double tol = 1e-12;
        answer = answer && GlobFitObjective::_is_all_data_provided();
        answer = answer && ( _curves.size() == 2 );
        return answer;
    }

// =============================================================
//                          SAME CURVATURE
// =============================================================

	void GlobFitObjective_SameCurvature::set_params(GlobFitCurveParametrization* curve1, const double t1, GlobFitCurveParametrization* curve2, const double t2) {
		_curves.clear();
		_curves.push_back(curve1);
		_curves.push_back(curve2);
		this->t1 = t1;
		this->t2 = t2;
	}

	int GlobFitObjective_SameCurvature::n_equations() const {
		return 1;
	}

	GlobFitObjectiveType GlobFitObjective_SameCurvature::get_type() const {
		return GLOBFIT_OBJECTIVE_SAME_CURVATURE;
	}

	void GlobFitObjective_SameCurvature::compute_objective_and_jacobian (
		Eigen::VectorXd& obj,
		Eigen::MatrixXd& dobj_dcurveparams ) const {

		assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );

		dobj_dcurveparams.resize ( n_equations(), _n_sum_curve_params() );
		obj.resize ( n_equations() );

		const GlobFitCurveParametrization& curve1 = *_curves.front();
		const GlobFitCurveParametrization& curve2 = *_curves.back();		

		//
		// Get curvatures
		//

		auto get_curvature_and_jacobian = [](const GlobFitCurveParametrization& curve, double t) -> std::pair<double, Eigen::MatrixXd>
		{
			Eigen::Vector2d dposdt = curve.get_curve()->dposdt(t);
			Eigen::Vector2d dposdtdt = curve.get_curve()->dposdtdt(t);
			Eigen::Matrix2Xd ddposdt_dparams = curve.dposdtdparams(t);
			Eigen::Matrix2Xd ddposdtdt_dparams = curve.dposdtdtdparams(t);

			double curvature;
			Eigen::RowVector2d dcurvature_dposdt;
			Eigen::RowVector2d dcurvature_ddposdtdt;
			SmoothCurveUtil::curvature_derivatives(dposdt, dposdtdt, curvature, dcurvature_dposdt, dcurvature_ddposdtdt);

			Eigen::MatrixXd dCurvInputdParams1(4, curve.n_params());
			dCurvInputdParams1.topRows<2>() = ddposdt_dparams;
			dCurvInputdParams1.bottomRows<2>() = ddposdtdt_dparams;

			Eigen::MatrixXd dCurvdInput1(1, 4);
			dCurvdInput1 << dcurvature_dposdt, dcurvature_ddposdtdt;

			auto dCurvdParams1 = dCurvdInput1 * dCurvInputdParams1;

			return { curvature, dCurvdParams1 };
		};

		double curvature1, curvature2;
		Eigen::MatrixXd dcurv1_dparams1, dcurv2_dparams2;
		std::tie(curvature1, dcurv1_dparams1) = get_curvature_and_jacobian(curve1, t1);
		std::tie(curvature2, dcurv2_dparams2) = get_curvature_and_jacobian(curve2, t2);

		double rec_curvature1 = 1.0 / curvature1;
		double rec_curvature2 = 1.0 / curvature2;
		Eigen::MatrixXd drec_curv1_dparams1 = dcurv1_dparams1, drec_curv2_dparams2 = dcurv2_dparams2;
		for (int i = 0; i < drec_curv1_dparams1.size(); ++i)
			drec_curv1_dparams1.coeffRef(i) = -drec_curv1_dparams1.coeff(i) / (curvature1 * curvature1);
		for (int i = 0; i < drec_curv2_dparams2.size(); ++i)
			drec_curv2_dparams2.coeffRef(i) = -drec_curv2_dparams2.coeff(i) / (curvature2 * curvature2);

		//
		// Get the objective
		//

		//if (std::abs(curvature1) < PF_EPS || std::abs(curvature2) < PF_EPS)
		//{
		//	// if one of the two curves is a line, we don't enforce G2 continuity there.
		//	obj(0) = 0;
		//	dobj_dcurveparams.row(0).setZero();
		//	return;
		//}

		const double eps = 1.e-3;
		double num = curvature1 - curvature2;
		double denom;
		
		Eigen::MatrixXd dnumdenom_dparams(2, curve1.n_params() + curve2.n_params());
		dnumdenom_dparams.row(0).leftCols(curve1.n_params()) = dcurv1_dparams1;
		dnumdenom_dparams.row(0).rightCols(curve2.n_params()) = -dcurv2_dparams2;		

#if MODIFIED_G2
		if (std::abs(curvature1) < std::abs(curvature2))
		{
			denom = std::abs(curvature2) + eps;
			dnumdenom_dparams.row(1).setZero();
			dnumdenom_dparams.row(1).rightCols(curve2.n_params()) = curvature2 > 0 ? dcurv2_dparams2 : -dcurv2_dparams2;
		}
		else
		{
			denom = std::abs(curvature1) + eps;
			dnumdenom_dparams.row(1).setZero();
			dnumdenom_dparams.row(1).leftCols(curve1.n_params()) = curvature1 > 0 ? dcurv1_dparams1 : -dcurv1_dparams1;
		}
#else
		denom = std::abs(curvature1) + std::abs(curvature2) + eps;		

		dnumdenom_dparams.row(1).leftCols(curve1.n_params()) = curvature1 > 0 ? dcurv1_dparams1 : -dcurv1_dparams1;
		dnumdenom_dparams.row(1).rightCols(curve2.n_params()) = curvature2 > 0 ? dcurv2_dparams2 : -dcurv2_dparams2;

#endif

		obj(0) = num / denom;

		Eigen::RowVector2d division_jacobian(1 / denom, -num / (denom * denom));
		dobj_dcurveparams = division_jacobian * dnumdenom_dparams;
	}

	bool GlobFitObjective_SameCurvature::_is_all_data_provided() const {
		bool answer = true;
		const double tol = PF_EPS;
		answer = answer && GlobFitObjective::_is_all_data_provided();
		answer = answer && ( _curves.size() == 2 );
		answer = answer && ( t1 > 0 - tol );
		answer = answer && ( t2 > 0 - tol );
		answer = answer && ( t1 < GlobFitCurve::t_end + tol );
		answer = answer && ( t2 < GlobFitCurve::t_end + tol );
		return answer;
	}

// =============================================================
//                BEZIER SMALL DEVIATION FROM INIT
// =============================================================

    void GlobFitObjective_BezierFairness::set_initial_bezier ( const Eigen::Matrix2Xd& control_points ) {
		_initial_bezier.set_control_points ( control_points );
        _is_initial_bezier_set = true;
    }

    void GlobFitObjective_BezierFairness::set_target_length ( const double length ) {
        _length_target = length;
    }

	void GlobFitObjective_BezierFairness::set_curve(GlobFitCurveParametrization * c)
	{
		auto bezier = dynamic_cast<BezierCurve*>(c->get_curve().get());
		assert_break(bezier != nullptr);
		_curves.clear();
		_curves.push_back(c);
	}

    void GlobFitObjective_BezierFairness::set_sub_weights (
        const double  weight_end_tangent_length_dev,
        const double  weight_end_tangent_angle_dev,
        const double  weight_end_points_dev,
        const double  weight_length,
        const double weight_prevent_control_point_overlap,
        const double thr_prevent_control_point_overlap,
        const double weight_not_too_far_control_points ) {

        // SHAYAN: This should be sqrt'd argh...
        // I forgot it. But now if I fix it, the results will change.
        // TODO: cooridate with ED.

        _weight_end_tangent_length_dev = weight_end_tangent_length_dev;
        _weight_end_tangent_angle_dev = weight_end_tangent_angle_dev;
        _weight_length = weight_length;
        _weight_end_points_dev = weight_end_points_dev;
        _weight_prevent_control_point_overlap = weight_prevent_control_point_overlap;
        _thr_prevent_control_point_overlap = thr_prevent_control_point_overlap;

        _weight_not_too_far_control_points_beg = 0;
        _weight_not_too_far_control_points_end = 0;

        if ( weight_not_too_far_control_points > 0 ) {
            assert_break ( _is_initial_bezier_set );

            const Eigen::Matrix2Xd control_points = _initial_bezier.get_control_points();
            double param_on_first_curve;
            double param_on_second_curve;
            bool do_intersect = LineUtils::intersect (
                                    control_points.col ( 0 ),
                                    control_points.col ( 1 ),
                                    control_points.col ( 3 ),
                                    control_points.col ( 2 ),
                                    param_on_first_curve,
                                    param_on_second_curve );

            if ( do_intersect && ( param_on_first_curve > 0 ) && ( param_on_second_curve > 0 ) ) { //
                const Eigen::Vector2d intersection = LineUtils::line_at (
                        control_points.col ( 0 ), control_points.col ( 1 ), param_on_first_curve );
                const double len1 = ( intersection - control_points.col ( 0 ) ).norm();
                const double len2 = ( intersection - control_points.col ( 3 ) ).norm();

                _weight_not_too_far_control_points_beg = weight_not_too_far_control_points / std::max ( len1, 0.1 );
                _weight_not_too_far_control_points_end = weight_not_too_far_control_points / std::max ( len2, 0.1 );
            } // if do_intersect
        } // if weight_not_too_far_control_points > 0
    } // all done


    int GlobFitObjective_BezierFairness::n_equations() const {
        return 13;
    }


    GlobFitObjectiveType GlobFitObjective_BezierFairness::get_type() const {
        return GLOBFIT_OBJECTIVE_BEZIER_FAIRNESS;
    }

    void GlobFitObjective_BezierFairness::compute_objective_and_jacobian (
        Eigen::VectorXd& obj,
        Eigen::MatrixXd& dobj_dcurveparams ) const {

        assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );

        dobj_dcurveparams.resize ( n_equations(), _n_sum_curve_params() );
        obj.resize ( n_equations() );

        GlobFitCurveParametrization& curve = *_curves.front();
		auto bezier = dynamic_cast<BezierCurve*>(curve.get_curve().get());
        assert_break (bezier != nullptr);        
        Eigen::VectorXd bezier_params = bezier->get_params();
        Eigen::VectorXd initial_params = _initial_bezier.get_params();


        // We must know the index mapping before hand.
        const int idx_l0 = 0;
        const int idx_l1 = 1;
        const int idx_pt0 = 2;
        const int idx_pt1 = 4;
        const int idx_angle0 = 6;
        const int idx_angle1 = 7;

        // Deviation of distance between control points
        obj ( 0 ) = _weight_end_tangent_length_dev * ( bezier_params ( idx_l0 ) - initial_params ( idx_l0 ) );
        obj ( 1 ) = _weight_end_tangent_length_dev * ( bezier_params ( idx_l1 ) - initial_params ( idx_l1 ) );
        // Deviation of begin and end angles
        obj ( 2 ) = _weight_end_tangent_angle_dev * ( bezier_params ( idx_angle0 ) - initial_params ( idx_angle0 ) );
        obj ( 3 ) = _weight_end_tangent_angle_dev * ( bezier_params ( idx_angle1 ) - initial_params ( idx_angle1 ) );
        // Deviation of begin and end control points
        obj.segment<2> ( 4 ) = _weight_end_points_dev
                               * ( bezier_params.segment<2> ( idx_pt0 ) - initial_params.segment<2> ( idx_pt0 ) );
        obj.segment<2> ( 6 ) = _weight_end_points_dev
                               * ( bezier_params.segment<2> ( idx_pt1 ) - initial_params.segment<2> ( idx_pt1 ) );
        // Regularizing length
        const double target_length = _length_target > 0. ? _length_target : _initial_bezier.length();
        obj ( 8 ) = _weight_length * ( bezier->length() - target_length  );

        // Don't let control points collapse on each other
        obj ( 9 ) = 0;
        obj ( 10 ) = 0;

        if ( bezier_params ( idx_l0 ) < _thr_prevent_control_point_overlap )    {
            obj ( 9 ) = _weight_prevent_control_point_overlap
                        * ( _thr_prevent_control_point_overlap - bezier_params ( idx_l0 ) );
        }

        if ( bezier_params ( idx_l1 ) < _thr_prevent_control_point_overlap ) {
            obj ( 10 ) = _weight_prevent_control_point_overlap
                         * ( _thr_prevent_control_point_overlap - bezier_params ( idx_l1 ) );
        }

        // Not too far control points
        obj ( 11 ) = _weight_not_too_far_control_points_beg * ( bezier_params ( idx_l0 )  );
        obj ( 12 ) = _weight_not_too_far_control_points_end * ( bezier_params ( idx_l1 )  );

        // Now derivatives
        dobj_dcurveparams.setZero();
        dobj_dcurveparams ( 0, idx_l0 ) = _weight_end_tangent_length_dev;
        dobj_dcurveparams ( 1, idx_l1 ) = _weight_end_tangent_length_dev;
        dobj_dcurveparams ( 2, idx_angle0 ) = _weight_end_tangent_angle_dev;
        dobj_dcurveparams ( 3, idx_angle1 ) = _weight_end_tangent_angle_dev;
        dobj_dcurveparams ( 4, idx_pt0 + 0 ) = _weight_end_points_dev;
        dobj_dcurveparams ( 5, idx_pt0 + 1 ) = _weight_end_points_dev;
        dobj_dcurveparams ( 6, idx_pt1 + 0 ) = _weight_end_points_dev;
        dobj_dcurveparams ( 7, idx_pt1 + 1 ) = _weight_end_points_dev;
        dobj_dcurveparams.row ( 8 ) = _weight_length * bezier->dlengthdparams();

        if ( bezier_params ( idx_l0 ) < _thr_prevent_control_point_overlap ) {
            dobj_dcurveparams ( 9, idx_l0 ) = -_weight_prevent_control_point_overlap;
        }

        if ( bezier_params ( idx_l1 ) < _thr_prevent_control_point_overlap ) {
            dobj_dcurveparams ( 10, idx_l1 ) = -_weight_prevent_control_point_overlap;
        }

        dobj_dcurveparams ( 11, idx_l0 ) = _weight_not_too_far_control_points_beg * ( bezier_params ( idx_l0 ) );
        dobj_dcurveparams ( 12, idx_l1 ) = _weight_not_too_far_control_points_end * ( bezier_params ( idx_l1 ) );
    }

    bool GlobFitObjective_BezierFairness::_is_all_data_provided() const {
        bool answer = true;
        answer = answer && GlobFitObjective::_is_all_data_provided();
        answer = answer && ( _weight_end_tangent_length_dev >= 0. );
        answer = answer && ( _weight_end_tangent_angle_dev >= 0. );
        answer = answer && ( _weight_length >= 0. );
        answer = answer && ( _weight_end_points_dev >= 0. );
        answer = answer && ( _weight_prevent_control_point_overlap >= 0. );
        answer = answer && ( _thr_prevent_control_point_overlap >= 0. );
        answer = answer && ( _weight_not_too_far_control_points_beg >= 0. );
        answer = answer && ( _weight_not_too_far_control_points_end >= 0. );
        answer = answer && ( _is_initial_bezier_set );
        return answer;
    }

// =============================================================
//                Point should lie on line
// =============================================================

	void GlobFitObjective_PointLieOnLine::set_params(GlobFitCurveParametrization * curve, const double t_curve, const Eigen::Vector2d & point, const Eigen::Vector2d & tangent)
	{
		_t_curve = t_curve;
		_line_point = point;
		_line_tangent = tangent.stableNormalized();
		_one_minus_ttT = Eigen::Matrix2d::Identity() - _line_tangent * _line_tangent.transpose();

		_curves.clear();
		_curves.push_back(curve);
	}

	int GlobFitObjective_PointLieOnLine::n_equations() const {
        return 2;
    }

    GlobFitObjectiveType GlobFitObjective_PointLieOnLine::get_type() const {
        return GLOBFIT_OBJECTIVE_POSITION_ON_CURVE_LIE_ON_LINE;
    }

    void GlobFitObjective_PointLieOnLine::compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const {

        assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );

        dobj_dcurveparams.resize ( n_equations(), _n_sum_curve_params() );
        obj.resize ( n_equations() );

        GlobFitCurveParametrization& curve0 = *_curves.front();

        obj =   _one_minus_ttT * ( curve0.get_curve()->pos ( _t_curve ) - _line_point ) ;
        dobj_dcurveparams = _one_minus_ttT * curve0.dposdparams ( _t_curve );
    }


    bool GlobFitObjective_PointLieOnLine::_is_all_data_provided() const {
        const double tol = 1e-12;
        bool answer = true;
        answer = answer && GlobFitObjective::_is_all_data_provided();
        answer = answer && ( _curves.size() == 1 );
        answer = answer && ( _t_curve > 0 - tol );
        answer = answer && ( _t_curve < GlobFitCurve::t_end + tol );
        answer = answer && ( _line_point.norm() < std::numeric_limits<double>::max() );
        answer = answer && ( _line_tangent.norm() < std::numeric_limits<double>::max() );
        return answer;
    }

// =============================================================
//                Curvature Variation Minimization
// =============================================================

	void GlobFitObjective_CurvatureVariation::set_params(GlobFitCurveParametrization* curve)
	{
		_curves.clear();
		_curves.push_back(curve);

		_M = 20;
		_h = 1.0 / (2 * _M);
	}

	int GlobFitObjective_CurvatureVariation::n_equations()  const
	{
		return 2 * _M + 1;
	}

	void GlobFitObjective_CurvatureVariation::compute_objective_and_jacobian(Eigen::VectorXd & obj, Eigen::MatrixXd & dobj_dcurveparams) const
	{
		dobj_dcurveparams.resize(n_equations(), _n_sum_curve_params());
		obj.resize(n_equations());

		double hOver3 = _h / 3.0;
		double sqrthOver3 = std::sqrt(hOver3);
		double sqrt2hOver3 = std::sqrt(2 * hOver3);
		double sqrt4hOver3 = std::sqrt(4 * hOver3);

		obj(0) = sqrthOver3 * calculate_curvature_variation(0, dobj_dcurveparams.row(0));
		dobj_dcurveparams.row(0) *= sqrthOver3;

		obj(1) = sqrthOver3 * calculate_curvature_variation(1, dobj_dcurveparams.row(1));
		dobj_dcurveparams.row(1) *= sqrthOver3;
		
		for (auto i = 1; i <= _M - 1; ++i)
		{
			double t = 2 * i * _h;
			int idx = 1 + i;
			obj(idx) = sqrt2hOver3 * calculate_curvature_variation(t, dobj_dcurveparams.row(idx));
			dobj_dcurveparams.row(idx) *= sqrt2hOver3;
		}

		for (auto i = 1; i <= _M; ++i)
		{
			double t = (2 * i - 1) * _h;
			int idx = _M + i;
			obj(idx) = sqrt4hOver3 * calculate_curvature_variation(t, dobj_dcurveparams.row(idx));
			dobj_dcurveparams.row(idx) *= sqrt4hOver3;
		}
	}

	double GlobFitObjective_CurvatureVariation::calculate_curvature_variation(double t, Eigen::Block<Eigen::MatrixXd, 1> gradient) const
	{
		auto curve = this->get_curves().front();
		auto d1 = curve->get_curve()->dposdt(t);
		auto d2 = curve->get_curve()->dposdtdt(t);
		auto d3 = curve->get_curve()->dposdtdtdt(t);

		auto d1Prime = curve->dposdtdparams(t);
		auto d2Prime = curve->dposdtdtdparams(t);
		auto d3Prime = curve->dposdtdtdtdparams(t);		

		//Planar cubic G 1 and quintic G 2 Hermite interpolations via curvature variation minimization
		//psi =   ||d1||^2   * (d1 cross  d3) - 3 * (d1 dot  d2) * (d1 cross  d2)
		//      '- term 1 -'   '-- term 2 --'       '- term 3 -'   '-- term 4 --'
		//      '--------- term 5 ----------'       '--------- term 6 ----------'

		double term1 = d1.squaredNorm();
		auto term1Prime = 2 * (d1.row(0) * d1Prime.row(0) + d1.row(1) * d1Prime.row(1));

		double term2 = d1.x() * d3.y() - d1.y() * d3.x();
		auto term2Prime = d1.x() * d3Prime.row(1) + d1Prime.row(0) * d3.y() - (d1.y() * d3Prime.row(0) + d1Prime.row(1) * d3.x());

		double term3 = d1.dot(d2);
		auto term3Prime = d1.x() * d2Prime.row(0) + d1Prime.row(0) * d2.x() + d1.y() * d2Prime.row(1) + d1Prime.row(1) * d2.y();

		double term4 = d1.x() * d2.y() - d1.y() * d2.x();
		auto term4Prime = d1.x() * d2Prime.row(1) + d1Prime.row(0) * d2.y() - (d1.y() * d2Prime.row(0) + d1Prime.row(1) * d2.x());
		
		double term5 = term1 * term2;
		auto term5Prime = term1 * term2Prime + term1Prime * term2;

		double term6 = term3 * term4;
		auto term6Prime = term3 * term4Prime + term3Prime * term4;

		double psi = term5 - 3 * term6;
		auto psiPrime = term5Prime - 3 * term6Prime;

		//kappaPrime = psi / ||d1||^5 = psi * term1^(-5/2) 
		double term1Pow = std::pow(term1, -2.5);
		double kappaPrime = psi * term1Pow;
		auto kappaPrimePrime = -2.5 * psi * term1Prime * std::pow(term1, -3.5) + psiPrime * term1Pow;

		gradient = kappaPrimePrime;

		return kappaPrime;
	}

// =============================================================
//        Regularize the scale and transition only Beziers
// =============================================================


    //void
    //GlobFitObjective_ScaleAndTranslationOnlyRegularizer::
    //set_sub_weights (
    //    const double  weight_for_scale,
    //    const double  weight_for_similarity,
    //    const double  weight_for_translation ) {
    //    _sqrt_weight_for_scale = sqrt ( weight_for_scale );
    //    _sqrt_weight_for_translation = sqrt ( weight_for_translation );
    //    _sqrt_weight_for_similarity = sqrt ( weight_for_similarity );
    //}

    //int
    //GlobFitObjective_ScaleAndTranslationOnlyRegularizer::
    //n_equations() {
    //    return 5; // two for translation, 2 for scale , one for similarity of scales
    //}

    //GlobFitObjectiveType
    //GlobFitObjective_ScaleAndTranslationOnlyRegularizer::
    //get_type() {
    //    return GLOBFIT_OBJECTIVE_SCALE_AND_TRANSITION_ONLY_REGULARIZER;
    //}

    //void
    //GlobFitObjective_ScaleAndTranslationOnlyRegularizer::
    //compute_objective_and_jacobian ( Eigen::VectorXd& obj, Eigen::MatrixXd& dobj_dcurveparams ) const {
    //    assert_break_msg ( _is_all_data_provided(), "Some data is not provided yet or is wrong" );

    //    dobj_dcurveparams.resize ( n_equations(), _n_sum_curve_params() );
    //    obj.resize ( n_equations() );

    //    GlobFitCurve& curve0 = *_cached_curves.front();
    //    Eigen::VectorXd params = curve0.get_params();

    //    // Must be consistent with  GlobFitCurve_BezierScaleAndTranslationOnly
    //    // If you change the param ordering there, this will break.
    //    const double sx = params ( 0 );
    //    const double sy = params ( 1 );
    //    const double tx = params ( 2 );
    //    const double ty = params ( 3 );

    //    obj ( 0 ) = _sqrt_weight_for_scale * ( sx - 1 );
    //    obj ( 1 ) = _sqrt_weight_for_scale * ( sy - 1 ) ;
    //    obj ( 2 ) = _sqrt_weight_for_translation * tx ;
    //    obj ( 3 ) = _sqrt_weight_for_translation * ty;
    //    obj ( 4 ) = _sqrt_weight_for_similarity * ( sx - sy );
    //    //
    //    dobj_dcurveparams.setZero();
    //    dobj_dcurveparams ( 0, 0 ) =  _sqrt_weight_for_scale;
    //    dobj_dcurveparams ( 1, 1 ) =  _sqrt_weight_for_scale;
    //    dobj_dcurveparams ( 2, 2 ) =  _sqrt_weight_for_translation;
    //    dobj_dcurveparams ( 3, 3 ) =  _sqrt_weight_for_translation;
    //    dobj_dcurveparams ( 4, 0 ) =  + _sqrt_weight_for_similarity;
    //    dobj_dcurveparams ( 4, 1 ) =  - _sqrt_weight_for_similarity;

    //    // force user to call set_curves() each time;
    //    _forget_curves();
    //}


    //bool
    //GlobFitObjective_ScaleAndTranslationOnlyRegularizer::
    //_is_all_data_provided() {
    //    const double tol = 1e-12;
    //    bool answer = true;
    //    answer = answer && GlobFitObjective::_is_all_data_provided();
    //    answer = answer && ( _cached_curves.size() == 1 );
    //    answer = answer && ( _sqrt_weight_for_scale >= -tol );
    //    answer = answer && ( _sqrt_weight_for_translation >= -tol );
    //    answer = answer && ( _sqrt_weight_for_similarity >= -tol );
    //    return answer;
    //}

// =============================================================
//                Regularization
// =============================================================

	void GlobFitObjective_Barrier::set_params(double soft_limit, double hard_limit)
	{
		this->soft_limit = soft_limit;
		this->hard_limit = hard_limit;

		double gap = (soft_limit - hard_limit);

		this->set_weight(1e8 / (gap * gap * gap * gap));
	}

	void GlobFitObjective_ParameterBound::set_params(GlobFitCurveParametrization* curve, int parameter, double soft_limit, double hard_limit)
	{
		_curves.clear();
		_curves.push_back(curve);

		this->parameter = parameter;
		GlobFitObjective_Barrier::set_params(soft_limit, hard_limit);
	}

	void GlobFitObjective_ParameterBound::compute_objective_and_jacobian(Eigen::VectorXd & obj, Eigen::MatrixXd & dobj_dcurveparams) const
	{
		obj.resize(n_equations());		

		const auto& params = _curves.front()->get_params();
		double value = params(parameter);
		Eigen::RowVectorXd value_jacobian(_n_sum_curve_params());
		value_jacobian.setZero();
		value_jacobian(parameter) = 1;
		
		obj(0) = this->compute_barrier(value, value_jacobian, dobj_dcurveparams);		
	}

	void GlobFitObjective_LineAngle::set_params(GlobFitLineParametrization * line, double angle)
	{
		_curves.clear();
		_curves.push_back(line);

		this->angle = angle;
	}

	void GlobFitObjective_LineAngle::compute_objective_and_jacobian(Eigen::VectorXd & obj, Eigen::MatrixXd & dobj_dcurveparams) const
	{
		auto curve = _curves.front();

		auto dir = curve->get_curve()->pos(1.0) - curve->get_curve()->pos(0.0);
		auto ddir_dparams = curve->dposdparams(1.0) - curve->dposdparams(0.0);

		auto angle = std::atan2(dir.y(), dir.x());
		Eigen::Matrix<double, 1, 2> atan2_jac;
		atan2_jac << -dir.y() / dir.squaredNorm(), dir.x() / dir.squaredNorm();
		auto dangle_dparams = atan2_jac * ddir_dparams;

		obj.resize(this->n_equations());
		obj(0) = angle - this->angle;

		dobj_dcurveparams.resize(this->n_equations(), curve->n_params());
		dobj_dcurveparams = dangle_dparams;
	}

	void GlobFitObjective_LineRegularization::set_params(GlobFitCurveParametrization * curve)
	{
		_curves.clear();
		_curves.push_back(curve);

		GlobFitObjective_Barrier::set_params(0.001, 0.0);
	}

	void GlobFitObjective_LineRegularization::compute_objective_and_jacobian(Eigen::VectorXd & obj, Eigen::MatrixXd & dobj_dcurveparams) const
	{
		obj.resize(n_equations());

		auto curve = _curves.front();

		auto dir = curve->get_curve()->pos(1.0) - curve->get_curve()->pos(0.0);
		auto ddir_dparams = curve->dposdparams(1.0) - curve->dposdparams(0.0);

		auto length_sq = dir.squaredNorm();
		Eigen::RowVectorXd length_sq_jac(2);
		length_sq_jac << 2 * dir.x(), 2 * dir.y();
		auto dlength_sq_dparams = length_sq_jac * ddir_dparams;
		
		obj(0) = this->compute_barrier(length_sq, dlength_sq_dparams, dobj_dcurveparams);
	}

} // end of polyvec