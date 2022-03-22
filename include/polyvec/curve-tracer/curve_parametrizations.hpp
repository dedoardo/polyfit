#pragma once

#include <polyvec/curve-tracer/curve_parametrization.hpp>

NAMESPACE_BEGIN(polyvec)

//Represents an orthogonal coordinate system
class CoordinateSystem
{
public:
	//Sets a new primary axis and updates the secondary one to be orthogonal
	void set_primary(const Eigen::Vector2d&);

	//Sets a new secondary axis and updates the primary one to be orthogonal
	void set_secondary(const Eigen::Vector2d&);

	Eigen::Matrix2d::ConstColXpr primary() const { return _matrix.col(0); }
	Eigen::Matrix2d::ConstColXpr secondary() const { return _matrix.col(1); }
	const Eigen::Matrix2d& matrix() const { return _matrix; }
private:
	Eigen::Matrix2d _matrix;
};

class GlobFitBezierAngleBasedParametrization : public GlobFitCurveParametrization {
public:
	GlobFitBezierAngleBasedParametrization(
        std::shared_ptr<BezierCurve> curve, 
        const bool invert_front_coordinate_system = false,
        const bool invert_back_coordinate_system = false);

	int n_params() const override { return 8; }

	void set_params(const Eigen::VectorXd& params) override;
	Eigen::VectorXd get_params() const override;
    void reverse()override;

	GlobFitCurveParametrization* clone() override;
	GlobFitCurveParametrization* create_for_curve(std::shared_ptr<GlobFitCurve> curve) const override;

	void reduce_degrees_of_freedom(DofOptions::Type options) override;

private:
	struct FittingParams {
		double beg_tangent_len;
		double end_tangent_len;
		//endpoints are represented in the respective local coordinate system
		Eigen::Vector2d beg_pt;
		Eigen::Vector2d end_pt;
		double beg_tangent_angle;
		double end_tangent_angle;
	};
	FittingParams params;

	//the initial endpoints, origins of the local coordinate systems
	Eigen::Vector2d startP, endP;

    CoordinateSystem front_system, back_system;
    bool front_system_inverted, back_system_inverted;

	void set_curve_parameters();
	void calculate_jacobian();
};

//This parametrization of a Bezier curve has reduced degrees of freedom. The endtangents
//are fixed and endpoints are only allowed to move orthogonally to endtangents.
class GlobFitBezierShiftBasedParametrization : public GlobFitCurveParametrization {
public:
	GlobFitBezierShiftBasedParametrization(std::shared_ptr<BezierCurve> curve);
	int n_params() const { return 4; }

	void set_params(const Eigen::VectorXd& params);
	Eigen::VectorXd get_params() const;

	GlobFitCurveParametrization* clone();
	GlobFitCurveParametrization* create_for_curve(std::shared_ptr<GlobFitCurve> curve) const;
private:
	struct FittingParams {
		double beg_tangent_len;
		double end_tangent_len;
		double beg_shift;
		double end_shift;
	};
	FittingParams params;

	//unit tangents, fixed
	Eigen::Vector2d startTang, endTang;

	//the initial endpoints; the endpoints are allowed to move orthogonally to their endtangents
	Eigen::Vector2d startP, endP;

	void set_curve_parameters();
	void calculate_jacobian();
};

//Represents a parametrization of a line where the endpoints are represented in a custom coordinate system. The default coordinate
//system has its primary axis orthogonal to the original line and its secondary axis parallel to it.
class GlobFitLineParametrization : public GlobFitCurveParametrization {
public:
    GlobFitLineParametrization(
        std::shared_ptr<GlobFitCurve_Line> curve,
        const bool invert_front_coordinate_system = false,
        const bool invert_back_coordinate_system = false);

    int n_params() const override { return 4; }

    void set_params(const Eigen::VectorXd& params) override;
    Eigen::VectorXd get_params() const override;

    GlobFitCurveParametrization* clone();
    GlobFitCurveParametrization* create_for_curve(std::shared_ptr<GlobFitCurve> curve) const override;

    //Sets a new custom coordinate system for the front point of the line. If rescale is set to true, the axis is re-scaled such that
    //its rejection from the original line is of unit-length.
    void set_front_primary_axis(const Eigen::Vector2d& d, bool rescale);

    //Sets a new custom coordinate system for the back point of the line. If rescale is set to true, the axis is re-scaled such that
    //its rejection from the original line is of unit-length.
    void set_back_primary_axis(const Eigen::Vector2d& d, bool rescale);

    const CoordinateSystem& get_front_coordinate_system() const { return front_system; }
    const CoordinateSystem& get_back_coordinate_system() const { return back_system; }

    //Replaces the endpoint of this line with the endpoint of the other line.
    void merge_with(GlobFitLineParametrization* other);

    void reduce_degrees_of_freedom(DofOptions::Type options) override;
    void reverse()override;

    bool is_front_reversed() const { return front_system_inverted; }
    bool is_back_reversed() const { return back_system_inverted; }
private:
	struct FittingParams {
		//Offsets are specified in the local coordinate systems of the front and back point.
		Eigen::Vector2d front;
		Eigen::Vector2d back;
	};
	FittingParams params;

	//The direction orthogonal to the original line
	Eigen::Vector2d original_orthogonal_dir;

	//The vectors that specify in what directions we can move the front and back point (i.e., local coordinate system) 
	CoordinateSystem front_system, back_system;
    bool front_system_inverted, back_system_inverted;

	//the initial endpoints and origins of the local coordinate systems
	Eigen::Vector2d startP, endP;

	void set_curve_parameters();
	void calculate_jacobian();
};

NAMESPACE_END(polyvec)
