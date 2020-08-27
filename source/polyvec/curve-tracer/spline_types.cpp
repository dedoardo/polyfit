#include <polyvec/curve-tracer/spline_types.hpp>

#include <polyvec/misc.hpp>
#include <polyvec/curve-tracer/find-fitting-points.hpp>
#include <polyvec/utils/directions.hpp>
#include <polyvec/utils/curve_sampling.hpp>

#define TANGENT_SAMPLES_PER_UNIT_ARCLENGTH 8

//Adds linearly interpolated tangent samples for the given Bezier curve to tangents/ts.
//This is the new scheme
void polyvec::FittingInfo::add_tangent_samples(GlobFitCurveParametrization* curve, double tMiddle) {
	const int T_SAMPLES = 1000;	
	polyvec::Sampling::DistanceSampling sampling(*curve, T_SAMPLES);	
	double distanceMiddle = sampling.get_distance(tMiddle);	

	auto startTangent = curve->get_curve()->dposdt(0.0).normalized();
	auto endTangent = curve->get_curve()->dposdt(1.0).normalized();

	auto startAngle = std::atan2(startTangent.y(), startTangent.x());
	auto endAngle = std::atan2(endTangent.y(), endTangent.x());

	auto diff = endAngle - startAngle;

	if (diff > M_PI) {
		diff -= 2 * M_PI;
	}

	if (diff < -M_PI) {
		diff += 2 * M_PI;
	}

	auto middleAngle = startAngle + 0.5 * diff;
	auto middleTangent = misc::lerp(startTangent, endTangent, 0.5);

	sparse_tangents.fit_tangents.clear();
	sparse_tangents.fit_tangents.push_back(startTangent);
	sparse_tangents.fit_tangent_ts.push_back(0.0);

	sparse_tangents.fit_tangents.push_back(Eigen::Vector2d(std::cos(middleAngle), std::sin(middleAngle)));
	sparse_tangents.fit_tangent_ts.push_back(tMiddle);

	sparse_tangents.fit_tangents.push_back(endTangent);
	sparse_tangents.fit_tangent_ts.push_back(1.0);

	const int tangent_samples = std::ceil(sampling.total_length() * 2 * TANGENT_SAMPLES_PER_UNIT_ARCLENGTH);
	const double arcLengthPerSample = sampling.total_length() / tangent_samples;
	dense_tangents.fit_tangents.clear();
	dense_tangents.fit_tangents.reserve(tangent_samples);
	dense_tangents.fit_tangent_ts.reserve(tangent_samples);

	sampling.sample_by_distances(polyvec::Sampling::UniformIterator(0, arcLengthPerSample), polyvec::Sampling::UniformIterator(tangent_samples, arcLengthPerSample),
		[&](double t, double distance, const Eigen::RowVectorXd& dt_dparams)
	{
		double fitAngle;

		if (t < tMiddle) {
			fitAngle = misc::lerp(startAngle, middleAngle, (distance - 0) / (distanceMiddle - 0));
		}
		else {
			fitAngle = misc::lerp(middleAngle, startAngle + diff, (distance - distanceMiddle) / (sampling.total_length() - distanceMiddle));
		}


		dense_tangents.fit_tangent_ts.push_back(t);
		dense_tangents.fit_tangents.push_back(Eigen::Vector2d(std::cos(fitAngle), std::sin(fitAngle)));
	});
}

void polyvec::FittingInfo::add_tangent_samples(GlobFitCurve_Line* curve) {
	/*Eigen::Matrix2Xd points = curve->get_points();
	Eigen::Vector2d tangent = (points.col(1) - points.col(0)).normalized();

	fit_tangents.push_back(tangent);
	fit_tangent_ts.push_back(0.5);*/
}

void polyvec::FittingInfo::add_midpoint_samples(const Eigen::Matrix2Xd& points, const Eigen::VectorXi& polygonV, const size_t firstEdge, const double firstEdgeT,

	const size_t lastEdge, const double lastEdgeT, bool circular) {
	polyfit::mat2xi fit;
	polyfit::mat2x p;
	polyfit::CurveTracer::find_pixel_centers_for_subpath(points, polygonV, firstEdge, firstEdgeT, lastEdge, lastEdgeT, fit, &p, circular);

	fit_midpoints.clear();
	fit_midpoint_normals.clear();

	for (int j = 0; j < fit.cols(); ++j) {
		fit_midpoints.emplace_back(p.col(j));
		fit_midpoint_normals.push_back(polyvec::util::normal_dir(points.col(fit(1, j)) - points.col(fit(0, j))));
	}
}

std::pair<polyvec::CurvePrimitive, polyvec::CurvePrimitive> polyvec::CurvePrimitive::split(double t_split) const
{
	std::pair<polyvec::CurvePrimitive, polyvec::CurvePrimitive> result;

	//split the underlying curves
	auto split_curves = curve->get_curve()->split(t_split);
	result.first.curve  = std::shared_ptr<polyvec::GlobFitCurveParametrization>(this->curve->create_for_curve(std::shared_ptr<polyvec::GlobFitCurve>(split_curves.first)));
	result.second.curve = std::shared_ptr<polyvec::GlobFitCurveParametrization>(this->curve->create_for_curve(std::shared_ptr<polyvec::GlobFitCurve>(split_curves.second)));
	
	//redistribute the midpoint constraints
	for (int i = 0; i < fitting_info.fit_midpoints.size(); ++i)
	{
		const double t = curve->get_curve()->project(fitting_info.fit_midpoints.at(i));

		if (t <= t_split)
		{
			result.first.fitting_info.fit_midpoints.push_back(fitting_info.fit_midpoints.at(i));
			result.first.fitting_info.fit_midpoint_normals.push_back(fitting_info.fit_midpoint_normals.at(i));
		}

		if (t >= t_split)
		{
			result.second.fitting_info.fit_midpoints.push_back(fitting_info.fit_midpoints.at(i));
			result.second.fitting_info.fit_midpoint_normals.push_back(fitting_info.fit_midpoint_normals.at(i));
		}
	}	

	auto redistribute_tangents = [&](const FittingInfo::Tangents& source_tangents, FittingInfo::Tangents& target_tangents_left, FittingInfo::Tangents& target_tangents_right)
	{
		for (int i = 0; i < source_tangents.fit_tangents.size(); ++i)
		{
			const double t = source_tangents.fit_tangent_ts.at(i);
			auto t_after_split = polyvec::GlobFitCurve::split_t(t, t_split);

			if (t_after_split.first <= 1.0)
			{
				target_tangents_left.fit_tangents.push_back(source_tangents.fit_tangents.at(i));
				target_tangents_left.fit_tangent_ts.push_back(t_after_split.first);
			}

			if (t_after_split.second >= 0.0)
			{
				target_tangents_right.fit_tangents.push_back(source_tangents.fit_tangents.at(i));
				target_tangents_right.fit_tangent_ts.push_back(t_after_split.second);
			}
		}
	};

	//redistribute the tangent constraints	
	redistribute_tangents(fitting_info.dense_tangents, result.first.fitting_info.dense_tangents, result.second.fitting_info.dense_tangents);
	redistribute_tangents(fitting_info.sparse_tangents, result.first.fitting_info.sparse_tangents, result.second.fitting_info.sparse_tangents);

	return result;
}

bool polyvec::CurvePrimitive::fits_same_corner(const CurvePrimitive & other, int polygon_corners) const
{
	if (corner == other.corner)
		return true;

	if (fits_next_corner && (corner + 1) % polygon_corners == other.corner)
		return true;

	if (other.fits_next_corner && (other.corner + 1) % polygon_corners == corner)
		return true;

	return false;
}
