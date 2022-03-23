#pragma once

#include <polyvec/core/macros.hpp>
#include <polyvec/core/types.hpp>
#include <polyvec/curve-tracer/curve.hpp>
#include <polyvec/curve-tracer/curve_parametrization.hpp>
#include <polyvec/misc.hpp>

NAMESPACE_BEGIN(polyvec)
NAMESPACE_BEGIN(Sampling)

class DistanceSampling
{
public:
	DistanceSampling(const polyvec::GlobFitCurveParametrization& curve, int n_samples)
		: curve(curve), distances(n_samples + 1), ddistances_dparams(n_samples + 1), h(1.0 / n_samples)
	{
		distances[0] = 0;
		ddistances_dparams[0].resize(curve.n_params());
		ddistances_dparams[0].setZero();

		for (int i = 1; i <= n_samples; ++i) {
			const Eigen::Vector2d diff = curve.get_curve()->pos(i * h) - curve.get_curve()->pos((i - 1) * h);
			double step = diff.norm();
			distances[i] = distances[i - 1] + step;

			const Eigen::Matrix2Xd ddiff_dparams = curve.get_curve()->dposdparams((i - 1) * h) - curve.get_curve()->dposdparams(i * h);
			ddistances_dparams[i] = ddistances_dparams[i - 1] + Eigen::RowVector2d(diff.x() / step, diff.y() / step) * ddiff_dparams;
		}
	}

	double get_distance(double t) const
	{
		int i = (int)std::floor(t / h);
		int j = (int)std::ceil(t / h);
		return polyvec::misc::lerp(distances[i], distances[j], (t - i * h) / h);
	}

	double total_length() const
	{
		return distances.back();
	}

	// Callback: void(double t, double distance, RowVector dt_dparams)
	template <typename DistanceIterator, typename Callback>
	void sample_by_distances(DistanceIterator distances_begin, DistanceIterator distances_end, Callback&& callback)
	{
		int iTSample = 0;
		for (auto it = distances_begin; it != distances_end; ++it) {
			auto desiredDistance = *it;

			// find the first sample that is behind the desired length
			while (iTSample < distances.size() && distances[iTSample] < desiredDistance) {
				++iTSample;
			}

			double t, distance;
			Eigen::RowVectorXd dt_dparams(1, curve.n_params());

			if (desiredDistance == 0) {
				t = 0;
				distance = 0;
				dt_dparams.setZero();
			}
			else {
				double lastDistance = distances[iTSample - 1];
				double thisDistance = distances[iTSample];
				double lastT = (iTSample - 1) * h;
				double thisT = iTSample * h;
				//solve (t - lastT) / (thisT - lastT) = (desiredDistance - lastDistance) / (thisDistance - lastDistance)
				distance = desiredDistance;
				auto diffDistance = thisDistance - lastDistance;
				t = (desiredDistance - lastDistance) / diffDistance * (thisT - lastT) + lastT;
				dt_dparams =
					(
						-ddistances_dparams[iTSample - 1] * diffDistance
						- (desiredDistance - lastDistance) * (ddistances_dparams[iTSample] - ddistances_dparams[iTSample - 1])
					) / (diffDistance * diffDistance) * (thisT - lastT);

				if (t > 1) {
					t = 1;
					dt_dparams.setZero();
				}
			}

			std::forward<Callback>(callback)(t, distance, dt_dparams);
		}
	}

private:
	const polyvec::GlobFitCurveParametrization& curve;
	std::vector<double> distances;
	std::vector<Eigen::RowVectorXd> ddistances_dparams;
	double h;
};

class UniformIterator
{
public:
	UniformIterator(int i, double step)
		: i(i), step(step)
	{ }

	double operator*() const { return i * step; }
	UniformIterator operator++() { ++i; return *this; }
	bool operator!=(const UniformIterator& other) const { return i != other.i || step != other.step; }

private:
	int i;
	double step;
};

NAMESPACE_END(Sampling)
NAMESPACE_END(polyvec)