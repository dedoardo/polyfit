#include <polyvec/curve-tracer/curve_bezier.hpp>

#include <random>
#include <iostream>

int main()
{
	std::mt19937 rnd;
	std::uniform_real_distribution<> c_dist(-10, 10);
	std::uniform_real_distribution<> t_dist(0.0, 1.0);

	for (int i = 0; i < 10; ++i)
	{
		//create a random Bezier curve
		Eigen::Matrix2Xd controls(2, 4);
		for (int j = 0; j < controls.size(); ++j)
			controls(j) = c_dist(rnd);

		polyvec::BezierCurve original_curve(controls);

		//split the curve
		double split_t = t_dist(rnd);
		auto split = original_curve.split(split_t);

		//check some points
		for (double t = 0; t <= 1.0; t += 0.01)
		{
			auto t_after_split = polyvec::GlobFitCurve::split_t(t, split_t);

			auto split_p = t_after_split.first < 1.0 ? split.first->pos(t_after_split.first) : split.second->pos(t_after_split.second);
			auto original_p = original_curve.pos(t);

			if ((split_p - original_p).squaredNorm() > 1e-8)
			{
				std::cout << "Incorrect reconstruction at t = " << t << " (left t = " << t_after_split.first << ", right t = " << t_after_split.second << "):" << std::endl
					<< original_p.transpose() << "  !=  " << split_p.transpose() << std::endl
					<< "split t = " << split_t << std::endl;
				return EXIT_FAILURE;
			}
		}

		delete split.first;
		delete split.second;
	}

	std::cout << "All tests succeeded." << std::endl;
	return EXIT_SUCCESS;
}