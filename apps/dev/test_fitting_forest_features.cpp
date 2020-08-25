// Polyvec
#include <polyvec/curve-tracer/fit_classifier.hpp>

// libc++
#include <cstdlib>
#include <random>

using namespace polyvec;

int main(int argc, char* argv[]) {
	std::mt19937 rnd;
	std::uniform_real_distribution<double> dist(0., 1.);

	FitClassifierRandomForest::Sample sample;
	for (int i = 0; i < 4; ++i) {
		sample.fits[i].distance_to_corner = dist(rnd);
		sample.fits[i].error.accuracy.e_pos = dist(rnd);
		sample.fits[i].error.accuracy.e_neg = dist(rnd);
		sample.fits[i].error.curvature.r_min = dist(rnd);
		sample.fits[i].error.curvature.r_max = dist(rnd);
	}

	sample.corner.accuracy.e_pos = dist(rnd);
	sample.corner.accuracy.e_neg = dist(rnd);
	sample.corner.dist_prev = dist(rnd);
	sample.corner.dist_next = dist(rnd);
	sample.corner.angle = 3.1415927 * dist(rnd);


	FitClassifierRandomForest::Sample sample_cp = sample;

	sample.corner.convexity[0] = 1;
	sample.label = 0;
	sample_cp.corner.convexity[0] = -1;
	sample_cp.label = 1;

	FitClassifierRandomForest::SampleSet samples;
	samples.samples.emplace_back(sample);
	samples.samples.emplace_back(sample_cp);
	samples.samples.emplace_back(sample);
	samples.samples.emplace_back(sample_cp);

	//TODO: Add polygon data

	FitClassifierRandomForest classifier;
	classifier.train(samples, "model.txt");

	return EXIT_SUCCESS;
}
