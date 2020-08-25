#pragma once

#include <polyvec/curve-tracer/fit_classifier.hpp>

// TPerSampleFunctor: void(const Sample&, const std::string& data_dir)
template <typename TPerSampleFunctor>
void load_samples_file(const std::string& path, FitClassifierRandomForest::SampleSet& samples, TPerSampleFunctor&& per_sample_functor)
{
	std::ifstream samples_file(path);
	if (!samples_file.good())
	{
		std::cerr << "Cannot read samples file " << path << std::endl;
		return;
	}

	const size_t STR_LEN = 255;
	char sample_classification_dir[STR_LEN];
	char sample_data_dir[STR_LEN];
	while (samples_file.getline(sample_classification_dir, STR_LEN) && samples_file.getline(sample_data_dir, STR_LEN))
	{
		int n_samples;
		samples_file >> n_samples;
		int samples_before = samples.samples.size();
		std::string sample_data_dir_with_slash = with_trailing_slash(sample_data_dir);
		for (int i = 0; i < n_samples; ++i)
		{
			std::string model;
			unsigned char label;
			samples_file >> model >> label;
			label -= '0';
			try {
				samples.samples.emplace_back(polyvec::FitClassifierRandomForest::Sample::read_from_data(sample_data_dir_with_slash, model, label));
			}
			catch (std::exception& e)
			{
				std::cerr << "Error reading model " << model << " from " << sample_data_dir_with_slash << ": " << e.what() << std::endl;
			}
		}
		auto samples_begin = samples.samples.begin() + samples_before;
		FitClassifierRandomForest::prepare_samples(samples, samples_begin, false, with_trailing_slash(sample_data_dir_with_slash));
		for (auto it = samples_begin; it != samples.samples.end(); ++it)
			std::forward<TPerSampleFunctor>(per_sample_functor)(*it, sample_data_dir_with_slash);
		samples_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	}
}