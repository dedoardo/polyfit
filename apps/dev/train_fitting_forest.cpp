// Polyvec
#include <polyvec/curve-tracer/fit_classifier.hpp>
#include "common_includes.hpp"
#include "training.hpp"

// libc++
#include <fstream>

using namespace polyfit;
using namespace polyvec;

const char* usage =
"Required parameters: <trained_model> <classification dir> <data dir> [<classification dir> <data dir>, ...] \n"
"<classification dir> must contain text files with corner classifications. Can also be \"--samples\" to denote reading from a samples file.\n"
"<data dir>           must contain the training samples that the classifications are based on. If \"--samples\" was specified, the path to the samples file.\n"
"<trained_model>      must contain the destination filename where the trained model will be written to.\n"
"Any other *pair* of data/label directories are also loaded. \n"
;

void noop(const polyvec::FitClassifierRandomForest::Sample& sample, const std::string& data_dir) { }

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr, usage);
		return EXIT_FAILURE;
	}

	const std::string model_uri = argv[1];

	FitClassifierRandomForest classifier;	
	FitClassifierRandomForest::SampleSet samples;

	const int n_training_sets = (argc - 2) / 2;

	if (n_training_sets < 1)
	{
		fprintf(stderr, usage);
		return EXIT_FAILURE;
	}

	std::cout << "Training with " << n_training_sets << " training sets." << std::endl;

	std::ofstream sampleFile(model_uri + ".samples.txt");
	size_t samples_loaded = 0;
	for (int i = 0; i < n_training_sets; ++i) {
		const std::string classification_dir = argv[2 + i * 2 + 0];
		const std::string data_dir = argv[2 + i * 2 + 1];

		if (classification_dir == "--samples")
		{
			std::cout << "Reading samples from sample file" << std::endl;
			load_samples_file(data_dir, samples, noop);
		}
		else
		{
			if (!FitClassifierRandomForest::read_samples(with_trailing_slash(data_dir), with_trailing_slash(classification_dir), samples)) {
				std::cerr << "No samples found in " << data_dir << ", " << classification_dir << std::endl;
			}
			else {
				sampleFile << classification_dir << std::endl << data_dir << std::endl << samples.samples.size() - samples_loaded << std::endl;
				for (size_t j = samples_loaded; j < samples.samples.size(); ++j)
					sampleFile << samples.samples[j].model << " " << (int)samples.samples[j].label << std::endl;
			}
		}		
		std::cout << "Read " << samples.samples.size() - samples_loaded << " samples from this training set." << std::endl;
		samples_loaded = samples.samples.size();
	}	

	try
	{
		classifier.train(samples, model_uri);
	}
	catch (std::exception& e)
	{
		std::cout << "Error training forest: " << e.what() << std::endl;
	}

	return EXIT_SUCCESS;
}