/*
	Utility functions for loading options from the command line and file.
	
	The command line requires the following positional arguments
	<input> <output> <log>

	If an 'argv' file is present in the working directory, one option per
	line is read with the following format.

	\s+option\s+=\s+value

	values are read as numbers and rounded if an integer is expected
*/
#pragma once

#define PVEC_OPT_ACCURACY_WEIGHT "accuracy-weight"
#define PVEC_OPT_

// map to polygon-tracer/options.hpp
#define OPT_SMOOTHNESS_LIMIT "smoothness-limit"
#define OPT_CONTINUITY_LIMIT "continuity-limit"
#define OPT_INFLECTION_LIMIT "inflection-limit"
#define OPT_INFLECTION_PENALTY "inflection-penalty"
#define OPT_SMOOTHNESS_WEIGHT "smoothness-weight"
#define OPT_ACCURACY_WEIGHT "accuracy-weight"
#define OPT_CONTINUITY_WEIGHT "continuity-weight"
#define OPT_INFLECTION_WEIGHT "inflection-weight"

#include <polyvec/core/macros.hpp>

NAMESPACE_BEGIN(pvec)

// reads one options per line in name=value format from an argv file 
// located in the current working directory
void  load_options(int argc, char** argv, bool silent = false);

// loads the default options for all the modules
void  load_options_default();

char* read_input(int argc, char** argv, bool silent = false);
char* read_output(int argc, char** argv, bool silent = false);

// redirects the logging to the specified directry
void init_logging(int argc, char** argv, bool silent = false);

NAMESPACE_END(pvec)
