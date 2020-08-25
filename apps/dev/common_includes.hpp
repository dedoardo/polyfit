// Polyvec includes
#include <polyvec/debug.hpp>
#include <polyvec/core/options.hpp>
#include <polyvec/utils/system.hpp>

// C++ STL
#include <fstream>
#include <string>

using namespace polyvec;

struct Arguments {
    str           verb;				 // Program to be executed
    str           verb_instructions; // extra instructions for verb, optional
    str           log_level;         // Log verbosity (verbose, info, warning, error)
    str           log_output;        // Output file address or std[out|err]
    str           read;              // Input file or directory (dependent on program)
    str           write;             // Output file or directory address
    VectorOptions opt;               // Vectorization options (see `dev/api.hpp`)

    bool          disable_develop_canvas; // Disable canvas even if in develop

    static Arguments& singleton() {
        static Arguments the_one_and_only;
        return the_one_and_only;
    }
};

static std::string with_trailing_slash(const std::string& str)
{
	if (str.back() == '/' || str.back() == '\\')
		return str;
	else
		return str + '/';
}

