// polyvec
#include <polyvec/utils/arguments.hpp>
#include <polyvec/core/options.hpp>
#include <polyvec/core/log.hpp>

// libc++
#include <string>
#include <fstream>

#define PVEC_READ_INPUT "D:/data/polyvec/svn/binary-perfect/dumbbell/32.png"
#define PVEC_READ_OUTPUT "D:/data/polyvec/polygon.svg"
#define PVEC_LOGGING stdout

using namespace polyvec;
using namespace std;

NAMESPACE_BEGIN ( pvec )

void load_options ( int argc, char** argv, bool silent ) {
    ifstream fs ( "argv" );

    if ( !fs.good() ) {
        if ( !silent ) {
            fprintf ( stderr, "./argv file not found, using default options\n" );
        }

        return;
    }

    string line;

    while ( std::getline ( fs, line ) ) {

    }
}

void load_options_default() {
    VectorOptions opt;
    opt.make_global();
    opt.make_thread_local();
}

char* read_input ( int argc, char** argv, bool silent ) {
#ifdef PVEC_READ_INPUT
    return PVEC_READ_INPUT;
#else

    if ( argc < 4 ) {
        if ( !silent ) {
            fprintf ( stderr, "expected <input> <output> <log> but only %d arguments.\n", argc );
        }

        return nullptr;
    }

    return argv[1];
#endif
}

char* read_output ( int argc, char** argv, bool silent ) {
#ifdef PVEC_READ_INPUT
    return PVEC_READ_OUTPUT;
#else

    if ( argc < 4 ) {
        if ( !silent ) {
            fprintf ( stderr, "expected <input> <output> <log> but only %d arguments.\n", argc );
        }

        return nullptr;
    }

    return argv[2];
#endif
}

void init_logging ( int argc, char** argv, bool silent ) {
#ifdef PVEC_LOGGING
    polyfit::Log::open ( ( char* ) PVEC_LOGGING );
#else

    if ( argc < 4 ) {
        if ( !silent ) {
            fprintf ( stderr, "expected <input> <output> <log> but only %d arguments.\n", argc );
        }

        return;
    }

    polyfit::Log::open ( argv[3] );
#endif
}

NAMESPACE_END ( pvec )