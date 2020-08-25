/*
    OS-specific utilities
*/
#ifndef polyvec_system_h_
#define polyvec_system_h_

// Polyvec
#include <polyvec/api.hpp>

namespace polyvec {

    // Platform specific code
    namespace os {
        // Misc routines
        std::string  user();                                         // os username
        std::string  time_str();                                     // formatted time string
        int  exec_shell  ( const std::string& cmd, const std::string& args ); // executes shell program
        void copy_file   ( const std::string& src, const std::string& dst );
        void make_dir    ( const std::string& path );                   // create directory if it doesn't exist
        bool file_remove ( const std::string& path );
        bool directory_exists ( const std::string& path );
        bool file_exists ( const std::string& addr );

        // non interprocess lock (critical section on windows)
        using Lock = void*;
        Lock cs_create();
        void cs_destroy ( Lock cs );
        void cs_lock ( Lock cs );
        void cs_unlock ( Lock cs );

        // High resolution clock. Usage
        // `start = timer_split()`
        // `end = timer_split()`
        // `milliseconds = timer_elapsed_ms(end - start)`
        using TimeSlice = uint64_t;
        TimeSlice timer_split();
        double    timer_elapsed_ms ( TimeSlice );

        void* aligned_malloc ( size_t bytes, size_t alignment );
        void  aligned_free ( void* ptr );

        void list_image_files_recursive ( std::vector<str>& addrs, const char* addr );
        void list_txt_files_recursive ( std::vector<str>& addrs, const char* addr );
    }
}

#endif // polyvec_system_h_