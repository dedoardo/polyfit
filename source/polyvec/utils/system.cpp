// Header
#include <polyvec/utils/system.hpp>

// Polyvec
#include <polyvec/api.hpp>
#include <polyvec/misc.hpp>
#include <polyvec/debug.hpp>

// libc
#include <ctime>
#include <iomanip>
#include <stdlib.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#include <Lmcons.h>
#include <shellapi.h>

namespace tronkko {
#pragma warning(push)
#pragma warning(disable : 4505)
#pragma warning(disable : 4267)
#include <dirent.h>
#pragma warning(pop)
}

namespace polyvec {

    namespace {
        thread_local WORD text_attr;
    }

    namespace os {
        std::string user() {
            char username[UNLEN + 1];
            DWORD username_len = UNLEN + 1;
            GetUserNameA ( username, &username_len );
            return username;
        }

        std::string time_str() {
            auto t = std::time ( nullptr );
            auto tm = *std::localtime ( &t );

            std::ostringstream oss;
            oss << std::put_time ( &tm, "%d-%m-%Y-at-%H-%M-%S" );
            return oss.str();
        }

        int exec_shell ( const std::string& cmd, const std::string& args ) {
            SHELLEXECUTEINFOA info = { 0 };
            info.cbSize = sizeof ( SHELLEXECUTEINFOA );
            info.fMask = SEE_MASK_NOCLOSEPROCESS;
            info.hwnd = NULL;
            info.lpVerb = NULL;
            info.lpFile = cmd.c_str();
            info.lpParameters = args.c_str();
            info.lpDirectory = NULL;
            info.nShow = SW_HIDE;
            info.hInstApp = NULL;

            dbg::info ( FMT ( "executing cmd %s args %s", cmd.c_str(), args.c_str() ) );

            if ( !ShellExecuteExA ( &info ) ) {
                return -1;
            }

            WaitForSingleObject ( info.hProcess, INFINITE );
            DWORD exit_code;
            GetExitCodeProcess ( info.hProcess, &exit_code );
            CloseHandle ( info.hProcess );
            return ( int ) exit_code;
        }

        void copy_file ( const std::string& src, const std::string& dst ) {
            CopyFileA ( src.c_str(), dst.c_str(), false );
        }

        void make_dir ( const std::string& path ) {
			static const std::string separators("\\/");

			DWORD fileAttributes = ::GetFileAttributesA(path.c_str());
			if (fileAttributes == INVALID_FILE_ATTRIBUTES) {

				// Recursively do it all again for the parent directory, if any
				std::size_t slashIndex = path.find_last_of(separators);
				if (slashIndex != std::string::npos) {
					make_dir(path.substr(0, slashIndex));
				}
			}

            CreateDirectoryA ( path.c_str(), nullptr );
        }

        bool file_remove ( const std::string& path ) {
            return DeleteFileA ( path.c_str() );
        }

        bool directory_exists ( const std::string& path ) {
            DWORD file_attr = GetFileAttributesA ( path.c_str() );
            return file_attr != INVALID_FILE_ATTRIBUTES &&
                   ( file_attr & FILE_ATTRIBUTE_DIRECTORY );
        }

        void debug_window_output ( const std::string& str ) {
            OutputDebugStringA ( "APPLICATION: " );
            OutputDebugStringA ( str.c_str() );
            OutputDebugStringA ( "\n" );
        }

        TimeSlice timer_split() {
            LARGE_INTEGER ts;
            QueryPerformanceCounter ( &ts );
            return ( TimeSlice ) ts.QuadPart;
        }

        double timer_elapsed_ms ( TimeSlice delta ) {
            thread_local static bool frequency_init = false;
            thread_local static LARGE_INTEGER frequency;

            if ( !frequency_init ) {
                QueryPerformanceFrequency ( &frequency );
            }

            return ( ( double ) delta / frequency.QuadPart ) * 1000;
        }

        Lock cs_create() {
            CRITICAL_SECTION* cs = new CRITICAL_SECTION;
            InitializeCriticalSection ( cs );
            return cs;
        }

        void cs_destroy ( Lock cs ) {
            delete cs;
        }

        void cs_lock ( Lock cs ) {
            assert_break ( cs );
            EnterCriticalSection ( ( LPCRITICAL_SECTION ) cs );
        }

        void cs_unlock ( Lock cs ) {
            assert_break ( cs );
            LeaveCriticalSection ( ( LPCRITICAL_SECTION ) cs );
        }

        void* aligned_malloc ( size_t bytes, size_t alignment ) {
            return _aligned_malloc ( bytes, alignment );
        }

        void  aligned_free ( void* ptr ) {
            _aligned_free ( ptr );
        }

        void list_image_files_recursive ( std::vector<str>& addrs, const char* addr ) {
#ifdef _WIN32
            using namespace tronkko;
#endif

            DIR* dir;
            struct dirent* entry;

            dir = opendir ( addr );

            if ( !dir ) {
                return;
            }

            while ( ( entry = readdir ( dir ) ) != NULL ) {
                if ( entry->d_type == DT_DIR ) {
                    char path[1024];

                    if ( strcmp ( entry->d_name, "." ) == 0 ||
                            strcmp ( entry->d_name, ".." ) == 0 ||
                            strcmp ( entry->d_name, "_" ) == 0 ) {
                        continue;
                    }

                    snprintf ( path, sizeof ( path ), "%s/%s", addr, entry->d_name );
                    list_image_files_recursive ( addrs, path );
                } else {
                    const char* exts[] = { "png", "jpg", "bmp" };
                    str img_addr = misc::sfmt ( "%s/%s", addr, entry->d_name );
                    char img_ext[3];
                    size_t img_addr_len = img_addr.size();

                    if ( img_addr_len > 3 ) {
                        memcpy ( img_ext, img_addr.data() + img_addr_len - 4, 3 );

                        for ( int iext = 0; iext < array_len ( exts ); ++iext ) {
                            if ( memcmp ( img_ext, exts[iext], 3 ) == 0 ) {
                                addrs.push_back ( img_addr );
                            }
                        }
                    }
                }
            }

            closedir ( dir );
        }

		void list_txt_files_recursive(std::vector<str>& addrs, const char* addr) {
#ifdef _WIN32
			using namespace tronkko;
#endif

			DIR* dir;
			struct dirent* entry;

			dir = opendir(addr);

			if (!dir) {
				return;
			}

			while ((entry = readdir(dir)) != NULL) {
				if (entry->d_type == DT_DIR) {
					char path[1024];

					if (strcmp(entry->d_name, ".") == 0 ||
						strcmp(entry->d_name, "..") == 0 ||
						strcmp(entry->d_name, "_") == 0) {
						continue;
					}

					snprintf(path, sizeof(path), "%s/%s", addr, entry->d_name);
					list_image_files_recursive(addrs, path);
				}
				else {
					const char* exts[] = { "txt" };
					str img_addr = misc::sfmt("%s/%s", addr, entry->d_name);
					char img_ext[3];
					size_t img_addr_len = img_addr.size();

					if (img_addr_len > 3) {
						memcpy(img_ext, img_addr.data() + img_addr_len - 4, 3);

						for (int iext = 0; iext < array_len(exts); ++iext) {
							if (memcmp(img_ext, exts[iext], 3) == 0) {
								addrs.push_back(img_addr);
							}
						}
					}
				}
			}

			closedir(dir);
		}
    }
}
#elif defined(polyvec_posix)

// linux
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>

// libc
#include <cstdio>

namespace polyvec {
    namespace os {

        std::string user() {
            return getenv ( "USER" );
        }

        std::string time_str() {
            assert_break ( false );
        }

        int exec_shell ( const std::string& cmd, const std::string& args ) {
            return system ( ( cmd + " " + args ).c_str() );
        }

        bool file_remove ( const std::string& path ) {
            return exec_shell ( "rm -f", path.c_str() );
        }


        void copy_file ( const std::string& src, const std::string& dst ) {
            FILE* fp_r = fopen ( src.c_str(), "rb" );
            FILE* fp_w = fopen ( dst.c_str(), "wb" );

            char buf[4096];
            size_t bytes = 0;

            while ( ( bytes = fread ( buf, 1, array_len ( buf ), fp_r ) ) ) {
                fwrite ( buf, 1, bytes, fp_w );
            }

            fclose ( fp_r );
            fclose ( fp_w );
        }

        void make_dir ( const std::string& path ) {
            constexpr int MAX_STRING_LEN = 4096;
            char tmp[MAX_STRING_LEN];
            char* p = NULL;
            size_t len;

            snprintf ( tmp, sizeof ( tmp ), "%s", path.c_str() );
            len = strlen ( tmp );

            if ( tmp[len - 1] == '/' ) {
                tmp[len - 1] = 0;
            }

            for ( p = tmp + 1; *p; p++ )
                if ( *p == '/' ) {
                    *p = 0;
                    mkdir ( tmp, S_IRWXU );
                    *p = '/';
                }

            mkdir ( tmp, S_IRWXU );

        }
        bool directory_exists ( const std::string& path ) {
            DIR* dir = opendir ( path.c_str() );

            if ( dir ) {
                closedir ( dir );
                return true;
            }

            assert_break ( errno == ENOENT );
            return false;
        }
        TimeSlice timer_split() {
            assert_break ( false );
        }

        double timer_elapsed_ms ( TimeSlice delta ) {
            assert_break ( false );
        }

        Lock cs_create() {
            assert_break ( false );
            return NULL;
        }

        void cs_destroy ( Lock cs ) {
            assert_break ( false );
        }

        void cs_lock ( Lock cs ) {
            assert_break ( false );
        }

        void cs_unlock ( Lock cs ) {
            assert_break ( false );
        }

        void* aligned_malloc ( size_t bytes, size_t alignment ) {
            void* ptr;
            posix_memalign ( &ptr, alignment, bytes );
            return ptr;
        }

        void aligned_free ( void* ptr ) {
            free ( ptr );
        }

        void list_image_files_recursive ( std::vector<str>& addrs, const char* addr ) {
            DIR* dir;
            struct dirent* entry;

            dir = opendir ( addr );

            if ( !dir ) {
                return;
            }

            while ( ( entry = readdir ( dir ) ) != NULL ) {
                if ( entry->d_type == DT_DIR ) {
                    char path[1024];

                    if ( strcmp ( entry->d_name, "." ) == 0 ||
                            strcmp ( entry->d_name, ".." ) == 0 ||
                            strcmp ( entry->d_name, "_" ) == 0 ) {
                        continue;
                    }

                    snprintf ( path, sizeof ( path ), "%s/%s", addr, entry->d_name );
                    list_image_files_recursive ( addrs, path );
                } else {
                    const char* exts[] = { "png", "jpg", "bmp" };
                    str img_addr = misc::sfmt ( "%s/%s", addr, entry->d_name );
                    char img_ext[3];
                    size_t img_addr_len = img_addr.size();

                    if ( img_addr_len > 3 ) {
                        memcpy ( img_ext, img_addr.data() + img_addr_len - 4, 3 );

                        for ( int iext = 0; iext < array_len ( exts ); ++iext ) {
                            if ( memcmp ( img_ext, exts[iext], 3 ) == 0 ) {
                                addrs.push_back ( img_addr );
                            }
                        }
                    }
                }
            }

            closedir ( dir );
        }

        void list_txt_files_recursive ( std::vector<str>& addrs, const char* addr ) {
            DIR* dir;
            struct dirent* entry;

            dir = opendir ( addr );

            if ( !dir ) {
                return;
            }

            while ( ( entry = readdir ( dir ) ) != NULL ) {
                if ( entry->d_type == DT_DIR ) {
                    char path[1024];

                    if ( strcmp ( entry->d_name, "." ) == 0 ||
                            strcmp ( entry->d_name, ".." ) == 0 ||
                            strcmp ( entry->d_name, "_" ) == 0 ) {
                        continue;
                    }

                    snprintf ( path, sizeof ( path ), "%s/%s", addr, entry->d_name );
                    list_image_files_recursive ( addrs, path );
                } else {
                    const char* exts[] = { "txt" };
                    str img_addr = misc::sfmt ( "%s/%s", addr, entry->d_name );
                    char img_ext[3];
                    size_t img_addr_len = img_addr.size();

                    if ( img_addr_len > 3 ) {
                        memcpy ( img_ext, img_addr.data() + img_addr_len - 4, 3 );

                        for ( int iext = 0; iext < array_len ( exts ); ++iext ) {
                            if ( memcmp ( img_ext, exts[iext], 3 ) == 0 ) {
                                addrs.push_back ( img_addr );
                            }
                        }
                    }
                }
            }

            closedir ( dir );
        }
    }
}
#endif

namespace polyvec {
    namespace os {
        bool file_exists ( const std::string& addr ) {
            FILE* fp_try = fopen ( addr.c_str(), "r" );

            if ( !fp_try ) {
                return false;
            }

            fclose ( fp_try );
            return true;
        }
    }
}