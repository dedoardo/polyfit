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
#include <dirent.h>
}
#else
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>

#include <cstdio>
#endif

namespace polyvec {
    namespace os {
#ifdef _WIN32
        void make_dir(const std::string& path) {
            static const std::string separators("\\/");

            DWORD fileAttributes = ::GetFileAttributesA(path.c_str());
            if (fileAttributes == INVALID_FILE_ATTRIBUTES) {

                // Recursively do it all again for the parent directory, if any
                std::size_t slashIndex = path.find_last_of(separators);
                if (slashIndex != std::string::npos) {
                    make_dir(path.substr(0, slashIndex));
                }
            }

            CreateDirectoryA(path.c_str(), nullptr);
        }

        void copy_file(const std::string& src, const std::string& dst) {
            CopyFileA(src.c_str(), dst.c_str(), false);
        }
#else
        void make_dir(const std::string& path) {
            constexpr int MAX_STRING_LEN = 4096;
            char tmp[MAX_STRING_LEN];
            char* p = NULL;
            size_t len;

            snprintf(tmp, sizeof(tmp), "%s", path.c_str());
            len = strlen(tmp);

            if (tmp[len - 1] == '/') {
                tmp[len - 1] = 0;
            }

            for (p = tmp + 1; *p; p++)
                if (*p == '/') {
                    *p = 0;
                    mkdir(tmp, S_IRWXU);
                    *p = '/';
                }

            mkdir(tmp, S_IRWXU);
        }


        void copy_file(const std::string& src, const std::string& dst) {
            FILE* fp_r = fopen(src.c_str(), "rb");
            FILE* fp_w = fopen(dst.c_str(), "wb");

            char buf[4096];
            size_t bytes = 0;

            while ((bytes = fread(buf, 1, array_len(buf), fp_r))) {
                fwrite(buf, 1, bytes, fp_w);
            }

            fclose(fp_r);
            fclose(fp_w);
        }
#endif

        void list_image_files_recursive(std::vector<str>& addrs, const char* addr) {
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
                    const char* exts[] = { "png", "jpg", "bmp" };
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

        bool file_exists(const std::string& addr) {
            FILE* fp_try = fopen(addr.c_str(), "r");

            if (!fp_try) {
                return false;
            }

            fclose(fp_try);
            return true;
        }
    }
}
