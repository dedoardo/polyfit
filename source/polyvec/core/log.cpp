#include <polyvec/core/log.hpp>

// libc++
#include <cstring>
#include <cstdarg>
#include <cstdio>

using namespace std;

namespace {
	struct Stream {
		~Stream() {
			if (fp && fp != stdout && fp != stderr) {
				fclose(fp);
			}
		}

		FILE* fp = nullptr;
	};

	Stream stream;
	int stream_flush = 0;
	int stream_channels = 0;
}

#define FP_FLUSH_COUNT 0

namespace polyfit {
	namespace Log {
		void open(const char* uri, const int channels) {
			if ((FILE*)uri == stdout || (FILE*)uri == stderr) {
				stream.fp = (FILE*)uri;
			} 
			else if (strcmp(uri, "stdout") == 0) {
				stream.fp = stdout;
			}
			else if (strcmp(uri, "stderr") == 0) {
				stream.fp = stderr;
			}
			else if (uri) {
				stream.fp = fopen(uri, "w");
			}

			stream_channels = channels;
		}

		void open(void* uri, const int channels) {
			open((char*)uri, channels);
		}

		void write(const char* fmt, ...) {
			if (stream.fp) {
				va_list args;
				va_start (args, fmt);
				vfprintf(stream.fp, fmt, args);
				if (stream_flush >= FP_FLUSH_COUNT) {
					fflush(stream.fp);
					stream_flush = 0;
				}
				else {
					++stream_flush;
				}

				va_end(args);
			}
		}

		void write_channel(const int channel, const char* fmt, ...) {
			if (!(stream_channels & channel)) {
				return;
			}

			if (stream.fp) {
				va_list args;
				va_start(args, fmt);
				vfprintf(stream.fp, fmt, args);
				if (stream_flush >= FP_FLUSH_COUNT) {
					fflush(stream.fp);
					stream_flush = 0;
				}
				else {
					++stream_flush;
				}

				va_end(args);
			}
		}
	}
}