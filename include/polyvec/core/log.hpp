#ifndef POLYFIT_LOG_H_
#define POLYFIT_LOG_H_

namespace polyfit {
	namespace Log {
		enum ChannelBits {
			CHANNEL_NONE = 0,
			CHANNEL_STATUS = 1 << 0,
			CHANNEL_VERBOSE = 1 << 1,
			CHANNEL_ERROR = 1 << 2,
            CHANNEL_DEV = 1 << 3
		};
#define PF_LOG_CHANNELS_DEFAULT (CHANNEL_ERROR)

		void open(const char* uri, const int channels = PF_LOG_CHANNELS_DEFAULT);
		void open(void* uri, const int channels = PF_LOG_CHANNELS_DEFAULT);
		void write(const char* fmt, ...);

		void write_channel(const int channel, const char* fmt, ...);
	}
}

#define PF_DEV_S(s) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_DEV, "%s\n", s);
#define PF_DEV_F(fmt, ...) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_DEV, fmt "\n", __VA_ARGS__);

#define PF_STATUS_S(s) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_STATUS, "%s\n", s);
#define PF_STATUS_F(fmt, ...) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_STATUS, fmt "\n", __VA_ARGS__);

#define PF_VERBOSE_S(s) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_VERBOSE, "%s\n", s);
#define PF_VERBOSE_F(fmt, ...) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_VERBOSE, fmt "\n", __VA_ARGS__);

#define PF_ERROR_S(s) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_ERROR, "%s\n", s);
#define PF_ERROR_F(fmt, ...) ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_ERROR, fmt "\n", __VA_ARGS__);

#define PF_LOGS(s) ::polyfit::Log::write(s "\n")
#define PF_LOGF(fmt, ...) ::polyfit::Log::write(fmt "\n", __VA_ARGS__)
#define PF_LOG_DIVIDER ::polyfit::Log::write_channel(::polyfit::Log::CHANNEL_VERBOSE, "----------------------------------------------\n");

#endif // POLYFIT_LOG_H_