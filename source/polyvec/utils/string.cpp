#include <polyvec/utils/string.hpp>

// libc
#include <cstdarg>
#include <algorithm>

using namespace std;

NAMESPACE_BEGIN(polyfit)
NAMESPACE_BEGIN(StringUtils)

std::string fmt(const char* fmt, ...) {
	va_list args, args_cp;
	va_start(args, fmt);
	va_copy(args_cp, args);
	int len = vsnprintf(nullptr, 0, fmt, args);
	va_end(args);
	string ret;
	ret.resize(len + 1);
	PF_ASSERT(vsnprintf(&ret[0], len + 1, fmt, args_cp));
	va_end(args_cp);
	return ret;
}

int replace_bounded(char* s, const char* token, const char* fmt, ...) {
	PF_ASSERT(s && token && fmt);

	const size_t token_len = strlen(token);

	va_list args, args_cp;
	va_start(args, fmt);
	va_copy(args_cp, args);
	const size_t token_new_len = min((size_t)vsnprintf(NULL, 0, fmt, args), token_len);
	va_end(args);
	std::vector<char> token_new;
	token_new.resize(token_new_len + 1);
	snprintf(token_new.data(), token_new_len, fmt, args_cp);
	va_end(args_cp);

	char* w = s;
	int replaced = 0;

	while (*w != '\0') {
		if (strcmp(w, token) == 0) {
			memcpy(w, token_new.data(), token_new_len);
			w += token_len;
			++replaced;
		}
	}

	return replaced;
}

std::string join_path(const char* p0, const char* p1, const char* ext) {
	PF_ASSERT(p0 && p1);
	
	std::string p = p0;
	if (p.back() != '/' && p.back() != '\\')
		p.append("/");

	if (p1[0] == '/' || p1[0] == '\\')
		++p1;

	p.append(p1);

	if (ext) {
		if (ext[0] != '.' && p.back() != '.') {
			p += '.';
		}
	
		p += ext;
	}

	return p;
}

NAMESPACE_END(StringUtils)
NAMESPACE_END(polyfit)