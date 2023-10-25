#pragma once

#if defined(_MSC_VER) && !defined(__clang__)
#define MSVC_WARNING_SUPPRESS(XXXX) __pragma(warning(suppress : XXXX))
#else
#define MSVC_WARNING_SUPPRESS(XXXX)
#endif
