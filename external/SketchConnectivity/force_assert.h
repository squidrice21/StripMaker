#pragma once

#ifdef NDEBUG

#include "detail/suppress_warning.h"

#include <cstdio>
#include <cstdlib>

#define force_assert(EXPR)                                                               \
  do {                                                                                   \
    MSVC_WARNING_SUPPRESS(4127);                                                         \
    if (!(EXPR)) {                                                                       \
      ::std::printf("Assertion %s failed in %s:%d\n", #EXPR, __FILE__, __LINE__);        \
      ::std::abort();                                                                    \
    }                                                                                    \
  } while (0)

#else // NDEBUG

#include <cassert>

#define force_assert assert

#endif // NDEBUG
