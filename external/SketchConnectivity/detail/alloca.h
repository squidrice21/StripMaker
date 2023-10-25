#pragma once

#define ALLOCA(TYPE, COUNT) (TYPE*)alloca((COUNT) * sizeof(TYPE))
/// Note ALLOCA_SPAN is an unsafe macro.
/// https://gcc.gnu.org/onlinedocs/cpp/Duplication-of-Side-Effects.html
#define ALLOCA_SPAN(TYPE, COUNT) span<TYPE>(ALLOCA(TYPE, COUNT), COUNT)
