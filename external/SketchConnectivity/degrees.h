#pragma once

#include "types.h"

namespace sketching {

/// Returns `degrees` converted to radians.
constexpr Float operator"" _deg(unsigned long long degrees) {
  return Float(degrees) * (M_PI / 180);
}

constexpr Float operator"" _deg(long double degrees) {
  return Float(degrees * (M_PI / 180));
}

} // namespace sketching
