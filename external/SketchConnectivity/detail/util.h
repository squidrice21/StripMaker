/**
 * Not to be included by any header!
 */
#pragma once

#include "../types.h"

namespace sketching {

constexpr auto infinity = std::numeric_limits<Float>::infinity();

constexpr Float square(const Float x) { return x * x; }

constexpr Float cube(const Float x) { return x * x * x; }

constexpr double clamped_lerp(double x, double y, double u) {
  const auto val = x + u * (y - x);
  if (x > y) {
    std::swap(x, y);
  }
  return std::clamp(val, x, y);
}

constexpr double fast_lerp(double x, double y, double u) { return x + u * (y - x); }

constexpr double unlerp(double a, double b, double t) { return (t - a) / (b - a); }

constexpr double linear_remap(double x, double x_min, double x_max, //
                              double out_min, double out_max) {
  return fast_lerp(out_min, out_max, unlerp(x_min, x_max, x));
}

inline Vec2 lerp(const Vec2& x, const Vec2& y, double u) { return x + u * (y - x); }

} // namespace sketching
