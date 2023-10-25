#pragma once

#include <Eigen/Core>
#include <nonstd/span.hpp>

namespace sketching {

using Float = double;
using Index = std::ptrdiff_t;
using Vec = Eigen::Matrix<Float, Eigen::Dynamic, 1>;
template <typename T = Float>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using CoordMat = Eigen::Matrix<Float, Eigen::Dynamic, 2, Eigen::RowMajor>;

template <typename T>
using span = ::nonstd::span<T>;

inline void operator*=(span<Float> array, Float scalar) {
  for (auto& x : array) {
    x *= scalar;
  }
}

struct Vec2 {
  Float x_;
  Float y_;

  constexpr Vec2()
    : x_(0)
    , y_(0){};

  constexpr Vec2(Float x_val, Float y_val)
    : x_(x_val)
    , y_(y_val){};

  /// Create a vector without initializing its components.
  static Vec2 Empty() { return Vec2(Uninitialized()); }

  struct Uninitialized {};
  explicit Vec2(Uninitialized){};

  // Eigen compatibility. Remove later.
  Float& x() { return x_; }
  Float x() const { return x_; }
  Float& y() { return y_; }
  Float y() const { return y_; }

  [[nodiscard]] Vec2 operator+(Vec2 other) const {
    return {x_ + other.x_, y_ + other.y_};
  }
  [[nodiscard]] Vec2 operator-(Vec2 other) const {
    return {x_ - other.x_, y_ - other.y_};
  }
  [[nodiscard]] Vec2 operator*(Float scalar) const { return {x_ * scalar, y_ * scalar}; }
  [[nodiscard]] Vec2 operator/(Float scalar) const { return {x_ / scalar, y_ / scalar}; }

  void operator+=(Vec2 other) {
    x_ += other.x_;
    y_ += other.y_;
  }
  void operator-=(Vec2 other) {
    x_ -= other.x_;
    y_ -= other.y_;
  }
  void operator*=(Float scalar) {
    x_ *= scalar;
    y_ *= scalar;
  }
  void operator/=(Float scalar) {
    x_ /= scalar;
    y_ /= scalar;
  }

  Vec2 operator-() const { return {-x_, -y_}; }

  Float norm() const { return std::sqrt(x_ * x_ + y_ * y_); }

  // Eigen compatibility. Rename to snake_case later.
  Float squaredNorm() const { return x_ * x_ + y_ * y_; }

  void normalize() {
    const auto n = norm();
    x_ /= n;
    y_ /= n;
  }

  [[nodiscard]] Vec2 normalized() const {
    const auto n = norm();
    return {x_ / n, y_ / n};
  }

  [[nodiscard]] Float dot(Vec2 other) const { return x_ * other.x_ + y_ * other.y_; }

  [[nodiscard]] Float cross(Vec2 other) const { return x_ * other.y_ - y_ * other.x_; }

  // Eigen compatibility. Rename to snake_case later.
  [[nodiscard]] bool hasNaN() const { return std::isnan(x_) || std::isnan(y_); }

  [[nodiscard]] bool operator==(Vec2 other) const {
    return x_ == other.x_ && y_ == other.y_;
  }
  [[nodiscard]] bool operator!=(Vec2 other) const { return !(*this == other); }

  bool isApprox(Vec2 other, Float prec = 1e-6) const {
    return std::abs(x_ - other.x_) < prec && std::abs(y_ - other.y_) < prec;
  }

  constexpr bool operator<(Vec2 other) const {
    if (x_ < other.x_) {
      return true;
    } else if (x_ > other.x_) {
      return false;
    }
    return y_ < other.y_;
  }
};

[[nodiscard]] inline Vec2 operator*(Float scalar, Vec2 v) { return v * scalar; }

struct Mat2 {
  /// Arranged like
  ///     a b
  ///     c d
  Float a, b, c, d;

  Mat2(Float a00, Float a01, Float a10, Float a11)
    : a(a00)
    , b(a01)
    , c(a10)
    , d(a11) {}

  void operator*=(Float scalar) {
    a *= scalar;
    b *= scalar;
    c *= scalar;
    d *= scalar;
  }
};
[[nodiscard]] inline Vec2 operator*(const Mat2& m, Vec2 v) {
  return {m.a * v.x_ + m.b * v.y_, m.c * v.x_ + m.d * v.y_};
}

} // namespace sketching
