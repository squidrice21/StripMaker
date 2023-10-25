#pragma once

#include "types.h"

namespace sketching {

struct Stroke;

namespace detail {

inline double bounding_box_mid(double min, double max) {
  double res = 0.5 * (min + max);

  // What did the above compute?
  //   * if min and max finite:                          the correct finite   mid-value
  //   * if min finite and max infinite (or vice-versa): the correct infinite mid-value
  //   * if min and max infinite, same sign:             the correct infinite mid-value
  //   * if min and max infinite, different signs:       NaN

  return std::isnan(res) ? 0.0 : res;
}

inline double bounding_box_distance(double a, double b) {
  if (a == b)
    return 0.0; // infinity case, avoid inf-inf = NaN
  return std::abs(b - a);
}

inline bool bounding_box_inverted(double min, double max) { return min > max + 1e-10; }

} // namespace detail

// Essentially copied from VPaint.
struct BoundingBox {
  /// Empty bounding box
  BoundingBox()
    : xMin_(std::numeric_limits<double>::infinity())
    , xMax_(-std::numeric_limits<double>::infinity())
    , yMin_(std::numeric_limits<double>::infinity())
    , yMax_(-std::numeric_limits<double>::infinity()) {}

  /// Single-point bounding box at position (x,y)
  BoundingBox(double x, double y)
    : xMin_(x)
    , xMax_(x)
    , yMin_(y)
    , yMax_(y) {}

  BoundingBox(double x, double y, double radius)
    : xMin_(x - radius)
    , xMax_(x + radius)
    , yMin_(y - radius)
    , yMax_(y + radius) {}

  /// Non-empty bounding box specified by its boundaries.
  /// It is safe to call this constructor with either x1==x2, x1<x2, or x2<x1.
  BoundingBox(double x1, double x2, double y1, double y2);

  bool isEmpty() const { return detail::bounding_box_inverted(xMin_, xMax_); }
  bool isDegenerate() const { return height() <= 1e-10 || width() <= 1e-10; }
  bool isInfinite() const {
    return height() == std::numeric_limits<double>::infinity() ||
           width() == std::numeric_limits<double>::infinity();
  }

  bool isProper() const { return !(isDegenerate() || isInfinite()); }

  double xMin() const { return xMin_; }
  double xMax() const { return xMax_; }
  double yMin() const { return yMin_; }
  double yMax() const { return yMax_; }

  /// Compute mid-points (0 if empty, or if min = -infinity and max = +infinity)
  double xMid() const { return detail::bounding_box_mid(xMin_, xMax_); }
  double yMid() const { return detail::bounding_box_mid(yMin_, yMax_); }

  double width() const {
    return isEmpty() ? 0.0 : detail::bounding_box_distance(xMin_, xMax_);
  }

  double height() const {
    return isEmpty() ? 0.0 : detail::bounding_box_distance(yMin_, yMax_);
  }

  double area() const { return isDegenerate() ? 0.0 : width() * height(); }

  /// Compute union.
  BoundingBox united(const BoundingBox& other) const;
  /// Compute intersection.
  BoundingBox intersected(const BoundingBox& other) const;

  /// In-place union.
  void unite(const BoundingBox& other);
  /// In-place intersection.
  void intersect(const BoundingBox& other);

  /// Returns whether the two bounding boxes intersect
  bool intersects(const BoundingBox& other) const {
    return !intersected(other).isEmpty();
  }

  /**
   * Returns whether the two bounding boxes intersect or touch within an epsilon fudge.
   * The overlap can have zero or negligible volume, which differs from `intersects`.
   */
  bool touches(const BoundingBox& other) const {
    return 2 * std::abs((xMin_ + 0.5 * width()) - (other.xMin_ + 0.5 * other.width())) <=
             (width() + other.width()) + 1e-10 &&
           2 * std::abs((yMin_ + 0.5 * height()) -
                        (other.yMin_ + 0.5 * other.height())) <=
             (height() + other.height()) + 1e-10;
  }

  bool contains(const Vec2 p) const {
    return p.x() >= xMin_ && p.x() <= xMax_ && p.y() >= yMin_ && p.y() <= yMax_;
  }

  bool contains(const BoundingBox& other) const {
    return (other.xMin_ >= xMin_ && other.xMax_ <= xMax_ && other.yMin_ >= yMin_ &&
            other.yMax_ <= yMax_);
  }

  double distanceTo(const Vec2 p) const { return distanceTo(p.x(), p.y()); }

  double distanceTo(const double x, const double y) const {
    const auto dx = std::max(0.0, std::max(xMin_ - x, x - xMax_));
    const auto dy = std::max(0.0, std::max(yMin_ - y, y - yMax_));
    return std::sqrt(dx * dx + dy * dy);
  }

  bool operator==(const BoundingBox& other) const noexcept {
    return xMin_ == other.xMin_ && xMax_ == other.xMax_ && yMin_ == other.yMin_ &&
           yMax_ == other.yMax_;
  }
  bool operator!=(const BoundingBox& other) const noexcept { return !(*this == other); }

  double largest_axis_length() const { return std::max(xMax_ - xMin_, yMax_ - yMin_); }

  double xMin_, xMax_, yMin_, yMax_;
};

BoundingBox bounds(const Stroke& s);

BoundingBox visual_bounds(const Stroke& s);

BoundingBox visual_bounds(span<const Stroke> strokes);

struct OrientedBoundingBox {
  /// Center point of bounding box.
  Vec2 center_;
  /// First local axis of the coordinate frame as a unit vector.
  /// The other axis is perpendicular, i.e. (-axis.y, axis.x).
  Vec2 axis_;
  /// Positive halfwidth extents along each axis.
  /// (Width/2, height/2) in the local coordinate frame, in other words.
  Vec2 extents_;

  [[nodiscard]] bool hasNaN() const {
    return center_.hasNaN() || axis_.hasNaN() || extents_.hasNaN();
  }

  Float distance_to(const Vec2& p) const { return std::sqrt(squared_distance_to(p)); }

  Float squared_distance_to(const Vec2& p) const {
    const auto v = p - center_;
    auto sq_dist = 0.0;
    {
      const auto d = v.dot(axis_);
      const auto e = extents_.x_;
      auto excess = Float(0);
      if (d < -e) {
        excess = d + e;
      } else if (d > e) {
        excess = d - e;
      }
      sq_dist += excess * excess;
    }
    {
      const auto axis2 = Vec2(-axis_.y_, axis_.x_);
      const auto d = v.dot(axis2);
      const auto e = extents_.y_;
      auto excess = Float(0);
      if (d < -e) {
        excess = d + e;
      } else if (d > e) {
        excess = d - e;
      }
      sq_dist += excess * excess;
    }
    return sq_dist;
  }
};

} // namespace sketching
