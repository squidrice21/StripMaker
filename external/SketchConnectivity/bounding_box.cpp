#include "bounding_box.h"

#include "sketching.h"

namespace sketching {

namespace {

constexpr auto infinity = std::numeric_limits<Float>::infinity();

}

BoundingBox::BoundingBox(double x1, double x2, double y1, double y2) {
  if (x1 < x2) {
    xMin_ = x1;
    xMax_ = x2;
  } else {
    xMin_ = x2;
    xMax_ = x1;
  }

  if (y1 < y2) {
    yMin_ = y1;
    yMax_ = y2;
  } else {
    yMin_ = y2;
    yMax_ = y1;
  }
}

BoundingBox BoundingBox::united(const BoundingBox& other) const {
  BoundingBox res(*this);
  res.unite(other);
  return res;
}

BoundingBox BoundingBox::intersected(const BoundingBox& other) const {
  BoundingBox res(*this);
  res.intersect(other);
  return res;
}

void BoundingBox::unite(const BoundingBox& other) {
  xMin_ = std::min(xMin_, other.xMin_);
  xMax_ = std::max(xMax_, other.xMax_);
  yMin_ = std::min(yMin_, other.yMin_);
  yMax_ = std::max(yMax_, other.yMax_);
}

void BoundingBox::intersect(const BoundingBox& other) {
  // Compute intersection
  xMin_ = std::max(xMin_, other.xMin_);
  xMax_ = std::min(xMax_, other.xMax_);
  yMin_ = std::max(yMin_, other.yMin_);
  yMax_ = std::min(yMax_, other.yMax_);

  // Handle empty special case
  if (detail::bounding_box_inverted(xMin_, xMax_) ||
      detail::bounding_box_inverted(yMin_, yMax_))
    *this = BoundingBox();
}

BoundingBox bounds(const Stroke& s) {
  if (s.size() == 0) {
    return BoundingBox();
  }
  Float xmin = infinity, ymin = infinity, xmax = -infinity, ymax = -infinity;
  for (Index i = 0; i < s.size(); ++i) {
    xmin = std::min(xmin, s.x(i));
    ymin = std::min(ymin, s.y(i));
    xmax = std::max(xmax, s.x(i));
    ymax = std::max(ymax, s.y(i));
  }
  return BoundingBox(xmin, xmax, ymin, ymax);
}

BoundingBox visual_bounds(const Stroke& s) {
  if (s.size() == 0) {
    return BoundingBox();
  }
  Float xmin = infinity, ymin = infinity, xmax = -infinity, ymax = -infinity;
  for (Index i = 0; i < s.size(); ++i) {
    const auto r = 0.5 * s.width(i);
    xmin = std::min(xmin, s.x(i) - r);
    ymin = std::min(ymin, s.y(i) - r);
    xmax = std::max(xmax, s.x(i) + r);
    ymax = std::max(ymax, s.y(i) + r);
  }
  return BoundingBox(xmin, xmax, ymin, ymax);
}

BoundingBox visual_bounds(const span<const Stroke> strokes) {
  auto out = BoundingBox();
  for (const auto& s : strokes) {
    out.unite(visual_bounds(s));
  }
  return out;
}

} // namespace sketching
