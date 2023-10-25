#pragma once

#include "sketching.h"

namespace sketching {

/**
 * Non-owning view of (part of) a `Stroke`.
 */
struct ConstStrokeView {
  ConstStrokeView(const Stroke& stroke)
    : stroke_(stroke)
    , start_(0)
    , end_(stroke.size()) {}

  /**
   * @param stroke The stroke to create a view of.
   * @param start Value in the interval [0, stroke.size()].
   * @param end Value in the interval [0, stroke.size()]. Must be greater than or equal to
   *            `start`.
   */
  ConstStrokeView(const Stroke& stroke, Index start, Index end)
    : stroke_(stroke)
    , start_(start)
    , end_(end) {}

  ConstStrokeView(const ConstStrokeView& view, Index start, Index end)
    : stroke_(view.stroke())
    , start_(start + view.start_)
    , end_(end + view.start_) {
    if (end_ > view.end_) {
      throw std::out_of_range("");
    }
  }

  Index size() const { return end_ - start_; }

  const Vec2 xy(Index i) const {
    assert(i < size() && "out of bounds");
    return stroke_.xy(start_ + i);
  }
  const Vec2 xy(Back_t) const {
    assert(size() > 0 && "out of bounds");
    return stroke_.xy(end_ - 1);
  }

  Float x(Index i) const {
    assert(i < size() && "out of bounds");
    return stroke_.x(start_ + i);
  }
  Float x(Back_t) const {
    assert(size() > 0 && "out of bounds");
    return stroke_.x(end_ - 1);
  }
  span<const Float> x() const {
    return {stroke_.x().data() + start_, size_t(end_ - start_)};
  }

  Float y(Index i) const {
    assert(i < size() && "out of bounds");
    return stroke_.y(start_ + i);
  }
  span<const Float> y() const {
    return {stroke_.y().data() + start_, size_t(end_ - start_)};
  }

  Float width(Index i) const {
    assert(i < size() && "out of bounds");
    return stroke_.width(start_ + i);
  }
  Float width(Back_t) const {
    assert(size() > 0 && "out of bounds");
    return stroke_.width(end_ - 1);
  }
  span<const Float> width() const {
    return {stroke_.width().data() + start_, size_t(end_ - start_)};
  }

  Float arclength(Index i) const {
    assert(i < size() && "out of bounds");
    return stroke_.arclength(start_ + i);
  }

  Vec2 pos(Float s) const { return stroke_.pos(s + stroke_.arclength(start_)); }

  Float length() const {
    stroke_.ensure_arclengths();
    return stroke_.arclength(end_ - 1) - stroke_.arclength(start_);
  }

  Stroke slice() const {
    auto start = start_;
    assert(start <= end_);
    const auto ht = stroke_.has_time();
    Stroke out(end_ - start, ht);
    for (Index i = 0; start < end_; ++start, ++i) {
      out.x(i) = stroke_.x(start);
      out.y(i) = stroke_.y(start);
      out.width(i) = stroke_.width(start);
      if (ht)
        out.time(i) = stroke_.time(start);
    }
    return out;
  }

  /** Return a reference to the underlying stroke. */
  const Stroke& stroke() const { return stroke_; }

private:
  const Stroke& stroke_;

public:
  Index start_;
  Index end_;
};

} // namespace sketching
