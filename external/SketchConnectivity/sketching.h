#pragma once

#include "types.h"

#include <memory>

namespace sketching {

/// Signify that you want the last element of something.
constexpr struct Back_t {
} Back;

struct Stroke {
  Stroke();

  Stroke(Index npoints, bool has_time);

  Index size() const { return size_; }

  span<Float> x() { return {x_, (size_t)size_}; }
  span<const Float> x() const { return {x_, (size_t)size_}; }

  span<Float> y() { return {y_, (size_t)size_}; }
  span<const Float> y() const { return {y_, (size_t)size_}; }

  span<Float> width() { return {width_, (size_t)size_}; }
  span<const Float> width() const { return {width_, (size_t)size_}; }

  Float &x(Index i) {
    assert(i >= 0 && i < size_ && "out of bounds");
    return x_[i];
  }
  Float x(Index i) const {
    assert(i >= 0 && i < size_ && "out of bounds");
    return x_[i];
  }
  Float &x(Back_t) {
    assert(size_ > 0 && "out of bounds");
    return x_[size_ - 1];
  }
  Float x(Back_t) const {
    assert(size_ > 0 && "out of bounds");
    return x_[size_ - 1];
  }

  Float &y(Index i) {
    assert(i >= 0 && i < size_ && "out of bounds");
    return y_[i];
  }
  Float y(Index i) const {
    assert(i >= 0 && i < size_ && "out of bounds");
    return y_[i];
  }
  Float &y(Back_t) {
    assert(size_ > 0 && "out of bounds");
    return y_[size_ - 1];
  }
  Float y(Back_t) const {
    assert(size_ > 0 && "out of bounds");
    return y_[size_ - 1];
  }

  const Vec2 xy(Index i) const {
    assert(i >= 0 && i < size_ && "out of bounds");
    return {x_[i], y_[i]};
  }
  const Vec2 xy(Back_t) const {
    assert(size_ > 0 && "out of bounds");
    return {x_[size() - 1], y_[size() - 1]};
  }

  const Vec2 pos(const Float s) const {
    assert(size_ > 0);
    if (size() == 1) {
      return {x_[0], y_[0]};
    }
    const auto v = fractional_index(s);
    return lerp_pos(v.first, v.second);
  }

  const Vec2 pos_norm(const Float normalized_arclen) const {
    return pos(normalized_arclen * length());
  }

  /**
   * Return the fractional index corresponding to given arc length.
   */
  std::pair<Index, Float> fractional_index(Float s) const;

  Float &width(Index i) {
    assert(i >= 0 && i < size_ && "out of bounds");
    return width_[i];
  }
  Float width(Index i) const {
    assert(i >= 0 && i < size_ && "out of bounds");
    return width_[i];
  }
  Float &width(Back_t) {
    assert(size_ > 0);
    return width_[size_ - 1];
  }
  Float width(Back_t) const {
    assert(size_ > 0);
    return width_[size_ - 1];
  }

  /// Get the width at arc length position s.
  Float width_at(Float s) const {
    const auto [i, u] = fractional_index(s);
    return lerp(width_, i, u);
  }

  bool has_time() const { return time_ != nullptr; }

  span<Float> time() {
    if (time_) {
      return {time_, (size_t)size_};
    }
    return span<Float>();
  }
  span<const Float> time() const {
    if (time_) {
      return {time_, (size_t)size_};
    }
    return span<Float>();
  }

  Float &time(Index i) {
    assert(i >= 0 && i < size_ && "out of bounds");
    return time_[i];
  }
  Float time(Index i) const {
    assert(i >= 0 && i < size_ && "out of bounds");
    return time_[i];
  }
  Float time(Back_t) const {
    assert(size_ > 0);
    return time_[size_ - 1];
  }

  Float length() const {
    assert(size_ > 0);
    ensure_arclengths();
    return arclength_[size() - 1];
  }

  /**
   * Return the maximum width along the stroke.
   */
  Float pen_width() const { return *std::max_element(width_, width_ + size_); }

  Stroke clone() const;

  /**
   * Make the `arclength` methods valid.  Needs to be manually called again if
   * there are any changes to `x` or `y`.  Warns about repeated vertices.
   */
  void compute_arclengths() const;

  /**
   * Compute arc lengths if they are not already computed.
   */
  void ensure_arclengths() const {
    if (!has_arclengths())
      compute_arclengths();
  }

  /**
   * Invalidate the `arclength` methods.  This should be called if you make any
   * changes to the `x` or `y`.
   */
  void invalidate_arclengths();

  /**
   * Return true iff arc lengths are available.  However, this doesn't guarantee
   * that they are up to date, especially if x or y have been updated since
   * `compute_arclengths` was last called.
   */
  bool has_arclengths() const { return arclength_ != nullptr; }

  span<Float> arclength() {
    assert(arclength_ && "call compute_arclengths first");
    return {arclength_.get(), (size_t)size_};
  }
  span<const Float> arclength() const {
    assert(arclength_ && "call compute_arclengths first");
    return {arclength_.get(), (size_t)size_};
  }

  const Float &arclength(Index i) {
    assert(arclength_ && "call compute_arclengths first");
    return arclength_[i];
  }
  const Float &arclength(Index i) const {
    assert(arclength_ && "call compute_arclengths first");
    return arclength_[i];
  }

  /**
   * Compute the mean distance between adjacent vertices in the stroke.
   */
  Float avg_sampling() const;

  void clear() {
    size_ = 0;
    invalidate_arclengths();
  }

  void insert(Index position, Float x, Float y, Float width, Float time = 0.0);

  void push_back(Float x, Float y, Float width, Float time = 0.0);

  void reserve(Index capacity);

  void resize(Index size) {
    reserve(size);
    size_ = size;
    invalidate_arclengths();
  }

  /**
   * Reverse the order of the vertices.
   */
  void reverse();

  /**
   * Shorten the stroke to just vertices [start, stop).
   */
  void trim(Index start_idx, Index stop_idx);

  void trim(Float start, Float stop);

  /**
   * Return a copy of part of this stroke.  This function will produce new
   * vertices at the cut points; to always cut exactly at vertices use the
   * `slice` function in `ConstStrokeView`.
   *
   * `stop` is inclusive, so passing `(0, stroke.length())` will get you a copy
   * of the original stroke.
   *
   * @param start Arc length value.
   * @param stop Arc length value.
   */
  [[nodiscard]] Stroke fractional_slice(Float start, Float stop) const;

  /**
   * Split (cut) a curve into many smaller curves.
   *
   * split_values must be ordered from smallest to largest. Values may not
   * repeat. If you want the first segment to start at the start of the curve,
   * the first split value must be 0. If you want the last segment to end at the
   * end of the curve, the last split value must equal length().
   *
   * The this pointer must _not_ be contained inside out_splits.
   *
   * Will never return length 0 strokes.
   */
  void split(span<const Float> split_values,
             std::vector<Stroke> &out_splits) const;

  std::string repr() const;

  /// Move constructors.
  Stroke(Stroke &&mE) noexcept = default;
  Stroke &operator=(Stroke &&mE) noexcept = default;

  // Use `clone` to do an explicit copy.
  Stroke(Stroke &mE) noexcept = delete;
  Stroke &operator=(Stroke &mE) noexcept = delete;

  /// Memory that this stroke owns.
  /// Can be nullptr for non-owning views of existing memory.
  std::unique_ptr<Float[]> memory_;

  /// Number of vertices in the polyline.
  Index size_;

  /// `memory_` will hold capacity for 3 x capacity_ Floats if `time_ ==
  /// nullptr` and 4 x capacity_ otherwise.
  Index capacity_;

  Float *x_;
  Float *y_;
  Float *width_;
  Float *time_;

  size_t index;

private:
  /// Not available until `compute_arclengths` is called.
  /// TODO: Support non-owning memory.
  mutable std::unique_ptr<Float[]> arclength_ = nullptr;

  Vec2 lerp_pos(const Index i, Float u) const {
    return Vec2((1 - u) * x(i) + u * x(i + 1), //
                (1 - u) * y(i) + u * y(i + 1));
  }

  static Float lerp(const Float arr[], const Index i, Float u) {
    return arr[i] + u * (arr[i + 1] - arr[i]);
  }
};

/// Copy vertices [src_start, src_end] or [src_end, src_start] (depending on if
/// `src_start < src_end`) from the src stroke to vertices [dst_start, dst_end]
/// or [dst_end, dst_start] of the dst stroke.  The src_* and dst_* ranges are
/// inclusive and must be the same size.  `dst` is allowed to point to the same
/// stroke as `src` -- in this case, the copy begins at `*_start` and proceeds
/// to `*_end`.  If there is overlap in the ranges, it is up to the caller to
/// make sure the order will result in the desired behaviour.
void copy_run(Stroke &dst, Index dst_start, Index dst_end, //
              const Stroke &src, Index src_start, Index src_end);

} // namespace sketching
