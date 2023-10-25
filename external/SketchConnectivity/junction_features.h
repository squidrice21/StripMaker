#pragma once

#include "types.h"

namespace sketching {

struct PolylineBVH;
struct Stroke;
struct StrokeGraph;

struct JunctionFeature {
  [[nodiscard]] virtual Float operator()( //
    const Stroke& s1, Float arclen1, std::size_t stroke_idx1, //
    const Stroke& s2, Float arclen2, std::size_t stroke_idx2) const = 0;

  /** Pre-compute caches and so on. */
  virtual void init(const PolylineBVH&) = 0;

  /**
   * Pre-compute caches and so on.
   *
   * Override this method if you want to compute features differently if connectivity
   * information is available.
   */
  virtual void init(const StrokeGraph&);

  virtual std::string description() const = 0;

  virtual void human_readable(span<Float> /*v*/) const {}

  virtual ~JunctionFeature() = default;

protected:
  JunctionFeature() = default;
};

} // namespace sketching
