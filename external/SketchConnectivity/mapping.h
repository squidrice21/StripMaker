#pragma once

#include "types.h"

#include <algorithm>
#include <vector>

namespace sketching {

struct StrokeMapping {
  /**
   * Evaluate f(arclen).  If `arclen` is outside the domain, then it will be clamped to
   * the domain.
   *
   * Preconditions: `is_valid()`.
   */
  [[nodiscard]] Float map_clamp(Float arclen) const;

  /**
   * Find the value x such that f(x) = arclen.  If `arclen` is outside the range, then one
   * of the domain boundaries will be returned depending on if the mapping is overall
   * increasing or overall decreasing.  For non-injective mappings, return the first such
   * match when beginning a search from the *start* of the mappings array.
   *
   * For injective mappings, both `inv_map_clamp_*` functions are equivalent.
   *
   * Preconditions: `is_valid()`.
   */
  [[nodiscard]] Float inv_map_clamp_from_start(Float arclen) const;

  /**
   * Find the value x such that f(x) = arclen.  If `arclen` is outside the range, then one
   * of the domain boundaries will be returned depending on if the mapping is overall
   * increasing or overall decreasing.  For non-injective mappings, return the first such
   * match when beginning a search from the *end* of the mappings array.
   *
   * For injective mappings, both `inv_map_clamp_*` functions are equivalent.
   *
   * Preconditions: `is_valid()`.
   */
  [[nodiscard]] Float inv_map_clamp_from_end(Float arclen) const;

  [[nodiscard]] std::pair<Float, Float> inv_map_shortest_interval(Float start,
                                                                  Float stop) const;

  [[nodiscard]] bool is_linear() const { return domain_arclens_.size() == 2; }

  [[nodiscard]] bool is_valid() const {
    return (domain_arclens_.size() == range_arclens_.size() &&
            domain_arclens_.size() >= 2 &&
            *std::min_element(domain_arclens_.begin(), domain_arclens_.end()) <
              *std::max_element(domain_arclens_.begin(), domain_arclens_.end()));
  }

  bool is_range_empty() const {
    return *std::min_element(range_arclens_.begin(), range_arclens_.end()) ==
           *std::max_element(range_arclens_.begin(), range_arclens_.end());
  }

  std::vector<Float> domain_arclens_;
  std::vector<Float> range_arclens_;
};

} // namespace sketching
