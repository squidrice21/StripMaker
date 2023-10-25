#include "mapping.h"

#include "force_assert.h"

namespace sketching {

Float StrokeMapping::map_clamp(const Float arclen) const {
  force_assert(!std::isnan(arclen));
  force_assert(is_linear()); // TODO: Implement for non-linear mappings.
  const auto front = domain_arclens_.front();
  const auto back = domain_arclens_.back();
  const auto u = std::clamp((arclen - front) / (back - front), 0.0, 1.0);
  return u * range_arclens_.back() + (1 - u) * range_arclens_[0];
}

Float StrokeMapping::inv_map_clamp_from_start(const Float arclen) const {
  force_assert(!std::isnan(arclen));
  const auto n = range_arclens_.size();
  const auto len = std::max(domain_arclens_[0], domain_arclens_.back());
  for (size_t i = 0; i < n - 1; ++i) {
    if ((range_arclens_[i] <= arclen && arclen <= range_arclens_[i + 1]) ||
        (range_arclens_[i] >= arclen && arclen >= range_arclens_[i + 1])) {
      const auto frac =
        (arclen - range_arclens_[i]) / (range_arclens_[i + 1] - range_arclens_[i]);
      if (!std::isnan(frac)) {
        assert(frac >= -1e-10);
        const auto from = domain_arclens_[i] * (1 - frac) + domain_arclens_[i + 1] * frac;
        return std::clamp(from, 0.0, len);
      }
    }
  }
  const auto range_increasing = range_arclens_[0] < range_arclens_.back();
  if (range_increasing) {
    if (arclen > range_arclens_.back()) {
      return domain_arclens_.back();
    } else {
      return domain_arclens_.front();
    }
  } else {
    if (arclen < range_arclens_.back()) {
      return domain_arclens_.back();
    } else {
      return domain_arclens_.front();
    }
  }
}

Float StrokeMapping::inv_map_clamp_from_end(const Float arclen) const {
  force_assert(!std::isnan(arclen));
  const auto n = range_arclens_.size();
  const auto len = std::max(domain_arclens_[0], domain_arclens_.back());
  for (int i = (int)n - 2; i >= 0; --i) {
    if ((range_arclens_[i] <= arclen && arclen <= range_arclens_[i + 1]) ||
        (range_arclens_[i] >= arclen && arclen >= range_arclens_[i + 1])) {
      const auto frac =
        (arclen - range_arclens_[i]) / (range_arclens_[i + 1] - range_arclens_[i]);
      if (!std::isnan(frac)) {
        assert(frac >= -1e-10);
        const auto from = domain_arclens_[i] * (1 - frac) + domain_arclens_[i + 1] * frac;
        return std::clamp(from, 0.0, len);
      }
    }
  }
  const auto range_increasing = range_arclens_[0] < range_arclens_.back();
  if (range_increasing) {
    if (arclen > range_arclens_.back()) {
      return domain_arclens_.back();
    } else {
      return domain_arclens_.front();
    }
  } else {
    if (arclen < range_arclens_.back()) {
      return domain_arclens_.back();
    } else {
      return domain_arclens_.front();
    }
  }
}

std::pair<Float, Float> StrokeMapping::inv_map_shortest_interval(const Float start,
                                                                 const Float stop) const {
  const auto sta1 = inv_map_clamp_from_start(start);
  const auto sta2 = inv_map_clamp_from_end(start);
  const auto sto1 = inv_map_clamp_from_start(stop);
  const auto sto2 = inv_map_clamp_from_end(stop);
  auto shortest = std::pair<Float, Float>(sta1, sto1);
  if (std::abs(shortest.second - shortest.first) > std::abs(sto2 - sta1)) {
    shortest = {sta1, sto2};
  }
  if (std::abs(shortest.second - shortest.first) > std::abs(sto1 - sta2)) {
    shortest = {sta2, sto1};
  }
  if (std::abs(shortest.second - shortest.first) > std::abs(sto2 - sta2)) {
    shortest = {sta2, sto2};
  }
  if (shortest.first > shortest.second) {
    std::swap(shortest.first, shortest.second);
  }
  return shortest;
}

} // namespace sketching
