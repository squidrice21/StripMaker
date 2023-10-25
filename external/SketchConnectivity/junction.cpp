#include "junction.h"

#include "busyness.h"
#include "closest.h"
#include "detail/util.h"
#include "intersect.h"
#include "stroke_graph.h"

#include <spdlog/spdlog.h>

#include <fstream>
#include <iomanip>
#include <map>
#include <unordered_set>

namespace sketching {

namespace junction {
namespace {

template <typename T>
void concat(std::vector<T>& a, const std::vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
}

template <typename T>
bool contains(const std::vector<T>& vector, const T& value) {
  return std::find(vector.cbegin(), vector.cend(), value) != vector.end();
}

} // namespace
} // namespace junction

Junction::Junction(std::vector<StrokeTime> pts, const std::string& t, const bool weak,
                   double pred_probability, bool skip_init)
  : points(std::move(pts))
  , is_weak(weak)
  , probability(pred_probability) {
  if (t == "r") {
    type = JunctionType::R;
  } else if (t == "t") {
    type = JunctionType::T;
  } else if (t == "x") {
    type = JunctionType::X;
  } else {
    throw std::domain_error("Unknown junction type");
  }
  if (!skip_init)
    init();
}

Junction::Junction(std::vector<StrokeTime> pts, const JunctionType::Type t, bool weak,
                   double pred_probability, bool skip_init)
  : points(std::move(pts))
  , type(t)
  , is_weak(weak)
  , probability(pred_probability) {
  if (!skip_init)
    init();
}

void Junction::throw_on_invalid() const {
  auto enforce_snapping = false;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& p = points[i];
    if (p.first < 0)
      throw std::domain_error(fmt::format("invalid index {}", p.first));
    if (type == JunctionType::X || (i == 0 && type == JunctionType::T)) {
      if (p.second < 0.0 || p.second > 1.0)
        throw std::domain_error(fmt::format("invalid arc length {}", p.second));
    } else {
      if (p.second != 0.0 && p.second != 1.0) {
        if (type != JunctionType::T || enforce_snapping) {
          throw std::domain_error(
            fmt::format("endpoints must be snapped, got arc length {}", p.second));
        }
      } else {
        // Now in endpoint territory.
        enforce_snapping = true;
      }
    }
  }
}

bool Junction::is_valid() const {
  auto enforce_snapping = false;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const auto& p = points[i];
    if (p.first < 0)
      return false;
    if (type == JunctionType::X || (i == 0 && type == JunctionType::T)) {
      if (p.second < 0.0 || p.second > 1.0)
        return false;
    } else {
      if (p.second != 0.0 && p.second != 1.0) {
        if (type != JunctionType::T || enforce_snapping) {
          return false;
        }
      } else {
        // Now in endpoint territory.
        enforce_snapping = true;
      }
    }
  }
  return true;
}

void Junction::init() {
  throw_on_invalid();
  sort_entries();
}

void Junction::sort_entries() {
  // Propagate non-endpoints to front.
  std::sort(points.begin(), points.end(), [](StrokeTime a, StrokeTime b) {
    if ((a.second == 0.0 || a.second == 1.0) && (b.second != 0.0 && b.second != 1.0)) {
      return false;
    } else if ((b.second == 0.0 || b.second == 1.0) &&
               (a.second != 0.0 && a.second != 1.0)) {
      return true;
    }
    return a < b;
  });
}

void Junction::remove_duplicate_points() {
  auto begin = points.begin();
  if (type == JunctionType::T) {
    ++begin;
  }
  std::sort(begin, points.end());
  points.erase(std::unique(begin, points.end()), points.end());
}

void Junction::add(StrokeTime p) {
  if (!::sketching::junction::contains(points, p)) {
    points.emplace_back(p);
    sort_entries();
  }
}

void Junction::merge(const Junction& other) {
  for (const auto& p : other.points) {
    if (type == JunctionType::R && other.type == JunctionType::T) {
      // Note this can only trigger on the first iteration of the loop.
      type = JunctionType::T;
      points.insert(points.begin(), p);
    } else if (type == JunctionType::R && other.type == JunctionType::X) {
      type = JunctionType::X;
      if (!::sketching::junction::contains(points, p)) {
        points.emplace_back(p);
      }
    } else if (!::sketching::junction::contains(points, p)) {
      points.emplace_back(p);
    }
  }
  sort_entries();
}

Junction Junction::merged(const Junction& other) const {
  Junction out = *this;
  out.merge(other);
  return out;
}

void Junction::normalize() {
  sort_entries();

  // Make the type as specific as possible.
  const auto first_is_endpoint = (points[0].second == 0.0 || points[0].second == 1.0);
  auto rest_has_non_endpoints = false;
  for (size_t i = 1; i < points.size(); ++i) {
    const auto arclen = points[i].second;
    if (arclen != 0.0 && arclen != 1.0) {
      rest_has_non_endpoints = true;
      break;
    }
  }
  if (first_is_endpoint && !rest_has_non_endpoints) {
    type = JunctionType::R;
  } else if (!first_is_endpoint && !rest_has_non_endpoints) {
    type = JunctionType::T;
  } else {
    type = JunctionType::X;
  }
}

void Junction::normalize(const span<const Stroke> strokes) {
  std::sort(points.begin(), points.end());
  if (!points.empty()) {
    auto [prev_stroke_idx, prev_arclen] = points[0];
    auto prev_pt_idx = size_t(0);
    for (size_t i = 1; i < points.size(); ++i) {
      const auto [idx, arclen] = points[i];
      if (idx == prev_stroke_idx) {
        const auto& stroke = strokes[idx];
        const auto pw = stroke.pen_width();
        const auto l = stroke.length();
        // Merge too-close points on the same stroke.
        // TODO: Is this the right thing to do if both are endpoints?
        if (l * std::abs(arclen - prev_arclen) < 5 * pw) {
          // Mark for deletion.
          if (arclen == 1.0) {
            points[prev_pt_idx].first = -1;
          } else {
            points[i].first = -1;
            auto& prev_arclen_ref = points[prev_pt_idx].second;
            if (prev_arclen_ref != 0.0 && prev_arclen_ref != 1.0) {
              prev_arclen_ref = 0.5 * (points[i].second + prev_arclen_ref);
            }
            continue; // Don't set prev_* to this one, because it's deleted.
          }
        }
      }
      std::tie(prev_stroke_idx, prev_arclen) = points[i];
      prev_pt_idx = i;
    }
  }
  points.erase(std::remove_if(points.begin(), points.end(),
                              [](const std::pair<int, Float>& p) { return p.first < 0; }),
               points.end());

  normalize();
}

unsigned char junction_type_to_char(JunctionType::Type t) {
  switch (t) {
    case JunctionType::T:
      return 't';
    case JunctionType::R:
      return 'r';
    case JunctionType::X:
      return 'x';
    default:
      throw std::invalid_argument("unknown junction type");
  }
}

JunctionType::Type junction_type_from_char(unsigned char t) {
  switch (t) {
    case 't':
      return JunctionType::T;
    case 'r':
      return JunctionType::R;
    case 'x':
      return JunctionType::X;
    default:
      throw std::invalid_argument("unknown junction type");
  }
}

void intersection_junctions(const PolylineBVH& bvh,
                            std::vector<Junction>& out_junctions) {
  const auto create_junction = //
    [&](const size_t stroke_i, const size_t stroke_j,
        const std::vector<Vec2>& intersections) {
      const auto len_i = bvh.nodes[stroke_i].geometry->length();
      const auto len_j = bvh.nodes[stroke_j].geometry->length();
      for (auto [arclen_i, arclen_j] : intersections) {
        arclen_i = std::clamp(arclen_i, 0.0, len_i);
        arclen_j = std::clamp(arclen_j, 0.0, len_j);
        out_junctions
          .emplace_back(std::vector<StrokeTime>{{(int)stroke_i, arclen_i / len_i},
                                                {(int)stroke_j, arclen_j / len_j}},
                        JunctionType::X)
          .sort_entries();
      }
    };

  const auto n = bvh.nodes.size();
  auto intersections = std::vector<Vec2>();
  intersections.reserve(4); // Arbitrary...
  for (size_t i = 0; i < n; ++i) {
    if (bvh.nodes[i].geometry->size() < 2)
      continue;
    intersect_self(bvh.nodes[i], intersections);
    create_junction(i, i, intersections);
    intersections.clear();
    for (size_t j = i + 1; j < n; ++j) {
      if (bvh.nodes[j].geometry->size() < 2)
        continue;
      intersect_different(bvh.nodes[i], bvh.nodes[j], intersections);
      create_junction(i, j, intersections);
      intersections.clear();
    }
  }
}

void intersection_junctions(const PolylineBVH& bvh,
                            std::vector<ExtendedJunction>& out_junctions) {
  assert(out_junctions.empty());

  const auto compute_priority = //
    [](const ExtendedJunction& a, const ExtendedJunction& b) -> int {
    const auto& r0_a = a.ranges_[0].second;
    const auto& r1_a = a.ranges_[1].second;
    const auto& r0_b = b.ranges_[0].second;
    const auto& r1_b = b.ranges_[1].second;
    if (r0_a.overlaps(r0_b)) {
      if (r0_a.start_ < 1e-5 || r0_b.start_ < 1e-5) {
        return (r0_a.mid_ < r0_b.mid_ ? 1 : -1);
      }
      if (r0_a.end_ > 1.0 - 1e-5 || r0_b.end_ > 1.0 - 1e-5) {
        return (r0_a.mid_ < r0_b.mid_ ? -1 : 1);
      }
    }
    if (r1_a.overlaps(r1_b)) {
      if (r1_a.start_ < 1e-5 || r1_b.start_ < 1e-5) {
        return (r1_a.mid_ < r1_b.mid_ ? 1 : -1);
      }
      if (r1_a.end_ > 1.0 - 1e-5 || r1_b.end_ > 1.0 - 1e-5) {
        return (r1_a.mid_ < r1_b.mid_ ? -1 : 1);
      }
    }
    return 0;
  };

  const auto create_junction = //
    [&](const size_t si, const size_t sj, const std::vector<Vec2>& intersections) {
      const auto& stroke_i = *bvh.nodes[si].geometry;
      const auto& stroke_j = *bvh.nodes[sj].geometry;
      const auto len_i = stroke_i.length();
      const auto len_j = stroke_j.length();
      const auto inv_len_i = 1.0 / len_i;
      const auto inv_len_j = 1.0 / len_j;
      const auto begin = out_junctions.size();
      for (const auto [arclen_i, arclen_j] : intersections) {
        const auto norm_arclen_i = arclen_i * inv_len_i;
        const auto norm_arclen_j = arclen_j * inv_len_j;
        const auto r_i = 0.5 * stroke_i.width_at(arclen_i);
        const auto r_j = 0.5 * stroke_j.width_at(arclen_j);
        // Rule for trimming: if the endpoint is covered by the other stroke
        // (approximate with the radius at the intersection point), then trim at the
        // intersection.
        const auto start_i = (arclen_i <= r_j ? 0.0 : norm_arclen_i - 1e-5);
        const auto start_j = (arclen_j <= r_i ? 0.0 : norm_arclen_j - 1e-5);
        const auto end_i = (len_i - arclen_i <= r_j ? 1.0 : norm_arclen_i + 1e-5);
        const auto end_j = (len_j - arclen_j <= r_i ? 1.0 : norm_arclen_j + 1e-5);
        const auto candidate = ExtendedJunction(
          // Range will clamp, so OK.
          {{(int)si, ExtendedJunction::Range(start_i, norm_arclen_i, end_i)},
           {(int)sj, ExtendedJunction::Range(start_j, norm_arclen_j, end_j)}});
        // Prioritize intersections that cut off more endpoint.
        for (auto it = begin; it != out_junctions.size(); ++it) {
          const auto priority = compute_priority(out_junctions[it], candidate);
          if (priority == 1) {
            out_junctions[it] = candidate;
            goto done;
          } else if (priority == -1) {
            goto done;
          }
        }
        out_junctions.emplace_back(candidate);
  done:;
      }
    };

  const auto n = bvh.nodes.size();
  auto intersections = std::vector<Vec2>();
  intersections.reserve(4); // Arbitrary...
  for (size_t i = 0; i < n; ++i) {
    const auto& geom_i = *bvh.nodes[i].geometry;
    if (skip_stroke(geom_i))
      continue;
    intersect_self(bvh.nodes[i], intersections);
    create_junction(i, i, intersections);
    intersections.clear();
    for (size_t j = i + 1; j < n; ++j) {
      const auto& geom_j = *bvh.nodes[j].geometry;
      if (skip_stroke(geom_j))
        continue;
      intersect_different(bvh.nodes[i], bvh.nodes[j], intersections);
      create_junction(i, j, intersections);
      intersections.clear();
    }
  }

  // We do not want to destroy legitimate intersection junctions, so we only trim off
  // what we can.
  auto head_trim = std::unordered_map<int, Float>();
  auto tail_trim = std::unordered_map<int, Float>();
  for (size_t i = 0; i < out_junctions.size(); ++i) {
    auto& junc = out_junctions[i];
    for (const auto& [stroke_idx, range] : junc.ranges_) {
      if (range.start_ == 0.0) {
        if (head_trim.find(stroke_idx) == head_trim.end()) {
          head_trim[stroke_idx] = range.mid_;
        } else {
          head_trim[stroke_idx] = std::min(head_trim[stroke_idx], range.mid_);
        }
      }
      if (range.end_ == 1.0) {
        if (tail_trim.find(stroke_idx) == tail_trim.end()) {
          tail_trim[stroke_idx] = range.mid_;
        } else {
          tail_trim[stroke_idx] = std::max(tail_trim[stroke_idx], range.mid_);
        }
      }
    }
  }
  for (size_t i = 0; i < out_junctions.size(); ++i) {
    auto& junc = out_junctions[i];
    for (auto& [stroke_idx, range] : junc.ranges_) {
      if (range.start_ == 0.0 && head_trim.find(stroke_idx) != head_trim.end() &&
          range.mid_ > head_trim[stroke_idx]) {
        range.start_ = range.mid_ - 1e-5;
      }
      if (range.end_ == 1.0 && tail_trim.find(stroke_idx) != tail_trim.end() &&
          range.mid_ < tail_trim[stroke_idx]) {
        range.end_ = range.mid_ + 1e-5;
      }
    }
  }
}

bool ExtendedJunction::Range::overlaps(const ExtendedJunction::Range& other) const {
  return (start_ <= other.end_ && other.start_ <= end_);
}

void ExtendedJunction::Range::merge(const ExtendedJunction::Range& other) {
  start_ = std::min(start_, other.start_);
  end_ = std::max(end_, other.end_);
  if (mid_ < 0.0) {
    mid_ = other.mid_;
  } else if (mid_ >= 0.0 && other.mid_ >= 0.0) {
    // Try to trim as much as we can.
    if (start_ < 1e-5) {
      mid_ = std::max(mid_, other.mid_);
    } else if (end_ > 1.0 - 1e-5) {
      mid_ = std::min(mid_, other.mid_);
    } else {
      // Interior; doesn't matter too much I think...
    }
  }
}

std::string ExtendedJunction::Range::repr() const {
  return fmt::format("ExtendedJunction::Range({}, {}, {})", start_, mid_, end_);
}

ExtendedJunction::ExtendedJunction(const Junction& junc) {
  ranges_.reserve(junc.size());
  for (const auto& p : junc.points) {
    ranges_.emplace_back(p.first, ExtendedJunction::Range(p.second, p.second));
  }
}

void ExtendedJunction::merge(const ExtendedJunction& other) {
  for (auto& [i, other_range] : other.ranges_) {
    for (auto& [j, range] : ranges_) {
      if (i == j && range.overlaps(other_range)) {
        range.merge(other_range);
        // Necessary if we care about the junction's ranges since it may overlap
        // multiple.
        resolve_overlaps();
        goto done;
      }
    }
    // Add new range.
    ranges_.emplace_back(i, other_range);
done:;
  }
}

void ExtendedJunction::resolve_overlaps() {
  std::sort(ranges_.begin(), ranges_.end(),
            [](const std::pair<int, ExtendedJunction::Range>& a,
               const std::pair<int, ExtendedJunction::Range>& b) {
              return std::pair(a.first, a.second.start_) <
                     std::pair(b.first, b.second.start_);
            });
  const auto n = ranges_.size();
  auto last = decltype(n){0};
  for (auto i = decltype(n){1}; i < n; ++i) {
    auto& [last_si, last_range] = ranges_[last];
    auto& [si, range] = ranges_[i];
    if (si == last_si && range.overlaps(last_range)) {
      last_range.merge(range);
      range.start_ = -1.0; // Mark for deletion.
    } else {
      last = i;
    }
  }
  ranges_.erase(std::remove_if(ranges_.begin(), ranges_.end(),
                               [](const std::pair<int, ExtendedJunction::Range>& ir) {
                                 return ir.second.start_ == -1.0;
                               }),
                ranges_.end());
}

Vec2 ExtendedJunction::position(span<const Stroke> strokes) const {
  // Highest priority: position of centerline intersection.
  for (const auto& [i, range] : ranges_) {
    if (range.mid_ >= 0.0 && range.start_ > 1e-5 && range.end_ < 1.0 - 1e-5) {
      return strokes[i].pos_norm(range.mid_);
    }
  }
  // Second choice: average of non-endpoints (endpoints are the most "move-able").
  {
    auto acc = Vec2(0, 0);
    auto denom = 0;
    for (const auto& [i, range] : ranges_) {
      if (range.start_ > 1e-5 || range.end_ < 1.0 - 1e-5) {
        acc += strokes[i].pos_norm(range.mid());
        denom++;
      }
    }
    if (denom != 0)
      return acc * (1.0 / denom);
  }
  // Last choice: average of endpoints.
  {
    auto acc = Vec2(0, 0);
    auto denom = 0;
    for (const auto& [i, range] : ranges_) {
      acc += strokes[i].pos_norm(range.mid());
      denom++;
    }
    return acc * (1.0 / denom);
  }
}

Float ExtendedJunction::junc_parameter(size_t s_idx) const {
  for (const auto& [_i, range] : ranges_) {
    if (_i != s_idx)
      continue;
    if (range.mid_ >= 0.0 && range.start_ > 1e-5 && range.end_ < 1.0 - 1e-5) {
      return range.mid_;
    } else if (range.start_ <= 1e-5)
      return 0;
    else if (range.end_ >= 1.0 - 1e-5)
      return 1;
  }
  return -1;
}

bool ExtendedJunction::is_grazing() const {
  for (const auto& [_i, range] : ranges_) {
    if (range.mid_ >= 0.0 || range.start_ < 1e-5 || range.end_ >= 1.0 - 1e-5) {
      return false;
    }
  }
  return true;
}

StrokeCoverage::StrokeCoverage(span<const ExtendedJunction> junctions)
  : split_junctions_(new ExtendedJunction[junctions.size()]) {
  for (size_t i = 0; i < junctions.size(); ++i) {
    split_junctions_[i] = junctions[i];
  }

  for (size_t i = 0; i < junctions.size(); ++i) {
    auto& junc = split_junctions_[i];
    auto to_merge = std::unordered_set<const ExtendedJunction*>();
    for (const auto& [stroke_idx, range] : junc.ranges_) {
      for (const auto& other_terr : coverage_[stroke_idx]) {
        if (other_terr.range_.overlaps(range)) {
          to_merge.insert(other_terr.junction_);
        }
      }
    }
    for (const auto* other : to_merge) {
      junc.merge(*other);
    }
    // Re-traverse junc to get transitivity.
    for (const auto& [stroke_idx, range] : junc.ranges_) {
      auto& territories = coverage_[stroke_idx];
      for (auto& territory : territories) {
        if (territory.range_.overlaps(range)) {
          territory.junction_ = &junc;
        }
      }
      territories.emplace_back(
        Territory{ExtendedJunction::Range(range.start_, range.mid_, range.end_), &junc});
      resolve_overlaps(territories);
    }
  }
}

void StrokeCoverage::resolve_overlaps(std::vector<Territory>& territories) {
  std::sort(territories.begin(), territories.end(),
            [](const Territory& a, const Territory& b) {
              return a.range_.start_ < b.range_.start_;
            });
  auto last = size_t(0);
  for (size_t i = 1; i < territories.size(); ++i) {
    if (territories[i].range_.overlaps(territories[last].range_)) {
      assert(territories[i].junction_ == territories[last].junction_);
      territories[last].range_.merge(territories[i].range_);
      territories[i].junction_ = nullptr; // Mark for deletion.
    } else {
      last = i;
    }
  }
  territories.erase(
    std::remove_if(territories.begin(), territories.end(),
                   [](const Territory& t) { return t.junction_ == nullptr; }),
    territories.end());
}

std::vector<const ExtendedJunction*> StrokeCoverage::junctions() const {
  auto already_added = std::unordered_set<const ExtendedJunction*>();
  auto out = std::vector<const ExtendedJunction*>();
  for (const auto& [si, terrs] : coverage_) {
    for (const auto& terr : terrs) {
      if (already_added.find(terr.junction_) == already_added.end()) {
        out.push_back(terr.junction_);
        already_added.insert(terr.junction_);
      }
    }
  }
  return out;
}

std::string StrokeCoverage::repr() const {
  auto ss = std::stringstream();
  ss << "StrokeCoverage(\n";
  for (const auto& [si, terrs] : coverage_) {
    ss << "  " << si << ": ";
    for (const auto& terr : terrs) {
      ss << "(" << terr.range_.start_ << ", " << terr.range_.mid_ << ", "
         << terr.range_.end_ << ", #" << (terr.junction_ - &split_junctions_[0]) << "), ";
    }
    ss << '\n';
  }
  ss << ")";
  return ss.str();
}

std::string StrokeCoverage::junctions_repr() const {
  auto ss = std::stringstream();
  ss << "StrokeCoverage.Junctions(\n";
  for (const auto* junc : junctions()) {
    ss << "  " << (junc - &split_junctions_[0]) << ": ";
    for (const auto& [si, range] : junc->ranges_) {
      ss << "(s" << si << ", " << range.start_ << ", " << range.mid_ << ", " << range.end_
         << "), ";
    }
    ss << '\n';
  }
  ss << ")";
  return ss.str();
}

std::string to_json(const Junction& junc) {
  auto ss = std::stringstream();
  // Enough precision to round-trip a double precision float.
  ss << std::fixed << std::setprecision(std::numeric_limits<double>::digits10 + 1);

  ss << '{';
  ss << "\"type\": \"" << ((junc.type == JunctionType::R) ? "r" : "t") << "\", ";
  ss << "\"probability\": " << junc.probability << ", ";
  ss << "\"is_weak\": " << ((junc.is_weak) ? "true" : "false") << ", ";

  std::set<size_t> search_str;
  ss << "\"points\": [";
  for (size_t i = 0; i < junc.points.size(); ++i) {
    search_str.emplace(junc.points[i].first);
    ss << "[" << junc.points[i].first << ", " << junc.points[i].second << "]";
    if (i + 1 < junc.points.size())
      ss << ", ";
  }
  ss << "], ";

  ss << "\"_searchstr\": \"";
  size_t i = 0;
  for (const auto sid : search_str) {
    ss << sid;
    if (i + 1 < search_str.size())
      ss << ",";
    i++;
  }
  if (search_str.size() == 1)
    ss << "," << *search_str.begin();
  ss << "\"";

  ss << "}";
  return ss.str();
}

void to_json(const std::vector<Junction>& junctions, const std::string& path) {
  auto file = std::ofstream();
  file.open(path);
  file << "{\"junctions\": [";

  for (size_t i = 0; i < junctions.size(); ++i) {
    const auto junc_str = to_json(junctions[i]);
    file << junc_str;
    if (i + 1 < junctions.size())
      file << ", ";
  }

  file << "]}\n";
  file.close();
}

} // namespace sketching
