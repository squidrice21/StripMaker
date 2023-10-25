#pragma once

#include "bvh.h"
#include "classifier.h"

#include <unordered_map>
#include <unordered_set>

namespace sketching {

struct Junction {
  // Stroke index and normalized (i.e. in [0, 1]) arc length pairs.
  // Non-endpoints come first, followed by endpoints.  Within these regions, they are
  // sorted using `operator<` of `StrokeTime`.
  std::vector<StrokeTime> points;
  std::vector<StrokeTime> point_tos;

  JunctionType::Type type;
  JunctionType::Type corner_type = JunctionType::X;

  bool is_weak = false;
  bool must_disconnect = false;
  bool to_reconsider = false;

  // Junction variable component
  int component_idx = -1;
  // Containing region
  int fi = -1;

  // The probability that this junction is connected. Since we only save positive
  // junctions for now. This value would be > 0.5 (unless we change the threshold).
  double probability = -1.0;
  double alt_probability = -1.0;

  // Used to record the junction distance based on the original strokes when the junction
  // is first created. Since we update the T-junction geometry, this distance would shift
  // from the correct one if not recorded.
  double orig_dist = -1.0;

  // The feature of this junction. Used for visualization.
  FeatureVector fea;

  std::string repr;

  explicit Junction(const JunctionType::Type t)
    : type(t) {}

  Junction(std::vector<StrokeTime> pts, const std::string& t, bool is_weak = false,
           double probability = 0.5, bool skip_init = false);
  Junction(std::vector<StrokeTime> pts, JunctionType::Type t, bool is_weak = false,
           double probability = 0.5, bool skip_init = false);

  /**
   * This method *must* be called after modifying points unless you are sure that you have
   * not violated the sort order.
   */
  void sort_entries();

  void remove_duplicate_points();

  /**
   * Note that this method sorts the vector of points, so it is more efficient to specify
   * all the points upfront in the constructor (if you can).
   */
  void add(StrokeTime p);

  void merge(const Junction& other);

  [[nodiscard]] Junction merged(const Junction& other) const;

  /**
   * Sort entries and ensure the type is as specific as possible.
   */
  void normalize();

  /**
   * Merge too-close points on the same stroke, sort entries, and ensure the type is as
   * specific as possible.
   */
  void normalize(span<const Stroke> strokes);

  [[nodiscard]] Junction normalized(const span<const Stroke> strokes) const {
    Junction out = *this;
    out.normalize(strokes);
    return out;
  }

  std::size_t size() const { return points.size(); }

  int stroke(const std::size_t i) const { return points[i].first; }

  double arclength(const std::size_t i) const { return points[i].second; }

  [[nodiscard]] bool empty() const { return points.empty(); }

  bool operator<(const Junction& other) const {
    return std::pair<JunctionType::Type, const decltype(points)&>(type, points) <
           std::pair<JunctionType::Type, const decltype(points)&>(other.type,
                                                                  other.points);
  }

  bool operator==(const Junction& other) const {
    return std::pair<JunctionType::Type, const decltype(points)&>(type, points) ==
           std::pair<JunctionType::Type, const decltype(points)&>(other.type,
                                                                  other.points);
  }

  void throw_on_invalid() const;

  bool is_valid() const;

private:
  void init();
};

unsigned char junction_type_to_char(JunctionType::Type t);

JunctionType::Type junction_type_from_char(unsigned char t);

[[nodiscard]] inline Vec2 centroid(const Junction& junc, const span<const Stroke> edges) {
  auto acc = Vec2(0, 0);
  for (const auto& [si, arclen] : junc.points) {
    const auto l = edges[si].length();
    acc += edges[si].pos(arclen * l);
  }
  return (1 / (Float)junc.points.size()) * acc;
}

[[nodiscard]] inline Vec2 centroid(const Junction& junc, const PolylineBVH& edges) {
  auto acc = Vec2(0, 0);
  for (const auto& [si, arclen] : junc.points) {
    const auto l = edges.nodes[si].geometry->length();
    acc += edges.nodes[si].geometry->pos(arclen * l);
  }
  return (1 / (Float)junc.points.size()) * acc;
}

struct ExtendedJunction {
  struct Range {
    Float start_;
    /// Is only defined for centerline intersections, where preserving the exact position
    /// is important.
    Float mid_;
    Float end_;

    Range(Float start, Float mid, Float end)
      : start_(std::clamp(start, 0.0, 1.0))
      , mid_(std::clamp(mid, 0.0, 1.0))
      , end_(std::clamp(end, 0.0, 1.0)) {}

    Range(Float start, Float end)
      : start_(std::clamp(start, 0.0, 1.0))
      , mid_(-1)
      , end_(std::clamp(end, 0.0, 1.0)) {}

    Float mid() const { return (mid_ >= 0.0 ? mid_ : (0.5 * (start_ + end_))); }

    bool overlaps(const Range& other) const;

    void merge(const Range& other);

    std::string repr() const;
  };

  ExtendedJunction() = default;

  explicit ExtendedJunction(std::vector<std::pair<int, Range>>&& ranges)
    : ranges_(std::move(ranges)) {}

  ExtendedJunction(const Junction&);

  // List of (stroke index, covered range).
  std::vector<std::pair<int, Range>> ranges_;

  void merge(const ExtendedJunction& other);

  void resolve_overlaps();

  Vec2 position(span<const Stroke> strokes) const;

  Float junc_parameter(size_t s_idx) const;

  bool is_grazing() const;
};

struct StrokeCoverage {
  explicit StrokeCoverage(span<const ExtendedJunction> junctions);

  struct Territory {
    ExtendedJunction::Range range_;
    const ExtendedJunction* junction_;
  };

  /// Maps stroke index to a series of regions covered by junctions.
  std::unordered_map<int, std::vector<Territory>> coverage_;

  std::vector<const ExtendedJunction*> junctions() const;

  std::string repr() const;

  std::string junctions_repr() const;

private:
  std::unique_ptr<ExtendedJunction[]> split_junctions_;
  static void resolve_overlaps(std::vector<Territory>& territories);
};

void intersection_junctions(const PolylineBVH& bvh, std::vector<Junction>& out_junctions);

void intersection_junctions(const PolylineBVH& bvh,
                            std::vector<ExtendedJunction>& out_junctions);

std::string to_json(const Junction& junc);
void to_json(const std::vector<Junction>& junctions, const std::string& json_path);

} // namespace sketching
