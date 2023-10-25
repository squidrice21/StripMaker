#include "resample.h"

#include "detail/util.h"
#include "force_assert.h"
#include "sketching.h"

#include <limits>

namespace sketching {

namespace {

/// Return the square of the shortest distance from p to the line segment defined by
/// start-end.
Float shortest_squared_dist(const Vec2& p, const Vec2& start, const Vec2& end) {
  const auto ab = end - start;
  const auto ac = p - start;
  const auto e = ac.dot(ab);
  if (e <= 0) {
    return ac.squaredNorm();
  }
  const auto f = ab.squaredNorm();
  if (e >= f) {
    return (p - end).squaredNorm();
  }
  return ac.squaredNorm() - e * e / f;
}

/// Variant on shortest_dist which keeps small hooks.
Float hook_preserving_squared_dist(const Vec2& p, const Vec2& start, const Vec2& end,
                                   Float& out_t) {
  const auto nsq = (end - start).squaredNorm();
  if (nsq > std::numeric_limits<Float>::epsilon()) {
    out_t = (p - start).dot(end - start) / nsq;
    if (0.0 <= out_t && out_t <= 1.0) {
      const auto proj = lerp(start, end, out_t);
      return (proj - p).squaredNorm();
    } else {
      return infinity; // We care more about t here.
    }
  }
  out_t = 0.0;
  return (start - p).squaredNorm();
}

// @optimize Use squared tolerance values everywhere.
void ramer_douglas_peucker(const Float x[], const Float y[], const size_t begin,
                           const size_t end, const double sq_tolerance, bool include[]) {
  force_assert(end - begin > 1);

  // Find the point with the maximum distance from line between start and end.
  auto dmax = 0.0;
  size_t index = 0;
  for (size_t i = begin + 1; i < end - 1; i++) {
    const auto d =
      shortest_squared_dist({x[i], y[i]}, {x[begin], y[begin]}, {x[end - 1], y[end - 1]});
    if (d > dmax) {
      index = i;
      dmax = d;
    }
  }

  // If max distance is greater than epsilon, recursively simplify.
  if (dmax > sq_tolerance) {
    // FIXME: We could overflow the stack here.
    ramer_douglas_peucker(x, y, begin, index + 1, sq_tolerance, include);
    ramer_douglas_peucker(x, y, index, end, sq_tolerance, include);
  } else {
    // Just include start and end points.
    include[begin] = true;
    include[end - 1] = true;
  }
}

void hook_preserving_ramer_douglas_peucker(const Stroke& points, const size_t begin,
                                           const size_t end, const double sq_tolerance,
                                           bool include[]) {
  force_assert(end - begin > 1);

  // Find the point with the maximum distance from line between start and end.
  auto dmax = Float(0.0);
  auto dtmax = Float(0.0);
  size_t index = 0;
  for (size_t i = begin + 1; i < end - 1; i++) {
    Float t;
    const auto d =
      hook_preserving_squared_dist(points.xy(i), points.xy(begin), points.xy(end - 1), t);
    if (d > dmax) {
      index = i;
      dmax = d;
    } else if (dmax == infinity && d == infinity) {
      // Compare t's.
      const auto dt = (t < 0.0 ? -t : t - 1.0);
      if (dt > dtmax) {
        index = i;
        dtmax = dt;
      }
    }
  }

  // If max distance is greater than epsilon, recursively simplify.
  if (dmax > sq_tolerance) {
    // FIXME: We could overflow the stack here.
    hook_preserving_ramer_douglas_peucker(points, begin, index + 1, sq_tolerance,
                                          include);
    hook_preserving_ramer_douglas_peucker(points, index, end, sq_tolerance, include);
  } else {
    // Just include start and end points.
    include[begin] = true;
    include[end - 1] = true;
  }
}

} // namespace

bool remove_duplicate_vertices(Stroke& s) {
  if (s.size() <= 1) {
    return false;
  }
  auto to_include = std::vector<Index>();
  to_include.reserve(s.size());
  to_include.push_back(0);
  auto last = to_include.back();
  for (Index i = 1; i < s.size(); ++i) {
    const auto r0r1 = 0.5 * std::abs(s.width(i) - s.width(last));
    if ((s.xy(i) - s.xy(last)).squaredNorm() >
        std::max(square(r0r1 * (1 + 1e-3)), 1e-5)) {
      to_include.push_back(i);
      last = i;
    } else if (s.width(i) > s.width(last)) {
      // Replace with the eclipsing vertex.
      to_include.back() = i;
      last = i;
    }
  }
  // Make a backward pass now
  // last = to_include.back(); // Not necessary.
  for (Index i = to_include.size() - 2; i >= 0; --i) {
    const auto cur = to_include[i];
    const auto r0r1 = 0.5 * std::abs(s.width(cur) - s.width(last));
    if ((s.xy(cur) - s.xy(last)).squaredNorm() <=
        std::max(square(r0r1 * (1 + 1e-3)), 1e-5)) {
      to_include[i] = -1; // Mark for deletion.
    } else {
      last = cur;
    }
  }
  to_include.erase(std::remove(to_include.begin(), to_include.end(), -1),
                   to_include.end());
  size_t i = 0;
  for (; i < to_include.size(); ++i) {
    const auto src = to_include[i];
    s.x(i) = s.x(src);
    s.y(i) = s.y(src);
    s.width(i) = s.width(src);
    if (s.has_time())
      s.time(i) = s.time(src);
  }
  if ((Index)i < s.size()) {
    s.resize(i);
    return true;
  }
  return false;
}

Stroke duplicate_vertices_removed(const Stroke& s) {
  if (s.size() <= 1) {
    return s.clone();
  }
  auto to_include = std::vector<Index>();
  to_include.reserve(s.size());
  to_include.push_back(0);
  auto last = to_include.back();
  for (Index i = 1; i < s.size(); ++i) {
    const auto r0r1 = 0.5 * std::abs(s.width(i) - s.width(last));
    if ((s.xy(i) - s.xy(last)).squaredNorm() >
        std::max(square(r0r1 * (1 + 1e-3)), 1e-5)) {
      to_include.push_back(i);
      last = i;
    } else if (s.width(i) > s.width(last)) {
      // Replace with the eclipsing vertex.
      to_include.back() = i;
      last = i;
    }
  }
  // Make a backward pass now
  // last = to_include.back(); // Not necessary.
  for (Index i = to_include.size() - 2; i >= 0; --i) {
    const auto cur = to_include[i];
    const auto r0r1 = 0.5 * std::abs(s.width(cur) - s.width(last));
    if ((s.xy(cur) - s.xy(last)).squaredNorm() <=
        std::max(square(r0r1 * (1 + 1e-3)), 1e-5)) {
      to_include[i] = -1; // Mark for deletion.
    } else {
      last = cur;
    }
  }
  to_include.erase(std::remove(to_include.begin(), to_include.end(), -1),
                   to_include.end());
  auto out_s = Stroke(to_include.size(), s.has_time());
  for (size_t i = 0; i < to_include.size(); ++i) {
    const auto src = to_include[i];
    out_s.x(i) = s.x(src);
    out_s.y(i) = s.y(src);
    out_s.width(i) = s.width(src);
    if (s.has_time())
      out_s.time(i) = s.time(src);
  }
  return out_s;
}

void ramer_douglas_peucker(const Stroke& s, const double tolerance, bool out_include[]) {
  const auto n = s.size();
  if (n < 3) {
    for (Index i = 0; i < n; ++i) {
      out_include[i] = true;
    }
    return;
  }
  for (auto i = 0; i < n; ++i) {
    out_include[i] = false;
  }
  ramer_douglas_peucker(s.x_, s.y_, 0, s.size(), /*sq_tolerance=*/tolerance * tolerance,
                        out_include);
  assert(out_include[0] && "endpoint missing");
  assert(out_include[n - 1] && "endpoint missing");
}

void ramer_douglas_peucker(const Float x[], const Float y[], const Index n,
                           const double tolerance, bool out_include[]) {
  if (n < 3) {
    for (Index i = 0; i < n; ++i) {
      out_include[i] = true;
    }
    return;
  }
  for (auto i = 0; i < n; ++i) {
    out_include[i] = false;
  }
  ramer_douglas_peucker(x, y, 0, n, /*sq_tolerance=*/tolerance * tolerance, out_include);
  assert(out_include[0] && "endpoint missing");
  assert(out_include[n - 1] && "endpoint missing");
}

void corners_and_hooks(const Stroke& s, double tolerance, bool out_include[]) {
  const auto n = s.size();
  if (n < 3) {
    for (Index i = 0; i < n; ++i) {
      out_include[i] = true;
    }
    return;
  }
  for (auto i = 0; i < n; ++i) {
    out_include[i] = false;
  }
  hook_preserving_ramer_douglas_peucker(
    s, 0, s.size(), /*sq_tolerance=*/tolerance * tolerance, out_include);
  assert(out_include[0] && "endpoint missing");
  assert(out_include[n - 1] && "endpoint missing");
}

void decimate(Stroke& stroke, const Float tolerance) {
  const auto n = stroke.size();
  if (n < 3) {
    return;
  }
  auto include = std::make_unique<bool[]>(n);
  ramer_douglas_peucker(stroke, /*tolerance=*/tolerance, include.get());

  auto out_idx = Index(1);
  for (Index i = 1; i < n; ++i) {
    if (include[i]) {
      stroke.x(out_idx) = stroke.x(i);
      stroke.y(out_idx) = stroke.y(i);
      stroke.width(out_idx) = stroke.width(i);
      if (stroke.has_time()) {
        stroke.time(out_idx) = stroke.time(i);
      }
      out_idx++;
    }
  }
  stroke.resize(out_idx);
}

void decimate(CoordMat& xy, const Float tolerance) {
  const auto n = xy.rows();
  if (n < 3) {
    return;
  }
  auto include = std::make_unique<bool[]>(n);
  auto xybuf = std::make_unique<Float[]>(2 * n);
  auto x = &xybuf[0];
  auto y = &xybuf[n];
  for (Index i = 0; i < n; ++i) {
    x[i] = xy(i, 0);
    y[i] = xy(i, 1);
  }

  ramer_douglas_peucker(x, y, n, tolerance, include.get());

  auto out_idx = Index(1);
  for (Index i = 1; i < n; ++i) {
    if (include[i]) {
      xy(out_idx, 0) = xy(i, 0);
      xy(out_idx, 1) = xy(i, 1);
      out_idx++;
    }
  }
  xy.conservativeResize(out_idx, Eigen::NoChange);
}

Stroke decimated(const Stroke& s, const Float tolerance) {
  const auto n = s.size();
  if (n < 3) {
    return s.clone();
  }
  auto include = std::make_unique<bool[]>(n);
  ramer_douglas_peucker(s, tolerance, include.get());

  auto new_size = Index(0);
  for (auto i = 0; i < n; ++i) {
    if (include[i])
      new_size++;
  }
  auto decimated = Stroke(new_size, s.has_time());
  auto out_idx = Index(0);
  for (auto i = 0; i < n; ++i) {
    if (include[i]) {
      decimated.x(out_idx) = s.x(i);
      decimated.y(out_idx) = s.y(i);
      decimated.width(out_idx) = s.width(i);
      if (s.has_time()) {
        decimated.time(out_idx) = s.time(i);
      }
      out_idx++;
    }
  }
  return decimated;
}

} // namespace sketching
