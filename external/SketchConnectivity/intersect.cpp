#include "intersect.h"

#include "bvh.h"
#include "detail/util.h"
#include "fitting.h"

#include <map>

namespace sketching {

namespace {

bool within_interval(Float x, Float min, Float max) {
  if (min > max)
    std::swap(min, max);
  return min <= x && x <= max;
}

void remove_duplicate_intersection_junctions(std::vector<Vec2> &arclens) {
  constexpr auto tolerance = 1e-6;
  for (auto i = 0u; i < arclens.size(); ++i) {
    for (auto j = i + 1; j < arclens.size(); ++j) {
      const auto diff = arclens[i] - arclens[j];
      if (std::abs(diff.x_) < tolerance && std::abs(diff.y_) < tolerance) {
        // Keep j but not i. This is a little arbitrary...
        arclens[i].x_ = infinity; // signal for deletion
      }
    }
  }
  arclens.erase(std::remove_if(arclens.begin(), arclens.end(),
                               [](Vec2 p) { return p.x_ == infinity; }),
                arclens.end());
}

/// Project point p onto segment ab.
double project_pt_on_segment(const Vec2 &p, const Vec2 &a, const Vec2 &b) {
  auto t =
    ((p.x() - a.x()) * (b.x() - a.x()) + (p.y() - a.y()) * (b.y() - a.y())) /
    (square(b.x() - a.x()) + square(b.y() - a.y()));
  return std::clamp(t, 0.0, 1.0);
}

/// Return arclength t on the stroke corresponding to the closest point if there
/// is an intersection, or -infinity if there is none.
double intersect_circle_with_envelope(const Vec2 &p, const double r,
                                      const EnvelopeBVHLeaf &env) {
  const auto p_bb = BoundingBox(p.x() - r, p.x() + r, p.y() - r, p.y() + r);
  if (!p_bb.touches(env.bb)) {
    return -infinity;
  }
  const auto &geom = *env.geometry;
  const auto n = geom.size();
  for (int i = 0; i < n - 1; ++i) {
    const auto a = geom.xy(i);
    const auto b = geom.xy(i + 1);
    const auto u = project_pt_on_segment(p, a, b);
    const auto q = lerp(a, b, u);
    const auto r_geom = 0.5 * fast_lerp(geom.width(i), geom.width(i + 1), u);
    const auto sqnorm = (p - q).squaredNorm();
    if (sqnorm <= square((r_geom + r) * (1 + 1e-3))) {
      return clamped_lerp(geom.arclength(i), geom.arclength(i + 1), u);
    }
  }
  return -infinity;
}

double intersect_circle_with_self_envelope(const int at, const Stroke &geom) {
  const auto n = geom.size();
  const auto p = geom.xy(at);
  const auto r = 0.5 * geom.width(at);
  // Check after.
  for (int i = at + 1; i < n - 1; ++i) {
    const auto a = geom.xy(i);
    const auto b = geom.xy(i + 1);
    const auto u = project_pt_on_segment(p, a, b);
    const auto q = lerp(a, b, u);
    const auto r_geom = 0.5 * fast_lerp(geom.width(i), geom.width(i + 1), u);
    const auto sqnorm = (p - q).squaredNorm();
    if (sqnorm <= square((r_geom + r) * (1 + 1e-3))) {
      const auto s = clamped_lerp(geom.arclength(i), geom.arclength(i + 1), u);
      // Check that it isn't spurious.
      if (3 * std::sqrt(sqnorm) < s - geom.arclength(at)) {
        return s;
      }
    }
  }
  if (at == n - 3) {
    // Check last vertex.
    const auto norm = (geom.xy(at) - geom.xy(n - 1)).norm();
    const auto r0r1 = 0.5 * (geom.width(at) + geom.width(n - 1));
    if (norm < r0r1 * (1 + 1e-3) &&
        3 * norm < geom.arclength(n - 1) - geom.arclength(at)) {
      return geom.arclength(n - 1);
    }
  }
  // Check before.
  for (int i = at - 2; i >= 0; --i) {
    const auto a = geom.xy(i);
    const auto b = geom.xy(i + 1);
    const auto u = project_pt_on_segment(p, a, b);
    const auto q = lerp(a, b, u);
    const auto r_geom = 0.5 * fast_lerp(geom.width(i), geom.width(i + 1), u);
    const auto sqnorm = (p - q).squaredNorm();
    if (sqnorm <= square((r_geom + r) * (1 + 1e-3))) {
      const auto s = clamped_lerp(geom.arclength(i), geom.arclength(i + 1), u);
      // Check that it isn't spurious.
      if (3 * std::sqrt(sqnorm) < geom.arclength(at) - s) {
        return s;
      }
    }
  }
  if (at == 2) {
    // Check first vertex.
    const auto norm = (geom.xy(at) - geom.xy(0)).norm();
    const auto r0r1 = 0.5 * (geom.width(at) + geom.width(0));
    if (norm < r0r1 * (1 + 1e-3) && 3 * norm < geom.arclength(at)) {
      return 0.0; // geom.arclength(0) in other words.
    }
  }
  return -infinity;
}

void intersect_segment_with_polyline(const Vec2 &a, const Vec2 &b,
                                     const EnvelopeBVHLeaf &polyline,
                                     std::vector<Vec2> &hits) {
  const auto ab_bb =
    BoundingBox(std::min(a.x(), b.x()), std::max(a.x(), b.x()),
                std::min(a.y(), b.y()), std::max(a.y(), b.y()));
  if (!ab_bb.touches(polyline.bb)) {
    return;
  }
  const auto &geom = *polyline.geometry;
  const auto n = geom.size();
  auto last = Vec2(infinity, infinity);
  for (int i = 0; i < n - 1; ++i) {
    const auto c = geom.xy(i);
    const auto d = geom.xy(i + 1);
    double u, v;
    const auto hit = segment_intersection(a, b, c, d, u, v);
    if (hit) {
      u /= (a - b).norm();
      if (u < 1.0 && v < (c - d).norm()) { // Do not include end of interval.
        const auto contender =
          Vec2(u, std::clamp(geom.arclength(i) + v, geom.arclength(i),
                             geom.arclength(i + 1)));
        // Avoid adding same intersection twice.
        if ((contender - last).squaredNorm() > 1e-6) {
          hits.emplace_back(contender.x_, contender.y_);
          last = contender;
        }
      }
    }
  }
}

bool segment_ray_intersection(const Vec2 &p, const Vec2 &q, const Vec2 &origin,
                              const Vec2 &dir, Float &out_s, Float &out_t) {
  const auto p0_x = p.x();
  const auto p0_y = p.y();
  const auto p1_x = q.x();
  const auto p1_y = q.y();
  const auto p2_x = origin.x();
  const auto p2_y = origin.y();
  const Vec2 origin_plus_dir = origin + dir;
  const auto p3_x = origin_plus_dir.x();
  const auto p3_y = origin_plus_dir.y();

  const Float s1_x = p1_x - p0_x;
  const Float s1_y = p1_y - p0_y;
  const Float s2_x = p3_x - p2_x;
  const Float s2_y = p3_y - p2_y;

  const Float s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) /
                  (-s2_x * s1_y + s1_x * s2_y);
  const Float t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) /
                  (-s2_x * s1_y + s1_x * s2_y);

  if (s >= 0.0 && t >= 0.0 && s <= 1.0 && t <= 1.0) {
    out_t = s;
    out_s = t;
    return true;
  }
  return false;
}

} // namespace

static bool intersect_segment_stroke_exclusive(const Vec2 p, const Vec2 q,
                                               const PolylineBVHLeaf &stroke,
                                               const Float min_s,
                                               const Float max_s, Float &out_s,
                                               Float &out_t) {
  // Note that the BoundingBox constructor does the appropriate comparisons to
  // make this a valid bounding box.
  const auto pq_bb = BoundingBox(p.x_, q.x_, p.y_, q.y_);
  if (stroke.bb.intersects(pq_bb)) {
    const auto pq_norm = (q - p).norm();
    const auto &geom = *stroke.geometry;
    const auto n = geom.size();
    for (auto i = 0; i < n - 1; ++i) {
      const auto pp = geom.xy(i);
      const auto qq = geom.xy(i + 1);
      Float s, t;
      const bool hit = segment_intersection(p, q, pp, qq, s, t);
      if (hit) {
        const auto sw = s / pq_norm;
        if (sw > min_s && sw < max_s) {
          geom.ensure_arclengths();
          out_s = sw;
          out_t = clamped_lerp(geom.arclength(i), geom.arclength(i + 1),
                               t / (p - q).norm());
          return true;
        }
      }
    }
  }
  return false;
}

static void intersect_segment_stroke_exclusive(const Vec2 p, const Vec2 q,
                                               const PolylineBVHLeaf &stroke,
                                               const Float min_s,
                                               const Float max_s,
                                               std::vector<Vec2> &out_arclens) {
  // Note that the BoundingBox constructor does the appropriate comparisons to
  // make this a valid bounding box.
  const auto pq_bb = BoundingBox(p.x_, q.x_, p.y_, q.y_);
  if (stroke.bb.intersects(pq_bb)) {
    const auto pq_norm = (q - p).norm();
    const auto &geom = *stroke.geometry;
    const auto n = geom.size();
    for (auto i = 0; i < n - 1; ++i) {
      const auto pp = geom.xy(i);
      const auto qq = geom.xy(i + 1);
      Float s, t;
      const bool hit = segment_intersection(p, q, pp, qq, s, t);
      if (hit) {
        const auto sw = s / pq_norm;
        if (sw > min_s && sw < max_s) {
          geom.ensure_arclengths();
          auto out_s = sw;
          auto out_t = clamped_lerp(geom.arclength(i), geom.arclength(i + 1),
                                    t / (p - q).norm());
          out_arclens.emplace_back(out_s, out_t);
        }
      }
    }
  }
}

bool intersect_segment_stroke_exclusive(const Vec2 p, const Vec2 q,
                                        const PolylineBVHLeaf &stroke,
                                        Float &out_s, Float &out_t) {
  const auto eps = Float(1e-5);
  return intersect_segment_stroke_exclusive(p, q, stroke, /*min_s=*/eps,
                                            /*max_s=*/1.0 - eps, out_s, out_t);
}

void intersect_segment_stroke_exclusive(Vec2 p, Vec2 q,
                                        const PolylineBVHLeaf &stroke,
                                        std::vector<Vec2> &out_arclens) {
  const auto eps = Float(1e-5);
  intersect_segment_stroke_exclusive(p, q, stroke, /*min_s=*/eps,
                                     /*max_s=*/1.0 - eps, out_arclens);
  remove_duplicate_intersection_junctions(out_arclens);
}

Bicapsule stroke_segment(const Stroke &s, const Index i) {
  return {s.xy(i), s.xy(i + 1), 0.5 * s.width(i), 0.5 * s.width(i + 1)};
}

// From https://www.shadertoy.com/view/4lcBWn (MIT).
Vec2 signed_dist_bicapsule(Vec2 p, Bicapsule bicap) {
  p.x_ -= bicap.p_a.x_;
  p.y_ -= bicap.p_a.y_;
  bicap.p_b.x_ -= bicap.p_a.x_;
  bicap.p_b.y_ -= bicap.p_a.y_;
  const auto h = bicap.p_b.squaredNorm();
  auto q =
    (1 / h) * Vec2(p.dot(Vec2(bicap.p_b.y_, -bicap.p_b.x_)), p.dot(bicap.p_b));

  q.x_ = std::abs(q.x_);

  const auto b = bicap.r_a - bicap.r_b;
  const auto c = Vec2(std::sqrt(h - b * b), b);

  const auto k = c.cross(q);
  const auto m = c.dot(q);
  const auto n = q.dot(q);

  if (k < 0.0) {
    // Closest is a.
    return {std::sqrt(h * (n)) - bicap.r_a, k};
  } else if (k > c.x_) {
    // Closest is b.
    return {std::sqrt(h * (n + 1.0 - 2.0 * q.y_)) - bicap.r_b, k};
  } else {
    // Closest is sides.
    return {m - bicap.r_a, k};
  }
}

Float signed_dist_stroke(const Vec2 p, const Stroke &stroke) {
  const auto n1 = stroke.size() - 1;
  auto sd = infinity;
  for (Index i = 0; i < n1; ++i) {
    const auto [env_dist, _u] =
      signed_dist_bicapsule(p, stroke_segment(stroke, i));
    sd = std::min(sd, env_dist);
  }
  return sd;
}

bool line_of_sight(Vec2 p, Vec2 q, span<const Stroke> strokes,
                   span<const BoundingBox> centerline_bbs) {
  assert(strokes.size() == centerline_bbs.size());
  Float s, t;
  const auto n = strokes.size();
  for (size_t i = 0; i < n; ++i) {
    if (intersect_segment_stroke_exclusive(
          p, q, {strokes[i], centerline_bbs[i]}, s, t)) {
      return false;
    }
  }
  return true;
}

bool line_of_sight(Vec2 p, Vec2 q, const PolylineBVH &bvh, const Float min_s,
                   const Float max_s) {
  assert(min_s <= 1.0 && "max_s should satisfy 0 <= min_s <= 1");
  assert(max_s <= 1.0 && "max_s should satisfy 0 <= max_s <= 1");

  Float s, t;
  for (const auto &node : bvh.nodes) {
    if (intersect_segment_stroke_exclusive(p, q, node, min_s, max_s, s, t))
      return false;
  }
  return true;
}

void intersect_self(const Stroke &geom, std::vector<Vec2> &out) {
  const auto n = geom.size();
  if (n < 3)
    return;
  for (int i = 0; i < n - 1; ++i) {
    const auto a = geom.xy(i);
    const auto b = geom.xy(i + 1);
    for (int j = i + 2; j < n - 1; ++j) {
      const auto c = geom.xy(j);
      const auto d = geom.xy(j + 1);
      Float s, t;
      if (segment_intersection(a, b, c, d, s, t)) {
        out.emplace_back(std::clamp(geom.arclength(i) + s, geom.arclength(i),
                                    geom.arclength(i + 1)),
                         std::clamp(geom.arclength(j) + t, geom.arclength(j),
                                    geom.arclength(j + 1)));
      }
    }
  }
  remove_duplicate_intersection_junctions(out);
}

void intersect_self(const PolylineBVHLeaf &node,
                    std::vector<Vec2> &out_arclens) {
  return intersect_self(*node.geometry, out_arclens);
}

void intersect_different(const PolylineBVHLeaf &node1,
                         const PolylineBVHLeaf &node2, std::vector<Vec2> &out) {
  assert(node1.geometry != node2.geometry);
  if (!node1.bb.touches(node2.bb)) {
    return;
  }
  const auto &geom1 = *node1.geometry;
  const auto &geom2 = *node2.geometry;
  const auto n1 = geom1.size();
  const auto n2 = geom2.size();
  for (int i = 0; i < n1 - 1; ++i) {
    const auto a = geom1.xy(i);
    const auto b = geom1.xy(i + 1);
    for (int j = 0; j < n2 - 1; ++j) {
      const auto c = geom2.xy(j);
      const auto d = geom2.xy(j + 1);
      Float s, t;
      if (segment_intersection(a, b, c, d, s, t)) {
        out.emplace_back(std::clamp(geom1.arclength(i) + s, geom1.arclength(i),
                                    geom1.arclength(i + 1)),
                         std::clamp(geom2.arclength(j) + t, geom2.arclength(j),
                                    geom2.arclength(j + 1)));
      }
    }
  }
  remove_duplicate_intersection_junctions(out);
}

void intersect_self(const EnvelopeBVHLeaf &node,
                    std::vector<IntersectionEvent> &out) {
  const auto &geom = *node.geometry;
  const auto n = geom.size();
  if (n < 4)
    return;
  for (int i = 0; i < n - 1; ++i) {
    const auto a = geom.xy(i);
    const auto b = geom.xy(i + 1);
    for (int j = i + 2; j < n - 1; ++j) {
      const auto c = geom.xy(j);
      const auto d = geom.xy(j + 1);
      Float u, v;
      if (segment_intersection(a, b, c, d, u, v)) {
        const auto s = std::clamp(geom.arclength(i) + u, geom.arclength(i),
                                  geom.arclength(i + 1));
        const auto t = std::clamp(geom.arclength(j) + v, geom.arclength(j),
                                  geom.arclength(j + 1));
        auto current_event = IntersectionEvent();
        current_event.s_start = s;
        current_event.t_start = t;
        current_event.s_mid = s;
        current_event.t_mid = t;
        current_event.s_end = s;
        current_event.t_end = t;
        out.emplace_back(current_event);
      }
    }
  }
  for (auto &x : out) {
    x.normalize();
  }

  auto current_event = IntersectionEvent();
  const auto pw = geom.pen_width();
  const auto try_to_add_current_event = [&]() {
    for (auto &event : out) {
      // Check for overlap.
      if ((event.s_start <= current_event.s_end &&
           current_event.s_start <= event.s_end) ||
          (event.t_start <= current_event.t_end &&
           current_event.t_start <= event.t_end)) {
        event.s_start = std::min(event.s_start, current_event.s_start);
        event.s_end = std::max(event.s_end, current_event.s_end);
        event.t_start = std::min(event.t_start, current_event.t_start);
        event.t_end = std::max(event.t_end, current_event.t_end);
        goto skip;
      }
      // Note that s and t lie on the same stroke in this case, so the event's s
      // interval might correspond to the current_event's t interval, or vice
      // versa.  Be careful here...
      if ((event.t_start <= current_event.s_end &&
           current_event.s_start <= event.t_end) ||
          (event.s_start <= current_event.t_end &&
           current_event.t_start <= event.s_end)) {
        event.t_start = std::min(event.t_start, current_event.s_start);
        event.t_end = std::max(event.t_end, current_event.s_end);
        event.s_start = std::min(event.s_start, current_event.t_start);
        event.s_end = std::max(event.s_end, current_event.t_end);
        goto skip;
      }
    }
    for (const auto &event : out) {
      // Check if too close to existing intersection.
      if ((event.s_start <= current_event.s_end + pw &&
           current_event.s_start <= event.s_end + pw) ||
          (event.t_start <= current_event.t_end + pw &&
           current_event.t_start <= event.t_end + pw) ||
          (event.t_start <= current_event.s_end + pw &&
           current_event.s_start <= event.t_end + pw) ||
          (event.s_start <= current_event.t_end + pw &&
           current_event.t_start <= event.s_end + pw)) {
        goto skip;
      }
    }
    out.emplace_back(current_event);
    out.back().normalize();
skip:
    current_event = IntersectionEvent();
  };

  // Find envelope intersections.
  for (int i = 0; i < n; ++i) {
    auto s = intersect_circle_with_self_envelope(i, geom);
    if (s >= 0.0) {
      if (!current_event) {
        current_event.s_start = s;
        current_event.t_start = geom.arclength(i);
      }
      current_event.s_end = s;
      current_event.t_end = geom.arclength(i);
    } else if (s < 0.0 && current_event) {
      try_to_add_current_event();
    }
  }
  if (current_event) {
    try_to_add_current_event();
  }
}

void intersect_different(const EnvelopeBVHLeaf &node1,
                         const EnvelopeBVHLeaf &node2,
                         std::vector<IntersectionEvent> &out) {
  assert(out.empty());
  assert(node1.geometry != node2.geometry);
  if (!node1.bb.touches(node2.bb)) {
    return;
  }

  // Intersect circles around vertices of geom1 with geom2.
  const auto &geom1 = *node1.geometry;
  const auto n1 = geom1.size();
  auto is_hit1_t = std::vector<Float>(n1, -infinity);
  for (int i = 0; i < n1; ++i) {
    is_hit1_t[i] =
      intersect_circle_with_envelope(geom1.xy(i), 0.5 * geom1.width(i), node2);
  }

  auto current_event = IntersectionEvent();

  // Intersect circles around vertices of geom2 with geom1 and record contiguous
  // runs.
  const auto &geom2 = *node2.geometry;
  const auto n2 = geom2.size();
  auto is_hit2 = std::vector<bool>(n2, false);
  auto reverse_events = std::vector<IntersectionEvent>();
  for (int i = 0; i < n2; ++i) {
    const auto s =
      intersect_circle_with_envelope(geom2.xy(i), 0.5 * geom2.width(i), node1);
    if (s >= 0.0) {
      is_hit2[i] = true;
      if (current_event) {
        // Check that this is actually the same event and that there wasn't a
        // discontinuous break in the connection with the other stroke.
        for (int j = 0; j < n1; ++j) {
          if (within_interval(geom1.arclength(j), current_event.s_start, s) &&
              is_hit1_t[j] == -infinity) {
            // There was a break.
            reverse_events.emplace_back(current_event);
            current_event = IntersectionEvent();
            break;
          }
          if (geom1.arclength(j) > std::max(current_event.s_start, s)) {
            break;
          }
        }
      }
      if (!current_event) {
        current_event.s_start = s;
        current_event.t_start = geom2.arclength(i);
      }
      current_event.s_end = s;
      current_event.t_end = geom2.arclength(i);
    } else if (s < 0.0 && current_event) {
      reverse_events.emplace_back(current_event);
      current_event = IntersectionEvent();
    }
  }
  if (current_event) {
    reverse_events.emplace_back(current_event);
    current_event = IntersectionEvent();
  }
  // Not sure if this normalization is necessary...
  for (auto &rev_event : reverse_events) {
    rev_event.normalize();
  }

  auto hits = std::vector<Vec2>();
  for (int i = 0; i < n1 - 1; ++i) {
    // Check if the start of the segment intersects with the other stroke.
    {
      const auto t = is_hit1_t[i];
      if (t >= 0.0) {
        if (current_event) {
          // Check that this is actually the same event and that there wasn't a
          // discontinuous break in the connection with the other stroke.
          for (int j = 0; j < n2; ++j) {
            if (within_interval(geom2.arclength(j), current_event.t_start, t) &&
                !is_hit2[j]) {
              // There was a break.
              out.emplace_back(current_event);
              current_event = IntersectionEvent();
              break;
            }
            if (geom2.arclength(j) > std::max(current_event.t_start, t)) {
              break;
            }
          }
        }
        if (!current_event) {
          current_event.s_start = geom1.arclength(i);
          current_event.t_start = t;
        }
        current_event.s_end = geom1.arclength(i);
        current_event.t_end = t;
      } else if (t < 0.0 && current_event) {
        out.emplace_back(current_event);
        current_event = IntersectionEvent();
      }
    }

    // Check if the interior of the segment intersects with the other stroke.
    {
      intersect_segment_with_polyline(geom1.xy(i), geom1.xy(i + 1), node2,
                                      hits);
      for (const auto &[u, t] : hits) {
        const auto s =
          clamped_lerp(geom1.arclength(i), geom1.arclength(i + 1), u);
        if (current_event.centerline_intersection() &&
            // The same intersection can be found twice across different
            // iterations due to numerical imprecision.
            std::abs(current_event.s_mid - s) >= 1e-3) {
          out.emplace_back(current_event);
          current_event = IntersectionEvent();
        }
        // Centerline intersections happen exactly at the intersection point.
        if (!current_event) {
          current_event.s_start = s;
          current_event.t_start = t;
        }
        current_event.s_mid = s;
        current_event.t_mid = t;
        current_event.s_end = s;
        current_event.t_end = t;
      }
      hits.clear();
      // Note that a miss does not necessarily mean that we have left the
      // intersection event, since we are intersecting with the *centerline*,
      // not the envelope.
    }

    if (i == n1 - 2) {
      // Check if the end of the segment intersects with the other stroke.
      // This is a special case for the last segment since for the other
      // segments we do not want to check vertex locations twice (one as the end
      // of a segment and once as the start).
      const auto t = is_hit1_t[i + 1];
      if (t >= 0.0) {
        if (current_event) {
          // Check that this is actually the same event and that there wasn't a
          // discontinuous break in the connection with the other stroke.
          for (int j = 0; j < n2; ++j) {
            if (within_interval(geom2.arclength(j), current_event.t_start, t) &&
                !is_hit2[j]) {
              // There was a break.
              out.emplace_back(current_event);
              current_event = IntersectionEvent();
              break;
            }
            if (geom2.arclength(j) > std::max(current_event.t_start, t)) {
              break;
            }
          }
        }
        if (!current_event) {
          current_event.s_start = geom1.arclength(i + 1);
          current_event.t_start = t;
        }
        current_event.s_end = geom1.arclength(i + 1);
        current_event.t_end = t;
      }
      if (current_event) {
        out.emplace_back(current_event);
        current_event = IntersectionEvent();
      }
    }
  }

  for (auto &x : out) {
    x.normalize();
  }

  // We need to make the t_start and t_end values accurate, and add possibly
  // missed envelope-envelope intersections.
  const auto pw = std::max(geom1.pen_width(), geom2.pen_width());
  for (const auto &rev_event : reverse_events) {
    for (auto &event : out) {
      // Check for overlap.
      if (event.s_start <= rev_event.s_end &&
          rev_event.s_start <= event.s_end) {
        // Try to update s, t, but do not break the junction.
        auto attempt = event;
        attempt.s_start = std::min(event.s_start, rev_event.s_start);
        attempt.s_end = std::max(event.s_end, rev_event.s_end);
        attempt.t_start = std::min(event.t_start, rev_event.t_start);
        attempt.t_end = std::max(event.t_end, rev_event.t_end);
        auto broken = false;
        for (int j = 0; j < n2; ++j) {
          if (within_interval(geom2.arclength(j), attempt.t_start,
                              attempt.t_end) &&
              !is_hit2[j]) {
            broken = true;
            break;
          }
          if (geom2.arclength(j) > std::max(attempt.t_start, attempt.t_end)) {
            break;
          }
        }
        if (!broken) {
          for (int j = 0; j < n1; ++j) {
            if (within_interval(geom1.arclength(j), attempt.s_start,
                                attempt.s_end) &&
                is_hit1_t[j] == -infinity) {
              broken = true;
              break;
            }
            if (geom1.arclength(j) > std::max(attempt.s_start, attempt.s_end)) {
              break;
            }
          }
        }
        if (!broken)
          event = attempt;
        goto skip; // Already added.
      }
    }
    for (auto &event : out) {
      if (event.s_start - pw <= rev_event.s_start &&
          rev_event.s_start <= event.s_end + pw) {
        goto skip; // Probably already added, but too close anyhow.
      }
    }
    out.emplace_back(rev_event);
skip:;
  }
}

std::string IntersectionEvent::repr() const {
  return fmt::format("IntersectionEvent(s_start={:.5g}, t_start={:.5g}, "
                     "s_mid={:.5g}, t_mid={:.5g}, s_end={:.5g}, t_end={:.5g})",
                     s_start, t_start, s_mid, t_mid, s_end, t_end);
}

Vec2 find_snap_point(const ConstStrokeView &to_snap, const bool head,
                     const EnvelopeBVHLeaf &other,
                     Float *out_closest_env_dist) {
  const auto p = (head ? to_snap.xy(0) : to_snap.xy(Back));
  auto radius = 0.5 * (head ? to_snap.width(0) : to_snap.width(Back));

  *out_closest_env_dist = INFINITY;
  auto best_dist = infinity;
  auto best_arclen = -infinity;
  if (BoundingBox(p.x_, p.y_, radius).touches(other.bb)) {
    const auto &other_stroke = *other.geometry;
    const auto n = other_stroke.size();
    Index start = 0, stop = n - 1;
    if (&to_snap.stroke() == &other_stroke) {
      if (head) {
        start = 1;
      } else {
        stop = n - 2;
      }
    }
    for (Index i = start; i < stop; ++i) {
      auto [env_dist, u] =
        signed_dist_bicapsule(p, stroke_segment(other_stroke, i));
      if (std::isnan(env_dist)) {
        // The bicapsule is degenerate (i.e. a circle).
        // Skipping is fine if this bicapsule occurs in the interior, because it
        // will be covered by the adjacent bicapsules, however this is a problem
        // if this is an endpoint with the end vertex having the larger radius.
        // I imagine that case is rare, so this probably okay...
        continue;
      }
      u = std::clamp(u, 0.0,
                     other_stroke.arclength(i + 1) - other_stroke.arclength(i));
      const auto arclen =
        std::min(other_stroke.arclength(i) + u, other_stroke.arclength(i + 1));
      u /= other_stroke.arclength(i + 1) - other_stroke.arclength(i);
      const auto q = lerp(other_stroke.xy(i), other_stroke.xy(i + 1), u);
      const auto connection = q - p;
      const auto cen_dist = connection.norm();
      if (&to_snap.stroke() == &other_stroke) {
        // For snapping to self, it needs to be non-spurious.
        const auto otherp_arclen = (head ? 0.0 : to_snap.length());
        const auto arc_dist = std::abs(arclen - otherp_arclen);
        // Keep this synced with others (grep for [spurious-condition]).
        if (!(std::max(3 * cen_dist, M_PI * to_snap.stroke().pen_width()) <
              arc_dist)) {
          continue; // Spurious.
        }
      }
      env_dist -= radius;
      *out_closest_env_dist = std::min(*out_closest_env_dist, env_dist);
      if (cen_dist >= best_dist) {
        continue;
      }
      if (env_dist <= 0) {
        best_dist = cen_dist;
        best_arclen = arclen;
      }
    }
  }
  return {best_dist, best_arclen};
}

Vec2 find_snap_point_legacy(const Stroke &to_snap, const bool head,
                            const EnvelopeBVHLeaf &other) {
  const auto p = (head ? to_snap.xy(0) : to_snap.xy(Back));
  constexpr auto fudge = 0.0;
  auto radius =
    (1 + fudge) * 0.5 * (head ? to_snap.width(0) : to_snap.width(Back));

  auto best_dist = infinity;
  auto best_arclen = -infinity;
  if (BoundingBox(p.x_, p.y_, radius).touches(other.bb)) {
    const auto &other_stroke = *other.geometry;
    const auto n = other_stroke.size();
    Index start = 0, stop = n - 1;
    if (&to_snap == &other_stroke) {
      if (head) {
        start = 1;
      } else {
        stop = n - 2;
      }
    }
    for (Index i = start; i < stop; ++i) {
      auto [env_dist, u] =
        signed_dist_bicapsule(p, stroke_segment(other_stroke, i));
      if (std::isnan(env_dist)) {
        // The bicapsule is degenerate (i.e. a circle).
        // Skipping is fine if this bicapsule occurs in the interior, because it
        // will be covered by the adjacent bicapsules, however this is a problem
        // if this is an endpoint with the end vertex having the larger radius.
        // I imagine that case is rare, so this probably okay...
        continue;
      }
      u = std::clamp(u, 0.0,
                     other_stroke.arclength(i + 1) - other_stroke.arclength(i));
      const auto arclen =
        std::min(other_stroke.arclength(i) + u, other_stroke.arclength(i + 1));
      u /= other_stroke.arclength(i + 1) - other_stroke.arclength(i);
      const auto q = lerp(other_stroke.xy(i), other_stroke.xy(i + 1), u);
      const auto connection = q - p;
      const auto cen_dist = connection.norm();
      if (cen_dist >= best_dist) {
        continue;
      }
      if (&to_snap == &other_stroke) {
        // For snapping to self, it needs to be non-spurious.
        const auto otherp_arclen = (head ? 0.0 : to_snap.length());
        const auto arc_dist = std::abs(arclen - otherp_arclen);
        if (4 * cen_dist >= arc_dist || arc_dist < 2 * to_snap.pen_width()) {
          continue; // Spurious.
        }
      }
      const auto other_w =
        0.5 * fast_lerp(other_stroke.width(i), other_stroke.width(i + 1), u);
      if (env_dist <= fudge * other_w) {
        // Endpoint falls within the other stroke's envelope (within fudge
        // factor).
        goto valid;
      }
      if (env_dist > radius) {
        // Envelopes don't even overlap.
        continue;
      }
valid:
      best_dist = cen_dist;
      best_arclen = arclen;
    }
  }
  return {best_dist, best_arclen};
}

std::vector<std::pair<size_t, Float>>
ray_intersections(const PolylineBVH &polylines, const Vec2 &p, const Vec2 &r) {
  assert(r.squaredNorm() > 1e-6);
  std::multimap<Float, std::pair<size_t, Float>> hits; // sort by distance
  for (size_t i = 0; i < polylines.nodes.size(); ++i) {
    const auto &node = polylines.nodes[i];
    // TODO: BB check for speed
    const auto &geom = *node.geometry;
    for (auto j = 0; j < geom.size() - 1; ++j) {
      const auto segstart = geom.xy(j);
      const auto segend = geom.xy(j + 1);
      Float s, t;
      if (segment_ray_intersection(segstart, segend, p, r, s, t)) {
        const auto narclen =
          ((1.0 - s) * geom.arclength(j) + s * geom.arclength(j + 1)) /
          geom.length();
        hits.insert({t, {i, narclen}});
      }
    }
  }

  std::vector<std::pair<size_t, Float>> out;
  out.reserve(hits.size());
  for (const auto &hit : hits) {
    out.push_back(hit.second);
  }

  return out;
}

} // namespace sketching
