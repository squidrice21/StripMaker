#include "closest.h"

#include "degrees.h"
#include "detail/util.h"
#include "endpoint.h"
#include "fitting.h"
#include "force_assert.h"
#include "intersect.h"
#include "junction.h"
#include "stroke_graph.h"

#include <map>

namespace sketching {

namespace {

Float acos(Float x) { return std::acos(std::min(x, 1.0)); }

} // namespace

Float closest_point(const Stroke& polyline, const Vec2 xy, Vec2& out_proj, Float& out_s) {
  polyline.ensure_arclengths();
  const auto npoints = polyline.size();
  Float best_sq_dist = infinity;
  for (auto i = 0; i < npoints - 1; ++i) {
    const auto p = polyline.xy(i);
    const auto q = polyline.xy(i + 1);
    const Float nsq = (q - p).squaredNorm();
    const auto t = std::clamp<Float>((xy - p).dot(q - p) / nsq, 0.0, 1.0);
    const Vec2 proj = lerp(p, q, t);
    const Float sq_dist = (proj - xy).squaredNorm();
    if (sq_dist < best_sq_dist) {
      best_sq_dist = sq_dist;
      out_proj = proj;
      out_s = clamped_lerp(polyline.arclength(i), polyline.arclength(i + 1), t);
    }
  }
  return std::sqrt(best_sq_dist);
}

Float closest_point_substroke( //
  const Stroke& polyline, const Float start_arclen, const Float end_arclen, const Vec2 xy,
  Vec2& out_proj, Float& out_s) {

  assert(start_arclen <= polyline.length() && "expected absolute arc length value");
  assert(end_arclen <= polyline.length() && "expected absolute arc length value");

  const auto start_frac_pair = polyline.fractional_index(start_arclen);
  const auto start_frac_idx = start_frac_pair.first + start_frac_pair.second;
  const auto end_frac_pair = polyline.fractional_index(end_arclen);
  const auto end_frac_idx = end_frac_pair.first + end_frac_pair.second;
  const auto end = (Index)std::ceil(end_frac_idx);

  Float best_sq_dist = infinity;
  for (auto i = (Index)start_frac_pair.first; i < end; ++i) {
    const auto p = polyline.xy(i);
    const auto q = polyline.xy(i + 1);
    const Float nsq = (q - p).squaredNorm();
    const auto lower = std::max(0.0, start_frac_idx - i);
    const auto upper = std::min(1.0, end_frac_idx - i);
    const auto t = std::clamp<Float>((xy - p).dot(q - p) / nsq, lower, upper);
    const Vec2 proj = lerp(p, q, t);
    const Float sq_dist = (proj - xy).squaredNorm();
    if (sq_dist < best_sq_dist) {
      best_sq_dist = sq_dist;
      out_proj = proj;
      out_s = clamped_lerp(polyline.arclength(i), polyline.arclength(i + 1), t);
    }
  }
  return std::sqrt(best_sq_dist);
}

static Float closest_point_to_own_stroke_interior_strict( //
  const Stroke& polyline, const Float arclen, Vec2& out_proj, Float& out_s) {

  const auto npoints = polyline.size();
  const auto xy = polyline.pos(arclen);
  const auto pen_width = polyline.pen_width();

  Float best_sq_dist = infinity;
  for (auto i = 0; i < npoints - 1; ++i) {
    const auto p = polyline.xy(i);
    const auto q = polyline.xy(i + 1);
    const Float nsq = (q - p).squaredNorm();
    const auto mid_t = std::clamp<Float>((xy - p).dot(q - p) / nsq, 0.0, 1.0);
    const auto t_values =
      (i == npoints - 2 ? std::initializer_list<Float>{0.0, mid_t, 1.0}
                        : std::initializer_list<Float>{0.0, mid_t});
    for (const auto t : t_values) {
      const Vec2 proj = lerp(p, q, t);
      const Float sq_dist = (proj - xy).squaredNorm();
      const auto s = clamped_lerp(polyline.arclength(i), polyline.arclength(i + 1), t);
      if (sq_dist < best_sq_dist &&
          // Check that it isn't spurious.
          // This isn't quite the same as the [spurious-condition].
          // If there are no endpoints involved, we should be more strict.
          4 * std::max(3 * std::sqrt(sq_dist), M_PI * pen_width) < std::abs(arclen - s)) {
        best_sq_dist = sq_dist;
        out_proj = proj;
        out_s = s;
      }
    }
  }
  return std::sqrt(best_sq_dist);
}

// TODO: Is the tangent check here really necessary?  Can we write this in terms of the
//       other overload instead?
Float closest_point_to_own_stroke(const Stroke& polyline, const bool head, Vec2& out_proj,
                                  Float& out_s) {
  const auto npoints = polyline.size();
  const auto xy = (head ? polyline.xy(0) : polyline.xy(Back));
  const auto tangent =
    (head ? head_tangent_lagrange(polyline) : tail_tangent_lagrange(polyline));
  const auto arclen = (head ? 0.0 : polyline.length());
  const auto pen_width = polyline.pen_width();

  Float best_sq_dist = infinity;
  for (auto i = 0; i < npoints - 1; ++i) {
    const auto p = polyline.xy(i);
    const auto q = polyline.xy(i + 1);
    const Float nsq = (q - p).squaredNorm();
    const auto t = std::clamp<Float>((xy - p).dot(q - p) / nsq, 0.0, 1.0);
    const Vec2 proj = lerp(p, q, t);
    const Float sq_dist = (proj - xy).squaredNorm();
    const auto angle = acos((proj - xy).normalized().dot(tangent));
    const auto s = clamped_lerp(polyline.arclength(i), polyline.arclength(i + 1), t);
    if (sq_dist < best_sq_dist && angle < 90_deg &&
        // Check that it isn't spurious.
        // Keep this synced with others (grep for [spurious-condition]).
        std::max(3 * std::sqrt(sq_dist), M_PI * pen_width) < std::abs(arclen - s)) {
      best_sq_dist = sq_dist;
      out_proj = proj;
      out_s = s;
    }
  }

  return std::sqrt(best_sq_dist);
}

Float closest_point_to_own_stroke(const Stroke& polyline, const Float arclen,
                                  Vec2& out_proj, Float& out_s) {
  const auto npoints = polyline.size();
  const auto xy = polyline.pos(arclen);
  const auto pen_width = polyline.pen_width();

  Float best_sq_dist = infinity;
  for (auto i = 0; i < npoints - 1; ++i) {
    const auto p = polyline.xy(i);
    const auto q = polyline.xy(i + 1);
    const Float nsq = (q - p).squaredNorm();
    const auto mid_t = std::clamp<Float>((xy - p).dot(q - p) / nsq, 0.0, 1.0);
    const auto t_values =
      (i == npoints - 2 ? std::initializer_list<Float>{0.0, mid_t, 1.0}
                        : std::initializer_list<Float>{0.0, mid_t});
    for (const auto t : t_values) {
      const Vec2 proj = lerp(p, q, t);
      const Float sq_dist = (proj - xy).squaredNorm();
      const auto s = clamped_lerp(polyline.arclength(i), polyline.arclength(i + 1), t);
      if (sq_dist < best_sq_dist &&
          // Check that it isn't spurious.
          // Keep this synced with others (grep for [spurious-condition]).
          std::max(3 * std::sqrt(sq_dist), M_PI * pen_width) < std::abs(arclen - s)) {
        best_sq_dist = sq_dist;
        out_proj = proj;
        out_s = s;
      }
    }
  }
  return std::sqrt(best_sq_dist);
}

Float closest_point_to_own_substroke( //
  const Stroke& polyline, const Float start_arclen, const Float end_arclen,
  const Float query_arclen, Vec2& out_proj, Float& out_s) {

  assert(start_arclen <= polyline.length() && "expected absolute arc length value");
  assert(end_arclen <= polyline.length() && "expected absolute arc length value");
  assert(query_arclen <= polyline.length() && "expected absolute arc length value");

  const auto xy = polyline.pos(query_arclen);
  const auto pen_width = polyline.pen_width();
  const auto start_frac_pair = polyline.fractional_index(start_arclen);
  const auto start_frac_idx = start_frac_pair.first + start_frac_pair.second;
  const auto end_frac_pair = polyline.fractional_index(end_arclen);
  const auto end_frac_idx = end_frac_pair.first + end_frac_pair.second;
  const auto end = (Index)std::ceil(end_frac_idx);

  Float best_sq_dist = infinity;
  for (auto i = (Index)start_frac_pair.first; i < end; ++i) {
    const auto p = polyline.xy(i);
    const auto q = polyline.xy(i + 1);
    const Float nsq = (q - p).squaredNorm();
    const auto lower = std::max(0.0, start_frac_idx - i);
    const auto upper = std::min(1.0, end_frac_idx - i);
    const auto t = std::clamp<Float>((xy - p).dot(q - p) / nsq, lower, upper);
    const Vec2 proj = lerp(p, q, t);
    const Float sq_dist = (proj - xy).squaredNorm();
    const auto s = clamped_lerp(polyline.arclength(i), polyline.arclength(i + 1), t);
    if (sq_dist < best_sq_dist &&
        // Check that it isn't spurious.
        // Keep this synced with others (grep for [spurious-condition]).
        std::max(3 * std::sqrt(sq_dist), M_PI * pen_width) < std::abs(query_arclen - s)) {
      best_sq_dist = sq_dist;
      out_proj = proj;
      out_s = s;
    }
  }
  return std::sqrt(best_sq_dist);
}

Float closest_points(const PolylineBVHLeaf& node1, const PolylineBVHLeaf& node2, //
                     Float& out_arclen1, Float& out_arclen2) {
  const auto& stroke1 = *node1.geometry;
  const auto& stroke2 = *node2.geometry;
  assert(stroke1.size() > 0);
  assert(stroke2.size() > 0);
  assert(&stroke1 != &stroke2);

  stroke1.ensure_arclengths();
  stroke2.ensure_arclengths();

  auto intersections = std::vector<Vec2>();
  intersect_different(node1, node2, intersections);
  if (!intersections.empty()) {
    out_arclen1 = intersections[0].x_;
    out_arclen2 = intersections[0].y_;
    return 0.0;
  }

  const auto n1 = stroke1.size();
  const auto n2 = stroke2.size();
  auto best_dist = Float(INFINITY);
  out_arclen1 = NAN;
  out_arclen2 = NAN;
  auto proj = Vec2(NAN, NAN);
  auto s = Float(NAN);
  // @optimize
  for (Index i = 0; i < n1; ++i) {
    const auto dist = closest_point(stroke2, stroke1.xy(i), proj, s);
    if (dist < best_dist) {
      best_dist = dist;
      out_arclen1 = stroke1.arclength(i);
      out_arclen2 = s;
    }
  }
  for (Index i = 0; i < n2; ++i) {
    const auto dist = closest_point(stroke1, stroke2.xy(i), proj, s);
    if (dist < best_dist) {
      best_dist = dist;
      out_arclen1 = s;
      out_arclen2 = stroke2.arclength(i);
    }
  }
  force_assert(std::isfinite(out_arclen1));
  force_assert(std::isfinite(out_arclen2));
  return (stroke1.pos(out_arclen1) - stroke2.pos(out_arclen2)).norm();
}

Float closest_points_on_self(const Stroke& stroke, Float& out_arclen1,
                             Float& out_arclen2) {
  stroke.ensure_arclengths();

  auto intersections = std::vector<Vec2>();
  intersect_self(stroke, intersections);
  if (!intersections.empty()) {
    out_arclen1 = intersections[0].x_;
    out_arclen2 = intersections[0].y_;
    return 0.0;
  }

  const auto n = stroke.size();
  auto best_dist = Float(INFINITY);
  out_arclen1 = NAN;
  out_arclen2 = NAN;
  auto proj = Vec2(NAN, NAN);
  auto s = Float(NAN);
  for (Index i = 0; i < n; ++i) {
    const auto dist =
      closest_point_to_own_stroke_interior_strict(stroke, stroke.arclength(i), proj, s);
    if (dist < best_dist) {
      best_dist = dist;
      out_arclen1 = stroke.arclength(i);
      out_arclen2 = s;
    }
  }
  if (std::isfinite(out_arclen1)) {
    assert(std::isfinite(out_arclen2));
    return (stroke.pos(out_arclen1) - stroke.pos(out_arclen2)).norm();
  } else {
    return INFINITY;
  }
}

std::vector<std::pair<std::size_t, Float>> pick(const PolylineBVH& polylines,
                                                const Vec2& xy, const Float tolerance) {
  std::multimap<Float, std::pair<std::size_t, Float>> hits; // sort by distance

  const auto nedges = polylines.nodes.size();
  for (auto i = decltype(nedges){0}; i < nedges; ++i) {
    if (polylines.nodes[i].bb.distanceTo(xy) <= tolerance) {
      const auto& geom = polylines.nodes[i];
      auto proj = Vec2::Empty();
      Float s;
      auto dist = closest_point(*geom.geometry, xy, proj, s);
      if (dist < tolerance) {
        hits.insert({dist, {i, s / geom.geometry->length()}});
      }
    }
  }

  std::vector<std::pair<std::size_t, Float>> out;
  out.reserve(hits.size());
  for (const auto& hit : hits) {
    out.push_back(hit.second);
  }

  return out;
}

std::vector<std::pair<std::size_t, Float>>
pick(const PolylineBVH& polylines, const Endpoint endpoint, const Float tolerance,
     const std::size_t max_results) {
  std::multimap<Float, std::pair<std::size_t, Float>> hits; // sort by distance

  const auto endpoint_idx = endpoint.stroke_idx();
  const auto& endpoint_stroke = *polylines.nodes[endpoint_idx].geometry;
  const auto xy = endpoint_stroke.xy(endpoint.is_head() ? 0 : endpoint_stroke.size() - 1);
  const auto nedges = polylines.nodes.size();
  for (auto i = decltype(nedges){0}; i < nedges; ++i) {
    if (i != endpoint_idx && polylines.nodes[i].bb.distanceTo(xy) <= tolerance) {
      const auto& geom = polylines.nodes[i];
      auto proj = Vec2::Empty();
      Float s;
      auto dist = closest_point(*geom.geometry, xy, proj, s);
      if (dist < tolerance) {
        hits.insert({dist, {i, s / geom.geometry->length()}});
      }
    }
  }

  // Special handling for own stroke.
  {
    auto proj = Vec2::Empty();
    Float s;
    const auto dist =
      closest_point_to_own_stroke(endpoint_stroke, endpoint.is_head(), proj, s);
    if (dist < tolerance) {
      hits.insert({dist, {endpoint_idx, s / endpoint_stroke.length()}});
    }
  }

  std::vector<std::pair<std::size_t, Float>> out;
  out.reserve(std::min(max_results, hits.size()));
  for (const auto& hit : hits) {
    if (out.size() >= max_results)
      break;
    out.push_back(hit.second);
  }

  return out;
}

std::vector<std::pair<size_t, Float>>
pick(const StrokeGraph& graph, const Endpoint endp, const size_t max_results,
     const std::unordered_set<size_t>* ignored_strokes) {
  auto hits = std::multimap<Float, std::pair<size_t, Float>>(); // sort by distance

  const auto nstrokes = graph.orig_strokes_.size();
  const auto& stroke = graph.orig_strokes_[endp.stroke_idx()];
  const auto xy = (endp.is_head() ? stroke.xy(0) : stroke.xy(Back));
  for (auto i = decltype(nstrokes){0}; i < nstrokes; ++i) {
    if (ignored_strokes && ignored_strokes->find(i) != ignored_strokes->end())
      continue;
    const auto& other_stroke = graph.orig_strokes_[i];
    if (other_stroke.size() == 0) {
      continue;
    }
    other_stroke.ensure_arclengths();
    auto proj = Vec2::Empty();
    Float s;
    auto dist = (&stroke == &other_stroke
                   ? closest_point_to_own_stroke(stroke, endp.is_head(), proj, s)
                   : closest_point(other_stroke, xy, proj, s));
    if (dist < infinity) {
      assert(s <= other_stroke.length());
      hits.insert({dist, {i, s / other_stroke.length()}});
    }
  }

  std::vector<std::pair<size_t, Float>> out;
  out.reserve(std::min(hits.size(), max_results));
  for (const auto& hit : hits) {
    if (out.size() >= max_results)
      break;
    out.push_back(hit.second);
  }
  return out;
}

std::vector<std::pair<std::size_t, Float>> pick_endpoints(const PolylineBVH& polylines,
                                                          const Vec2 xy) {
  std::multimap<Float, std::pair<std::size_t, Float>> hits; // sort by distance

  const auto nedges = polylines.nodes.size();
  for (auto i = decltype(nedges){0}; i < nedges; ++i) {
    const auto& stroke = *polylines.nodes[i].geometry;
    for (auto qi : {Index(0), stroke.size() - 1}) {
      if (stroke.size() < 1)
        continue;
      const auto q = stroke.xy(qi);
      const auto dist = Float{(xy - q).norm()};
      hits.insert({dist, {i, (qi == 0 ? 0.0 : 1.0)}});
    }
  }

  std::vector<std::pair<std::size_t, Float>> out;
  out.reserve(hits.size());
  for (const auto& hit : hits) {
    out.push_back(hit.second);
  }

  return out;
}

std::vector<Endpoint> pick_endpoints(const StrokeGraph& graph, const Endpoint endp,
                                     const size_t max_results,
                                     const std::unordered_set<size_t>* ignored_strokes) {
  auto hits = std::multimap<Float, Endpoint::IdType>(); // sort by distance

  const auto& stroke = graph.orig_strokes_[endp.stroke_idx()];
  const auto p = (endp.is_head() ? stroke.xy(0) : stroke.xy(Back));

  for (size_t si = 0; si < graph.orig_strokes_.size(); ++si) {
    const auto& other_stroke = graph.orig_strokes_[si];
    if (skip_stroke(other_stroke) ||
        (ignored_strokes && ignored_strokes->find(si) != ignored_strokes->end())) {
      continue;
    }
    for (auto is_head : {true, false}) {
      const auto other_endp = Endpoint(si, is_head);
      if (endp != other_endp) {
        const auto it = graph.endpoint2vertex_.find(other_endp);
        // No need to check dangling since any endpoint of an original stroke is allowed
        // here
        if (it != graph.endpoint2vertex_.end()) {
          const auto other_vertex = graph.vertex(it->second);
          if (other_vertex.is_active() && other_vertex.is_dangling()) {
            const auto q = other_vertex.pos();
            const auto dist = (q - p).norm();
            if (endp.stroke_idx() == si && 3 * dist >= stroke.length()) {
              continue;
            }
            hits.insert({dist, other_endp.as_int()});
          }
        }
      }
    }
  }

  auto out = std::vector<Endpoint>();
  out.reserve(std::min(hits.size(), max_results));
  for (const auto& hit : hits) {
    if (out.size() >= max_results)
      break;
    out.emplace_back(hit.second);
  }
  return out;
}

std::vector<std::size_t> pick_junctions(const std::vector<Junction>& junctions,
                                        const PolylineBVH& strokes, const Vec2 xy,
                                        const Float tolerance, const size_t max_results) {
  const auto sq_tolerance = tolerance * tolerance;
  const auto njunctions = junctions.size();
  std::multimap<Float, std::size_t> hits; // sort by distance
  for (auto i = decltype(njunctions){0}; i < njunctions; ++i) {
    const auto c = centroid(junctions[i], strokes);
    const Float sq_dist = (c - xy).squaredNorm();
    if (sq_dist < sq_tolerance) {
      hits.insert({sq_dist, i});
    }
  }
  std::vector<std::size_t> out;
  out.reserve(std::min(max_results, hits.size()));
  for (const auto& hit : hits) {
    if (out.size() >= max_results)
      break;
    out.push_back(hit.second);
  }
  return out;
}

} // namespace sketching
