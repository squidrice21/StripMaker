#include "junction_type.h"

#include "force_assert.h"
#include "intersect.h"
#include "stroke_graph_extra.h"

namespace sketching {

/**
 * Assess whether a connection crosses and diverges, e.g. something like
 *
 *               │
 *     ──────────┼────
 *               │
 *               │
 *               │
 */
static bool is_crossing_end_end_diverging(const Stroke& stroke1, Float arclen1,
                                          const Stroke& stroke2, Float arclen2) {
  force_assert(arclen1 == 0.0 || arclen1 == stroke1.length());
  force_assert(arclen2 == 0.0 || arclen2 == stroke2.length());

  if (&stroke1 == &stroke2 && arclen1 > arclen2) {
    // Get consistent behaviour for the returned intersections.
    std::swap(arclen1, arclen2);
  }

  // Find intersection point.
  stroke1.ensure_arclengths();
  stroke2.ensure_arclengths();
  auto arclens = std::vector<Vec2>();
  if (&stroke1 == &stroke2) {
    intersect_self(PolylineBVHLeaf(stroke1, bounds(stroke1)), arclens);
  } else {
    intersect_different(PolylineBVHLeaf(stroke1, bounds(stroke1)),
                        PolylineBVHLeaf(stroke2, bounds(stroke2)), arclens);
  }
  if (arclens.empty()) {
    return false;
  }
  if (&stroke1 == &stroke2) {
    for (const auto& intersection : arclens) {
      // We need this for the below arc length distances to be correct.
      force_assert(intersection.x_ <= intersection.y_);
    }
  }

  // Determine if they share the same closest intersection point (because they could
  // intersect each other multiple times).
  auto closest_intersection_idx1 = size_t(0);
  for (size_t i = 1; i < arclens.size(); ++i) {
    if (arclen1 == 0.0) {
      if (arclens[i].x_ < arclens[closest_intersection_idx1].x_) {
        closest_intersection_idx1 = i;
      }
    } else {
      if (arclens[i].x_ > arclens[closest_intersection_idx1].x_) {
        closest_intersection_idx1 = i;
      }
    }
  }
  auto closest_intersection_idx2 = size_t(0);
  for (size_t i = 1; i < arclens.size(); ++i) {
    if (arclen2 == 0.0) {
      if (arclens[i].y_ < arclens[closest_intersection_idx2].y_) {
        closest_intersection_idx2 = i;
      }
    } else {
      if (arclens[i].y_ > arclens[closest_intersection_idx2].y_) {
        closest_intersection_idx2 = i;
      }
    }
  }
  if (closest_intersection_idx1 != closest_intersection_idx2) {
    return false;
  }

  const auto& intersection = arclens[closest_intersection_idx1];
  const auto arc_dist1 = std::abs(arclen1 - intersection.x_);
  const auto arc_dist2 = std::abs(arclen2 - intersection.y_);
  const auto arc_dist = arc_dist1 + arc_dist2;
  const auto cen_dist = (stroke1.pos(arclen1) - stroke2.pos(arclen2)).norm();
  const auto pen_width = std::max(stroke1.pen_width(), stroke2.pen_width());

  // Keep this synced with others (grep for [spurious-condition]).
  return !(std::max(3 * cen_dist, M_PI * pen_width) < arc_dist);
}

bool intersecting_and_diverging_t(const Stroke& endpoint, const Float arclen1,
                                  const Stroke& interior, const Float arclen2) {
  force_assert(arclen1 == 0.0 || arclen1 == endpoint.length());
  force_assert(0.0 <= arclen2 || arclen2 <= interior.length());

  // Find intersection point.
  interior.ensure_arclengths();
  endpoint.ensure_arclengths();
  auto arclens = std::vector<Vec2>();
  if (&interior == &endpoint) {
    intersect_self(PolylineBVHLeaf(interior, bounds(interior)), arclens);
  } else {
    intersect_different(PolylineBVHLeaf(endpoint, bounds(endpoint)),
                        PolylineBVHLeaf(interior, bounds(interior)), arclens);
  }
  if (arclens.empty()) {
    return false;
  }
  if (&interior == &endpoint) {
    for (const auto& intersection : arclens) {
      // We need this for the below arc length distances to be correct.
      force_assert(intersection.x_ <= intersection.y_);
    }
  }

  // Find the closest intersection to the endpoint.
  auto closest_intersection_idx1 = size_t(0);
  auto arc_dist1 = std::abs(arclens[closest_intersection_idx1].x_ - arclen1);
  for (size_t i = 1; i < arclens.size(); ++i) {
    const auto dist = std::abs(arclens[i].x_ - arclen1);
    if (dist < arc_dist1) {
      closest_intersection_idx1 = i;
      arc_dist1 = dist;
    }
    if (&interior == &endpoint) {
      // Both arclens in the intersection are from this stroke so we must check both.
      const auto dist2 = std::abs(arclens[i].y_ - arclen1);
      if (dist2 < arc_dist1) {
        closest_intersection_idx1 = i;
        arc_dist1 = dist;
      }
    }
  }

  // See if we have another intersection between the intersection and arclen2.
  Float intersection_arclen2 = 0.0;
  Float arc_dist2 = 0.0;
  if (&interior == &endpoint &&
      std::abs(arclen2 - arclens[closest_intersection_idx1].x_) <
        std::abs(arclen2 - arclens[closest_intersection_idx1].y_)) {
    intersection_arclen2 = arclens[closest_intersection_idx1].x_;
    arc_dist2 = std::abs(arclen2 - arclens[closest_intersection_idx1].x_);
  } else {
    intersection_arclen2 = arclens[closest_intersection_idx1].y_;
    arc_dist2 = std::abs(arclen2 - arclens[closest_intersection_idx1].y_);
  }
  for (const auto& intersection : arclens) {
    if ((intersection_arclen2 < intersection.y_ && intersection.y_ < arclen2) ||
        (intersection_arclen2 > intersection.y_ && intersection.y_ > arclen2)) {
      return false;
    }
    if (&interior == &endpoint) {
      if ((intersection_arclen2 < intersection.x_ && intersection.x_ < arclen2) ||
          (intersection_arclen2 > intersection.x_ && intersection.x_ > arclen2)) {
        return false;
      }
    }
  }

  const auto arc_dist = arc_dist1 + arc_dist2;
  const auto cen_dist = (endpoint.pos(arclen1) - interior.pos(arclen2)).norm();
  const auto pen_width = std::max(endpoint.pen_width(), interior.pen_width());

  // Keep this synced with others (grep for [spurious-condition]).
  return !(std::max(3 * cen_dist, M_PI * pen_width) < arc_dist);
}

std::uint8_t end_end_junction_type(const Stroke& s1, const Float arclen1,
                                   const Stroke& s2, const Float arclen2) {
  // Arc lengths should be absolute.
  force_assert(arclen1 == 0.0 || arclen1 == s1.length());
  force_assert(arclen2 == 0.0 || arclen2 == s2.length());

  if (is_crossing_end_end_diverging(s1, arclen1, s2, arclen2)) {
    // Override this case to be type 2 (diverging) regardless of projections.
    return 2;
  }

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  Vec2 proj_onto_s2 = Vec2::Empty();
  Float s1_onto_s2;
  smart_projection_dist(s1, arclen1, s2, proj_onto_s2, s1_onto_s2);
  Vec2 proj_onto_s1 = Vec2::Empty();
  Float s2_onto_s1;
  smart_projection_dist(s2, arclen2, s1, proj_onto_s1, s2_onto_s1);

  const auto margin1 = 0.5 * s1.width_at(arclen1);
  const auto margin2 = 0.5 * s2.width_at(arclen2);
  const auto s1_projects_onto_end =
    (s1_onto_s2 <= margin2 || s1_onto_s2 >= s2.length() - margin2);
  // Should be std::abs(arclen2 - s1_onto_s2) < margin2;
  const auto s2_projects_onto_end =
    (s2_onto_s1 <= margin1 || s2_onto_s1 >= s1.length() - margin1);
  // Should be std::abs(arclen1 - s2_onto_s1) < margin1;

  if (s1_projects_onto_end && s2_projects_onto_end) {
    return 0;
  } else if (s1_projects_onto_end != s2_projects_onto_end) {
    return 1;
  } else {
    double s = -INFINITY, t = -INFINITY;
    const auto intersecting =
      segment_intersection(p, proj_onto_s2, q, proj_onto_s1, s, t);
    const auto norm1 = (p - proj_onto_s2).norm();
    const auto norm2 = (q - proj_onto_s1).norm();
    constexpr auto eps = 1e-5;
    if (intersecting && s > eps && s < norm1 - eps && t > eps && t < norm2 - eps) {
      // Projections intersect.
      assert(std::isfinite(s));
      assert(std::isfinite(t));
      return 2;
    } else {
      return 3;
    }
  }
}

std::uint8_t end_stroke_junction_type(const Stroke& stroke, const Float arclen) {
  assert(arclen <= stroke.length() && "arclen should be absolute");

  const auto margin = 0.5 * stroke.width_at(arclen);
  const auto s1_projects_onto_end =
    (arclen <= margin || arclen >= stroke.length() - margin);
  if (s1_projects_onto_end) {
    return 0; // Not a T-junction.
  } else {
    return 1;
  }
}

std::uint8_t end_stroke_junction_type(const Stroke& stroke, Float end_arclen,
                                      Float arclen, Float ratio) {
  assert(arclen <= stroke.length() && "arclen should be absolute");

  // const auto margin = ratio * stroke.width_at(arclen);
  const auto margin = ratio * stroke.pen_width();
  const auto s1_projects_onto_end = std::abs(end_arclen - arclen) < margin;
  if (s1_projects_onto_end) {
    return 0; // Not a T-junction.
  } else {
    return 1;
  }
}

std::uint8_t end_stroke_junction_type(StrokeGraph::HedgeView edge, const Float arclen) {
  if (!edge.forward()) {
    edge = edge.twin();
  }
  const auto& stroke = edge.stroke();
  assert(arclen <= stroke.length() && "arclen should be absolute");
  const auto margin = 0.5 * stroke.width_at(arclen);
  if (arclen < 0.5 * stroke.length()) {
    const auto tol = (edge.origin().is_dangling() ? margin : 0.0);
    if (arclen <= tol) {
      return 0;
    }
  } else {
    const auto tol = (edge.dest().is_dangling() ? margin : 0.0);
    if (arclen >= stroke.length() - tol) {
      return 0;
    }
  }
  return 1;
}

} // namespace sketching
