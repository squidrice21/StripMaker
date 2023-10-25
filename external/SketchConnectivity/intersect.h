#pragma once

#include "types.h"

#include "3rdparty/Triangle/Triangle/triangle.h"

#include <vector>

namespace sketching {

struct Stroke;
struct PolylineBVH;
struct PolylineBVHLeaf;
struct BoundingBox;
struct EnvelopeBVHLeaf;
struct ConstStrokeView;

/**
 * Positive iff triangle pqr winds counter-clockwise in a Y-up system.
 */
inline Float signed_area(Vec2 p, Vec2 q, Vec2 r) {
  return 0.5 * (p.x() * q.y() + q.x() * r.y() + r.x() * p.y() //
                - r.x() * q.y() - p.x() * r.y() - q.x() * p.y());
}

/**
 * Positive iff triangle pqr winds counter-clockwise in a Y-up system.  Falls
 * back to robust arithmetic if necessary.
 */
inline Float counter_clockwise(Vec2 p, Vec2 q, Vec2 r) {
  constexpr auto epsilon = Float(0.5) * std::numeric_limits<Float>::epsilon();
  constexpr auto ccwerrboundA = Float((3.0 + 16.0 * epsilon) * epsilon);

  const auto detleft = (p.x_ - r.x_) * (q.y_ - r.y_);
  const auto detright = (p.y_ - r.y_) * (q.x_ - r.x_);
  const auto det = detleft - detright;
  Float detsum;
  if (detleft > 0.0) {
    if (detright <= 0.0) {
      return det;
    } else {
      detsum = detleft + detright;
    }
  } else if (detleft < 0.0) {
    if (detright >= 0.0) {
      return det;
    } else {
      detsum = -detleft - detright;
    }
  } else {
    return det;
  }
  const auto errbound = ccwerrboundA * detsum;
  if ((det >= errbound) || (-det >= errbound)) {
    return det;
  }
  return counterclockwiseadapt(&p.x_, &q.x_, &r.x_, detsum);
}

// Based off of VPaint's VectorAnimationComplex/Intersection.cpp.
//
// returns if there is an intersection or not
// in addition, return in s and t the position of intersection:
//   I = A + s*(B-A)/||B-A|| = C + t*(D-C)/||D-C||
//
//   In a math world, 0 <= s <= ||B-A|| and 0 <= t <= ||D-C||
//
//   Here, sometimes returns true with -epsilon < s < 0
//   or ||B-A|| < s < ||B-A|| + epsilon
//
//   This is  to be sure we  don't miss any  intersection due to
//   rounding   error,  but   caution  then,   maybe   the  same
//   intersection could be found twice.
[[nodiscard]] inline bool segment_intersection( //
  const Vec2 &a, const Vec2 &b, const Vec2 &c, const Vec2 &d, double &s,
  double &t) {

  // avoid useless computation if too far
  Vec2 u = (b - a);
  const Float normU = u.norm();
  Vec2 v = (d - c);
  const Float normV = v.norm();
  if ((c - a).norm() < normU + normV) {
    if (normU > 0)
      u *= (1 / normU);
    if (normV > 0)
      v *= (1 / normV);

    const Float detA = v.x() * u.y() - v.y() * u.x();
    const Float epsilon = 1e-10;
    // reason of epsilon:
    //   don't intersect if nearly parallel
    if (detA < -epsilon || epsilon < detA) {
      Mat2 invA(-v.y(), v.x(), -u.y(), u.x());
      invA *= (1 / detA);
      Vec2 st = invA * (c - a);

      // reason of eps:
      //   can't miss an intersection due
      //   to rounding errors, but eventually
      //   find the same twice.
      const Float eps = 1e-10;
      if (-eps <= st.x() && st.x() < normU + eps && -eps <= st.y() &&
          st.y() < normV + eps) {
        // here, we know they intersect
        s = st.x();
        t = st.y();
        return true;
      }
    }
  }

  return false;
}

/**
 * Excludes p and q from the intersections.  Useful for line of sight
 * calculations.
 *
 * @param p Start of line segment.
 * @param q End of line segment.
 * @param stroke Stroke to intersect with.
 * @param out_s Value in [0, 1] indicating where on segment pq intersection was
 * found.
 * @param out_t Arc length value indicating where on stroke intersection was
 * found.
 * @return true if and only if an intersection was found.
 */
[[nodiscard]] bool
intersect_segment_stroke_exclusive(Vec2 p, Vec2 q,
                                   const PolylineBVHLeaf &stroke, Float &out_s,
                                   Float &out_t);

void intersect_segment_stroke_exclusive(Vec2 p, Vec2 q,
                                        const PolylineBVHLeaf &stroke,
                                        std::vector<Vec2> &out_arclens);

bool line_of_sight(Vec2 p, Vec2 q, span<const Stroke> strokes,
                   span<const BoundingBox> centerline_bbs);

bool line_of_sight(Vec2 p, Vec2 q, const PolylineBVH &bvh, //
                   Float min_s = 1e-5, Float max_s = 1.0 - 1e-5);

/// Returns list of (s, t) arc lengths.
void intersect_self(const Stroke &node, std::vector<Vec2> &out_arclens);

void intersect_self(const PolylineBVHLeaf &node,
                    std::vector<Vec2> &out_arclens);

/// Returns list of (s, t) arc lengths.
void intersect_different(const PolylineBVHLeaf &node1,
                         const PolylineBVHLeaf &node2,
                         std::vector<Vec2> &out_arclens);

struct IntersectionEvent {
  double s_start = -1;
  double t_start = -1;
  double s_mid = -1;
  double t_mid = -1;
  double s_end = -1;
  double t_end = -1;

  explicit operator bool() const { return s_start >= 0.0; }

  /// Only centerline intersections have a well-defined midpoint.
  bool centerline_intersection() const { return s_mid >= 0.0; }

  void normalize() {
    if (s_start > s_end) {
      std::swap(s_start, s_end);
    }
    if (t_start > t_end) {
      std::swap(t_start, t_end);
    }
  }

  std::string repr() const;
};

void intersect_self(const EnvelopeBVHLeaf &node,
                    std::vector<IntersectionEvent> &out_intersections);

void intersect_different(const EnvelopeBVHLeaf &node1,
                         const EnvelopeBVHLeaf &node2,
                         std::vector<IntersectionEvent> &out_intersections);

/**
 * An uneven capsule with different radii at each end.  Tapered strokes are made
 * of these.
 */
struct Bicapsule {
  Vec2 p_a;
  Vec2 p_b;
  Float r_a;
  Float r_b;
};

Bicapsule stroke_segment(const Stroke &s, const Index i);
/**
 * Return (signed distance, projection amount) to a bicapsule.  Projection
 * amount is not clamped to [0, 1].
 */
Vec2 signed_dist_bicapsule(Vec2 p, Bicapsule bicap);

Float signed_dist_stroke(Vec2 p, const Stroke &stroke);

/**
 * Returns (distance from p, arc length on stroke).  Distance will be infinity
 * if no good snap point is found.  Distance is the centerline distance.
 *
 * @param to_snap Stroke containing endpoint to snap.
 * @param head True iff you are trying to snap the head, otherwise it is assumed
 * that you are trying to snap the tail.
 * @param other BVH node containing the stroke to consider snapping to.
 * @param out_closest_env_dist Will be set to the smallest envelope distance
 * from the endpoint to the stroke.
 */
Vec2 find_snap_point(const ConstStrokeView &to_snap, bool head,
                     const EnvelopeBVHLeaf &other, Float *out_closest_env_dist);

/// Deprecated.
Vec2 find_snap_point_legacy(const Stroke &to_snap, bool head,
                            const EnvelopeBVHLeaf &other);

/**
 * Cast a ray in direction `r` from point `p`, and return all points which
 * intersect the ray.  The results are sorted by distance from `r`.
 */
std::vector<std::pair<std::size_t, Float>>
ray_intersections(const PolylineBVH &polylines, const Vec2 &p, const Vec2 &r);

} // namespace sketching
