#pragma once

#include "stroke_graph.h"

namespace sketching {

struct Stroke;

inline Float smart_projection_dist(const Stroke& s1, const Float arclen1,
                                   const Stroke& s2, Vec2& out_proj, Float& out_s) {
  if (&s1 == &s2) {
    const auto dist = closest_point_to_own_stroke(s1, arclen1, out_proj, out_s);
    return dist;
  } else {
    return closest_point(s2, s1.pos(arclen1), out_proj, out_s);
  }
}

inline Float smart_projection_dist(const Stroke& s1, const Float arclen1,
                                   const Stroke& s2, Float& out_s) {
  auto proj = Vec2::Empty();
  return smart_projection_dist(s1, arclen1, s2, proj, out_s);
}

inline Float smart_projection_dist(const Stroke& s1, const Float arclen1,
                                   const Stroke& s2) {
  auto proj = Vec2::Empty();
  Float s;
  return smart_projection_dist(s1, arclen1, s2, proj, s);
}

// TODO: Create an enum for the return type so we can refer to types by name.
std::uint8_t end_end_junction_type(const Stroke& stroke1, Float arclen1,
                                   const Stroke& stroke2, Float arclen2);

/// This type is different from the broad classification into end-end vs end-stroke, but
/// instead describes the configuration which is determined by the endpoint's projection
/// to an endpoint or the interior.
std::uint8_t end_stroke_junction_type(const Stroke& stroke, Float arclen);
std::uint8_t end_stroke_junction_type(const Stroke& stroke, Float end_arclen,
                                      Float arclen, Float ratio);

std::uint8_t end_stroke_junction_type(StrokeGraph::HedgeView edge, Float arclen);

[[nodiscard]] bool intersecting_and_diverging_t( //
  const Stroke& endpoint_stroke, Float arclen1, const Stroke& interior, Float arclen2);

} // namespace sketching
