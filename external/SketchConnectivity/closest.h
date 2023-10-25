#pragma once

#include "bvh.h"

#include <unordered_set>

namespace sketching {

struct Endpoint;
struct Junction;
struct StrokeGraph;

/**
 * Polyline must not have any segments of length 0 (i.e. repeated vertices).
 */
Float closest_point(const Stroke& polyline, Vec2 xy, Vec2& out_proj, Float& out_s);

/**
 * Return the distance to the closest point to `xy` on `polyline`.
 *
 * If the return value is infinity, then there were no points within `max_dist` distance
 * (although it is possible for the return value to be greater than `max_dist`).
 *
 * Polyline must not have any segments of length 0 (i.e. repeated vertices).
 */
inline Float closest_point(const PolylineBVHLeaf& polyline, Vec2 xy, Float max_dist,
                           Vec2& out_proj, Float& out_s) {
  if (polyline.bb.distanceTo(xy) > max_dist)
    return std::numeric_limits<Float>::infinity();
  return closest_point(*polyline.geometry, xy, out_proj, out_s);
}

/**
 * Return the closest point on the substroke starting from `start_arclen` and ending at
 * `end_arclen`
 *
 * Polyline must not have any segments of length 0 (i.e. repeated vertices).
 */
Float closest_point_substroke( //
  const Stroke& polyline, Float start_arclen, Float end_arclen, Vec2 xy, //
  Vec2& out_proj, Float& out_s);

Float closest_point_to_own_stroke(const Stroke& polyline, bool head, Vec2& out_proj,
                                  Float& out_s);
Float closest_point_to_own_stroke(const Stroke& polyline, Float arclen, Vec2& out_proj,
                                  Float& out_s);

/**
 * Substroke version of `closest_point_to_own_stroke`.
 *
 * See also `closest_point_substroke`.
 */
Float closest_point_to_own_substroke( //
  const Stroke& polyline, Float start_arclen, Float end_arclen, Float query_arclen, //
  Vec2& out_proj, Float& out_s);

Float closest_points(const PolylineBVHLeaf& stroke1, const PolylineBVHLeaf& stroke2, //
                     Float& out_arclen1, Float& out_arclen2);

Float closest_points_on_self(const Stroke& stroke, Float& out_arclen1,
                             Float& out_arclen2);

/**
 * Return the closest points on all strokes in `polylines` within `tolerance` distance of
 * (x, y).  The results are sorted by distance.
 */
std::vector<std::pair<std::size_t, Float>> pick(const PolylineBVH& polylines,
                                                const Vec2& xy, Float tolerance);

/**
 * Return the closest points on all strokes in `polylines` within `tolerance` distance of
 * (x, y).  The results are sorted by distance.
 */
inline std::vector<std::pair<std::size_t, Float>>
pick(const PolylineBVH& polylines, Float x, Float y, Float tolerance) {
  return pick(polylines, Vec2(x, y), tolerance);
}

std::vector<std::pair<std::size_t, Float>> pick(const PolylineBVH& polylines, Endpoint,
                                                Float tolerance, std::size_t max_results);

std::vector<std::pair<std::size_t, Float>>
pick(const StrokeGraph& graph, Endpoint, size_t max_results,
     const std::unordered_set<size_t>* ignored_strokes = nullptr);

/**
 * Return all endpoints in `polylines` within `tolerance` distance of (x, y).  The results
 * are sorted by distance.
 */
std::vector<std::pair<std::size_t, Float>> pick_endpoints(const PolylineBVH& polylines,
                                                          Vec2 xy);

/**
 * Return all endpoints in `polylines` within `tolerance` distance of (x, y).  The results
 * are sorted by distance.
 */
inline std::vector<std::pair<std::size_t, Float>>
pick_endpoints(const PolylineBVH& polylines, Float x, Float y) {
  return pick_endpoints(polylines, Vec2(x, y));
}

/**
 * Returns endpoints with indices corresponding to strokes in `graph.orig_strokes_`.
 */
std::vector<Endpoint>
pick_endpoints(const StrokeGraph& graph, Endpoint endp,
               size_t max_results = std::numeric_limits<size_t>::max(),
               const std::unordered_set<size_t>* ignored_strokes = nullptr);

std::vector<std::size_t>
pick_junctions(const std::vector<Junction>& junctions, const PolylineBVH& strokes,
               Vec2 xy, Float tolerance,
               size_t max_results = std::numeric_limits<size_t>::max());

inline std::vector<std::size_t>
pick_junctions(const std::vector<Junction>& junctions, const PolylineBVH& strokes,
               Float x, Float y, Float tolerance,
               size_t max_results = std::numeric_limits<size_t>::max()) {
  return pick_junctions(junctions, strokes, Vec2(x, y), tolerance, max_results);
}

} // namespace sketching
