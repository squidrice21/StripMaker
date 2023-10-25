#include "incremental.h"

#include "classifier.h"
#include "deform.h"
#include "degrees.h"
#include "detail/alloca.h"
#include "detail/suppress_warning.h"
#include "detail/util.h"
#include "fairing.h"
#include "features/junction_features_impl.h"
#include "fitting.h"
#include "force_assert.h"
#include "incremental_param.h"
#include "incremental_region_util.h"
#include "intersect.h"
#include "junction_type.h"
#include "stroke_graph_extra.h"

namespace sketching {

FeatureType prediction_feature_type = FeatureType::OrigStroke;

namespace {

using HedgeRecord = StrokeGraph::HedgeRecord;
using VertexRecord = StrokeGraph::VertexRecord;
using FaceView = StrokeGraph::FaceView;
using HedgeView = StrokeGraph::HedgeView;
using VertexView = StrokeGraph::VertexView;

constexpr auto invalid = StrokeGraph::invalid;
constexpr auto n_candidates = 3;

Float acos(Float x) { return std::acos(std::min(x, 1.0)); }

void delete_stroke_and_update_bvh(StrokeGraph &graph, const size_t si) {
  graph.bvh_.clear(si);
  graph.strokes2orig_[si].clear();
  graph.hedges_[2 * si].deactivate();
  graph.hedges_[2 * si + 1].deactivate();
}

std::pair<size_t, size_t> make_twins(std::vector<HedgeRecord> &hedges) {
  hedges.emplace_back(); // Do not use this reference, it will become invalid if
                         // the vector resizes on the next line.
  hedges.emplace_back();
  const auto new_edge_idx = hedges.size() - 2;
  const auto new_twin_idx = new_edge_idx + 1;
  return {new_edge_idx, new_twin_idx};
}

void assign_prev_next(StrokeGraph &graph, const size_t prev,
                      const size_t next) {
  graph.hedges_[prev].next_ = next;
  graph.hedges_[next].prev_ = prev;
}

size_t snap_if_overlapping(StrokeGraph &graph, const VertexView vert,
                           IncrementalCache &cache) {
  // Find snap location if it exists.
  const auto vert_pos = vert.pos();
  auto best_dist = infinity;
  auto stroke_idx = invalid;
  auto best_arclen = -infinity;
  for (size_t si = 0; si < graph.strokes_.size(); ++si) {
    const auto &other_stroke = graph.strokes_[si];
    const auto other_he = graph.hedge(2 * si);
    if (other_stroke.size() > 1) {
      other_stroke.ensure_arclengths();
      const auto &stroke_bb = graph.bvh_.envelope_bbs()[si];
      // @optimize Sort bounding boxes by distance, then early-out.
      if (stroke_bb.distanceTo(vert_pos) < best_dist) {
        const auto &this_stroke = vert.hedge().stroke();
        auto env_dist = Float(infinity);
        auto [cen_dist, arclen] = find_snap_point(
          this_stroke, /*head=*/vert.hedge().forward(),
          EnvelopeBVHLeaf(other_he.stroke(), stroke_bb), &env_dist);
        force_assert(!std::isfinite(cen_dist) || env_dist <= 0);
        if (other_he.origin() == vert || other_he.dest() == vert) {
          // For snapping to self, it needs to be non-spurious.
          const auto otherp_arclen =
            (other_he.origin() == vert ? 0.0 : other_stroke.length());
          const auto arc_dist = std::abs(arclen - otherp_arclen);
          // Keep this synced with others (grep for [spurious-condition]).
          if (!(std::max(3 * cen_dist, M_PI * other_stroke.pen_width()) <
                arc_dist)) {
            continue; // Spurious.
          }
        }

        if (cen_dist < best_dist) {
          best_dist = cen_dist;
          stroke_idx = si;
          best_arclen = arclen;
        }
      }
    }
  }

  // Do the snap.
  if (std::isfinite(best_dist)) {
    auto new_snap_vertex_idx = invalid;
    if (would_form_sliver_vert_interior(
          graph, vert.index_,
          {(int)stroke_idx,
           best_arclen / graph.strokes_[stroke_idx].length()})) {
      if (vert.is_dangling()) {
        const auto he = vert.hedge();
        new_snap_vertex_idx = he.dest().index_;
        const auto si = he.stroke_idx();
        drop_vertex_and_neighbors(graph, vert.index_, cache.candidates_);
        delete_stroke_and_update_bvh(graph, si);
      }
    } else {
      // Check line of sight.  If we are blocked by another centerline, connect
      // to that instead.
      Float s = -1, t = -1;
      auto snap_pos = graph.strokes_[stroke_idx].pos(best_arclen);
      const auto bvh = graph.bvh_.polyline_bvh();
      for (size_t i = 0; i < graph.strokes_.size(); ++i) {
        if (intersect_segment_stroke_exclusive(vert_pos, snap_pos, bvh.nodes[i],
                                               s, t)) {
          stroke_idx = i;
          best_arclen = t;
          snap_pos = graph.strokes_[i].pos(best_arclen);
          best_dist = (snap_pos - vert_pos).norm();
        }
      }

      const auto snap_vertex = graph.snap_endpoint_to_edge(
        vert, stroke_idx, best_arclen, graph.snapping_method_type_);
      if (snap_vertex.is_valid() && snap_vertex != vert) {
        new_snap_vertex_idx = snap_vertex.index_;
      }
    }
    return new_snap_vertex_idx;
  }

  return invalid;
}

/**
 * @param endpoint A dangling endpoint.
 * @param out_candidates Output array of (stroke index, normalized arc length)
 * pairs.  The input size determines the maximum number of candidates to return.
 *                       If there are fewer candidates than this, then the size
 * of the span will be modified to reflect the actual number of returned
 *                       candidates.
 */
void find_end_end_candidates(const VertexView vertex,
                             span<StrokeTime> &out_candidates) {
  if (out_candidates.empty()) {
    return;
  }

  const auto &graph = *vertex.graph_;
  const auto centerline_bbs = graph.bvh_.centerline_bbs();
  auto hits = std::multimap<Float, size_t>(); // Sort by distance.

  const auto n_vertices = graph.vertices_.size();
  for (size_t i = 0; i < n_vertices; ++i) {
    const auto other_vertex = graph.vertex(i);
    if (other_vertex.is_active() &&
        is_snap_valid(graph, vertex.index_, other_vertex.index_)) {
      auto has_diverging_type = false;
      const auto he1 = vertex.hedge();
      const auto he2 = other_vertex.hedge();
      auto it1 = he1;
      auto iterations = 0;
      do {
        auto it2 = he2;
        do {
          // TODO: When we switch to computing features on original strokes,
          // this needs to
          //       change. Grep for [featonorig].
          const auto &s1 = it1.stroke();
          const auto &s2 = it2.stroke();
          const auto arclen1 = (it1.forward() ? 0.0 : s1.length());
          const auto arclen2 = (it2.forward() ? 0.0 : s2.length());
          if (end_end_junction_type(s1, arclen1, s2, arclen2) == 2) {
            has_diverging_type = true;
            goto jump_out; // Jump out of nested do-whiles.
          }

          it2 = it2.twin().next();
          iterations++;
          force_assert(iterations < 1024 && "likely infinite loop found");
        } while (it2 != he2);
        it1 = it1.twin().next();
      } while (it1 != he1);

jump_out:
      if (!has_diverging_type) {
        // This candidate passes the pre-filter.
        const auto dist = (vertex.pos() - other_vertex.pos()).norm();
        hits.insert({dist, other_vertex.index_});
      }
    }
  }

  size_t i = 0;
  for (const auto &hit : hits) {
    const auto other_vertex = graph.vertex(hit.second);
    if (line_of_sight(vertex.pos(), other_vertex.pos(), graph.strokes_,
                      centerline_bbs)) {
      const auto endp = vertex_to_endpoint(other_vertex);
      out_candidates[i++] =
        StrokeTime((int)endp.stroke_idx(), endp.is_head() ? 0.0 : 1.0);
      if (i == out_candidates.size()) {
        break;
      }
    }
  }
  out_candidates = out_candidates.first(i);
}

} // namespace

static bool is_parallel_slow(const Stroke &stroke1, const Float arclen1, //
                             const Stroke &stroke2, const Float arclen2) {
#ifdef HAS_GUROBI
  if (&stroke1 == &stroke2) {
    return false;
  }

  const auto ee_type =
    end_end_junction_type(stroke1, arclen1, stroke2, arclen2);
  if (ee_type != 0 && ee_type != 1) {
    return false;
  }

  // Copied from paired_parallel_vertices in adjacent_junction.cpp.
  const auto angle_aligned = [](const Float angle) {
    return std::abs(90_deg - angle) <= end_alignment_angle_offset;
  };
  const auto end_tangent_alignment = [](const Stroke &s1, Float arclen1,
                                        const Stroke &s2,
                                        Float arclen2) -> bool {
    const Vec2 p = s1.pos(arclen1);
    const Vec2 q = s2.pos(arclen2);
    if ((p - q).squaredNorm() < 1e-10) {
      return false;
    } else {
      auto v = q - p;
      const auto centerline_dist = v.norm();
      // Min step away size
      const auto centerline_dist1 =
        std::max(centerline_dist, 2 * s1.pen_width());
      const auto centerline_dist2 =
        std::max(centerline_dist, 2 * s2.pen_width());
      const auto tangent1 = features::stepaway_tangent(
        s1, arclen1 < 0.5 * s1.length(), centerline_dist1);
      const auto tangent2 = features::stepaway_tangent(
        s2, arclen2 < 0.5 * s2.length(), centerline_dist2);
      const auto angle =
        std::acos(std::clamp(tangent1.dot(tangent2), -1.0, 1.0));
      return angle <= end_alignment_angle_offset;
    }
  };

  const auto cen_dist = (stroke1.pos(arclen1) - stroke2.pos(arclen2)).norm();

  features::StepawayTangentAngle1 connection_end_angle;
  const auto angle1 =
    connection_end_angle(stroke1, arclen1, 0, stroke2, arclen2, 1);
  const auto angle2 =
    connection_end_angle(stroke2, arclen2, 1, stroke1, arclen1, 0);
  if (!angle_aligned(angle1) || !angle_aligned(angle2) ||
      !end_tangent_alignment(stroke1, arclen1, stroke2, arclen2)) {
    return false;
  }

  constexpr Float accuracy = 1;
  constexpr bool cut_strokes = true;
  constexpr Float sampling_px_step = 2;
  constexpr size_t min_samples = 5;
  Float overlapping_ratio = -1;
  const auto max_dist = strokewise_measure(
    stroke1, stroke2, pointwise_euclidean_difference, accuracy, cut_strokes,
    sampling_px_step, min_samples, &overlapping_ratio, /*find_max=*/true);
  // @optimize This finds the common parameterization twice!
  const auto avg_angle_diff = strokewise_measure(
    stroke1, stroke2, pointwise_angular_difference, accuracy, cut_strokes,
    sampling_px_step, min_samples, &overlapping_ratio);
  if (overlapping_ratio > 0 && max_dist / cen_dist <= dist_change_ratio &&
      avg_angle_diff <= stroke_alignment_angle) {
    return true;
  }
#endif
  return false;
}

static bool is_parallel(const StrokeGraph &graph, const Stroke &stroke1,
                        const Float arclen1, //
                        const Stroke &stroke2, const Float arclen2) {
  assert(arclen1 == 0.0 || arclen1 == stroke1.length());
  assert(arclen2 == 0.0 || arclen2 == stroke2.length());
  const auto end1 =
    Endpoint(&stroke1 - graph.orig_strokes_.data(), arclen1 == 0.0);
  force_assert(end1.stroke_idx() < graph.orig_strokes_.size());
  const auto end2 =
    Endpoint(&stroke2 - graph.orig_strokes_.data(), arclen2 == 0.0);
  force_assert(end2.stroke_idx() < graph.orig_strokes_.size());

  auto key = std::make_pair(end1.as_int(), end2.as_int());
  if (key.first > key.second) {
    std::swap(key.first, key.second);
  }
  const auto hash =
    uint64_t(key.first) << 32 | (uint64_t(uint32_t(key.second)));
  auto it = graph.parallel_endpoints_cache_.find(hash);
  if (it != graph.parallel_endpoints_cache_.end()) {
    return it->second;
  } else {
    const auto result = is_parallel_slow(stroke1, arclen1, stroke2, arclen2);
    graph.parallel_endpoints_cache_[hash] = result;
    return result;
  }
}

bool end_end_pre_filter(VertexView vertex, VertexView other_vertex,
                        bool train_time) {
  const auto &graph = *vertex.graph_;
  if (other_vertex.is_active() &&
      is_snap_valid(graph, vertex.index_, other_vertex.index_)) {
    auto has_original_endpoints = false;
    auto too_far = false;
    auto has_diverging_type = false;
    auto has_parallel = false;
    auto cen_dist = std::numeric_limits<Float>::quiet_NaN();
    const auto he1 = vertex.hedge();
    const auto he2 = other_vertex.hedge();
    auto it1 = he1;
    auto iterations = 0;
    do {
      auto it2 = he2;
      do {
        auto end1 =
          StrokeTime((int)it1.stroke_idx(), it1.forward() ? 0.0 : 1.0);
        auto end2 =
          StrokeTime((int)it2.stroke_idx(), it2.forward() ? 0.0 : 1.0);
        if (convert_strokes2orig(graph, end1) &&
            convert_strokes2orig(graph, end2) &&
            (end1.second == 0 || end1.second == 1) &&
            (end2.second == 0 || end2.second == 1)) {
          has_original_endpoints = true;
          const auto &s1 = graph.orig_strokes_[end1.first];
          const auto &s2 = graph.orig_strokes_[end2.first];
          if (endpoint_too_far(s1, end1, s2, end2, &cen_dist)) {
            too_far = true;
            goto jump_out; // Jump out of nested do-whiles.
          }
          if (end_end_junction_type(s1, end1.second * s1.length(), s2,
                                    end2.second * s2.length()) == 2) {
            has_diverging_type = true;
            goto jump_out; // Jump out of nested do-whiles.
          }
          if (!train_time && vertex.is_originally_dangling() &&
              other_vertex.is_originally_dangling() &&
              is_parallel(graph, s1, end1.second * s1.length(), s2,
                          end2.second * s2.length())) {
            has_parallel = true;
            goto jump_out; // Jump out of nested do-whiles.
          }
        }

        it2 = it2.twin().next();
        iterations++;
        force_assert(iterations < 1024 && "likely infinite loop found");
      } while (it2 != he2);
      it1 = it1.twin().next();
    } while (it1 != he1);

jump_out:
    if (has_original_endpoints && !too_far && !has_diverging_type &&
        !has_parallel) {
      force_assert(std::isfinite(cen_dist));
      return true;
    }
  }
  return false;
}

namespace {

void find_end_end_candidates_stroke(const VertexView vertex,
                                    span<StrokeTime> &out_candidates,
                                    bool train_time) {
  if (out_candidates.empty()) {
    return;
  }

  const auto &graph = *vertex.graph_;
  const auto centerline_bbs = graph.bvh_.centerline_bbs();
  auto hits = std::multimap<Float, size_t>(); // Sort by distance.

  const auto n_vertices = graph.vertices_.size();
  for (size_t i = 0; i < n_vertices; ++i) {
    const auto other_vertex = graph.vertex(i);
    if (other_vertex.is_active()) {
      // Check if endpoint.
      const auto other_he = other_vertex.hedge();
      const auto cen_dist = (vertex.pos() - other_vertex.pos()).norm();
      auto it = other_he;
      do {
        auto end = StrokeTime((int)it.stroke_idx(), it.forward() ? 0.0 : 1.0);
        if (convert_strokes2orig(graph, end) &&
            (end.second == 0.0 || end.second == 1.0)) {
          hits.insert({cen_dist, other_vertex.index_});
          break;
        }
        it = it.twin().next();
      } while (it != other_he);
    }
  }

  size_t i = 0;
  auto n_hits_looked_at = 0;
  for (const auto &hit : hits) {
    const auto other_vertex = graph.vertex(hit.second);
    if (end_end_pre_filter(vertex, other_vertex, train_time) &&
        line_of_sight(vertex.pos(), other_vertex.pos(), graph.strokes_,
                      centerline_bbs)) {
      const auto he = other_vertex.hedge();
      auto it = he;
      do {
        // Find a StrokeTime for `other_vertex` that corresponds to an original
        // endpoint. This makes future conversions to `Junction` using original
        // indexing easier.
        const auto edge_st =
          StrokeTime((int)it.stroke_idx(), it.forward() ? 0.0 : 1.0);
        auto st = edge_st;
        if (convert_strokes2orig(graph, st) &&
            (st.second == 0.0 || st.second == 1.0)) {
          out_candidates[i++] = edge_st;
          break;
        }
        it = it.twin().next();
      } while (it != he);
    }
    n_hits_looked_at++;
    if (n_hits_looked_at >= out_candidates.size()) {
      break;
    }
  }
  out_candidates = out_candidates.first(i);
}

void tjunc_classifier_candidates(
  const VertexView vertex,
  const std::unordered_set<std::size_t> &ignored_strokes, const int n_closest,
  std::vector<StrokeTime> &out_candidates, bool filter_close_to_end = true) {
  out_candidates.clear();
  if (n_closest < 1)
    throw std::invalid_argument("n_closest must be at least 1");

  assert(ignored_strokes.empty() && "ignored_strokes no longer supported");

  if (!vertex.is_active() || !vertex.is_dangling()) {
    return;
  }
  const StrokeGraph &graph = *vertex.graph_;
  const auto si = vertex.hedge().orig_stroke_idx();
  const auto &s = graph.orig_strokes_[si];
  const auto endp = Endpoint(si, vertex.hedge().forward());
  const auto hits = pick(graph, endp, n_closest, &ignored_strokes);

  auto end1 = StrokeTime((int)vertex.hedge().stroke_idx(),
                         vertex.hedge().forward() ? 0.0 : 1.0);
  const auto ok = convert_strokes2orig(graph, end1);
  force_assert(ok && "couldn't map from strokes to orig");
  if (!((endp.is_head() && end1.second == 0.0) ||
        (!endp.is_head() && end1.second == 1.0))) {
    return; // Cannot form a T-junction starting from an original stroke
            // interior.
  }
  assert((size_t)end1.first == si);

  for (const auto &[other_si, other_norm_arclen] : hits) {
    const auto &other_stroke = graph.orig_strokes_[other_si];
    const auto other_arclen = other_norm_arclen * other_stroke.length();
    const auto p = (endp.is_head() ? s.xy(0) : s.xy(Back));
    const auto q = other_stroke.pos(other_arclen);
    const auto cen_dist = (p - q).norm();
    const auto env_dist = cen_dist - 0.5 * other_stroke.width_at(other_arclen) -
                          vertex_radius(vertex);
    if ((env_dist < env_distance_over_pen_width_hard_threshold * 0.5 *
                      (s.pen_width() + other_stroke.pen_width()) ||
         env_dist < std::min(s.length(), other_stroke.length())) &&
        (!filter_close_to_end ||
         end_stroke_junction_type(other_stroke, other_arclen) != 0) &&
        // !intersecting_and_diverging_t(s, end1.second * s.length(),
        // other_stroke,
        //                               other_arclen) &&
        line_of_sight(p, q, *graph.orig_bvh_)) {

#if 1
      auto end = StrokeTime(-1, -1.0);
      auto best_dist = Float(INFINITY);
      for (const auto &[edge_idx, map] : graph.orig2strokes_[other_si]) {
        Vec2 proj;
        Float arclen = 0;
        graph.strokes_[edge_idx].ensure_arclengths();
        const auto dist = //
          (edge_idx == vertex.hedge().stroke_idx()
             ? closest_point_to_own_stroke(graph.strokes_[edge_idx],
                                           /*head=*/vertex.hedge().forward(),
                                           proj, arclen)
             : closest_point(graph.strokes_[edge_idx], p, proj, arclen));
        if (dist < best_dist) {
          best_dist = dist;
          end = StrokeTime((int)edge_idx,
                           arclen / graph.strokes_[edge_idx].length());
        }
      }
      assert(end.first != -1);
      assert(end.second != -1.0);
#else
      StrokeTime end(other_si, other_norm_arclen);
      const auto ok2 = convert_orig2strokes(graph, end);
      force_assert(ok2 && "couldn't map from orig to strokes");
#endif

      out_candidates.emplace_back(end);
    }
  }
}

void find_end_stroke_candidates(const VertexView vertex,
                                span<StrokeTime> &out_candidates) {
  if (out_candidates.empty()) {
    return;
  }

  const auto &graph = *vertex.graph_;
  auto hits = std::multimap<Float, StrokeTime>(); // Sort by distance.

  const auto nstrokes = graph.strokes_.size();
  const auto he = vertex.hedge();
  const auto xy = vertex.pos();
  const auto tangent = (he.forward() ? head_tangent_lagrange(he.stroke())
                                     : tail_tangent_lagrange(he.stroke()));
  for (auto i = decltype(nstrokes){0}; i < nstrokes; ++i) {
    const auto &other_stroke = graph.strokes_[i];
    if (other_stroke.size() > 1 &&
        !(graph.hedges_[2 * i].flags_ & StrokeGraph::HedgeRecord::Bridge)) {
      auto proj = Vec2::Empty();
      Float s;
      const auto cen_dist = find_projection(vertex, i, proj, s);

      auto has_original_endpoints = false;
      auto too_far = false;
      auto at_endpoint = false;
      auto is_original_spurious = false;

      auto end1 = StrokeTime((int)he.stroke_idx(), he.forward() ? 0.0 : 1.0);
      auto end2 = StrokeTime((int)i, s / other_stroke.length());
      if (convert_strokes2orig(graph, end1) &&
          convert_strokes2orig(graph, end2)) {
        has_original_endpoints = true;
        const auto &s1 = graph.orig_strokes_[end1.first];
        const auto &s2 = graph.orig_strokes_[end2.first];
        s2.ensure_arclengths();

        const auto p = s1.pos(end1.second * s1.length());
        auto orig_s2 = Float(-1);
        Vec2 _proj;
        if (end1.first == end2.first) {
          closest_point_to_own_stroke(s2, /*head=*/end1.first == 0, _proj,
                                      orig_s2);
        } else {
          closest_point(s2, p, _proj, orig_s2);
        }
        if (orig_s2 == -1) {
          is_original_spurious = true;
        } else {
          if (end_stroke_junction_type(s2, orig_s2) == 0) {
            at_endpoint = true;
          }

          const auto env_dist =
            cen_dist - 0.5 * (s1.width_at(end1.second * s1.length()) +
                              s2.width_at(orig_s2));
          if (env_dist > env_distance_over_pen_width_hard_threshold * 0.5 *
                           (s1.pen_width() + s2.pen_width()) &&
              env_dist > std::min(s1.length(), s2.length())) {
            too_far = true;
          }
        }
      }

      constexpr auto eps = 1e-4;
      if (cen_dist < infinity && //
          has_original_endpoints && //
          !is_original_spurious && //
          !too_far && //
          !at_endpoint && //
          end_stroke_junction_type(graph.hedge(2 * i), s) != 0 && //
          // For numerical reasons, let's not create candidates to close to
          // vertices.
          eps < s && s < other_stroke.length() - eps && //
          is_snap_valid(graph, vertex.index_,
                        {(int)i, s / other_stroke.length()})) {
        const auto angle =
          acos(tangent.dot((other_stroke.pos(s) - xy).normalized()));
        if (angle <= 180_deg - 45_deg) {
          hits.insert({cen_dist, {(int)i, s / other_stroke.length()}});
        }
      }
    }
  }

  const auto centerline_bbs = graph.bvh_.centerline_bbs();
  size_t i = 0;
  for (const auto &[_dist, cand] : hits) {
    const auto [other_si, other_narclen] = cand;
    const auto &other_stroke = graph.strokes_[other_si];
    if (line_of_sight(vertex.pos(), other_stroke.pos_norm(other_narclen),
                      graph.strokes_, centerline_bbs)) {
      out_candidates[i++] = cand;
      if (i == out_candidates.size()) {
        break;
      }
    }
  }
  out_candidates = out_candidates.first(i);
}

/**
 * Output valid vertex-vertex candidates that are also in the input
 * `intersection_set` (could be `enforced_connections_` or
 * `allowed_connections_`).
 */
void find_filtered_end_end_candidates(const VertexView vertex,
                                      const IncrementalCache &cache,
                                      const Connections &intersection_set,
                                      FeatureType feature_type,
                                      span<StrokeTime> &out_candidates,
                                      bool train_time) {
  if (out_candidates.empty()) {
    return;
  }

  const auto &graph = *vertex.graph_;
  const auto centerline_bbs = graph.bvh_.centerline_bbs();

  if (intersection_set.empty()) {
    if (feature_type == FeatureType::Graph)
      find_end_end_candidates(vertex, out_candidates);
    else if (feature_type == FeatureType::OrigStroke)
      find_end_end_candidates_stroke(vertex, out_candidates, train_time);
    else
      abort();
  } else {
    auto orig_endp = hedge_to_endpoint(vertex.hedge()).as_pair();
    {
      const auto ok = convert_strokes2orig(graph, orig_endp);
      force_assert(ok && "couldn't map from strokes to orig");
    }
    const auto last_orig_si = (int)graph.orig_strokes_.size() - 1;
    auto hits = std::multimap<Float, StrokeTime>(); // Sort by distance.
    for (const auto &conn :
         intersection_set.associated_connections(orig_endp)) {
      auto other = conn.second;
      if ((size_t)other.first >= graph.orig2strokes_.size() ||
          ((other.second == 0 || other.second == 1) &&
           cache.enforced_connections_.largest_connected_stroke_idx(other) >
             last_orig_si)) {
        continue; // Wait for future strokes to arrive...
      }
      if (!convert_orig2strokes(graph, other)) { // Note the side effect.
        continue; // This original stroke was dropped.
      }
      assert(last_orig_si < cache.n_orig_strokes_);
      const auto other_v =
        endpoint_to_vertex(graph, Endpoint(other.first, other.second == 0));
      if ((other.second == 0 || other.second == 1) &&
          is_snap_valid(graph, vertex.index_, other_v.index_) &&
          line_of_sight(vertex.pos(), other_v.pos(), graph.strokes_,
                        centerline_bbs) &&
          line_of_sight(vertex.pos(), other_v.pos(),
                        {&cache.orig_strokes_[last_orig_si + 1],
                         cache.n_orig_strokes_ - last_orig_si - 1},
                        {&cache.orig_strokes_bb_[last_orig_si + 1],
                         cache.n_orig_strokes_ - last_orig_si - 1})) {
        const auto dist = (vertex.pos() - other_v.pos()).norm();
        hits.insert({dist, other});
      }
    }

    size_t i = 0;
    for (const auto &[_dist, other] : hits) {
      out_candidates[i++] = other;
      if (i >= out_candidates.size())
        break;
    }
    out_candidates = out_candidates.first(i);
  }
}

/**
 * Return valid vertex-stroke candidates after applying the constraints in
 * the input `intersection_set` (could be `enforced_connections_` or
 * `allowed_connections_`).
 */
void find_filtered_end_stroke_candidates(const VertexView vertex,
                                         const IncrementalCache &cache,
                                         const Connections &intersection_set,
                                         span<StrokeTime> &out_candidates) {
  if (out_candidates.empty()) {
    return;
  }

  const auto &graph = *vertex.graph_;

  if (intersection_set.empty()) {
    find_end_stroke_candidates(vertex, out_candidates);
  } else {
    auto orig_endp = hedge_to_endpoint(vertex.hedge()).as_pair();
    {
      const auto ok = convert_strokes2orig(graph, orig_endp);
      force_assert(ok && "couldn't map from strokes to orig");
    }
    const auto last_orig_si = (int)graph.orig_strokes_.size() - 1;
    auto hits = std::multimap<Float, StrokeTime>(); // Sort by distance.
    for (const auto &conn :
         intersection_set.associated_connections(orig_endp)) {
      auto other = conn.second;
      if ((size_t)other.first >= graph.orig2strokes_.size() ||
          ((other.second == 0 || other.second == 1) &&
           cache.enforced_connections_.largest_connected_stroke_idx(other) >
             last_orig_si)) {
        continue; // Wait for future strokes to arrive...
      }
      if (!convert_orig2strokes(graph, other)) { // Note the side effect.
        continue; // This original stroke was dropped.
      }
      if (other.second != 0 && other.second != 1) {
        const auto orig_si = int{conn.second.first};
        // Find the closest edge that originates from original stroke orig_si.
        // Note this assumes we only need to generate at most one candidate for
        // all the edges associated with orig_si.
        const auto maybe_cand2 = reprojected_vertex_stroke_candidate(
          vertex, /*orig_si=*/orig_si,
          /*last_orig_si=*/(int)graph.orig_strokes_.size() - 1, cache);
        const auto cand2 = StrokeTime{maybe_cand2.value_or(StrokeTime(0, 0))};
        const auto &other_stroke = graph.strokes_[cand2.first];
        if (maybe_cand2.has_value() &&
            end_stroke_junction_type(
              other_stroke, cand2.second * other_stroke.length()) != 0) {
          const auto dist =
            (vertex.pos() - other_stroke.pos_norm(cand2.second)).norm();
          hits.insert({dist, cand2});
        }
      }
    }

    size_t i = 0;
    for (const auto &[_dist, other] : hits) {
      out_candidates[i++] = other;
      if (i >= out_candidates.size())
        break;
    }
    out_candidates = out_candidates.first(i);
  }
}

SnapInfo snap_with_classifier(
  StrokeGraph &graph, const size_t vi, IncrementalCache *const cache,
  FeatureType feature_type, const Float accept_threshold, bool train_time,
  std::vector<ClassifierPrediction> *out_predictions = nullptr) {
  auto snap_info = SnapInfo();
  const auto vertex = graph.vertex(vi);
  assert(vertex && "vertex was invalidated");

  // Note we can snap *to* a vertex with a continuous edge, but we cannot snap
  // *from* such a vertex. (The difference is in the position of the resulting
  // vertex after snapping.) We also shouldn't snap from vertices that weren't
  // originally dangling.
  if (has_continuous_edge(vertex) || !vertex.is_originally_dangling()) {
    snap_info.status = SnapInfo::Skip;
    return snap_info;
  }

  graph.orig_bvh_ = std::make_unique<PolylineBVH>(graph.orig_strokes_);

  auto orig_endp = hedge_to_endpoint(vertex.hedge()).as_pair();
  const auto ok = convert_strokes2orig(graph, orig_endp);
  force_assert(ok && "couldn't map from strokes to orig");
  const auto last_orig_si = (int)graph.orig_strokes_.size() - 1;
  force_assert(last_orig_si >= 0);
  const auto enforced_connections =
    cache->enforced_connections_.associated_connections(orig_endp);
  for (const auto &conn : enforced_connections) {
    const auto other_orig_si = int{conn.second.first};
    if (other_orig_si > last_orig_si) {
      // The endpoint needs to connect to a later stroke first.
      snap_info.status = SnapInfo::ReservedDelay;
      return snap_info;
    }
  }

  const auto cand1_val =
    Endpoint(vertex.hedge().stroke_idx(), vertex.hedge().forward());

  StrokeTime candidates_buf[n_candidates];
  auto cand2 = span<StrokeTime>(candidates_buf, 3);
  find_filtered_end_end_candidates(vertex, *cache, cache->allowed_connections_,
                                   feature_type, cand2, train_time);
  if (!cand2.empty()) {
    auto cand1 = ALLOCA_SPAN(Endpoint, cand2.size());
    for (size_t i = 0; i < cand2.size(); ++i) {
      cand1[i] = cand1_val;
    }
    auto prob = ALLOCA_SPAN(Float, cand2.size());
    compute_junction_probabilities(graph, cand1, cand2, prediction_feature_type,
                                   prob, out_predictions);

    // Index in order of decreasing probability.
    auto indices = ALLOCA_SPAN(size_t, cand2.size());
    for (size_t i = 0; i < cand2.size(); ++i) {
      indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(),
              [=](size_t a, size_t b) { return prob[a] > prob[b]; });

    for (const auto i : indices) {
      const auto [si2, narclen2] = cand2[i];
      const auto other_endp = Endpoint(si2, narclen2 == 0);
      const auto other_vi = endpoint_to_vertex(graph, other_endp).index_;
      snap_info.max_prob = prob[i];
      snap_info.prediction_type = JunctionType::R;
      snap_info.other_vertex = other_vi;

      if (prob[i] <= accept_threshold) {
        snap_info.status = SnapInfo::PredictionDelay;
        break;
      }

      // Try to snap the edge.
      if (graph.snap_vertices(vi, other_vi, graph.snapping_method_type_)) {
        // Snapping succeeded.
        snap_info.status = SnapInfo::PredictionPos;
        return snap_info;
      }
    }
  }

  {
    // Try finding a T-junction.
    MSVC_WARNING_SUPPRESS(4456); // declaration hides previous local declaration
    StrokeTime candidates_buf[3]{{0, 0}, {0, 0}, {0, 0}};
    MSVC_WARNING_SUPPRESS(4456);
    auto cand2 = span<StrokeTime>(candidates_buf, 3);
    find_filtered_end_stroke_candidates(vertex, *cache,
                                        cache->allowed_connections_, cand2);

    if (!cand2.empty()) {
      auto cand1 = ALLOCA_SPAN(Endpoint, cand2.size());
      for (size_t i = 0; i < cand2.size(); ++i) {
        cand1[i] = cand1_val;
      }
      auto prob = ALLOCA_SPAN(Float, cand2.size());
      compute_junction_probabilities(
        graph, cand1, cand2, prediction_feature_type, prob, out_predictions);

      // Index in order of decreasing probability.
      auto indices = ALLOCA_SPAN(size_t, cand2.size());
      for (size_t i = 0; i < cand2.size(); ++i) {
        indices[i] = i;
      }
      std::sort(indices.begin(), indices.end(),
                [=](size_t a, size_t b) { return prob[a] > prob[b]; });

      for (const auto i : indices) {
        if (prob[i] <= accept_threshold) {
          snap_info.status = SnapInfo::PredictionDelay;
          if (prob[i] > snap_info.max_prob) {
            snap_info.max_prob = prob[i];
            snap_info.prediction_type = JunctionType::T;
          }
          break;
        }

        // Try to snap the edge.
        const auto [other_si, other_narclen] = cand2[i];
        const auto snap_vertex = graph.snap_endpoint_to_edge(
          vertex, other_si, other_narclen * graph.strokes_[other_si].length(),
          graph.snapping_method_type_);
        if (snap_vertex.is_valid() && snap_vertex != vertex) {
          // Snap succeeded.
          assert(!graph.vertices_[vertex.index_].is_active());
          snap_info.status = SnapInfo::PredictionPos;
          snap_info.max_prob = prob[i];
          snap_info.prediction_type = JunctionType::T;
          snap_info.other_vertex = snap_vertex.index_;
          return snap_info;
        }
      }
    }
  }
  return snap_info;
}

} // namespace

bool endpoint_too_far(const Stroke &s1, const StrokeTime &end1,
                      const Stroke &s2, const StrokeTime &end2,
                      Float *const out_cen_dist) {
  const auto p1 = s1.pos_norm(end1.second);
  const auto p2 = s2.pos_norm(end2.second);
  const auto cen_dist = (p1 - p2).norm();
  if (out_cen_dist) {
    *out_cen_dist = cen_dist;
  }
  const auto env_dist =
    cen_dist - 0.5 * (s1.width_at(end1.second * s1.length()) +
                      s2.width_at(end2.second * s2.length()));
  return (env_dist > env_distance_over_pen_width_hard_threshold * 0.5 *
                       (s1.pen_width() + s2.pen_width()) &&
          env_dist > std::min(s1.length(), s2.length()));
}

static std::pair<size_t, size_t>
find_original_endpoint_vertices(const StrokeGraph &graph,
                                const Index orig_idx) {
  auto head_vertex_idx = invalid;
  auto tail_vertex_idx = invalid;
  auto head_position = StrokeTime(orig_idx, 0);
  if (convert_orig2strokes(graph, head_position)) {
    const auto he = graph.hedge(2 * head_position.first);
    if (head_position.second == 0) {
      head_vertex_idx = he.origin().index_;
    } else if (head_position.second == 1) {
      head_vertex_idx = he.dest().index_;
    }
  }
  auto tail_position = StrokeTime(orig_idx, 1);
  if (convert_orig2strokes(graph, tail_position)) {
    const auto he = graph.hedge(2 * tail_position.first);
    if (tail_position.second == 0) {
      tail_vertex_idx = he.origin().index_;
    } else if (tail_position.second == 1) {
      tail_vertex_idx = he.dest().index_;
    }
  }
  return {head_vertex_idx, tail_vertex_idx};
}

Float find_projection(const VertexView vertex, const size_t si, Vec2 &proj,
                      Float &s) {
  const auto &graph = *vertex.graph_;
  const auto vertex_is_origin =
    (graph.hedges_[2 * si].origin_ == vertex.index_);
  const auto vertex_is_dest =
    (graph.hedges_[2 * si + 1].origin_ == vertex.index_);
  const auto &other_stroke = graph.strokes_[si];
  other_stroke.ensure_arclengths();
  return (
    (vertex_is_origin || vertex_is_dest)
      ? closest_point_to_own_stroke(other_stroke, vertex_is_origin, proj, s)
      : closest_point(other_stroke, vertex.pos(), proj, s));
}

std::optional<StrokeTime>
reprojected_vertex_stroke_candidate(const VertexView v, const int orig_si,
                                    const int last_orig_si,
                                    const IncrementalCache &cache) {
  const auto &graph = *v.graph_;
  auto best_dist = infinity;
  auto best_si = 0;
  auto best_arclen = Float(0);
  for (const auto &[si, _] : graph.orig2strokes_[orig_si]) {
    Vec2 proj;
    Float arclen;
    const auto dist = find_projection(v, si, proj, arclen);
    if (dist < best_dist) {
      best_dist = dist;
      best_si = (int)si;
      best_arclen = arclen;
    }
  }
  if (best_dist < infinity) {
    const auto cand =
      StrokeTime{best_si, best_arclen / graph.strokes_[best_si].length()};
    const auto p = v.pos();
    const auto q = graph.strokes_[best_si].pos(best_arclen);
    if (is_snap_valid(graph, v.index_, cand) &&
        line_of_sight(p, q, graph.strokes_, graph.bvh_.centerline_bbs()) &&
        line_of_sight(p, q,
                      {&cache.orig_strokes_[last_orig_si + 1],
                       cache.n_orig_strokes_ - last_orig_si - 1},
                      {&cache.orig_strokes_bb_[last_orig_si + 1],
                       cache.n_orig_strokes_ - last_orig_si - 1})) {
      return cand;
    }
  }
  return std::nullopt;
}

bool should_dissolve_vertex(const VertexView v) {
  if (v.is_active() && v.valence() == 2) {
    const auto he1 = v.hedge();
    const auto he2 = he1.twin().next();
    const Stroke *stroke1 = &he1.stroke();
    const Stroke *stroke2 = &he2.stroke();
    auto sample_head1 = he1.forward();
    auto sample_head2 = he2.forward();
    // If we have a real centerline intersection, we can ignore the overshot
    // region and look at the edge tangents directly.  But if it's only an
    // envelope overlap, then our deformation will have messed up the tangents,
    // so we should use the original tangents.
    if (!(v.flags() & VertexRecord::CenterlineIntersection)) {
      auto end1 = StrokeTime((int)he1.stroke_idx(), he1.forward() ? 0.0 : 1.0);
      if (!convert_strokes2orig(*v.graph_, end1)) {
        return false;
      }
      auto end2 = StrokeTime((int)he2.stroke_idx(), he2.forward() ? 0.0 : 1.0);
      if (!convert_strokes2orig(*v.graph_, end2)) {
        return false;
      }
      if (end1.second != 0.0 && end1.second != 1.0) {
        return false;
      }
      if (end2.second != 0.0 && end2.second != 1.0) {
        return false;
      }
      stroke1 = &v.graph_->orig_strokes_[end1.first];
      stroke2 = &v.graph_->orig_strokes_[end2.first];
      sample_head1 = (end1.second == 0.0);
      sample_head2 = (end2.second == 0.0);
    }

    for (const auto factor : {1.0, 2.0}) {
      const auto v1 = (sample_head1 ? head_tangent_lagrange(*stroke1, factor)
                                    : tail_tangent_lagrange(*stroke1, factor));
      const auto v2 = -(sample_head2 ? head_tangent_lagrange(*stroke2, factor)
                                     : tail_tangent_lagrange(*stroke2, factor));
      const auto turn_angle = std::acos(std::clamp(v1.dot(v2), -1.0, 1.0));
      if (turn_angle < 20_deg) {
        return true;
      }
    }
  }
  return false;
}

std::unique_ptr<IncrementalCache>
make_incremental_cache(const Float accept_threshold) {
  auto cache = std::make_unique<IncrementalCache>();
  if (accept_threshold < 0.0) {
    cache->accept_threshold_ = 0.8;
  } else {
    cache->accept_threshold_ = accept_threshold;
  }
  cache->intersections_.reserve(8);
  return cache;
}

/**
 * Enforce the connections from the input strokes.  Additional connections will
 * also be allowed, but the algorithm will prioritize these first.
 *
 * Note that these constraints are not hard ones.  The algorithm is free to drop
 * constraints if they would result in bad topology.
 */
void set_enforced_connections(IncrementalCache *cache,
                              const EnvelopeBVH &strokes) {
  cache->enforced_connections_ = Connections(strokes);
}

/**
 * Allow all potential connections from the input graph (and no others, besides
 * overlapping edges and the ones in `enforced_connections_`).
 */
static void add_allowed_connections(IncrementalCache *const cache,
                                    const StrokeGraph &graph,
                                    const size_t num_candidates,
                                    FeatureType feature_type, bool train_time) {
  auto all_cand1 = std::vector<StrokeTime>(); // Uses original indexing.
  auto all_cand2 = std::vector<StrokeTime>(); // Uses original indexing.
  auto cand2_buf = ALLOCA(StrokeTime, num_candidates);
  for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
    const auto v = graph.vertex(vi);
    if (v.is_active()) {
      // Add end-end candidates.
      auto cand2 = span<StrokeTime>(cand2_buf, num_candidates);
      if (feature_type == FeatureType::Graph)
        find_end_end_candidates(v, cand2);
      else if (feature_type == FeatureType::OrigStroke)
        find_end_end_candidates_stroke(v, cand2, train_time);
      else
        abort();

      if (!cand2.empty()) {
        const auto he = v.hedge();
        auto it = he;
        do {
          auto c1 = hedge_to_endpoint(it).as_pair();
          if (convert_strokes2orig(graph, c1)) {
            if (c1.second == 0 || c1.second == 1) {
              for (const auto &pre_expansion_c2 : cand2) {
                const auto he2 = endpoint_to_hedge(
                  graph, Endpoint(pre_expansion_c2.first,
                                  pre_expansion_c2.second == 0));
                auto it2 = he2;
                do {
                  auto c2 = hedge_to_endpoint(it2).as_pair();
                  if (convert_strokes2orig(graph, c2)) {
                    if ((c2.second != 0 && c2.second != 1) ||
                        c1 < c2) { // Deduplicate.
                      all_cand1.push_back(c1);
                      all_cand2.push_back(c2);
                    }
                    // Create a corresponding end-stroke candidate.
                    // This is because the vertex may get dissolved due to
                    // deformation.
                    all_cand1.push_back(c1);
                    all_cand2.emplace_back(c2.first, 0.5);
                  }

                  it2 = it2.twin().next();
                } while (it2 != he2);
              }
            }
          }
          it = it.twin().next();
        } while (it != he);
      }

      // Add end-stroke candidates.
      cand2 = span<StrokeTime>(cand2_buf, num_candidates);
      find_end_stroke_candidates(v, cand2);
      if (!cand2.empty()) {
        const auto he = v.hedge();
        auto it = he;
        do {
          auto c1 = hedge_to_endpoint(it).as_pair();
          if (convert_strokes2orig(graph, c1)) {
            if (c1.second == 0 || c1.second == 1) {
              for (auto c2 : cand2) {
                if (convert_strokes2orig(graph, c2)) {
                  all_cand1.push_back(c1);
                  all_cand2.push_back(c2);
                }
              }
            }
          }
          it = it.twin().next();
        } while (it != he);
      }
    }
  }

  cache->allowed_connections_.insert(all_cand1, all_cand2);
}

void set_future_constraints(IncrementalCache *cache,
                            const span<const Stroke> strokes,
                            const size_t num_candidates, bool train_time) {
  cache->n_orig_strokes_ = strokes.size();
  cache->orig_strokes_.reset(new Stroke[cache->n_orig_strokes_]);
  cache->orig_strokes_bb_.reset(new BoundingBox[cache->n_orig_strokes_]);
  for (size_t i = 0; i < strokes.size(); ++i) {
    cache->orig_strokes_[i] = strokes[i].clone();
    cache->orig_strokes_bb_[i] = bounds(cache->orig_strokes_[i]);
  }
  set_enforced_connections(cache, EnvelopeBVH(strokes));

  {
    auto final_cache = IncrementalCache();
    final_cache.enforced_connections_ = cache->enforced_connections_;
    auto final_graph = StrokeGraph(StrokeGraph::SnappingType::Connection);
    for (size_t i = 0; i < strokes.size(); ++i) {
      add_stroke_incremental_topological(final_graph, &final_cache, strokes[i],
                                         false);
    }
    add_allowed_connections(cache, final_graph, num_candidates,
                            prediction_feature_type, train_time);
  }
  {
    auto final_cache = IncrementalCache();
    final_cache.enforced_connections_ = cache->enforced_connections_;
    auto final_graph = StrokeGraph(StrokeGraph::SnappingType::Deformation);
    for (size_t i = 0; i < strokes.size(); ++i) {
      add_stroke_incremental_topological(final_graph, &final_cache, strokes[i],
                                         false);
    }
    // Sometimes, candidates will be blocked in the straight-line connection
    // version of the graph, but not when deforming.  We allow a connection if
    // it is valid in either.
    add_allowed_connections(cache, final_graph, num_candidates,
                            prediction_feature_type, train_time);
  }
}

StrokeSnapInfo
add_stroke_incremental_topological(StrokeGraph &graph, IncrementalCache *cache,
                                   const Stroke &new_stroke, bool to_dissolve) {
  auto info = StrokeSnapInfo();
  graph.orig_strokes_.emplace_back(new_stroke.clone());
  if (!graph.orig_bvh_)
    graph.orig_bvh_ = std::make_unique<PolylineBVH>();
  graph.orig_bvh_->nodes.emplace_back(graph.orig_strokes_.back(),
                                      bounds(graph.orig_strokes_.back()));
  for (size_t i = 0; i < graph.orig_strokes_.size(); ++i) {
    // Is there a better way to deal with pointer invalidation in std::vector?
    graph.orig_bvh_->nodes[i].geometry = &graph.orig_strokes_[i];
  }

  if (new_stroke.size() <= 1) {
    graph.orig2strokes_.emplace_back();
    assert(graph.orig_strokes_.size() == graph.orig2strokes_.size());
    info.head.status = info.tail.status = SnapInfo::Skip;
    return info;
  }

  new_stroke.ensure_arclengths();

  const auto new_bb = bounds(new_stroke);
  const auto new_node = PolylineBVHLeaf(new_stroke, new_bb);

  std::vector<CutPoint> cut_points;

  auto &intersections = cache->intersections_;
  const auto n_strokes = graph.strokes_.size();
  for (size_t si = 0; si < n_strokes; ++si) {
    // Check for intersections.
    cache->intersections_.clear();
    graph.strokes_[si].ensure_arclengths();
    intersect_different(new_node, graph.bvh_.polyline_bvh_leaf(si),
                        intersections);
    if (!cache->intersections_.empty()) {
      auto &old_stroke = graph.strokes_[si];

      // Prepare to cut pre-existing edge.
      std::sort(intersections.begin(), intersections.end(),
                [](const Vec2 &a, const Vec2 &b) { return a.y_ < b.y_; });
      bool trim_start = false, trim_end = false;
      if (cache->trim_overshoots_) {
        trim_start =
          intersections[0].y_ <=
          (graph.hedge(2 * si).origin().is_dangling()
             ? 0.5 * std::max(new_stroke.width_at(intersections[0].x_),
                              old_stroke.width_at(intersections[0].y_))
             : 1e-10);
        trim_end =
          intersections.back().y_ >=
          (graph.hedge(2 * si + 1).origin().is_dangling()
             ? old_stroke.length() -
                 0.5 * std::max(new_stroke.width_at(intersections.back().x_),
                                old_stroke.width_at(intersections.back().y_))
             : old_stroke.length() - 1e-10);
      } else {
        trim_start = (intersections[0].y_ <= 1e-10);
        trim_end = (intersections.back().y_ >= old_stroke.length() - 1e-10);
      }
      if (intersections.size() == 1 && trim_start && trim_end) {
        const auto s = intersections[0].y_;
        trim_start = trim_end = false;
        if (s <= 1e-10) {
          trim_start = true;
        } else if (s >= old_stroke.length() - 1e-10) {
          trim_end = true;
        }
      }
      // @optimize Create a bump allocator for small allocations like these.
      auto split_values_unique = std::make_unique<Float[]>(
        intersections.size() + (trim_start ? 0 : 1) + (trim_end ? 0 : 1));
      auto split_values = span<Float>(
        split_values_unique.get(),
        intersections.size() + (trim_start ? 0 : 1) + (trim_end ? 0 : 1));
      {
        size_t i = 0;
        if (!trim_start) {
          split_values[i++] = 0.0;
        }
        for (const auto &intersection : intersections) {
          split_values[i++] = intersection.y_;
        }
        if (!trim_end) {
          split_values[i] = old_stroke.length();
        }
      }

      if (split_values.size() == 2) {
        // No need to split; only trim.
        const auto new_domain =
          std::pair<Float, Float>(split_values[0] / old_stroke.length(),
                                  split_values[1] / old_stroke.length());
        old_stroke.trim(split_values[0], split_values[1]);
        clamp_edge_mapping(graph, si, new_domain);
        const auto he = graph.hedge(2 * si);
        if (trim_start) {
          const auto orig_vi = he.origin().index_;
          cut_points.push_back({intersections[0].x_, orig_vi, false});
          graph.vertices_[orig_vi].p_ = old_stroke.xy(0);
          graph.vertices_[orig_vi].flags_ |=
            VertexRecord::CenterlineIntersection;
        }
        if (trim_end) {
          const auto dest_vi = he.dest().index_;
          cut_points.push_back({intersections.back().x_, dest_vi, false});
          graph.vertices_[dest_vi].p_ = old_stroke.xy(Back);
          graph.vertices_[dest_vi].flags_ |=
            VertexRecord::CenterlineIntersection;
        }
        graph.bvh_.full_update(si);
      } else {
        const auto start_index = graph.strokes_.size();
        MSVC_WARNING_SUPPRESS(
          4456); // declaration hides previous local declaration
        const auto old_stroke = std::move(graph.strokes_[si]);
        graph.bvh_.clear(si);

        Float old_stroke_length = old_stroke.length();
        std::vector<std::pair<Float, Float>> domain_intervals;
        domain_intervals.reserve(split_values.size() - 1);
        for (size_t i = 0; i + 1 < split_values.size(); i++) {
          domain_intervals.emplace_back(split_values[i] / old_stroke_length,
                                        split_values[i + 1] /
                                          old_stroke_length);
        }

        auto split_strokes = std::vector<Stroke>();
        split_strokes.reserve(split_values.size() - 1);
        // TODO: There can be duplicates here.
        old_stroke.split(split_values, split_strokes);
        for (auto &split_stroke : split_strokes) {
          assert(split_stroke.size() > 1);
          split_stroke.compute_arclengths();
          graph.bvh_.add(std::move(split_stroke));
          graph.strokes2vid_.emplace_back();
        }
        split_strokes.clear();

        // @optimize Create a bump allocator for small allocations like these.
        auto t_vertices_unique =
          std::make_unique<size_t[]>(split_values.size());
        auto t_vertices =
          span<size_t>(t_vertices_unique.get(), split_values.size());
        t_vertices[0] = graph.hedges_[2 * si].origin_;
        if (trim_start) {
          cut_points.push_back({intersections[0].x_, t_vertices[0], false});
          // Update position (in case of trimming).
          graph.vertices_[t_vertices[0]].p_ = graph.strokes_[start_index].xy(0);
          graph.vertices_[t_vertices[0]].flags_ |=
            VertexRecord::CenterlineIntersection;
        }
        for (size_t i = 1; i < split_values.size() - 1; ++i) {
          // Create vertex at split point.
          t_vertices[i] = graph.vertices_.size();
          cut_points.push_back({intersections[i - (trim_start ? 0 : 1)].x_,
                                graph.vertices_.size(), false});
          auto &vertex_rec =
            graph.vertices_.emplace_back(graph.strokes_[start_index + i].xy(0));
          vertex_rec.flags_ |= VertexRecord::CenterlineIntersection;
        }
        t_vertices[split_values.size() - 1] = graph.hedges_[2 * si + 1].origin_;
        if (trim_end) {
          cut_points.push_back({intersections.back().x_,
                                t_vertices[split_values.size() - 1], false});
          // Update position (in case of trimming).
          graph.vertices_[t_vertices[split_values.size() - 1]].p_ =
            graph.strokes_.back().xy(Back);
          graph.vertices_[t_vertices[split_values.size() - 1]].flags_ |=
            VertexRecord::CenterlineIntersection;
        }

        size_t to_split_si = si;
        Float split_end = domain_intervals.back().second;
        Float split_v =
          (domain_intervals[0].second - domain_intervals[0].first) /
          (split_end - domain_intervals[0].first);
        for (size_t i = 0; i + 1 < split_values.size(); ++i) {
          const auto [new_edge_idx, new_twin_idx] = make_twins(graph.hedges_);
          auto &edge_rec = graph.hedges_[new_edge_idx];
          auto &twin_rec = graph.hedges_[new_twin_idx];
          edge_rec.origin_ = t_vertices[i];
          graph.vertices_[t_vertices[i]].hedge_ = new_edge_idx;
          twin_rec.origin_ = t_vertices[i + 1];
          if (i + 2 == split_values.size()) {
            graph.vertices_[t_vertices[i + 1]].hedge_ = new_twin_idx;
          }

          if (i == 0) {
            if (graph.hedges_[2 * si].prev_ == 2 * si + 1) {
              // Dangling edge
              edge_rec.prev_ = new_twin_idx;
              twin_rec.next_ = new_edge_idx;
            } else {
              assign_prev_next(graph, graph.hedges_[2 * si].prev_,
                               new_edge_idx);
              assign_prev_next(graph, new_twin_idx,
                               graph.hedges_[2 * si + 1].next_);

              if (graph.hedges_[2 * si].continuity_ != StrokeGraph::invalid) {
                graph.hedges_[new_edge_idx].continuity_ =
                  graph.hedges_[2 * si].continuity_;
                graph.hedges_[graph.hedges_[2 * si].continuity_].continuity_ =
                  new_edge_idx;
              }
            }
          } else {
            edge_rec.prev_ = new_edge_idx - 2;
            twin_rec.next_ = new_twin_idx - 2;

            graph.hedges_[new_edge_idx].continuity_ = new_edge_idx - 1;
            assert(graph.hedges_[new_edge_idx - 1].continuity_ == new_edge_idx);
          }
          if (i + 2 == split_values.size()) {
            if (graph.hedges_[2 * si].next_ == 2 * si + 1) {
              // Dangling edge
              edge_rec.next_ = new_twin_idx;
              twin_rec.prev_ = new_edge_idx;
            } else {
              assign_prev_next(graph, new_edge_idx,
                               graph.hedges_[2 * si].next_);
              assign_prev_next(graph, graph.hedges_[2 * si + 1].prev_,
                               new_twin_idx);

              if (graph.hedges_[2 * si + 1].continuity_ !=
                  StrokeGraph::invalid) {
                graph.hedges_[new_twin_idx].continuity_ =
                  graph.hedges_[2 * si + 1].continuity_;
                graph.hedges_[graph.hedges_[2 * si + 1].continuity_]
                  .continuity_ = new_twin_idx;
              }
            }
          } else {
            edge_rec.next_ = new_edge_idx + 2;
            twin_rec.prev_ = new_twin_idx + 2;

            graph.hedges_[new_twin_idx].continuity_ = new_edge_idx + 2;
          }

          edge_rec.face_ = graph.hedges_[2 * si].face_;
          twin_rec.face_ = graph.hedges_[2 * si + 1].face_;

          if (i + 2 < split_values.size()) {
            size_t new_stroke_index =
              (i == 0) ? graph.strokes2orig_.size() : to_split_si;
            split_stroke_mapping(graph, to_split_si, split_v, new_stroke_index,
                                 new_stroke_index + 1);
            to_split_si = new_stroke_index + 1;

            split_v =
              (domain_intervals[i + 1].second - domain_intervals[i + 1].first) /
              (split_end - domain_intervals[i + 1].first);
          }
        } // end of for (size_t i = 0; i + 1 < split_values.size(); ++i)
      } // end of if (split_values.size() == 2) { ... } else
    } // end of if (!cache->intersections_.empty())
  }

  // Deactivate strokes marked for deletion.
  for (size_t si = 0; si < graph.strokes_.size(); ++si) {
    if (graph.strokes_[si].size() == 0 &&
        (graph.hedge(2 * si).is_active() ||
         graph.hedge(2 * si + 1).is_active())) {
      delete_stroke_and_update_bvh(graph, si);
    }
  }

  intersections.clear();
  intersect_self(PolylineBVHLeaf(new_stroke, bounds(new_stroke)),
                 intersections);
  for (const auto &hit : intersections) {
    const auto vi = graph.vertices_.size();
    auto &vertex_rec = graph.vertices_.emplace_back(
      new_stroke.pos(hit.x_)); // We will set hedge_ later.
    vertex_rec.flags_ |= VertexRecord::CenterlineIntersection;
    cut_points.push_back({hit.x_, vi, true});
    cut_points.push_back({hit.y_, vi, true});
  }

  std::sort(cut_points.begin(), cut_points.end(),
            [](const CutPoint &a, const CutPoint &b) { return a.s < b.s; });
  // Remove duplicate and near-duplicate cut points.
  for (size_t i = 1; i < cut_points.size(); ++i) {
    if (cut_points[i].s - cut_points[i - 1].s < 1e-10) {
      cut_points[i - 1].s = -infinity;
    }
  }
  cut_points.erase(
    std::remove_if(cut_points.begin(), cut_points.end(),
                   [](const CutPoint &cut) { return cut.s == -infinity; }),
    cut_points.end());

  // Prepare to cut new_stroke.
  bool trim_start = false, trim_end = false;
  if (cache->trim_overshoots_) {
    trim_start = //
      (!cut_points.empty() &&
       cut_points[0].s <=
         std::max(0.5 * new_stroke.width_at(cut_points[0].s),
                  vertex_radius(graph.vertex(cut_points[0].vertex))));
    trim_end = //
      (!cut_points.empty() &&
       cut_points.back().s >=
         new_stroke.length() -
           std::max(0.5 * new_stroke.width_at(cut_points.back().s),
                    vertex_radius(graph.vertex(cut_points.back().vertex))));
  } else {
    trim_start = (!cut_points.empty() && cut_points[0].s <= 1e-10);
    trim_end = (!cut_points.empty() &&
                cut_points.back().s >= new_stroke.length() - 1e-10);
  }
  if (cut_points.size() == 1 && trim_start && trim_end) {
    trim_start = trim_end = false;
    const auto s = cut_points[0].s;
    if (s <= 1e-10) {
      trim_start = true;
    } else if (s >= new_stroke.length() - 1e-10) {
      trim_end = true;
    }
  }
  if (!trim_start) {
    const auto vi = graph.vertices_.size();
    cut_points.insert(cut_points.begin(), CutPoint{0.0, vi, false});
    graph.vertices_.emplace_back(new_stroke.xy(0)); // We will set hedge_ later.
    const auto endp = Endpoint(graph.orig_strokes_.size() - 1, true);
    graph.endpoint2vertex_.insert({endp, vi});
    graph.vertex2endpoint_.insert({vi, endp});
  }
  if (!trim_end) {
    const auto vi = graph.vertices_.size();
    cut_points.push_back({new_stroke.length(), vi, false});
    graph.vertices_.emplace_back(
      new_stroke.xy(Back)); // We will set hedge_ later.
    const auto endp = Endpoint(graph.orig_strokes_.size() - 1, false);
    graph.endpoint2vertex_.insert({endp, vi});
    graph.vertex2endpoint_.insert({vi, endp});
  }
  auto split_values_unique = std::make_unique<Float[]>(cut_points.size());
  auto split_values = span<Float>(split_values_unique.get(), cut_points.size());
  for (size_t i = 0; i < cut_points.size(); ++i) {
    split_values[i] = cut_points[i].s;
  }
  Float new_stroke_length = new_stroke.length();
  std::vector<std::pair<Float, Float>> domain_intervals;
  domain_intervals.reserve(split_values.size() - 1);
  for (size_t i = 0; i + 1 < split_values.size(); i++) {
    domain_intervals.emplace_back(split_values[i] / new_stroke_length,
                                  split_values[i + 1] / new_stroke_length);
  }

  auto split_strokes = std::vector<Stroke>();
  split_strokes.reserve(split_values.size() - 1);
  new_stroke.split(split_values, split_strokes);
  const auto start_index =
    graph.strokes_.size(); // First stroke index of splits.
  for (auto &split_stroke : split_strokes) {
    assert(split_stroke.size() > 1);
    graph.bvh_.add(std::move(split_stroke));
  }
  split_strokes.clear();

  domain_intervals[0].first = 0.0;
  domain_intervals.back().second = 1.0;
  auto &o2s = graph.orig2strokes_.emplace_back();
  for (size_t si = start_index; si < graph.strokes_.size(); ++si) {
    const auto &stroke = graph.strokes_[si];
    stroke.compute_arclengths();
    graph.strokes2orig_.emplace_back(
      std::vector<size_t>{graph.orig_strokes_.size() - 1});
    graph.strokes2vid_.emplace_back();

    StrokeMapping mapping;
    mapping.domain_arclens_.push_back(domain_intervals[si - start_index].first);
    mapping.domain_arclens_.push_back(
      domain_intervals[si - start_index].second);
    mapping.range_arclens_.push_back(0);
    mapping.range_arclens_.push_back(1);
    o2s.emplace_back(si, mapping);

    const auto [new_edge_idx, new_twin_idx] = make_twins(graph.hedges_);
    auto &edge_rec = graph.hedges_[new_edge_idx];
    auto &twin_rec = graph.hedges_[new_twin_idx];

    edge_rec.origin_ = cut_points[si - start_index].vertex;
    twin_rec.origin_ = cut_points[si - start_index + 1].vertex;

    if (si == start_index) {
      // Dangling edge
      edge_rec.prev_ = new_twin_idx;
      twin_rec.next_ = new_edge_idx;
    } else {
      edge_rec.prev_ = new_edge_idx - 2;
      twin_rec.next_ = new_twin_idx - 2;

      graph.hedges_[new_edge_idx].continuity_ = new_edge_idx - 1;
      assert(graph.hedges_[new_edge_idx - 1].continuity_ == new_edge_idx);
    }
    if (si == graph.strokes_.size() - 1) {
      // Dangling edge
      edge_rec.next_ = new_twin_idx;
      twin_rec.prev_ = new_edge_idx;
    } else {
      edge_rec.next_ = new_edge_idx + 2;
      twin_rec.prev_ = new_twin_idx + 2;

      graph.hedges_[new_twin_idx].continuity_ = new_edge_idx + 2;
    }

    if (si != start_index || trim_start) {
    } else {
      graph.vertices_[edge_rec.origin_].hedge_ = new_edge_idx;
    }
    if (si != graph.strokes_.size() - 1 || trim_end) {
    } else {
      graph.vertices_[twin_rec.origin_].hedge_ = new_twin_idx;
    }
    if (cut_points[si - start_index].is_self_intersection) {
      graph.vertices_[edge_rec.origin_].hedge_ = new_edge_idx;
    }

    if (graph.faces_.empty()) {
      graph.faces_.emplace_back();
    }
    edge_rec.face_ = graph.boundary_face_;
    twin_rec.face_ = graph.boundary_face_;
  }
  graph.snap_endpoints();
  // Orient vertex stars.
  for (size_t si = start_index; si < graph.strokes_.size(); ++si) {
    const auto new_edge_idx = 2 * si;
    const auto new_twin_idx = 2 * si + 1;
    auto &edge_rec = graph.hedges_[new_edge_idx];
    auto &twin_rec = graph.hedges_[new_twin_idx];
    if (si != start_index || trim_start) {
      insert_into_star(graph, new_edge_idx, edge_rec.origin_);
    }
    if (si != graph.strokes_.size() - 1 || trim_end) {
      insert_into_star(graph, new_twin_idx, twin_rec.origin_);
    }
  }

  // @optimize Make face finding incremental too.
  graph.faces_.clear();
  for (auto &edge_rec : graph.hedges_) {
    edge_rec.face_ = invalid;
  }

  // Snap existing vertices if they overlap.
  for (const auto ci : {(size_t)0, cut_points.size() - 1}) {
    const auto vi = cut_points[ci].vertex;
    auto vert = graph.vertex(vi);
    if (vert.is_active() && vert.is_dangling()) {
      const auto snap_vi = snap_if_overlapping(graph, vert, *cache);
      if (snap_vi != invalid) {
        cut_points[ci].vertex = snap_vi;
      }
    }
  }
  for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
    auto vert = graph.vertex(vi);
    if (vert.is_active() && vert.is_dangling()) {
      snap_if_overlapping(graph, vert, *cache);
    }
  }
  for (const auto ci : {(size_t)0, cut_points.size() - 1}) {
    const auto vert = graph.vertex(cut_points[ci].vertex);
    if (!vert || !vert.is_dangling()) {
      (ci == 0 ? info.head : info.tail).status = SnapInfo::Overlap;
    }
  }

  if (to_dissolve) {
    for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
      if (should_dissolve_vertex(graph.vertex(vi))) {
        dissolve_vertex(graph, vi);
      }
    }
  }

  // Determine and save junction types; assign unique vertex IDs
  for (size_t i = 0; i < graph.vertices_.size(); ++i) {
    auto &v = graph.vertices_[i];
    if (!v.ids_.empty())
      continue;
    // Vertex id
    v.ids_.emplace_back(
      StrokeGraph::VertexID({StrokeGraph::VertexID::Initialization, i}));
  }

  {
    graph.faces_.clear();
    for (size_t j = 0; j < graph.hedges_.size(); ++j) {
      graph.hedges_[j].face_ = StrokeGraph::invalid;
    }
    construct_faces(graph);
  }

  const auto [head_vertex_idx, tail_vertex_idx] =
    find_original_endpoint_vertices(graph, graph.orig_strokes_.size() - 1);
  if (head_vertex_idx != invalid &&
      cache->enforced_connections_
        .associated_connections(
          StrokeTime((int)graph.orig_strokes_.size() - 1, 0.0))
        .empty()) {
    graph.vertices_[head_vertex_idx].flags_ |= VertexRecord::OriginallyDangling;
  }
  if (tail_vertex_idx != invalid &&
      cache->enforced_connections_
        .associated_connections(
          StrokeTime((int)graph.orig_strokes_.size() - 1, 1.0))
        .empty()) {
    graph.vertices_[tail_vertex_idx].flags_ |= VertexRecord::OriginallyDangling;
  }

#ifndef NDEBUG
  check_consistency(graph);
  check_continuity(graph);
  check_mapping(graph);
  // check_intersection(graph);
#endif

  return info;
}

StrokeSnapInfo
add_stroke_incremental(StrokeGraph &graph, IncrementalCache *const cache,
                       const Stroke &new_stroke, bool to_dissolve) {
  StrokeSnapInfo info =
    add_stroke_incremental_topological(graph, cache, new_stroke, to_dissolve);

  size_t endpoint_vertex_ids[2] = {invalid, invalid};
  auto &head_vertex_idx = endpoint_vertex_ids[0];
  auto &tail_vertex_idx = endpoint_vertex_ids[1];
  std::tie(head_vertex_idx, tail_vertex_idx) =
    find_original_endpoint_vertices(graph, graph.orig_strokes_.size() - 1);

  // Enforce connections.
  // This can be needed if the original strokes overlap, but because of
  // deformation, they no longer overlap in the current graph.
  for (const auto &conn : cache->enforced_connections_.connections()) {
    const auto current_orig_si = (int)graph.orig_strokes_.size() - 1;
    StrokeTime cand1 = conn.first;
    StrokeTime cand2 = conn.second;
    if (!(cand1.first == current_orig_si || cand2.first == current_orig_si) ||
        cand1.first > current_orig_si || cand2.first > current_orig_si) {
      continue; // Already applied or can't be applied yet.
    }
    const auto orig_si2 = cand2.first;
    if (!convert_orig2strokes(graph, cand1) ||
        !convert_orig2strokes(graph, cand2)) {
      continue;
    }
    if (cand1.second != 0.0 && cand1.second != 1.0) {
      continue; // Can't snap from stroke interior.
    }
    const auto v1 =
      endpoint_to_vertex(graph, Endpoint(cand1.first, cand1.second == 0));
    if (cand2.second == 0 || cand2.second == 1) { // Vertex-vertex.
      const auto v2 =
        endpoint_to_vertex(graph, Endpoint(cand2.first, cand2.second == 0));
      if (is_snap_valid(graph, v1.index_, v2.index_)) {
        const auto bvh = graph.bvh_.polyline_bvh();
        if (line_of_sight(v1.pos(), v2.pos(), bvh) &&
            graph.snap_vertices(v1.index_, v2.index_)) {
          if (v1.index_ == head_vertex_idx) {
            info.head.status = SnapInfo::ReservedPos;
          } else if (v1.index_ == tail_vertex_idx) {
            info.tail.status = SnapInfo::ReservedPos;
          }
          if (v2.index_ == head_vertex_idx) {
            info.head.status = SnapInfo::ReservedPos;
          } else if (v2.index_ == tail_vertex_idx) {
            info.tail.status = SnapInfo::ReservedPos;
          }
        }
      }
    } else { // Vertex-stroke interior
      auto already_connected = false;
      const auto he = v1.hedge();
      auto it = he;
      do {
        const auto &s2o = graph.strokes2orig_[it.stroke_idx()];
        if (std::find(s2o.begin(), s2o.end(), orig_si2) != s2o.end()) {
          // We could also look up at the endpoint position specifically,
          // however that might cause spurious connections...
          already_connected = true;
          break;
        }
        it = it.twin().next();
      } while (it != he);
      if (already_connected) {
        continue;
      }
      const auto maybe_cand2 = reprojected_vertex_stroke_candidate(
        v1, /*orig_si=*/orig_si2,
        /*last_orig_si=*/(int)graph.orig_strokes_.size() - 1, *cache);
      if (maybe_cand2.has_value()) {
        const auto [si, norm_arclen] = maybe_cand2.value();
        graph.snap_endpoint_to_edge(v1, si,
                                    norm_arclen * graph.strokes_[si].length(),
                                    graph.snapping_method_type_);
        if (v1.index_ == head_vertex_idx) {
          info.head.status = SnapInfo::ReservedPos;
        } else if (v1.index_ == tail_vertex_idx) {
          info.tail.status = SnapInfo::ReservedPos;
        }
      }
    }
  }

  if (cache->accept_threshold_ < 1.0) {
    // Look for snap points using the classifier.
    for (const auto ci : {0, 1}) {
      if (endpoint_vertex_ids[ci] == invalid ||
          (ci == 0 && info.head.status != SnapInfo::Dangling) ||
          (ci != 0 && info.tail.status != SnapInfo::Dangling)) {
        continue;
      }
      auto &snap_info = (ci == 0 ? info.head : info.tail);
      snap_info = snap_with_classifier(
        graph, endpoint_vertex_ids[ci], cache, prediction_feature_type,
        cache->accept_threshold_, &cache->predictions);
      if (snap_info.is_connected()) {
        endpoint_vertex_ids[ci] = snap_info.other_vertex;
        if (to_dissolve &&
            should_dissolve_vertex(graph.vertex(snap_info.other_vertex))) {
          dissolve_vertex(graph, snap_info.other_vertex);
        }
        if (ci == 0 && snap_info.other_vertex == tail_vertex_idx) {
          assert(tail_vertex_idx != invalid);
          info.tail =
            snap_info; // We form a closed stroke and head info is same as tail.
          tail_vertex_idx = snap_info.other_vertex;
        }
      }
    }
    auto did_snap = false;
    auto predictions = std::vector<ClassifierPrediction>();
    do {
      did_snap = false;
      for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
        const auto vert = graph.vertex(vi);
        if (vert.is_active()) {
          predictions.clear();
          const auto snap_info =
            snap_with_classifier(graph, vi, cache, prediction_feature_type,
                                 cache->accept_threshold_, &predictions);
          if (snap_info.is_connected()) {
            did_snap = true;
            if (to_dissolve &&
                should_dissolve_vertex(graph.vertex(snap_info.other_vertex))) {
              dissolve_vertex(graph, snap_info.other_vertex);
            }
          }
          for (const auto &pred : predictions) {
            if (pred.prob > 0.4) {
              cache->predictions.emplace_back(pred);
            }
          }
        }
      }
    } while (did_snap);

    for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
      const auto v = graph.vertex(vi);
      if (to_dissolve && should_dissolve_vertex(v)) {
        dissolve_vertex(graph, v.index_);
      }
    }
  }

  // @optimize Make face finding incremental too.
  {
    graph.faces_.clear();
    for (size_t j = 0; j < graph.hedges_.size(); ++j) {
      graph.hedges_[j].face_ = StrokeGraph::invalid;
    }
    construct_faces(graph);
  }

#ifndef NDEBUG
  check_consistency(graph);
  check_continuity(graph);
  check_mapping(graph);
  // check_intersection(graph);
#endif

  info.predictions = std::move(cache->predictions);
  std::sort(info.predictions.begin(), info.predictions.end());
  cache->predictions.clear();

  return info;
}

void snap_candidates(const StrokeGraph &graph, const size_t vi,
                     const IncrementalCache *const cache,
                     const Connections &intersection_set,
                     FeatureType feature_type,
                     std::vector<std::vector<SnapInfo>> &candidate_sets,
                     size_t cand_count, bool to_expand, bool to_intersect,
                     bool train_time) {
  if (intersection_set.empty() && to_intersect)
    return;

  const auto vertex = graph.vertex(vi);
  assert(vertex && "vertex was invalidated");

  if (!vertex.is_originally_dangling()) {
    return;
  }

  // To avoid assertion triggered
  if (to_expand && has_continuous_edge(vertex))
    return;

  {
    auto cand2 = ALLOCA_SPAN(StrokeTime, cand_count);
    find_filtered_end_end_candidates(vertex, *cache, intersection_set,
                                     feature_type, cand2, train_time);

    if (!to_expand) {
      candidate_sets.emplace_back();
      for (size_t i = 0; i < cand2.size(); ++i) {
        SnapInfo snap_info;
        const auto &other_si = cand2[i].first;
        snap_info.status = SnapInfo::PredictionDelay;
        snap_info.max_prob = -1;
        snap_info.prediction_type = JunctionType::R;
        snap_info.stroke_pos1 = StrokeTime((int)vertex.hedge().stroke_idx(),
                                           (vertex.hedge().forward()) ? 0 : 1);
        snap_info.stroke_pos2 = StrokeTime(other_si, cand2[i].second);
        candidate_sets.back().emplace_back(snap_info);
      }
    } else {
      // Expand high-valence junction candidates
      auto n_evaluations = std::vector<Index>(cand2.size());
      size_t n_evaluations_total = 0;
      auto c1 = vertex_to_endpoint(vertex);
      for (size_t i = 0; i < cand2.size(); ++i) {
        const auto &c2 = cand2[i];
        size_t n = n_classifier_evaluations_vv(
          graph, c1, Endpoint{(size_t)c2.first, c2.second == 0.0});
        n_evaluations[i] = n;
        n_evaluations_total += n;
      }

      auto expanded_cand1 =
        std::vector<Endpoint>(n_evaluations_total, Endpoint(0));
      auto expanded_cand2 = std::vector<StrokeTime>(n_evaluations_total);
      for (size_t i = 0, k = 0; i < cand2.size(); ++i) {
        size_t n = n_evaluations[i];
        const auto &c2 = cand2[i];
        expand_candidates_vv(graph, c1,
                             Endpoint{(size_t)c2.first, c2.second == 0.0},
                             {&expanded_cand1[k], n}, {&expanded_cand2[k], n});

        candidate_sets.emplace_back();
        for (size_t j = k; j < k + n; ++j) {
          SnapInfo snap_info;
          const auto &other_si = expanded_cand2[j].first;
          snap_info.status = SnapInfo::PredictionDelay;
          snap_info.max_prob = -1;
          snap_info.prediction_type = JunctionType::R;
          snap_info.stroke_pos1 =
            StrokeTime{(int)expanded_cand1[j].stroke_idx(),
                       (expanded_cand1[j].is_head()) ? 0 : 1};
          snap_info.stroke_pos2 =
            StrokeTime{(int)other_si, expanded_cand2[j].second};
          candidate_sets.back().emplace_back(snap_info);
        }

        k += n;
      }
    }
  }

  // Note we can snap *to* a vertex with a continuous edge, but we cannot snap
  // *from* such a vertex. (The difference is in the position of the resulting
  // vertex after snapping.) We also shouldn't snap from vertices that weren't
  // originally dangling.
  if (has_continuous_edge(vertex))
    return;

  {
    // Try finding a T-junction.
    std::vector<StrokeTime> cand2;
    tjunc_classifier_candidates(vertex, std::unordered_set<std::size_t>(),
                                (int)cand_count, cand2);

    if (cand2.empty())
      return;

    if (!to_expand) {
      candidate_sets.emplace_back();
      for (size_t i = 0; i < cand2.size(); ++i) {
        SnapInfo snap_info;
        const auto [other_si, other_arclen] = cand2[i];
        snap_info.status = SnapInfo::PredictionDelay;
        snap_info.max_prob = -1;
        snap_info.prediction_type = JunctionType::T;
        snap_info.stroke_pos1 = StrokeTime((int)vertex.hedge().stroke_idx(),
                                           (vertex.hedge().forward()) ? 0 : 1);
        snap_info.stroke_pos2 = StrokeTime(other_si, other_arclen);
        candidate_sets.back().emplace_back(snap_info);
      }
    } else {
      // Expand high-valence junction candidates
      size_t n_evaluations_total = 0;
      auto c1 = vertex_to_endpoint(vertex);
      n_evaluations_total = n_classifier_evaluations_vs(graph, c1);

      auto expanded_cand1 =
        std::vector<Endpoint>(n_evaluations_total, Endpoint(0));
      expand_candidates_vs(graph, c1,
                           {&expanded_cand1[0], n_evaluations_total});

      for (size_t i = 0; i < cand2.size(); ++i) {
        candidate_sets.emplace_back();
        const auto [other_si, other_arclen] = cand2[i];
        for (size_t j = 0; j < expanded_cand1.size(); ++j) {
          SnapInfo snap_info;
          snap_info.status = SnapInfo::PredictionDelay;
          snap_info.max_prob = -1;
          snap_info.prediction_type = JunctionType::T;
          snap_info.stroke_pos1 =
            StrokeTime{(int)expanded_cand1[j].stroke_idx(),
                       (expanded_cand1[j].is_head()) ? 0 : 1};
          snap_info.stroke_pos2 = StrokeTime{(int)other_si, other_arclen};
          candidate_sets.back().emplace_back(snap_info);
        }
      }
    }
  }
}

void snap_candidates(const StrokeGraph &graph, size_t vi,
                     const IncrementalCache *const cache,
                     const Connections &intersection_set,
                     FeatureType feature_type,
                     std::vector<SnapInfo> &candidates, size_t cand_count,
                     bool to_expand, bool train_time) {
  std::vector<std::vector<SnapInfo>> candidate_sets;
  snap_candidates(graph, vi, cache, intersection_set, feature_type,
                  candidate_sets, cand_count, to_expand, true, train_time);
  for (auto &cand : candidate_sets) {
    for (auto &snap : cand) {
      candidates.emplace_back(std::move(snap));
    }
  }
}

void snap_candidates(const StrokeGraph &graph, size_t vi,
                     const IncrementalCache *const cache,
                     FeatureType feature_type,
                     std::vector<SnapInfo> &candidates, size_t cand_count,
                     bool to_expand, bool train_time) {
  std::vector<std::vector<SnapInfo>> candidate_sets;
  snap_candidates(graph, vi, cache, Connections(), feature_type, candidate_sets,
                  cand_count, to_expand, false, train_time);
  for (auto &cand : candidate_sets) {
    for (auto &snap : cand) {
      candidates.emplace_back(std::move(snap));
    }
  }
}

void snap_corner_candidates(const StrokeGraph &graph, size_t vi,
                            const IncrementalCache *const cache,
                            FeatureType /*feature_type*/,
                            std::vector<SnapInfo> &candidates,
                            size_t cand_count) {
  const auto vertex = graph.vertex(vi);
  assert(vertex && "vertex was invalidated");

  // Note we can snap *to* a vertex with a continuous edge, but we cannot snap
  // *from* such a vertex. (The difference is in the position of the resulting
  // vertex after snapping.) We also shouldn't snap from vertices that weren't
  // originally dangling.
  if (has_continuous_edge(vertex) || !vertex.is_originally_dangling()) {
    return;
  }

  auto orig_endp = hedge_to_endpoint(vertex.hedge()).as_pair();
  const auto ok = convert_strokes2orig(graph, orig_endp);
  force_assert(ok && "couldn't map from strokes to orig");
  const auto last_orig_si = (int)graph.orig_strokes_.size() - 1;
  force_assert(last_orig_si >= 0);
  const auto enforced_connections =
    cache->enforced_connections_.associated_connections(orig_endp);
  for (const auto &conn : enforced_connections) {
    const auto other_orig_si = int{conn.second.first};
    if (other_orig_si > last_orig_si) {
      // The endpoint needs to connect to a later stroke first.
      return;
    }
  }

  // Try finding a T-junction.
  std::vector<StrokeTime> candidates_buf;
  candidates_buf.resize(cand_count, StrokeTime{0, false});

  // Disable too close to stroke endpoint check
  std::vector<StrokeTime> cand2;
  tjunc_classifier_candidates(vertex, std::unordered_set<std::size_t>(),
                              (int)cand_count, cand2, false);

  if (cand2.empty())
    return;

  for (size_t i = 0; i < cand2.size(); ++i) {
    SnapInfo snap_info;
    const auto [other_si, other_arclen] = cand2[i];
    snap_info.status = SnapInfo::PredictionDelay;
    snap_info.max_prob = -1;
    snap_info.prediction_type = JunctionType::T;
    snap_info.stroke_pos1 = StrokeTime((int)vertex.hedge().stroke_idx(),
                                       (vertex.hedge().forward()) ? 0 : 1);
    snap_info.stroke_pos2 = StrokeTime(other_si, other_arclen);
    candidates.emplace_back(snap_info);
  }
}

void all_snap_candidates(const StrokeGraph &graph, size_t orig_si,
                         const IncrementalCache *const cache,
                         const Connections &intersection_set,
                         FeatureType feature_type,
                         std::vector<std::vector<Junction>> &junction_sets,
                         size_t cand_count, bool train_time) {
  if (intersection_set.empty())
    return;

  // Find all vertices on this original stroke
  std::vector<std::vector<SnapInfo>> candidate_sets;
  std::unordered_set<size_t> new_ends;
  for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
    auto v = graph.vertex(vi);
    if (!v.is_valid() || !v.is_active())
      continue;
    const auto he = v.hedge();
    auto it = he;
    do {
      StrokeTime st((int)it.stroke_idx(), (Float)!it.forward());
      auto ok = convert_strokes2orig(graph, st);
      force_assert(ok && "couldn't map from strokes to orig");

      if (st.first == orig_si) {
        new_ends.emplace(vi);
        break;
      }
      it = it.twin().next();
    } while (it != he);
  }
  std::unordered_set<size_t> new_vv;
  for (auto vi : new_ends) {
    const auto v = graph.vertex(vi);
    if (!v.is_valid() || !v.is_active())
      continue;
    snap_candidates(graph, v.index_, cache, intersection_set, feature_type,
                    candidate_sets, cand_count, true, true, train_time);
    new_vv.emplace(v.index_);
  }

  junction_sets.reserve(candidate_sets.size());
  bool seen_self_connection = false;
  for (const auto &candidates : candidate_sets) {
    junction_sets.emplace_back();
    junction_sets.back().reserve(candidates.size());
    for (const auto &snap : candidates) {
      Junction junc(
        std::vector<std::pair<int, double>>{
          std::make_pair((int)snap.stroke_pos2.first, snap.stroke_pos2.second),
          std::make_pair((int)snap.stroke_pos1.first, snap.stroke_pos1.second)},
        snap.prediction_type, false, 0.5, true);

      const auto ok = stroke_to_original_stroke_indexing(graph, junc);
      force_assert(ok && "couldn't map from strokes to orig");

      // Avoid two copies of self regular junctions
      if (junc.type == JunctionType::Type::R &&
          junc.points[0].first == junc.points[1].first) {
        if (seen_self_connection)
          continue;
        else
          seen_self_connection = true;
      }

      // Write the initial junction distance based on the original strokes. This
      // value should never be changed.
      junc.orig_dist = junction_distance_init(graph, junc);

      junction_sets.back().emplace_back(std::move(junc));
    }
  }

  if (graph.orig_strokes_.empty())
    return;

  // From existing vertices to the new stroke
  std::vector<std::vector<SnapInfo>> to_new_candidate_sets;
  size_t new_orig_si = graph.orig_strokes_.size() - 1;

  // Filter out allowed/enforced set given the newest original stroke index
  Connections new_intersection_set =
    intersection_set.subset_involving(new_orig_si);

  for (size_t i = 0;
       !new_intersection_set.empty() && i < graph.vertices_.size(); ++i) {
    StrokeGraph::VertexView v = graph.vertex(i);
    if (!v.is_valid() || !v.is_active() || new_vv.count(i))
      continue;

    to_new_candidate_sets.clear();
    snap_candidates(graph, v.index_, cache, new_intersection_set, feature_type,
                    to_new_candidate_sets, cand_count, true, true, train_time);
    for (const auto &to_new_stroke_candidates : to_new_candidate_sets) {
      if (junction_sets.empty() || !junction_sets.back().empty())
        junction_sets.emplace_back();
      junction_sets.back().reserve(to_new_stroke_candidates.size());
      for (const auto &can : to_new_stroke_candidates) {
        // Regular ones are duplicates
        if (can.prediction_type == JunctionType::Type::R)
          continue;

        Junction junc(
          std::vector<std::pair<int, double>>{
            std::make_pair((int)can.stroke_pos2.first, can.stroke_pos2.second),
            std::make_pair((int)can.stroke_pos1.first, can.stroke_pos1.second)},
          can.prediction_type, false, 0.5, true);

        const auto ok = stroke_to_original_stroke_indexing(graph, junc);
        force_assert(ok && "couldn't map from strokes to orig");

        force_assert(junc.points[1].first != orig_si);
        if (junc.points[0].first == orig_si) {
          // Write the initial junction distance based on the original strokes.
          // This value should never be changed.
          junc.orig_dist = junction_distance_init(graph, junc);
          junction_sets.back().emplace_back(std::move(junc));
        }
      }

      if (junction_sets.back().empty())
        junction_sets.pop_back();
    }
  }
}

std::vector<ClassifierPrediction>
finalize_incremental(StrokeGraph &graph, IncrementalCache *const cache,
                     bool to_dissolve) {
  // @optimize Make face finding incremental too.
  graph.faces_.clear();
  for (auto &edge_rec : graph.hedges_) {
    edge_rec.face_ = invalid;
  }

  auto predictions = std::vector<ClassifierPrediction>();
  auto did_snap = false;
  do {
    did_snap = false;
    for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
      const auto vert = graph.vertex(vi);
      if (vert.is_active()) {
        const auto snap_info = snap_with_classifier(
          graph, vi, cache, prediction_feature_type, 0.5, &predictions);
        if (snap_info.is_connected()) {
          did_snap = true;
          if (to_dissolve &&
              should_dissolve_vertex(graph.vertex(snap_info.other_vertex))) {
            dissolve_vertex(graph, snap_info.other_vertex);
          }
#ifndef NDEBUG
          check_mapping(graph);
#endif
        }
      }
    }
  } while (did_snap);

  // @optimize Make face finding incremental too.
  construct_faces(graph);

  return predictions;
}

} // namespace sketching
