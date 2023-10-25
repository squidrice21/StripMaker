#include "incremental_util.h"

#include "classifier.h"
#include "detail/alloca.h"
#include "features/junction_features_impl.h"
#include "fitting.h"
#include "force_assert.h"
#include "incremental_param.h"
#include "incremental_region_util.h"
#include "incremental_util.h"
#include "intersect.h"
#include "junction_type.h"
#include "stroke_graph_extra.h"
#include "stroke_view.h"

#include <fstream>
#include <iomanip>

#if defined(_MSC_VER) && !defined(__clang__)
#pragma warning(disable : 4267) // conversion, possible loss of data
#endif

namespace sketching {
JunctionType::Type recosider_type = JunctionType::T;
bool include_ts = true;
bool include_corners = false;

Float corner_projection_ratio = 4;
namespace {

Float get_final_prediction(const Junction &junc,
                           const std::vector<Junction> &final_predictions) {
  Float prob = -1;

  // Match the index and junction type
  Float time_diff = 1;
  Junction closest_junc(JunctionType::X);
  Junction sorted_junc = junc;
  sorted_junc.sort_entries();
  for (const auto &j : final_predictions) {
    if (j.type != junc.type)
      continue;

    Junction sorted_j = j;
    sorted_j.sort_entries();
    if ((sorted_junc.type == JunctionType::R &&
         sorted_junc.points[0] == sorted_j.points[0] &&
         sorted_junc.points[1] == sorted_j.points[1])) {
      closest_junc = sorted_j;
      break;
    } else if (sorted_junc.type == JunctionType::T &&
               sorted_junc.points[0].first == sorted_j.points[0].first &&
               sorted_junc.points[1] == sorted_j.points[1]) {
      Float diff =
        std::abs(sorted_junc.points[0].second - sorted_j.points[0].second);
      if (diff < time_diff) {
        time_diff = diff;
        closest_junc = sorted_j;
      }
      if (time_diff == 0)
        break;
    }
  }
  if (closest_junc.type != JunctionType::X)
    return closest_junc.probability;

  return prob;
}

bool should_delay(const StrokeGraph &final_stroke_graph, const Junction &junc,
                  size_t cur_stroke_step) {
  // Check the end containing original stroke if any snapped vertices is close
  // to this end
  auto is_vertex_close = [&final_stroke_graph](const StrokeTime &end,
                                               size_t exclude_si) {
    StrokeTime s_end = end;
    auto ok = convert_orig2strokes(final_stroke_graph, s_end);
    // This original stroke is entirely delete in the final graph
    if (!ok)
      return false;
    if (s_end.second != 0.0 && s_end.second != 1.0)
      return false;
    const auto v = endpoint_to_vertex(
      final_stroke_graph, Endpoint(s_end.first, s_end.second == 0.0));
    assert(v.is_valid() && v.is_active());
    if (v.is_dangling())
      return false;

    // Check if the stroke is snapped to a stroke that is not the one in this
    // proposed candidate
    bool seen_others = false;
    const auto he = v.hedge();
    auto it = he;
    do {
      size_t orig_si = it.orig_stroke_idx();
      if (orig_si > exclude_si && orig_si > end.first) {
        seen_others = true;
        break;
      }
      it = it.twin().next();
    } while (it != he);

    return seen_others;
  };

  if (junc.type == JunctionType::T)
    return is_vertex_close(junc.points[1], cur_stroke_step);
  else if (junc.type == JunctionType::R)
    return is_vertex_close(junc.points[0], cur_stroke_step) ||
           is_vertex_close(junc.points[1], cur_stroke_step);
  std::abort();
}

bool is_visible_in_final(span<const Stroke> strokes, const Vec2 p,
                         const Vec2 q) {
  Float ss, t;
  for (const auto &s : strokes) {
    if (intersect_segment_stroke_exclusive(p, q, {s, bounds(s)}, ss, t)) {
      return false;
    }
  }

  return true;
}

bool is_final_prediction_acceptable(
  const Junction &junc, const std::vector<Junction> &final_predictions) {
  Float prob = get_final_prediction(junc, final_predictions);
  if (prob >= 0 && prob < min_final_prob)
    return false;

  return true;
}

void get_stroke_idx_round(
  const StrokeGraph &graph, const std::pair<int, double> end,
  std::vector<std::unordered_map<size_t, int>> &junc_si) {
  size_t orig_si = end.first;
  Float orig_s_time = end.second;
  Float dist = 1;
  std::vector<Endpoint> best_ends;
  for (const auto &si_mapping : graph.orig2strokes_[orig_si]) {
    Float d1 = std::abs(orig_s_time - si_mapping.second.domain_arclens_[0]);
    Float d2 = std::abs(orig_s_time - si_mapping.second.domain_arclens_[1]);
    if (d1 < dist) {
      best_ends.clear();
      dist = d1;
      best_ends.emplace_back(
        Endpoint{si_mapping.first, si_mapping.second.range_arclens_[0] < 0.5});
    } else if (d1 == dist && d1 < 1e-10) { // So we only match both continuous
                                           // edges if at the vertex
      best_ends.emplace_back(
        Endpoint{si_mapping.first, si_mapping.second.range_arclens_[0] < 0.5});
    }

    if (d2 < dist) {
      best_ends.clear();
      dist = d2;
      best_ends.emplace_back(
        Endpoint{si_mapping.first, si_mapping.second.range_arclens_[1] < 0.5});
    } else if (d2 == dist && d2 < 1e-10) { // So we only match both continuous
                                           // edges if at the vertex
      best_ends.emplace_back(
        Endpoint{si_mapping.first, si_mapping.second.range_arclens_[1] < 0.5});
    }
  }
  for (const auto &best_end : best_ends)
    junc_si.back().emplace(best_end.stroke_idx() * 10 + best_end.is_tail(), 0);
}

bool is_high_valence_junction_on_cycle(const StrokeGraph &graph,
                                       const span<const Junction> candidates,
                                       size_t fi,
                                       const StrokeGraph::VertexView v,
                                       const std::string &vid_str,
                                       bool boundary_only = true) {
  // 1. Find the junction to get the adjacent strokes
  int cand_idx = -1;
  Junction strokes_junc = candidates.front();
  std::vector<std::unordered_map<size_t, int>> junc_si;
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (candidates[i].repr != vid_str)
      continue;
    strokes_junc = candidates[i];
    cand_idx = (int)i;

    junc_si.emplace_back();
    get_stroke_idx_round(graph, strokes_junc.points[0], junc_si);
    junc_si.emplace_back();
    get_stroke_idx_round(graph, strokes_junc.points[1], junc_si);

    break;
  }

  // It's possible if a junction is removed since one of its two corresponding
  // original stroke is removed
  if (cand_idx < 0)
    return false;

  assert(junc_si.size() == 2);

  // 2. Circulate the face cycle
  const auto &f = graph.face(fi);
  std::vector<std::vector<size_t>> face_edge_ends;
  face_edge_ends.reserve(f.cycles().size());
  for (const auto hi : f.cycles()) {
    SPDLOG_DEBUG("============");
    size_t iterations = 0;
    face_edge_ends.emplace_back();
    auto &s_indices = face_edge_ends.back();
    const auto he = f.graph_->hedge(hi);
    auto it = he;
    do {
      size_t si = it.stroke_idx();
      size_t si_hash = si * 10 + !it.forward();
      size_t si_hash_dest = si * 10 + !it.twin().forward();
      SPDLOG_DEBUG("In cycle: {} -> {}", si_hash, si_hash_dest);

      // Ignore the strokes within the region
      if (boundary_only && !s_indices.empty() && s_indices.back() == si_hash) {
        s_indices.pop_back();
        s_indices.pop_back();
      } else {
        s_indices.emplace_back(si_hash);
        s_indices.emplace_back(si_hash_dest);
      }

      it = it.next();
      iterations++;
      assert(iterations < 1024 && "likely infinite loop found");
    } while (it != he);
    if (face_edge_ends.back().empty())
      face_edge_ends.pop_back();
  }

  // 3. See if the two sides are adjacent in the cycle
  size_t passing_count = 0;
  // Use seen_pairs to handle the special case of a single edge connected to
  // itself forming a cycle. In this case, we would check the same edge pair
  // twice.
  std::set<std::pair<size_t, size_t>> seen_pairs;
  for (const auto &c : face_edge_ends) {
    seen_pairs.clear();
    for (size_t i = 0; i < c.size(); ++i) {
      auto si_hash = c[i];
      auto si_hash_next = c[(i + 1) % c.size()];
      auto p = std::make_pair(std::min(i, (i + 1) % c.size()),
                              std::max(i, (i + 1) % c.size()));
      if (seen_pairs.count(p))
        continue;
      seen_pairs.emplace(p);
      if ((junc_si[0].count(si_hash) && junc_si[1].count(si_hash_next)) ||
          (junc_si[1].count(si_hash) && junc_si[0].count(si_hash_next)))
        passing_count++;
    }
  }

  if (boundary_only) {
    return passing_count == 1;
  } else {
    return passing_count > 0;
  }
}

// When would we have a false from this?
bool can_expand(const StrokeGraph &graph, const Endpoint e1,
                bool continuity_check = true) {
  auto n_1st_endp_entries = Index(0);
  const auto he = endpoint_to_hedge(graph, e1);
  auto it = he;
  do {
    if (continuity_check && it.continuity() != StrokeGraph::invalid)
      return false;
    n_1st_endp_entries++;
    it = it.twin().next();
    if (n_1st_endp_entries >= 1024)
      return false;
  } while (it != he);
  return true;
}

void update_high_valence_prediction(StrokeGraph &varying_graph,
                                    const std::vector<size_t> &combine_he_idx,
                                    Junction &hv_junc) {
  auto vid = to_v_idx(varying_graph, hv_junc.points[0]);

  assert(varying_graph.vertex(vid).valence() > 1);

  // Mark other strokes then recompute the predictions
  std::unordered_set<int> junc_sid{hv_junc.points[0].first,
                                   hv_junc.points[1].first};
  const auto he = varying_graph.vertex(vid).hedge();
  auto it = he;
  do {
    auto he_sid = it.stroke_idx();
    StrokeTime end((int)he_sid, (it.forward()) ? 0 : 1);
    const auto ok = convert_strokes2orig(varying_graph, end);
    force_assert(ok && "couldn't map from strokes to orig");

    // Find the chain in the init graph
    StrokeTime end_init = end;
    const auto ok2 = convert_orig2strokes(varying_graph, end_init);
    force_assert(ok2 && "couldn't map from orig to strokes");

    std::vector<std::pair<size_t, bool>> orig_indices;
    size_t he_init = end_init.first * 2;
    if (end_init.second == 0.0 || end_init.second == 1.0) {
      if (end_init.second == 1.0)
        he_init++;
      get_forward_chain_original_indices(varying_graph, he_init, orig_indices);
    } else {
      get_forward_chain_original_indices(varying_graph, he_init, orig_indices);
      get_forward_chain_original_indices(varying_graph, he_init + 1,
                                         orig_indices);
    }

    for (const auto chain_orig_sid : orig_indices)
      if (!junc_sid.count(chain_orig_sid.first)) {
        varying_graph.orig_bvh_->masked_nodes.emplace_back(
          chain_orig_sid.first);
      }

    it = it.twin().next();
  } while (it != he);

  // Predict on the same graph used for the initial predictions. Only difference
  // is the masked original strokes.
  std::vector<Junction> junc{hv_junc};
  assert(hv_junc.corner_type == JunctionType::T ||
         hv_junc.corner_type == JunctionType::R);
  // junc.front().type = hv_junc.corner_type;
  // TODO: support the high-valence R
  junc.front().type = JunctionType::T;
  if (junc.front().to_reconsider)
    junc.front().type = recosider_type;
  junc.front().fea.type_ = FeatureVector::EndStroke;
  varying_graph.orig_bvh_->combine_he_indices = combine_he_idx;
  update_disconnected_junction_predictions(
    varying_graph, junc, prediction_feature_type, false, true);

  hv_junc.fea = junc.front().fea;
  hv_junc.probability = junc.front().probability;

  // Reset mask
  varying_graph.orig_bvh_->masked_nodes.clear();
  varying_graph.orig_bvh_->combine_he_indices.clear();
}

bool graph_visible(const StrokeGraph &stroke_graph, const Junction &in_junc) {
  Junction junc = in_junc;
  if (!original_stroke_to_stroke_indexing(stroke_graph, junc)) {
    return false;
  }
  Vec2 p =
    stroke_graph.strokes_[junc.points[0].first].pos_norm(junc.points[0].second);
  Vec2 q =
    stroke_graph.strokes_[junc.points[1].first].pos_norm(junc.points[1].second);
  return line_of_sight(p, q, stroke_graph.strokes_,
                       stroke_graph.bvh_.centerline_bbs());
}

bool envelope_visible(const StrokeGraph &stroke_graph, const Junction &junc) {
  Vec2 p = stroke_graph.orig_strokes_[junc.points[0].first].pos_norm(
    junc.points[0].second);
  Vec2 q = stroke_graph.orig_strokes_[junc.points[1].first].pos_norm(
    junc.points[1].second);

  size_t orig_sid = junc.points[1].first;

  bool end_inside = false;
  for (size_t i = 0; i < stroke_graph.orig_strokes_.size(); ++i) {
    if (i == junc.points[0].first || i == junc.points[1].first)
      continue;
    if (stroke_graph.orig_bvh_->nodes[i].bb.contains(p)) {
      Float dist = signed_dist_stroke(p, stroke_graph.orig_strokes_[i]);
      if (dist <= 0) {
        end_inside = true;
        break;
      }
    }
  }
  if (!end_inside)
    return true;

  // Take a step along the connection line direction
  const auto &s1 = stroke_graph.orig_strokes_[junc.points[0].first];
  s1.ensure_arclengths();
  Float end_width = s1.width_at(s1.length() * junc.points[0].second);
  Vec2 p2 = end_width * (q - p).normalized() + p;

  if (!line_of_sight(p2, q, *stroke_graph.orig_bvh_))
    return false;

  for (size_t i = 0; i < stroke_graph.orig_strokes_.size(); ++i) {
    if (i == junc.points[0].first || i == junc.points[1].first)
      continue;
    if (stroke_graph.orig_bvh_->nodes[i].bb.contains(p2)) {
      Float dist = signed_dist_stroke(p2, stroke_graph.orig_strokes_[i]);
      if (dist <= 0)
        return false;
    }
  }

  return true;
}

bool corner_projection_check(
  const StrokeGraph &stroke_graph, const Junction &in_junc,
  bool to_check_t_visible = true,
  const StrokeGraph::VertexView corner_v = StrokeGraph::VertexView()) {
  Junction junc = in_junc;

  // Find the two adjacent strokes
  std::vector<Junction> bind_junc;
  expand_adjacent_junction(stroke_graph, junc, bind_junc);
  std::unordered_set<size_t> adj_orig_sid;
  for (const auto &j : bind_junc) {
    adj_orig_sid.emplace(j.points[0].first);
    adj_orig_sid.emplace(j.points[1].first);
  }

  if (to_check_t_visible && !(adj_orig_sid.count(junc.points[0].first) &&
                              adj_orig_sid.count(junc.points[1].first)))
    return false;

  // Project and find distances
  bool close_to_all = true;
  std::unordered_map<size_t, bool> orig_proj_close;
  for (const auto &proj_junc : bind_junc) {
    const auto orig_si = proj_junc.points[0].first;
    const auto dangling_si = proj_junc.points[1].first;

    auto _proj = Vec2::Empty();
    Float s;
    Vec2 proj;
    Float arclen = 0;
    stroke_graph.orig_strokes_[orig_si].ensure_arclengths();
    stroke_graph.orig_strokes_[dangling_si].ensure_arclengths();
    const auto dist = //
      (proj_junc.points[0].first == proj_junc.points[1].first)
        ? closest_point_to_own_stroke(stroke_graph.orig_strokes_[orig_si],
                                      proj_junc.points[1].first < 0.5, proj,
                                      arclen)
        : closest_point(stroke_graph.orig_strokes_[orig_si],
                        stroke_graph.orig_strokes_[dangling_si].pos_norm(
                          proj_junc.points[1].second),
                        proj, arclen);

    Float end_arclen = (proj_junc.points[0].second < 0.5)
                         ? 0
                         : stroke_graph.orig_strokes_[orig_si].length();

    if (!corner_v.is_valid()) {
      close_to_all &= end_stroke_junction_type(
                        stroke_graph.orig_strokes_[orig_si], end_arclen, arclen,
                        corner_projection_ratio) == 0;
      // If the rounded vertex is close to the endpoint
      close_to_all &= end_stroke_junction_type(
                        stroke_graph.orig_strokes_[orig_si], end_arclen,
                        proj_junc.points[0].second *
                          stroke_graph.orig_strokes_[orig_si].length(),
                        corner_projection_ratio) == 0;
    } else {
      // If the rounded vertex is close to the projection
      close_to_all &=
        end_stroke_junction_type(stroke_graph.orig_strokes_[orig_si], arclen,
                                 proj_junc.points[0].second *
                                   stroke_graph.orig_strokes_[orig_si].length(),
                                 corner_projection_ratio) == 0;
    }

    // To handle the self intersection case of closed strokes
    orig_proj_close[orig_si] |= close_to_all;
  }

  for (auto const [orig_si, close] : orig_proj_close) {
    if (!close)
      return false;
  }

  return true;
}

bool trial_add_corner_junction(const StrokeGraph &end_stroke_graph,
                               const StrokeGraph &stroke_graph,
                               const Junction &junc_round, const Junction &junc,
                               const Junction &s_junc,
                               std::vector<Junction> &corner_predictions) {
  bool added = false;
  Vec2 p = stroke_graph.orig_strokes_[junc_round.points[0].first].pos_norm(
    junc_round.points[0].second);
  Vec2 q = stroke_graph.orig_strokes_[junc_round.points[1].first].pos_norm(
    junc_round.points[1].second);

  // Check distance to the nearest vertex
  const auto &s = stroke_graph.strokes_[s_junc.points[0].first];
  const auto &orig_s = stroke_graph.orig_strokes_[junc.points[0].first];
  size_t orig_sid = junc.points[1].first;
  orig_s.ensure_arclengths();
  s.ensure_arclengths();
  if ( // end_stroke_junction_type(s, s.length() * s_junc.points[0].second,
       // corner_projection_ratio) == 0
    end_stroke_junction_type(
      orig_s, (junc_round.points[0].second < 0.5) ? 0 : orig_s.length(),
      orig_s.length() * junc.points[0].second, corner_projection_ratio) == 0 &&
    end_stroke_junction_type(
      orig_s, (junc_round.points[0].second < 0.5) ? 0 : orig_s.length(),
      orig_s.length() * junc_round.points[0].second,
      corner_projection_ratio) == 0 &&
    std::find_if(corner_predictions.begin(), corner_predictions.end(),
                 [&junc_round](const Junction &j) {
                   return (j.type == junc_round.type) &&
                          (j.points[0].first == junc_round.points[0].first) &&
                          (j.points[0].second == junc_round.points[0].second) &&
                          (j.points[1].first == junc_round.points[1].first) &&
                          (j.points[1].second == junc_round.points[1].second);
                 }) == corner_predictions.end() &&
    graph_visible(end_stroke_graph, junc_round)
    // line_of_sight(p, q, *stroke_graph.orig_bvh_)
  ) {
    int end_d0 = end_degree(stroke_graph, junc_round.points[0]),
        end_d1 = end_degree(stroke_graph, junc_round.points[1]);
    if ((end_d0 > 1 || end_d1 > 1) && (end_d0 == 1 || end_d1 == 1)) {
      bool close_to_all = corner_projection_check(stroke_graph, junc_round);
      if (close_to_all) {
        corner_predictions.emplace_back(junc_round);
        corner_predictions.back().orig_dist =
          junction_distance_init(stroke_graph, junc_round);
        corner_predictions.back().corner_type = JunctionType::T;
        added = true;
      }
    }
  }

  return added;
}

bool high_valence_dangling_visible(const StrokeGraph &stroke_graph,
                                   const Junction &junc) {
  bool visible = false;

  Vec2 q = stroke_graph.orig_strokes_[junc.points[1].first].pos_norm(
    junc.points[1].second);

  Junction s_junc = junc;
  auto ok = original_stroke_to_stroke_indexing(stroke_graph, s_junc);
  force_assert(ok && "couldn't map from orig to strokes");
  const auto v =
    endpoint_to_vertex(stroke_graph, Endpoint(s_junc.points[0].first,
                                              s_junc.points[0].second == 0.0));

  auto he = v.hedge();
  auto hit = he;
  do {
    StrokeTime s_end(hit.stroke_idx(), (hit.forward()) ? 0.0 : 1.0);
    bool ok2 = convert_strokes2orig(stroke_graph, s_end);
    force_assert(ok2 && "couldn't map from strokes to orig");

    // Skip interior points
    if (s_end.second == 0.0 || s_end.second == 1.0) {
      Vec2 p = stroke_graph.orig_strokes_[s_end.first].pos_norm(s_end.second);
      if (line_of_sight(p, q, *stroke_graph.orig_bvh_)
          // && envelope_visible(stroke_graph, junc)
      ) {
        visible = true;
        break;
      }
    }

    hit = hit.twin().next();
  } while (hit != he);

  return visible;
}

void add_e2e_corner_junction(const StrokeGraph &end_stroke_graph,
                             const StrokeGraph &stroke_graph, Junction junc,
                             std::vector<Junction> &corner_predictions) {
  if (high_valence_dangling_visible(end_stroke_graph, junc)) {
    Junction s_junc = junc;
    auto ok = original_stroke_to_stroke_indexing(stroke_graph, s_junc);
    force_assert(ok && "couldn't map from orig to strokes");

    // Round the end to the closest vertex to check the valence
    Junction junc_round = s_junc;
    junc_round.points[0].second =
      (junc_round.points[0].second < 0.5) ? 0.0 : 1.0;
    ok = stroke_to_original_stroke_indexing(stroke_graph, junc_round);
    force_assert(ok && "couldn't map from strokes to orig");

    // Swap if sorting doesn't work
    {
      int end_d0 = end_degree(stroke_graph, junc_round.points[0]),
          end_d1 = end_degree(stroke_graph, junc_round.points[1]);
      if (end_d1 != 1) {
        std::swap(junc.points[0], junc.points[1]);
        std::swap(s_junc.points[0], s_junc.points[1]);
        std::swap(junc_round.points[0], junc_round.points[1]);
      }
    }

    const auto v = endpoint_to_vertex(
      stroke_graph,
      Endpoint(s_junc.points[0].first, s_junc.points[0].second == 0.0));
    assert(v.is_valid());

    size_t orig_sid = junc.points[1].first;

    auto corner_close =
      corner_projection_check(stroke_graph, junc_round, false, v);
    // If we are generating corners for reconsideration with end-end classifier
    // predictions, keep the projection test.
    if ((&end_stroke_graph == &stroke_graph &&
         recosider_type == JunctionType::R) ||
        corner_close) {
      junc.orig_dist = junction_distance_init(stroke_graph, junc);
      junc.corner_type = JunctionType::R;
      corner_predictions.emplace_back(junc);
    }
  }
}

bool connection_corner_visible(const StrokeGraph &stroke_graph,
                               const Junction &junc) {
  // Find the two adjacent strokes
  std::vector<Junction> bind_junc;
  expand_adjacent_junction(stroke_graph, junc, bind_junc);

  // Check end visibility
  bool seen_endpoint = false;
  bool blocked = false;
  bool seen_interior = false;
  for (const auto &j : bind_junc) {
    for (const auto &end : j.points) {
      if (end != junc.points[1]) {
        if (end.second == 0.0 || end.second == 1.0)
          seen_endpoint = true;
        else
          seen_interior = true;
        Vec2 p = stroke_graph.orig_strokes_[j.points[0].first].pos_norm(
          j.points[0].second);
        Vec2 q = stroke_graph.orig_strokes_[j.points[1].first].pos_norm(
          j.points[1].second);
        if (!line_of_sight(p, q, *stroke_graph.orig_bvh_) ||
            !graph_visible(stroke_graph, j)) {
          blocked = true;
          break;
        }
      }
    }
  }

  if (blocked || !seen_endpoint)
    return false;
  if (recosider_type == JunctionType::R && seen_interior)
    return false;

  return true;
}

} // namespace

int to_v_idx(const StrokeGraph &stroke_graph, const StrokeTime &end) {
  StrokeTime s_end = end;
  auto ok = convert_orig2strokes(stroke_graph, s_end);
  // This original stroke is entirely deleted in the final graph
  assert(ok);
  assert(s_end.second == 0.0 || s_end.second == 1.0);
  const auto v = endpoint_to_vertex(stroke_graph,
                                    Endpoint(s_end.first, s_end.second == 0.0));
  assert(v.is_valid() && v.is_active());
  return v.index_;
}

void print_junctions(const StrokeGraph &stroke_graph,
                     const std::vector<Junction> &candidates) {
  SPDLOG_DEBUG("===========================");
  std::vector<size_t> non_corner_indices;
  std::map<std::pair<int, int>, std::vector<size_t>> high_valence_junctions;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto &junc = candidates[i];
    if (!is_corner_candidate(stroke_graph, junc)) {
      non_corner_indices.emplace_back(i);
    } else {
      const auto junc_p =
        std::make_pair(std::min(to_v_idx(stroke_graph, junc.points[0]),
                                to_v_idx(stroke_graph, junc.points[1])),
                       std::max(to_v_idx(stroke_graph, junc.points[0]),
                                to_v_idx(stroke_graph, junc.points[1])));

      if (high_valence_junctions.count(junc_p) == 0) {
        high_valence_junctions[junc_p] = std::vector<size_t>();
      }

      high_valence_junctions[junc_p].emplace_back(i);
    }
  }

  for (const auto &[p, juncs] : high_valence_junctions) {
    for (const auto j_idx : juncs) {
      const auto &junc = candidates[j_idx];
      SPDLOG_DEBUG(
        "Junc: {}, {} - {}, {}:\n\t{} comp. {}; {}; {}; {}",
        junc.points[0].first, junc.points[0].second, junc.points[1].first,
        junc.points[1].second, junc.component_idx, junc.probability,
        (junc.type == JunctionType::R) ? "R" : "T",
        (junc.repr.empty()) ? "Dis" : "Conn", (junc.is_weak) ? "WEAK" : "");
    }
  }
  if (!non_corner_indices.empty()) {
    SPDLOG_DEBUG("= Non corners =");
    for (const auto j_idx : non_corner_indices) {
      const auto &junc = candidates[j_idx];
      SPDLOG_DEBUG(
        "Junc: {}, {} - {}, {}:\n\t{} comp. {}; {}; {}; {}",
        junc.points[0].first, junc.points[0].second, junc.points[1].first,
        junc.points[1].second, junc.component_idx, junc.probability,
        (junc.type == JunctionType::R) ? "R" : "T",
        (junc.repr.empty()) ? "Dis" : "Conn", (junc.is_weak) ? "WEAK" : "");
    }
  }
}

void print_junctions(const std::vector<Junction> &candidates) {
  SPDLOG_INFO("===========================");
  for (const auto &junc : candidates)
    SPDLOG_INFO("Junc: {}, {} - {}, {}:\n\t{} comp. {}; {}; {}; {}",
                junc.points[0].first, junc.points[0].second,
                junc.points[1].first, junc.points[1].second, junc.component_idx,
                junc.probability, (junc.type == JunctionType::R) ? "R" : "T",
                (junc.repr.empty()) ? "Dis" : "Conn",
                (junc.is_weak) ? "WEAK" : "");
}

int end_degree(const StrokeGraph &stroke_graph, const StrokeTime &end) {
  StrokeTime s_end = end;
  auto ok = convert_orig2strokes(stroke_graph, s_end);
  // This original stroke is entirely deleted in the final graph
  if (!ok)
    return -1;
  if (s_end.second != 0.0 && s_end.second != 1.0)
    return 0;
  const auto v = endpoint_to_vertex(stroke_graph,
                                    Endpoint(s_end.first, s_end.second == 0.0));
  assert(v.is_valid() && v.is_active());
  return (int)v.valence();
}

bool has_bridge_vertex(const StrokeGraph &stroke_graph, const Junction &junc) {
  return false;
}

bool is_snap_valid(const StrokeGraph &stroke_graph, const Junction &junc) {
  auto vert = orig2endpoint(stroke_graph, junc.points[1]);
  // Invalid junction
  if (!vert.is_valid() || !vert.is_active())
    return false;

  auto vert2 = orig2endpoint(stroke_graph, junc.points[0]);

  // Special case of snapping two copies of the same vertices. This is
  // considered valid since this can emerge naturally in sequenctially snapping
  // high-valence vertices.
  if (vert == vert2)
    return true;

  if (junc.type == JunctionType::T) {
    auto pos = junc.points[0];
    auto ok = convert_orig2strokes(stroke_graph, pos);
    force_assert(ok && "couldn't map from orig to strokes");
    return is_snap_valid(stroke_graph, vert.index_, pos);
  } else {
    // Invalid junction
    if (!vert2.is_valid() || !vert2.is_active())
      return false;

    return is_snap_valid(stroke_graph, vert.index_, vert2.index_);
  }
}

bool is_high_valence_junction_on_boundary_cycle(
  const StrokeGraph &graph, const span<const Junction> candidates, size_t fi,
  const StrokeGraph::VertexView v, const std::string &vid_str) {
  return is_high_valence_junction_on_cycle(graph, candidates, fi, v, vid_str,
                                           true);
}

bool is_high_valence_junction_in_region_cycle(
  const StrokeGraph &graph, const span<const Junction> candidates, size_t fi,
  const StrokeGraph::VertexView v, const std::string &vid_str) {
  return is_high_valence_junction_on_cycle(graph, candidates, fi, v, vid_str,
                                           false);
}

bool is_candidate_valid(const StrokeGraph &stroke_graph,
                        span<const Stroke> strokes,
                        const StrokeGraph &final_stroke_graph,
                        size_t cur_stroke_step, const Junction &junc,
                        const std::vector<Junction> &final_predictions) {
  // Use the stroke based indexing for the position since the original strokes
  // may have been deformed
  Junction s_junc = junc;
  if (!original_stroke_to_stroke_indexing(stroke_graph, s_junc)) {
    return false;
  }
  Junction s_junc_final = junc;
  if (!original_stroke_to_stroke_indexing(final_stroke_graph, s_junc_final)) {
    return false;
  }

  Vec2 p = stroke_graph.strokes_[s_junc.points[0].first].pos_norm(
    s_junc.points[0].second);
  Vec2 q = stroke_graph.strokes_[s_junc.points[1].first].pos_norm(
    s_junc.points[1].second);

  // Use the current strokes in the graph and all the future strokes for the
  // intersection test
  assert(cur_stroke_step + 1 <= strokes.size());

  // We also drop prob=0 junctions
  if (junc.probability > std::numeric_limits<double>::epsilon() &&
      (final_predictions.empty() ||
       is_final_prediction_acceptable(junc, final_predictions)) &&
      is_snap_valid(stroke_graph, junc) &&
      !should_delay(final_stroke_graph, junc, cur_stroke_step) &&
      line_of_sight(p, q, stroke_graph.strokes_,
                    stroke_graph.bvh_.centerline_bbs()) &&
      is_visible_in_final(
        span<const Stroke>(strokes.begin() + cur_stroke_step + 1,
                           strokes.size() - (cur_stroke_step + 1)),
        p, q)) {
    return true;
  }
  return false;
}

bool is_weak_candidate_valid(const StrokeGraph &stroke_graph,
                             span<const Stroke> strokes,
                             const StrokeGraph &final_stroke_graph,
                             size_t cur_stroke_step, const Junction &junc) {
  // Use the stroke based indexing for the position since the original strokes
  // may have been deformed
  Junction s_junc = junc;
  if (!original_stroke_to_stroke_indexing(stroke_graph, s_junc)) {
    return false;
  }
  Junction s_junc_final = junc;
  if (!original_stroke_to_stroke_indexing(final_stroke_graph, s_junc_final)) {
    return false;
  }

  Vec2 p = stroke_graph.strokes_[s_junc.points[0].first].pos_norm(
    s_junc.points[0].second);
  Vec2 q = stroke_graph.strokes_[s_junc.points[1].first].pos_norm(
    s_junc.points[1].second);

  // Use the current strokes in the graph and all the future strokes for the
  // intersection test
  assert(cur_stroke_step + 1 <= strokes.size());

  // We also drop prob=0 junctions
  if (junc.probability > std::numeric_limits<double>::epsilon() &&
      is_snap_valid(stroke_graph, junc) &&
      !should_delay(final_stroke_graph, junc, cur_stroke_step) &&
      line_of_sight(p, q, stroke_graph.strokes_,
                    stroke_graph.bvh_.centerline_bbs()) &&
      is_visible_in_final(
        span<const Stroke>(strokes.begin() + cur_stroke_step + 1,
                           strokes.size() - (cur_stroke_step + 1)),
        p, q)) {
    return true;
  }
  return false;
}

void expand_junction(const StrokeGraph &graph, const Junction &in_junc,
                     std::vector<Junction> &bind_junc) {
  Junction junc = in_junc;
  if (!original_stroke_to_stroke_indexing(graph, junc)) {
    return;
  }
  Endpoint c1(junc.points[1].first, junc.points[1].second < 0.5);
  std::pair<size_t, Float> cand2 = junc.points[0];
  if (junc.type == JunctionType::R) {
    // Expand high-valence junction candidates
    size_t n_evaluations;
    Endpoint c2{cand2.first, cand2.second == 0.0};

    // c1 can't be continuous
    if (!can_expand(graph, c1) || !can_expand(graph, c2, false)) {
      if (!can_expand(graph, c2) || !can_expand(graph, c1, false)) {
        return;
      }
      std::swap(c1, c2);
    }

    size_t n = n_classifier_evaluations_vv(graph, c1, c2);
    n_evaluations = n;

    auto expanded_cand1 = std::vector<Endpoint>(n_evaluations, Endpoint(0));
    auto expanded_cand2 = std::vector<StrokeTime>(n_evaluations);
    expand_candidates_vv(graph, c1, c2, {&expanded_cand1[0], n},
                         {&expanded_cand2[0], n});

    for (size_t j = 0; j < n; ++j) {
      SnapInfo snap_info;
      const auto &other_si = expanded_cand2[j].first;
      snap_info.status = SnapInfo::PredictionDelay;
      snap_info.max_prob = -1;
      snap_info.prediction_type = JunctionType::R;
      snap_info.stroke_pos1 = StrokeTime{(int)expanded_cand1[j].stroke_idx(),
                                         (expanded_cand1[j].is_head()) ? 0 : 1};
      snap_info.stroke_pos2 =
        StrokeTime{(int)other_si, expanded_cand2[j].second};
      Junction new_junc(
        std::vector<std::pair<int, double>>{
          std::make_pair((int)snap_info.stroke_pos2.first,
                         snap_info.stroke_pos2.second),
          std::make_pair((int)snap_info.stroke_pos1.first,
                         snap_info.stroke_pos1.second)},
        snap_info.prediction_type, false, 0.5, true);
      const auto ok = stroke_to_original_stroke_indexing(graph, new_junc);
      force_assert(ok && "couldn't map from strokes to orig");
      bind_junc.emplace_back(new_junc);
    }
  } else if (junc.type == JunctionType::T) {
    // Expand high-valence junction candidates
    size_t n_evaluations_total = 0;
    if (!can_expand(graph, c1)) {
      return;
    }

    n_evaluations_total = n_classifier_evaluations_vs(graph, c1);

    auto expanded_cand1 =
      std::vector<Endpoint>(n_evaluations_total, Endpoint(0));
    expand_candidates_vs(graph, c1, {&expanded_cand1[0], n_evaluations_total});

    const auto [other_si, other_arclen] = cand2;
    for (size_t j = 0; j < expanded_cand1.size(); ++j) {
      SnapInfo snap_info;
      snap_info.status = SnapInfo::PredictionDelay;
      snap_info.max_prob = -1;
      snap_info.prediction_type = JunctionType::T;
      snap_info.stroke_pos1 = StrokeTime{(int)expanded_cand1[j].stroke_idx(),
                                         (expanded_cand1[j].is_head()) ? 0 : 1};
      snap_info.stroke_pos2 = StrokeTime{(int)other_si, other_arclen};
      Junction new_junc(
        std::vector<std::pair<int, double>>{
          std::make_pair((int)snap_info.stroke_pos2.first,
                         snap_info.stroke_pos2.second),
          std::make_pair((int)snap_info.stroke_pos1.first,
                         snap_info.stroke_pos1.second)},
        snap_info.prediction_type, false, 0.5, true);
      const auto ok = stroke_to_original_stroke_indexing(graph, new_junc);
      force_assert(ok && "couldn't map from strokes to orig");
      bind_junc.emplace_back(new_junc);
    }
  } else
    abort();
}

void expand_adjacent_junction(const StrokeGraph &graph, const Junction &junc_,
                              std::vector<Junction> &bind_junc) {
  bind_junc.clear();

  Junction in_junc = junc_;
  Junction junc = in_junc;
  if (!original_stroke_to_stroke_indexing(graph, junc) ||
      (junc.points[1].second != 0.0 && junc.points[1].second != 1.0)) {
    return;
  }
  if (junc.type == JunctionType::T &&
      (junc.points[0].second != 0.0 && junc.points[0].second != 1.0)) {
    bind_junc.emplace_back(in_junc);
    return;
  }

  // Swap so the first element is non-dangling
  {
    const auto v = endpoint_to_vertex(
      graph, Endpoint(junc.points[0].first, junc.points[0].second == 0.0));
    const auto dangling_v = endpoint_to_vertex(
      graph, Endpoint(junc.points[1].first, junc.points[1].second == 0.0));
    if (v.valence() == 1) {
      std::swap(in_junc.points[0], in_junc.points[1]);
      std::swap(junc.points[0], junc.points[1]);
    }
  }

  // Get the two adjacent strokes
  const auto v = endpoint_to_vertex(
    graph, Endpoint(junc.points[0].first, junc.points[0].second == 0.0));
  assert(v.is_valid() && v.is_active());
  if (v.valence() == 1) {
    bind_junc.emplace_back(in_junc);
    return;
  }

  auto dangling_v = endpoint_to_vertex(
    graph, Endpoint(junc.points[1].first, junc.points[1].second == 0.0));
  assert(dangling_v.is_valid() && dangling_v.is_active());
  if (dangling_v.valence() != 1) {
    return;
  }

  StrokeGraph varying_graph = graph.clone();
  if (!is_snap_valid(varying_graph, v.index_, dangling_v.index_))
    return;
  bool snapped = varying_graph.snap_vertices(
    v.index_, dangling_v.index_, varying_graph.snapping_method_type_);
  if (!snapped)
    return;

  dangling_v = varying_graph.vertex(dangling_v.index_);

  // Circulate
  std::vector<StrokeTime> star_he;
  std::vector<size_t> star_he_idx;
  std::vector<size_t> out_he_idx;
  auto he = dangling_v.hedge();
  auto hit = he;
  do {
    StrokeTime s_end(hit.stroke_idx(), (hit.forward()) ? 0.0 : 1.0);
    bool ok = convert_strokes2orig(varying_graph, s_end);
    force_assert(ok && "couldn't map from strokes to orig");
    star_he.emplace_back(s_end);
    star_he_idx.emplace_back(hit.index_);
    hit = hit.twin().next();
  } while (hit != he);
  for (size_t i = 0; i < star_he.size(); ++i) {
    size_t prev = (i + star_he.size() - 1) % star_he.size();
    size_t next = (i + 1) % star_he.size();
    if ((star_he[i].first == in_junc.points[1].first) &&
        (star_he[i].second == in_junc.points[1].second)) {
      Junction hv_junc = in_junc;
      hv_junc.type = JunctionType::R;
      hv_junc.points[0] = star_he[prev];
      hit = varying_graph.hedge(star_he_idx[prev]);
      StrokeTime s_end1(hit.stroke_idx(), (hit.forward()) ? 0.0 : 1.0);
      StrokeTime to_s_end1(hit.stroke_idx(), (hit.forward()) ? 1.0 : 0.0);
      bool ok = convert_strokes2orig(varying_graph, s_end1);
      force_assert(ok && "couldn't map from strokes to orig");
      ok = convert_strokes2orig(varying_graph, to_s_end1);
      force_assert(ok && "couldn't map from strokes to orig");
      hv_junc.point_tos.clear();
      hv_junc.point_tos.emplace_back(to_s_end1);
      hv_junc.corner_type = in_junc.corner_type;
      bind_junc.emplace_back(hv_junc);
      out_he_idx.emplace_back(hit.index_);
      // Assuming no dissolving
      assert(to_s_end1.first == hv_junc.points[0].first);

      hv_junc.points[0] = star_he[next];
      hit = varying_graph.hedge(star_he_idx[next]);
      StrokeTime s_end2(hit.stroke_idx(), (hit.forward()) ? 0.0 : 1.0);
      StrokeTime to_s_end2(hit.stroke_idx(), (hit.forward()) ? 1.0 : 0.0);
      ok = convert_strokes2orig(varying_graph, s_end2);
      force_assert(ok && "couldn't map from strokes to orig");
      ok = convert_strokes2orig(varying_graph, to_s_end2);
      force_assert(ok && "couldn't map from strokes to orig");
      hv_junc.point_tos.clear();
      hv_junc.point_tos.emplace_back(to_s_end2);
      hv_junc.corner_type = in_junc.corner_type;
      bind_junc.emplace_back(hv_junc);
      out_he_idx.emplace_back(hit.index_);
      // Assuming no dissolving
      assert(to_s_end2.first == hv_junc.points[0].first);
    }
  }

  // Convert the dangling-continous pairs to a single T-junction (so the region
  // code can handle it)
  if (bind_junc.size() == 2 &&
      varying_graph.hedge(out_he_idx.front()).continuity() ==
        out_he_idx.back()) {
    std::vector<Junction> tmp_bind_junc;
    tmp_bind_junc.emplace_back(bind_junc.front());
    tmp_bind_junc.back().type = JunctionType::T;
    bind_junc = std::move(tmp_bind_junc);
  }
}

bool update_disconnected_junction_predictions(
  const StrokeGraph &stroke_graph, std::vector<Junction> &varying_candidates,
  FeatureType feature_type, bool to_expand, bool use_input_junc_type) {
  prediction_feature_type = feature_type;

  if (varying_candidates.empty())
    return true;

  auto cand1 = std::vector<Endpoint>();
  cand1.reserve(varying_candidates.size());
  auto cand2 = std::vector<StrokeTime>();
  cand2.reserve(varying_candidates.size());
  // Because not all candidates will be valid, we need to map from
  // `varying_candidates` indices into a possibly reduced size cand1 and cand2.
  auto indices = std::vector<size_t>();
  indices.reserve(varying_candidates.size());

  std::vector<JunctionType::Type> junc_types;
  bool all_valid = true;
  for (size_t i = 0; i < varying_candidates.size(); ++i) {
    auto &junc = varying_candidates[i];

    auto s_end = junc.points[1];
    bool ok = convert_orig2strokes(stroke_graph, s_end);
    force_assert(ok && "couldn't map from orig to strokes");

    // This is an enforced junction
    if (junc.probability == 1.0) {
      indices.push_back((size_t)-2);
      continue;
    }

    // This is an invalid junction only used for high-valence check
    if (junc.is_weak) {
      if (!(s_end.second < 1e-10 || std::abs(1.0 - s_end.second) < 1e-10)) {
        junc.must_disconnect = true;
        indices.push_back((size_t)-1);
        continue;
      }
    }

    assert(s_end.second < 1e-10 || std::abs(1.0 - s_end.second) < 1e-10);
    Junction pred_junc = junc;
    const auto could_map =
      original_stroke_to_stroke_indexing(stroke_graph, pred_junc);
    if (!could_map || !pred_junc.is_valid() ||
        !is_snap_valid(stroke_graph, junc) ||
        has_continuous_edge(endpoint_to_vertex(
          stroke_graph, Endpoint(pred_junc.points[1].first,
                                 pred_junc.points[1].second < 0.5)))) {
      if (!junc.is_weak)
        SPDLOG_WARN("Incorrect junction: {}, {} - {}, {} => {}, {} - {}, {}",
                    junc.points[0].first, junc.points[0].second,
                    junc.points[1].first, junc.points[1].second,
                    pred_junc.points[0].first, pred_junc.points[0].second,
                    pred_junc.points[1].first, pred_junc.points[1].second);
      all_valid = false;
      indices.push_back((size_t)-1);
      if (junc.is_weak)
        junc.must_disconnect = true;
    } else {
      cand1.emplace_back(pred_junc.points[1].first,
                         pred_junc.points[1].second < 0.5);
      // Stroke interiors come first in a T-junction Junction.
      cand2.emplace_back(pred_junc.points[0].first, pred_junc.points[0].second);
      indices.push_back(cand1.size() - 1);
      if (use_input_junc_type) {
        junc_types.emplace_back(pred_junc.type);
      }
    }
  }

  auto prob = std::vector<Float>(cand1.size());
  std::vector<ClassifierPrediction> out_predictions;
  compute_junction_probabilities(stroke_graph, cand1, cand2, feature_type, prob,
                                 &out_predictions, to_expand, junc_types);

  // Assign new probabilities.
  assert(indices.size() == varying_candidates.size());
  for (size_t i = 0; i < varying_candidates.size(); ++i) {
    if (indices[i] == (size_t)-1) {
      varying_candidates[i].probability = 0;
    } else if (indices[i] != (size_t)-2) {
      varying_candidates[i].probability = prob[indices[i]];
      varying_candidates[i].fea = out_predictions[indices[i]].fea;
    }
  }

  return all_valid;
}

bool update_high_valence_junction_predictions(
  const StrokeGraph &stroke_graph, std::vector<Junction> &varying_candidates,
  FeatureType feature_type, bool to_dedup) {
  // 1. Predict in the old way for all junctions
  update_disconnected_junction_predictions(stroke_graph, varying_candidates,
                                           feature_type, false);
  // SPDLOG_INFO("Regular predictions:");
  // print_junctions(varying_candidates);

  // 2. Gather the high-valence junctions
  std::vector<size_t> reorder_indices;
  std::map<std::pair<int, int>, std::vector<size_t>> high_valence_junctions;
  for (size_t i = 0; i < varying_candidates.size(); ++i) {
    const auto &junc = varying_candidates[i];
    const auto junc_p =
      std::make_pair(std::min(to_v_idx(stroke_graph, junc.points[0]),
                              to_v_idx(stroke_graph, junc.points[1])),
                     std::max(to_v_idx(stroke_graph, junc.points[0]),
                              to_v_idx(stroke_graph, junc.points[1])));

    if (high_valence_junctions.count(junc_p) == 0) {
      high_valence_junctions[junc_p] = std::vector<size_t>();
    }

    high_valence_junctions[junc_p].emplace_back(i);
  }

  // 3. Predict individually
  StrokeGraph varying_stroke_graph = stroke_graph.clone();
  for (auto &[p, juncs] : high_valence_junctions) {
    if (juncs.size() == 1 &&
        varying_candidates[juncs.front()].point_tos.empty()) {
      reorder_indices.emplace_back(juncs.front());
      continue;
    }
    // For now, set all junctions to have the same prediction
    std::vector<size_t> combine_he_idx;
    for (const auto i : juncs) {
      const auto &junc = varying_candidates[i];
      const auto &v = varying_stroke_graph.vertex(
        to_v_idx(varying_stroke_graph, junc.points[0]));
      assert(!junc.point_tos.empty());

      bool found_he = false;
      const auto he = v.hedge();
      auto it = he;
      do {
        StrokeTime point_end(it.stroke_idx(), it.forward());
        auto ok = convert_strokes2orig(varying_stroke_graph, point_end);
        force_assert(ok && "couldn't map from strokes to orig");
        if (junc.point_tos.front() == point_end) {
          found_he = true;
          combine_he_idx.emplace_back(it.index_);
          break;
        }
        it = it.twin().next();
      } while (it != he);
      assert(found_he);
    }

    // For continuous corner
    bool found_he = true;
    if (juncs.size() == 1) {
      assert(combine_he_idx.size() == 1);
      const auto &junc = varying_candidates[juncs.front()];
      const auto &v = varying_stroke_graph.vertex(
        to_v_idx(varying_stroke_graph, junc.points[0]));
      assert(!junc.point_tos.empty());

      found_he = false;
      const auto he = v.hedge();
      auto it = he;
      do {
        if (it.continuity() == combine_he_idx.front()) {
          found_he = true;
          combine_he_idx.emplace_back(it.index_);
          break;
        }
        it = it.twin().next();
      } while (it != he);

      if (!found_he) {
        varying_candidates[juncs.front()].probability = 0;
      }
    }

    if (found_he) {
      assert(combine_he_idx.size() == 2);

      std::set<std::pair<Float, size_t>> corner_probabilities;
      for (const auto i : juncs) {
        update_high_valence_prediction(varying_stroke_graph, combine_he_idx,
                                       varying_candidates[i]);
        corner_probabilities.emplace(-varying_candidates[i].probability, i);
      }
      size_t junc_idx = corner_probabilities.begin()->second;
      for (const auto i : juncs) {
        auto &junc = varying_candidates[i];
        junc.probability = varying_candidates[junc_idx].probability;
        junc.fea = varying_candidates[junc_idx].fea;
      }

      // Reorder based on the probability so the visualization is consistent
      std::vector<size_t> reorder_juncs;
      for (const auto [neg_prob, i] : corner_probabilities) {
        reorder_indices.emplace_back(i);
        reorder_juncs.emplace_back(i);
      }
      force_assert(reorder_juncs.size() == juncs.size());
      juncs = std::move(reorder_juncs);
    } else {
      reorder_indices.emplace_back(juncs.front());
    }
  }

  if (to_dedup) {
    std::vector<Junction> dedup_candidates;
    for (const auto &[p, juncs] : high_valence_junctions) {
      if (juncs.size() == 1 &&
          varying_candidates[juncs.front()].point_tos.empty()) {
        dedup_candidates.emplace_back(varying_candidates[juncs.front()]);
        continue;
      }

      for (const auto i : juncs) {
        auto &junc = varying_candidates[i];
        dedup_candidates.emplace_back(junc);
        break;
      }
    }
    varying_candidates = std::move(dedup_candidates);
  } else // Reorder
  {
    force_assert(reorder_indices.size() == varying_candidates.size());
    std::vector<Junction> reorder_candidates;
    for (const auto i : reorder_indices) {
      reorder_candidates.emplace_back(std::move(varying_candidates[i]));
    }
    varying_candidates = std::move(reorder_candidates);
  }

  // SPDLOG_INFO("Corner predictions:");
  // print_junctions(varying_candidates);

  return true;
}

void complete_candidate_graph(const StrokeGraph &graph,
                              const std::vector<Junction> &candidates,
                              std::vector<Junction> &out_candidates,
                              const std::vector<bool> &connectivity) {
  // 1. Build an endpoint graph
  std::unordered_map<size_t, std::vector<double>> v2ends;
  std::vector<double> endpoint_times;
  endpoint_times.reserve(candidates.size() * 2);

  auto find_endpoint_time =
    [&endpoint_times](const double junc_hash) -> size_t {
    auto junc_itr = std::find_if(
      endpoint_times.begin(), endpoint_times.end(),
      [&junc_hash](const double t) { return std::abs(t - junc_hash) < 1e-10; });
    return junc_itr - endpoint_times.begin();
  };

  std::unordered_map<size_t, std::pair<int, double>> endpoint_indexing;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto &junc = candidates[i];
    for (const auto &end : junc.points) {
      StrokeTime s_end = end;
      auto ok = convert_orig2strokes(graph, s_end);
      force_assert(ok && "couldn't map from orig to strokes");

      double junc_hash = end.first * 10 + end.second;
      if (s_end.second == 0.0 || s_end.second == 1.0) {
        auto v =
          endpoint_to_vertex(graph, Endpoint(s_end.first, s_end.second == 0.0));
        if (!v2ends.count(v.index_))
          v2ends[v.index_] = std::vector<double>();
        if (std::find(v2ends[v.index_].begin(), v2ends[v.index_].end(),
                      junc_hash) == v2ends[v.index_].end())
          v2ends[v.index_].emplace_back(junc_hash);
      }
      size_t idx = endpoint_times.size();

      size_t end_i = find_endpoint_time(junc_hash);

      if (end_i == endpoint_times.size()) {
        endpoint_times.emplace_back(junc_hash);
        endpoint_indexing[idx] = end;
      }
    }
  }

  // 2. Label connected components
  std::vector<int> components(endpoint_times.size(), -1);
  int max_color = 0;

  auto merge_components = [&max_color, find_endpoint_time,
                           &components](double junc_hash1, double junc_hash2) {
    size_t end_idx1 = find_endpoint_time(junc_hash1);
    size_t end_idx2 = find_endpoint_time(junc_hash2);
    if (components[end_idx1] == -1 && components[end_idx2] == -1) {
      components[end_idx1] = components[end_idx2] = max_color++;
    } else if (components[end_idx1] != components[end_idx2]) {
      int min_c = std::min(components[end_idx1], components[end_idx2]);
      int max_c = std::max(components[end_idx1], components[end_idx2]);
      if (min_c == -1)
        components[end_idx1] = components[end_idx2] = max_c;
      else {
        for (auto &c : components)
          if (c == min_c)
            c = max_c;
      }
    }
  };

  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto &junc = candidates[i];
    double junc_hash1 = junc.points[0].first * 10 + junc.points[0].second;
    double junc_hash2 = junc.points[1].first * 10 + junc.points[1].second;
    merge_components(junc_hash1, junc_hash2);
  }

  // Now color using the snapped ends
  for (const auto &[v, ends] : v2ends) {
    for (size_t i = 0; i + 1 < ends.size(); ++i) {
      for (size_t j = i + 1; j < ends.size(); ++j) {
        double junc_hash1 = ends[i];
        double junc_hash2 = ends[j];
        merge_components(junc_hash1, junc_hash2);
      }
    }
  }

  // 3. Make connected components complete
  std::unordered_set<int> colors;
  colors.insert(components.begin(), components.end());
  for (auto c : colors) {
    for (size_t i = 0; c >= 0 && i + 1 < components.size(); ++i) {
      if (components[i] != c)
        continue;
      for (size_t j = i + 1; j < components.size(); ++j) {
        if (components[j] != c)
          continue;
        Junction junc(std::vector<std::pair<int, double>>{endpoint_indexing[i],
                                                          endpoint_indexing[j]},
                      JunctionType::R, true, 0.5, true);
        junc.sort_entries();
        auto s_end = junc.points[0];
        bool ok = convert_orig2strokes(graph, s_end);
        force_assert(ok && "couldn't map from orig to strokes");
        junc.type =
          (s_end.second < 1e-10 || std::abs(1.0 - s_end.second) < 1e-10)
            ? JunctionType::R
            : JunctionType::T;
        auto s_end2 = junc.points[1];
        ok = convert_orig2strokes(graph, s_end2);
        force_assert(ok && "couldn't map from orig to strokes");
        // if (!(s_end2.second < 1e-10 || std::abs(1.0 - s_end2.second) <
        // 1e-10))
        //  junc.must_disconnect = true;

        // Avoid creating interior-interior junctions since it's discarded even
        // before graph modification
        if (junc.type == JunctionType::T &&
            !(s_end2.second < 1e-10 || std::abs(1.0 - s_end2.second) < 1e-10))
          continue;

        // Avoid creating high-valence-high-valence junctions
        if (junc.type == JunctionType::R) {
          const auto v1 = endpoint_to_vertex(
            graph, Endpoint(s_end.first, s_end.second < 0.5));
          const auto v2 = endpoint_to_vertex(
            graph, Endpoint(s_end2.first, s_end2.second < 0.5));
          if (!v1.is_dangling() && !v2.is_dangling())
            continue;
        }

        // Avoid creating high-valence-interior junctions
        if (junc.type == JunctionType::T) {
          const auto v2 = endpoint_to_vertex(
            graph, Endpoint(s_end2.first, s_end2.second < 0.5));
          if (!v2.is_dangling())
            continue;
        }

        auto junc_itr = std::find_if(
          candidates.begin(), candidates.end(), [&junc](const Junction &junc_) {
            return (junc_.type == junc.type) &&
                   (((junc_.points[0].first == junc.points[0].first) &&
                     (std::abs(junc_.points[0].second - junc.points[0].second) <
                      1e-10) &&
                     (junc_.points[1].first == junc.points[1].first) &&
                     (std::abs(junc_.points[1].second - junc.points[1].second) <
                      1e-10)) ||
                    ((junc_.points[0].first == junc.points[1].first) &&
                     (std::abs(junc_.points[0].second - junc.points[1].second) <
                      1e-10) &&
                     (junc_.points[1].first == junc.points[0].first) &&
                     (std::abs(junc_.points[1].second - junc.points[0].second) <
                      1e-10)));
          });
        // Add to the result set
        if (junc_itr == candidates.end()) {
          auto v1 = endpoint_to_vertex(
            graph, Endpoint{(size_t)s_end.first, s_end.second < 0.5});
          auto v2 = endpoint_to_vertex(
            graph, Endpoint{(size_t)s_end2.first, s_end2.second < 0.5});
          // Ends are not snapped
          if (v1 != v2) {
            junc.orig_dist = junction_distance_init(graph, junc);
            out_candidates.emplace_back(junc);
          }
        }
      }
    }
  }
}

bool reproject_t_junc(const StrokeGraph &graph, Junction &junc) {
  // This would only happen if both ends get dissolved
  auto v = orig2endpoint(graph, junc.points[1]);
  if (!v.is_valid() || !v.is_active()) {
    return false;
  }

  // Reproject and redo the line-of-sight check to account for deformed
  // geometry.
  auto cand = Junction{junc}; // Copy.

  auto ok = original_stroke_to_stroke_indexing(graph, cand);
  if (!ok)
    return false;

  const auto si = cand.points[0].first;
  const auto &strokes = graph.strokes_;
  auto _proj = Vec2::Empty();
  Float s;
  const auto dist = find_projection(v, si, _proj, s);
  if (dist == std::numeric_limits<Float>::infinity()) {
    return false;
  }
  cand.points[0].second = s / strokes[si].length();
  if (!line_of_sight(
        strokes[cand.points[1].first].pos_norm(cand.points[1].second),
        strokes[si].pos(s), graph.strokes_, graph.bvh_.centerline_bbs())) {
    return false;
  }

  ok = stroke_to_original_stroke_indexing(graph, cand);
  if (!ok)
    return false;

  junc = cand;

  return true;
}

bool reproject_t_junc_orig(const StrokeGraph &graph, Junction &junc) {
  Vec2 v;
  bool is_head;
  size_t si, v_si, write_pi;
  if (junc.points[1].second == 0.0 || junc.points[1].second == 1.0) {
    v =
      graph.orig_strokes_[junc.points[1].first].pos_norm(junc.points[1].second);
    is_head = junc.points[1].second < 0.5;
    si = 0;
    v_si = 1;
  } else {
    v =
      graph.orig_strokes_[junc.points[0].first].pos_norm(junc.points[0].second);
    is_head = junc.points[0].second < 0.5;
    si = 1;
    v_si = 0;
  }
  write_pi = si;
  si = junc.points[si].first;
  Float v_s = junc.points[v_si].second;
  v_si = junc.points[v_si].first;
  const auto &other_stroke = graph.orig_strokes_[si];
  other_stroke.ensure_arclengths();
  Vec2 _proj;
  Float s;
  const auto dist =
    ((junc.points[0].first == junc.points[1].first)
       ? closest_point_to_own_stroke(other_stroke, is_head, _proj, s)
       : closest_point(other_stroke, v, _proj, s));

  // Reproject and redo the line-of-sight check to account for deformed
  // geometry.
  auto cand = Junction{junc}; // Copy.
  const auto &strokes = graph.orig_strokes_;
  if (dist == std::numeric_limits<Float>::infinity()) {
    return false;
  }
  cand.points[write_pi].second = s / strokes[si].length();
  if (!line_of_sight(strokes[v_si].pos_norm(v_s), strokes[si].pos(s),
                     *graph.orig_bvh_)) {
    return false;
  }

  junc = cand;

  return true;
}

bool replay_deformations(span<const Stroke> strokes,
                         span<const Junction> candidates,
                         StrokeGraph &stroke_graph) {
  auto cache = IncrementalCache();
  stroke_graph = StrokeGraph(StrokeGraph::SnappingType::LinearSolve);
  std::vector<std::pair<size_t, size_t>> adj_faces;
  std::vector<Float> junc_distances;
  std::vector<StrokeGraph::VertexID> junc_vertex_ids;
  for (size_t i = 0; i < strokes.size(); ++i) {
    add_stroke_incremental_topological(stroke_graph, &cache, strokes[i], false);
    // Dissolve
    for (size_t vi = 0; vi < stroke_graph.vertices_.size(); ++vi) {
      const auto v = stroke_graph.vertex(vi);
      if (should_dissolve_vertex(v)) {
        dissolve_vertex(stroke_graph, v.index_);
      }
    }

    // Connect
    for (size_t j = 0; j < candidates.size(); ++j) {
      int max_sid =
        std::max(candidates[j].points[0].first, candidates[j].points[1].first);
      if (max_sid == i) {
        // Reproject
        auto junc = candidates[j];
        {
          auto test_junc = candidates[j];
          auto ok = original_stroke_to_stroke_indexing(stroke_graph, test_junc);
          if (!ok)
            continue;
        }
        // Don't return false directly if reproject_t_junc fails. Let
        // modify_graph determine that.
        if (junc.type == JunctionType::T)
          reproject_t_junc(stroke_graph, junc);
        std::vector<Junction> connected_junc{junc};
        auto modified_graph =
          modify_graph(stroke_graph, connected_junc, std::vector<bool>{true},
                       adj_faces, junc_distances, junc_vertex_ids);
        if (!modified_graph)
          return false;
        stroke_graph = std::move(*modified_graph);
      } else if (max_sid > i) {
        break;
      }
    }
  }

  return true;
}

bool replay_connection(const std::vector<Junction> &candidates,
                       const std::vector<bool> &junction_connected,
                       StrokeGraph &stroke_graph,
                       std::vector<Junction> &out_junctions) {
  out_junctions = candidates;
  std::vector<size_t> connected_indices, component_labels;
  group_connected_junctions(stroke_graph, candidates, junction_connected,
                            connected_indices, component_labels);

  std::map<size_t, std::vector<size_t>> group_indices;
  for (size_t i = 0; i < connected_indices.size(); ++i) {
    if (!group_indices.count(component_labels[i]))
      group_indices[component_labels[i]] = std::vector<size_t>();
    group_indices[component_labels[i]].emplace_back(connected_indices[i]);
  }

  // Connect by groups. Reproject T-junctions.
  for (const auto &[c, indices] : group_indices) {
    std::map<StrokeTime, StrokeTime> reprojections;
    for (const auto i : indices) {
      // Ignore weak junctions. They are added to complete the component.
      if (candidates[i].is_weak || candidates[i].type == JunctionType::R)
        continue;

      if (reprojections.count(candidates[i].points[0]))
        continue;

      Junction reproj_junc = candidates[i];
      bool ok = reproject_t_junc(stroke_graph, reproj_junc);
      force_assert(ok && "couldn't reproject T-junction");
      reprojections[candidates[i].points[0]] = reproj_junc.points[0];
    }

    // Update component
    std::vector<Junction> connected_junc;
    for (const auto i : indices) {
      connected_junc.emplace_back(candidates[i]);
      if (candidates[i].type == JunctionType::R)
        continue;
      assert(reprojections.count(connected_junc.back().points[0]));
      connected_junc.back().points[0] =
        reprojections[connected_junc.back().points[0]];
    }

    // Connect
    std::vector<bool> connected_flags;
    connected_flags.resize(connected_junc.size(), true);
    std::vector<std::pair<size_t, size_t>> adj_faces;
    std::vector<Float> junc_distances;
    std::vector<StrokeGraph::VertexID> junc_vertex_ids;
    auto modified_graph =
      modify_graph(stroke_graph, connected_junc, connected_flags, adj_faces,
                   junc_distances, junc_vertex_ids);
    if (!modified_graph)
      return false;
    stroke_graph = std::move(*modified_graph);

    for (size_t i = 0; i < indices.size(); ++i) {
      size_t cand_i = indices[i];
      out_junctions[cand_i].repr = connected_junc[i].repr;
    }
  }
  return true;
}

void color_endpoint_graph(const StrokeGraph &stroke_graph,
                          const span<Junction> candidates,
                          std::vector<int> &junc_components,
                          const Float min_prob) {
  // 1. Build an endpoint graph
  std::unordered_map<size_t, std::vector<double>> v2ends;
  std::vector<double> endpoint_times;
  endpoint_times.reserve(candidates.size() * 2);

  auto find_endpoint_time =
    [&endpoint_times](const double junc_hash) -> size_t {
    auto junc_itr = std::find_if(
      endpoint_times.begin(), endpoint_times.end(),
      [&junc_hash](const double t) { return std::abs(t - junc_hash) < 1e-10; });
    return junc_itr - endpoint_times.begin();
  };

  std::unordered_map<size_t, std::pair<int, double>> endpoint_indexing;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto &junc = candidates[i];
    for (const auto &end : junc.points) {
      StrokeTime s_end = end;
      auto ok = convert_orig2strokes(stroke_graph, s_end);
      force_assert(ok && "couldn't map from orig to strokes");

      double junc_hash = end.first * 10 + end.second;
      if (s_end.second == 0.0 || s_end.second == 1.0) {
        auto v = endpoint_to_vertex(stroke_graph,
                                    Endpoint(s_end.first, s_end.second == 0.0));
        if (!v2ends.count(v.index_))
          v2ends[v.index_] = std::vector<double>();
        if (std::find(v2ends[v.index_].begin(), v2ends[v.index_].end(),
                      junc_hash) == v2ends[v.index_].end())
          v2ends[v.index_].emplace_back(junc_hash);
      }
      size_t idx = endpoint_times.size();

      size_t end_i = find_endpoint_time(junc_hash);

      if (end_i == endpoint_times.size()) {
        endpoint_times.emplace_back(junc_hash);
        endpoint_indexing[idx] = end;
      }
    }
  }

  // 2. Label connected components
  std::vector<int> components(endpoint_times.size(), -1);
  int max_color = 0;

  auto merge_components = [&max_color, find_endpoint_time,
                           &components](double junc_hash1, double junc_hash2) {
    size_t end_idx1 = find_endpoint_time(junc_hash1);
    size_t end_idx2 = find_endpoint_time(junc_hash2);
    if (components[end_idx1] == -1 && components[end_idx2] == -1) {
      components[end_idx1] = components[end_idx2] = max_color++;
    } else if (components[end_idx1] != components[end_idx2]) {
      int min_c = std::min(components[end_idx1], components[end_idx2]);
      int max_c = std::max(components[end_idx1], components[end_idx2]);
      if (min_c == -1)
        components[end_idx1] = components[end_idx2] = max_c;
      else {
        for (auto &c : components)
          if (c == min_c)
            c = max_c;
      }
    }
  };

  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto &junc = candidates[i];

    // Ignore low probability if given a threshold
    if (junc.probability < min_prob)
      continue;

    // If given a threshold, also ignore the weak junctions (the ones used to
    // complete the component graph)
    if (min_prob > 0.0 && junc.is_weak)
      continue;

    double junc_hash1 = junc.points[0].first * 10 + junc.points[0].second;
    double junc_hash2 = junc.points[1].first * 10 + junc.points[1].second;
    merge_components(junc_hash1, junc_hash2);
  }

  // Now color using the snapped ends
  for (const auto &[v, ends] : v2ends) {
    for (size_t i = 0; i + 1 < ends.size(); ++i) {
      for (size_t j = i + 1; j < ends.size(); ++j) {
        double junc_hash1 = ends[i];
        double junc_hash2 = ends[j];
        merge_components(junc_hash1, junc_hash2);
      }
    }
  }

  junc_components.resize(candidates.size());
  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto &junc = candidates[i];
    double junc_hash1 = junc.points[0].first * 10 + junc.points[0].second;
    double junc_hash2 = junc.points[1].first * 10 + junc.points[1].second;
    size_t end_idx1 = find_endpoint_time(junc_hash1);
    size_t end_idx2 = find_endpoint_time(junc_hash2);
    assert(min_prob != 0.0 || components[end_idx1] == components[end_idx2]);
    if (min_prob != 0.0) {
      if (components[end_idx1] == components[end_idx2])
        junc_components[i] = components[end_idx1];
      else
        junc_components[i] = -1;
    } else
      junc_components[i] = components[end_idx1];
  }
}

void trim_overshoots(StrokeGraph &stroke_graph) {
  // Detect overshoots.
  stroke_graph.faces_.clear();
  for (auto &edge_rec : stroke_graph.hedges_) {
    edge_rec.face_ = (size_t)-1;
  }
  for (size_t vi = 0; vi < stroke_graph.vertices_.size(); ++vi) {
    auto v = stroke_graph.vertex(vi);
    if (v.is_active() && v.is_dangling()) {
      auto he = v.hedge();
      auto edge_length_to_intersection = he.stroke().length();
      while (he.dest().valence() == 2) {
        if (he.twin().continuity_edge().is_valid()) {
          he = he.twin().continuity_edge();
          edge_length_to_intersection += he.stroke().length();
        } else {
          break;
        }
      }
      const auto dest = he.dest();
      if (dest.is_dangling()) {
        continue; // If it's floating then it's not an overshoot.
      }

      auto end =
        Endpoint(v.hedge().stroke_idx(), v.hedge().forward()).as_pair();
      const auto ok = convert_strokes2orig(stroke_graph, end);
      force_assert(ok && "couldn't map from strokes to orig");

      const auto &orig_stroke = stroke_graph.orig_strokes_[end.first];
      const auto d2 = (v.pos() - dest.pos()).norm() + vertex_radius(v);

      if (d2 < 1.5 * 2 * vertex_radius(dest) &&
          edge_length_to_intersection <
            0.15 * std::max(orig_stroke.length(), 2 * vertex_radius(dest))) {

        // Trim the edge.
        drop_vertex_and_neighbors(stroke_graph, vi, span<Junction>());

        // Or if you want to mark as overshot vertex.
        // stroke_graph.vertices_[vi].flags_ |=
        // StrokeGraph::VertexRecord::Overshoot;
      }
    }
  }
  construct_faces(stroke_graph);
}

void mark_missed_overlapping_endpoints(StrokeGraph &stroke_graph) {
  auto orig_env_bvh = EnvelopeBVH(stroke_graph.orig_strokes_);
  for (size_t vi = 0; vi < stroke_graph.vertices_.size(); ++vi) {
    const auto vert = stroke_graph.vertex(vi);
    if (vert.is_active() && vert.is_dangling()) {
      auto orig_pos = StrokeTime((int)vert.hedge().stroke_idx(),
                                 vert.hedge().forward() ? 0.0 : 1.0);
      const auto ok = convert_strokes2orig(stroke_graph, orig_pos);
      force_assert(ok && "couldn't map from strokes to orig");
      force_assert(orig_pos.second == 0.0 || orig_pos.second == 1.0);

      for (size_t orig_si = 0; orig_si < stroke_graph.orig_strokes_.size();
           ++orig_si) {
        Float env_dist = INFINITY;
        const auto &other_orig_stroke = stroke_graph.orig_strokes_[orig_si];
        other_orig_stroke.ensure_arclengths();
        const auto &this_orig_stroke =
          stroke_graph.orig_strokes_[orig_pos.first];
        const auto [cen_dist, arclen] =
          find_snap_point(stroke_graph.orig_strokes_[orig_pos.first],
                          /*head=*/orig_pos.second == 0.0,
                          orig_env_bvh.nodes[orig_si], &env_dist);
        if (cen_dist < INFINITY) {
          assert(env_dist <= 0);

          // Find intersection point with the other.
          auto shares_intersection = false;
          auto he = vert.hedge();
          auto edge_length_to_intersection = he.stroke().length();
          auto edge_pen_width = he.stroke().pen_width();
          while (he.twin().continuity_edge().is_valid()) {
            const auto twin = he.twin();
            auto it = twin;
            auto valence = 0;
            do {
              if (it.orig_stroke_idx() == orig_si) {
                shares_intersection = true;
                break;
              }
              valence++;
              assert(valence < 1024 && "likely infinite loop found");
              it = it.twin().next();
            } while (it != twin);
            if (shares_intersection)
              break;

            he = he.twin().continuity_edge();
            edge_length_to_intersection += he.stroke().length();
            edge_pen_width = std::max(edge_pen_width, he.stroke().pen_width());
          }
          {
            const auto twin = he.twin();
            auto it = twin;
            auto valence = 0;
            do {
              if (it.orig_stroke_idx() == orig_si) {
                shares_intersection = true;
                break;
              }
              valence++;
              assert(valence < 1024 && "likely infinite loop found");
              it = it.twin().next();
            } while (it != twin);
          }

          // If snap point is close to other end, must be non-spurious.
          if (!shares_intersection ||
              // Keep this synced with others (grep for [spurious-condition]).
              std::max(3 * cen_dist, M_PI * this_orig_stroke.pen_width()) <
                edge_length_to_intersection) {

            // This vertex failed to snap, though it should.  Mark it.
            stroke_graph.vertices_[vi].flags_ |=
              StrokeGraph::VertexRecord::Overlapping;
            break;
          }
        } else {
          // Check for self-loop.
          const auto loop_cen_dist =
            (this_orig_stroke.xy(0) - this_orig_stroke.xy(Back)).norm();
          const auto loop_env_dist =
            loop_cen_dist -
            0.5 * (this_orig_stroke.width(0) + this_orig_stroke.width(Back));
          if (loop_env_dist <= 0 &&
              // This is a bit weaker than [spurious-condition].
              3 * loop_cen_dist < this_orig_stroke.length()) {
            // This vertex failed to snap, though it should.  Mark it.
            stroke_graph.vertices_[vi].flags_ |=
              StrokeGraph::VertexRecord::Overlapping;
            break;
          }
        }
      }
    }
  }
}

void fix_up_corner_positions(StrokeGraph &stroke_graph) {
  // Other snapping method types not supported.
  force_assert(stroke_graph.snapping_method_type_ == StrokeGraph::Connection);

  for (size_t vi = 0; vi < stroke_graph.vertices_.size(); ++vi) {
    const auto v = stroke_graph.vertex(vi);
    if (!v.is_active()) {
      continue;
    }

    if (v.valence() == 2) {
      // Try to find this snap in the snap history.
      auto *snap_rec = (const ClassifierPrediction *)nullptr;
      for (const auto &snap : stroke_graph.snap_history_) {
        assert(snap.key.cand1 !=
               (int)v.index_); // v.index_ is always going to be cand2.
        if (snap.key.type == JunctionType::R &&
            snap.key.cand2 == (int)v.index_) {
          snap_rec = &snap;
          break;
        }
      }
      if (snap_rec &&
          (snap_rec->orig_a.second == 0.0 || snap_rec->orig_a.second == 1.0) &&
          (snap_rec->orig_b.second == 0.0 || snap_rec->orig_b.second == 1.0)) {
        const auto &orig_stroke1 =
          stroke_graph.orig_strokes_[snap_rec->orig_a.first];
        const auto &orig_stroke2 =
          stroke_graph.orig_strokes_[snap_rec->orig_b.first];
        const auto is_head1 = (snap_rec->orig_a.second == 0.0);
        const auto is_head2 = (snap_rec->orig_b.second == 0.0);
        const auto tangent1 = (is_head1 ? head_tangent_lagrange(orig_stroke1)
                                        : tail_tangent_lagrange(orig_stroke1));
        const auto tangent2 = (is_head2 ? head_tangent_lagrange(orig_stroke2)
                                        : tail_tangent_lagrange(orig_stroke2));
        const auto connection = (snap_rec->p_b - snap_rec->p_a).normalized();
        // Turn angles.
        const auto angle1 =
          std::acos(std::clamp(tangent1.dot(connection), 0.0, 1.0));
        const auto angle2 =
          std::acos(std::clamp(tangent2.dot(-connection), 0.0, 1.0));
        if (angle2 < angle1) {
          assert(
            (stroke_graph.vertices_[vi].p_ - snap_rec->p_b).squaredNorm() <=
            (stroke_graph.vertices_[vi].p_ - snap_rec->p_a).squaredNorm());
          auto end1 = snap_rec->orig_a;
          auto ok = convert_orig2strokes(stroke_graph, end1);
          force_assert(ok && "couldn't map from orig to strokes");
          auto end2 = snap_rec->orig_b;
          ok = convert_orig2strokes(stroke_graph, end2);
          force_assert(ok && "couldn't map from orig to strokes");
          auto &edge1 = stroke_graph.strokes_[end1.first];
          auto &edge2 = stroke_graph.strokes_[end2.first];
          // If edge1's size is <= 2, then some other connection got made in the
          // middle of the connection. In this case, we can't fix up the vertex
          // location without majorly reworking the topology.
          if (edge1.size() > 2) {
            // We need to switch the vertex location.
            stroke_graph.vertices_[vi].p_ = snap_rec->p_a;
            if (end1.second == 0.0) {
              copy_run(edge1, 0, edge1.size_ - 2, edge1, 1, edge1.size_ - 1);
            }
            edge1.size_--;
            if (end2.second == 0.0) {
              edge2.insert(0, v.pos().x_, v.pos().y_, edge2.width(0));
            } else {
              edge2.push_back(v.pos().x_, v.pos().y_, edge2.width(Back));
            }
            stroke_graph.bvh_.full_update(end1.first);
            stroke_graph.bvh_.full_update(end2.first);
          }
        }
      }
    }
  }
}

void build_plane_graph(span<const Stroke> strokes, StrokeGraph &stroke_graph) {
  stroke_graph = StrokeGraph(StrokeGraph::SnappingType::Connection);
  IncrementalCache cache;
  cache.trim_overshoots_ = false;
  for (size_t i = 0; i < strokes.size(); ++i) {
    add_stroke_incremental_topological(stroke_graph, &cache, strokes[i], false);
  }
  stroke_graph.orig_bvh_ =
    std::make_unique<PolylineBVH>(stroke_graph.orig_strokes_);

  trim_overshoots(stroke_graph);
  mark_missed_overlapping_endpoints(stroke_graph);
  fix_up_corner_positions(stroke_graph);

#ifndef NDEBUG
  check_consistency(stroke_graph);
#endif
}

StrokeSnapInfo increment_strokes(StrokeGraph &graph, span<const Stroke> strokes,
                                 size_t start_idx, size_t end_idx) {
  force_assert(end_idx <= strokes.size());
  IncrementalCache cache;
  cache.trim_overshoots_ = false;
  StrokeSnapInfo snap_info;
  for (size_t i = start_idx; i < end_idx; ++i) {
    snap_info =
      add_stroke_incremental_topological(graph, &cache, strokes[i], false);
  }
  graph.orig_bvh_ = std::make_unique<PolylineBVH>(graph.orig_strokes_);
  trim_overshoots(graph);
  mark_missed_overlapping_endpoints(graph);
  fix_up_corner_positions(graph);

#ifndef NDEBUG
  check_consistency(graph);
#endif

  return snap_info;
}

static void vanilla_candidates(const StrokeGraph &end_stroke_graph,
                               const StrokeGraph &stroke_graph,
                               const FeatureType feature_type,
                               std::vector<Junction> &dangling_predictions,
                               std::vector<Junction> &corner_predictions,
                               bool train_time, bool generate_corners = false,
                               int num_cand = -1) {
  dangling_predictions.clear();
  corner_predictions.clear();

  size_t van_num_cand = (num_cand < 0) ? snap_candidate_count : num_cand;

  IncrementalCache cache;
  std::vector<Junction> tmp_corner_predictions;
  // 1. Generate corner (high-valence) candidates
  if (generate_corners) {
    std::vector<SnapInfo> stroke_pos_candidates;
    for (size_t vi = 0; vi < stroke_graph.vertices_.size(); ++vi) {
      if (!stroke_graph.vertex(vi).is_active() ||
          stroke_graph.vertex(vi).valence() != 1)
        continue;

      size_t orig_sid = stroke_graph.vertex(vi).hedge().orig_stroke_idx();

      // Check if this vertex is dangling in the end_stroke_graph
      bool is_end_dangling = false;
      {
        auto end = vertex_to_endpoint(stroke_graph.vertex(vi));
        StrokeTime st(end.stroke_idx(), end.is_tail());
        if (convert_strokes2orig(stroke_graph, st)) {
          if (convert_orig2strokes(end_stroke_graph, st) &&
              (std::abs(st.second) < 1e-10 ||
               std::abs(st.second - 1.0) < 1e-10)) {
            auto end_v = endpoint_to_vertex(
              end_stroke_graph, Endpoint(st.first, st.second < 0.5));
            if (end_v.is_active() && end_v.is_dangling()) {
              is_end_dangling = true;
            }
          }
        }
      }

      if (!is_end_dangling)
        continue;

      snap_corner_candidates(stroke_graph, vi, &cache, feature_type,
                             stroke_pos_candidates, van_num_cand);
    }
    tmp_corner_predictions.reserve(stroke_pos_candidates.size());
    for (size_t i = 0; i < stroke_pos_candidates.size(); ++i) {
      const auto &snap = stroke_pos_candidates[i];
      if (snap.prediction_type != JunctionType::T)
        continue;
      Junction junc(
        {StrokeTime((int)snap.stroke_pos2.first, snap.stroke_pos2.second),
         StrokeTime((int)snap.stroke_pos1.first, snap.stroke_pos1.second)},
        snap.prediction_type, false, 0.5, true);
      auto ok = stroke_to_original_stroke_indexing(stroke_graph, junc);
      force_assert(ok && "couldn't map from strokes to orig");

      junc.sort_entries();

      Junction s_junc = junc;
      ok = original_stroke_to_stroke_indexing(stroke_graph, s_junc);
      force_assert(ok && "couldn't map from orig to strokes");

      // Round the end to the closest vertex to check the valence
      Junction junc_round = s_junc;
      junc_round.points[0].second =
        (junc_round.points[0].second < 0.5) ? 0.0 : 1.0;
      ok = stroke_to_original_stroke_indexing(stroke_graph, junc_round);
      force_assert(ok && "couldn't map from strokes to orig");

      // Swap if sorting doesn't work
      {
        int end_d0 = end_degree(stroke_graph, junc_round.points[0]),
            end_d1 = end_degree(stroke_graph, junc_round.points[1]);
        if (end_d1 != 1) {
          std::swap(junc.points[0], junc.points[1]);
          std::swap(s_junc.points[0], s_junc.points[1]);
          std::swap(junc_round.points[0], junc_round.points[1]);
        }
      }

      // Try two different rounding
      auto &s = stroke_graph.strokes_[s_junc.points[0].first];
      Junction alt_junc_round = s_junc;
      alt_junc_round.points[0].second =
        (alt_junc_round.points[0].second < 0.5) ? 1.0 : 0.0;
      ok = stroke_to_original_stroke_indexing(stroke_graph, alt_junc_round);

      Float edge_length = s.length();
      Float orig_length =
        stroke_graph.orig_strokes_[junc_round.points[0].first].length() *
        std::abs(alt_junc_round.points[0].second - junc_round.points[0].second);
      if (!trial_add_corner_junction(end_stroke_graph, stroke_graph, junc_round,
                                     junc, s_junc, tmp_corner_predictions) &&
          edge_length > 1.5 * orig_length) {
        junc_round = s_junc;
        junc_round.points[0].second =
          (junc_round.points[0].second < 0.5) ? 1.0 : 0.0;
        ok = stroke_to_original_stroke_indexing(stroke_graph, junc_round);
        trial_add_corner_junction(end_stroke_graph, stroke_graph, junc_round,
                                  junc, s_junc, tmp_corner_predictions);
      }
    }
  }

  // 2. Generate other candidates
  {
    std::vector<SnapInfo> stroke_pos_candidates;
    for (size_t vi = 0; vi < stroke_graph.vertices_.size(); ++vi) {
      if (!stroke_graph.vertex(vi).is_active())
        continue;

      // Check if a stroke endpoint is in its neighborhood
      bool seen_endpoint = stroke_graph.vertex(vi).valence() == 1;
      if (!seen_endpoint) {
        const auto he = stroke_graph.vertex(vi).hedge();
        auto it = he;
        do {
          auto end = StrokeTime((int)it.stroke_idx(), it.forward() ? 0.0 : 1.0);
          if (convert_strokes2orig(stroke_graph, end) &&
              (end.second == 0.0 || end.second == 1.0)) {
            seen_endpoint = true;
            break;
          }
          it = it.twin().next();
        } while (it != he);
      }

      if (!seen_endpoint)
        continue;

      // Check if this vertex is dangling in the end_stroke_graph
      bool is_end_dangling = false;
      {
        auto end = vertex_to_endpoint(stroke_graph.vertex(vi));
        StrokeTime st(end.stroke_idx(), end.is_tail());
        if (convert_strokes2orig(stroke_graph, st)) {
          if (convert_orig2strokes(end_stroke_graph, st) &&
              (std::abs(st.second) < 1e-10 ||
               std::abs(st.second - 1.0) < 1e-10)) {
            auto end_v = endpoint_to_vertex(
              end_stroke_graph, Endpoint(st.first, st.second < 0.5));
            if (end_v.is_active() && end_v.is_dangling()) {
              is_end_dangling = true;
            }
          }
        }
      }

      if (!is_end_dangling)
        continue;

      std::vector<SnapInfo> tmp_stroke_pos_candidates;
      snap_candidates(stroke_graph, vi, &cache, feature_type,
                      tmp_stroke_pos_candidates, van_num_cand, false,
                      train_time);
      for (const auto &junc : tmp_stroke_pos_candidates) {
        if (stroke_graph.vertex(vi).valence() == 1 ||
            junc.prediction_type == JunctionType::R)
          stroke_pos_candidates.emplace_back(junc);
      }
    }

    dangling_predictions.reserve(stroke_pos_candidates.size());
    for (const auto &snap : stroke_pos_candidates) {
      Junction junc(
        {StrokeTime((int)snap.stroke_pos2.first, snap.stroke_pos2.second),
         StrokeTime((int)snap.stroke_pos1.first, snap.stroke_pos1.second)},
        snap.prediction_type, false, 0.5, true);
      const auto v1 =
        endpoint_to_vertex(stroke_graph, Endpoint(junc.points[0].first,
                                                  junc.points[0].second == 0));
      const auto v2 =
        endpoint_to_vertex(stroke_graph, Endpoint(junc.points[1].first,
                                                  junc.points[1].second == 0));
      const auto ok = stroke_to_original_stroke_indexing(stroke_graph, junc);
      force_assert(ok && "couldn't map from strokes to orig");

      junc.sort_entries();
      // Disable high valence
      int end_d0 = end_degree(stroke_graph, junc.points[0]),
          end_d1 = end_degree(stroke_graph, junc.points[1]);
      size_t orig_sid =
        (end_d1 == 1) ? junc.points[1].first : junc.points[0].first;

      // Skip the T-junction that can potentially at a bridge vertex
      if (junc.type == JunctionType::T && std::max(end_d0, end_d1) > 1)
        continue;

      // Skip the overlapping endpoints
      if ((v1.flags() & StrokeGraph::VertexRecord::Overlapping) ||
          (v2.flags() & StrokeGraph::VertexRecord::Overlapping))
        continue;

      size_t dangling_predictions_size = dangling_predictions.size();

      bool unseen_junc =
        std::find_if(dangling_predictions.begin(), dangling_predictions.end(),
                     [&junc](const Junction &j) {
                       return (j.type == junc.type) &&
                              (j.points[0].first == junc.points[0].first) &&
                              (j.points[0].second == junc.points[0].second) &&
                              (j.points[1].first == junc.points[1].first) &&
                              (j.points[1].second == junc.points[1].second);
                     }) == dangling_predictions.end() &&
        std::find_if(tmp_corner_predictions.begin(),
                     tmp_corner_predictions.end(), [&junc](const Junction &j) {
                       return (j.type == junc.type) &&
                              (j.points[0].first == junc.points[0].first) &&
                              (j.points[0].second == junc.points[0].second) &&
                              (j.points[1].first == junc.points[1].first) &&
                              (j.points[1].second == junc.points[1].second);
                     }) == tmp_corner_predictions.end();
      if (unseen_junc && (end_d0 >= 0 && end_d0 <= 1) &&
          (end_d1 >= 0 && end_d1 <= 1) && (end_d0 == 1 || end_d1 == 1) &&
          graph_visible(end_stroke_graph, junc)) {
        junc.orig_dist = junction_distance_init(stroke_graph, junc);
        dangling_predictions.emplace_back(junc);
      } else if (generate_corners && unseen_junc &&
                 (end_d0 > 1 || end_d1 > 1) && (end_d0 == 1 || end_d1 == 1) &&
                 graph_visible(end_stroke_graph, junc)) {
        if (end_d1 > 1)
          std::swap(junc.points[0], junc.points[1]);

        add_e2e_corner_junction(end_stroke_graph, stroke_graph, junc,
                                tmp_corner_predictions);
      }
    }
  }

  // 3. Move the corner junction position once deduplication is done
  std::set<std::pair<int, int>> seen_end_pairs;
  for (auto &junc : tmp_corner_predictions) {
    if (junc.type == JunctionType::R)
      continue;

    auto ok = original_stroke_to_stroke_indexing(stroke_graph, junc);
    force_assert(ok && "couldn't map from orig to strokes");
    junc.type = JunctionType::R;
    ok = stroke_to_original_stroke_indexing(stroke_graph, junc);
    force_assert(ok && "couldn't map from strokes to orig");
  }
  auto to_v_idx = [&stroke_graph](const StrokeTime &end) -> int {
    StrokeTime s_end = end;
    auto ok = convert_orig2strokes(stroke_graph, s_end);
    // This original stroke is entirely deleted in the final graph
    assert(ok);
    assert(s_end.second == 0.0 || s_end.second == 1.0);
    const auto v = endpoint_to_vertex(
      stroke_graph, Endpoint(s_end.first, s_end.second == 0.0));
    assert(v.is_valid() && v.is_active());
    return v.index_;
  };
  for (auto &junc : tmp_corner_predictions) {
    size_t vid1 = to_v_idx(junc.points[0]);
    size_t vid2 = to_v_idx(junc.points[1]);
    const auto junc_p =
      std::make_pair(std::min(vid1, vid2), std::max(vid1, vid2));
    if (seen_end_pairs.count(junc_p) != 0)
      continue;

    auto v1 = stroke_graph.vertex(vid1);
    auto v2 = stroke_graph.vertex(vid2);
    if (v1.valence() == 1)
      std::swap(v1, v2);

    // Check if it's a closed stroke
    if (v1.valence() == 2) {
      std::unordered_set<size_t> adj_sid;
      const auto he = v1.hedge();
      auto it = he;
      do {
        adj_sid.emplace(it.stroke_idx());
        it = it.twin().next();
      } while (it != he);

      if (adj_sid.size() == 1)
        continue;
    }

    if (!is_snap_spurious(stroke_graph, v1.index_, v2.index_)) {
      corner_predictions.emplace_back(junc);
      seen_end_pairs.emplace(junc_p);
    }
  }
}

void vanilla_candidates(const StrokeGraph &stroke_graph,
                        std::vector<Junction> &dangling_predictions,
                        bool train_time, int num_cand) {
  std::vector<Junction> corner_predictions;
  vanilla_candidates(stroke_graph, stroke_graph, FeatureType::OrigStroke,
                     dangling_predictions, corner_predictions, train_time,
                     false, num_cand);
}

void vanilla_candidates(const StrokeGraph &end_stroke_graph,
                        const StrokeGraph &stroke_graph,
                        std::vector<Junction> &dangling_predictions,
                        std::vector<Junction> &corner_predictions,
                        bool train_time) {
  vanilla_candidates(end_stroke_graph, stroke_graph, FeatureType::OrigStroke,
                     dangling_predictions, corner_predictions, train_time,
                     false, -1);
}

void predicted_corner_candidates(const StrokeGraph &end_stroke_graph,
                                 const StrokeGraph &stroke_graph,
                                 const FeatureType feature_type,
                                 std::vector<Junction> &corner_predictions,
                                 bool to_dedup, bool include_prev_connections) {
  prediction_feature_type = feature_type;

  // 1. Generate the prototype corner candidates
  std::vector<Junction> dangling_predictions, proto_corner_predictions;
  vanilla_candidates(end_stroke_graph, stroke_graph, FeatureType::OrigStroke,
                     dangling_predictions, proto_corner_predictions, false,
                     true, -1);

  // Consider the connections in the intermediate graph for the corner
  // generation
  if (include_prev_connections) {
    std::set<std::pair<int, int>> high_valence_junction_dedup;
    for (size_t i = 0; i < proto_corner_predictions.size(); ++i) {
      const auto &junc = proto_corner_predictions[i];
      const auto junc_p =
        std::make_pair(std::min(to_v_idx(end_stroke_graph, junc.points[0]),
                                to_v_idx(end_stroke_graph, junc.points[1])),
                       std::max(to_v_idx(end_stroke_graph, junc.points[0]),
                                to_v_idx(end_stroke_graph, junc.points[1])));

      if (high_valence_junction_dedup.count(junc_p) == 0) {
        high_valence_junction_dedup.emplace(junc_p);
      }
    }

    std::vector<Junction> conn_proto_corner_predictions;
    dangling_predictions.clear();
    vanilla_candidates(end_stroke_graph, end_stroke_graph,
                       FeatureType::OrigStroke, dangling_predictions,
                       conn_proto_corner_predictions, false, true, -1);

    // Add to the result set
    std::vector<Junction> connection_corner_candidates;
    for (auto &junc : conn_proto_corner_predictions) {
      const auto junc_p =
        std::make_pair(std::min(to_v_idx(end_stroke_graph, junc.points[0]),
                                to_v_idx(end_stroke_graph, junc.points[1])),
                       std::max(to_v_idx(end_stroke_graph, junc.points[0]),
                                to_v_idx(end_stroke_graph, junc.points[1])));
      if (high_valence_junction_dedup.count(junc_p))
        continue;
      junc.to_reconsider = true;
      if (std::find(proto_corner_predictions.begin(),
                    proto_corner_predictions.end(),
                    junc) == proto_corner_predictions.end() &&
          std::find(connection_corner_candidates.begin(),
                    connection_corner_candidates.end(),
                    junc) == connection_corner_candidates.end()) {
        bool connection_check =
          connection_corner_visible(end_stroke_graph, junc);
        if (connection_check)
          connection_corner_candidates.emplace_back(junc);
      }
    }

    // Debug:
    proto_corner_predictions.insert(proto_corner_predictions.end(),
                                    connection_corner_candidates.begin(),
                                    connection_corner_candidates.end());
    // proto_corner_predictions = std::move(connection_corner_candidates);
    //
  }

  // 2. Expand
  expand_corner_candidates((include_prev_connections) ? end_stroke_graph
                                                      : stroke_graph,
                           proto_corner_predictions, corner_predictions);

  // 3. Predict
  update_high_valence_junction_predictions(
    (include_prev_connections) ? end_stroke_graph : stroke_graph,
    corner_predictions, prediction_feature_type, to_dedup);
}

void complete_predicted_corner_candidates(
  const StrokeGraph &end_stroke_graph, const StrokeGraph &stroke_graph,
  const FeatureType feature_type, std::vector<Junction> &corner_predictions,
  Float min_probability, bool include_prev_connections) {
  predicted_corner_candidates(end_stroke_graph, stroke_graph, feature_type,
                              corner_predictions, false,
                              include_prev_connections);

  // 4. Filter
  if (min_probability >= 0) {
    std::vector<Junction> filtered_corner_predictions;
    for (const auto &junc : corner_predictions) {
      if (junc.probability <= min_probability)
        continue;
      filtered_corner_predictions.emplace_back(junc);
    }
    corner_predictions = std::move(filtered_corner_predictions);
  }

  print_junctions(end_stroke_graph, corner_predictions);

  // 5. Complete
  {
    std::vector<Junction> complementary_new_candidates;
    std::vector<Junction> to_complete_candidates;
    to_complete_candidates.reserve(corner_predictions.size());
    to_complete_candidates.insert(to_complete_candidates.end(),
                                  corner_predictions.begin(),
                                  corner_predictions.end());
    complete_graph_candidates(
      end_stroke_graph, stroke_graph.orig_strokes_, stroke_graph,
      stroke_graph.orig_strokes_.size() - 1, to_complete_candidates,
      complementary_new_candidates);
    update_disconnected_junction_predictions(
      end_stroke_graph, complementary_new_candidates, feature_type, false);
    for (auto &junc : complementary_new_candidates) {
      if (junc.probability == 0.0)
        junc.must_disconnect = true;
      bool unseen_junc =
        std::find_if(corner_predictions.begin(), corner_predictions.end(),
                     [&junc](const Junction &j) {
                       return (j.points[0].first == junc.points[0].first) &&
                              (j.points[0].second == junc.points[0].second) &&
                              (j.points[1].first == junc.points[1].first) &&
                              (j.points[1].second == junc.points[1].second);
                     }) == corner_predictions.end();
      if (unseen_junc)
        corner_predictions.emplace_back(junc);
    }
  }

  print_junctions(end_stroke_graph, corner_predictions);
}

bool is_corner_candidate(const StrokeGraph &stroke_graph,
                         const Junction &junc) {
  StrokeTime s_end = junc.points[0];
  auto ok = convert_orig2strokes(stroke_graph, s_end);
  // This original stroke is entirely deleted in the final graph
  assert(ok);
  if (std::abs(s_end.second) >= 1e-10 && std::abs(1.0 - s_end.second) >= 1e-10)
    return false;
  s_end.second = (s_end.second < 0.5) ? 0.0 : 1.0;

  auto to_v_idx = [&stroke_graph](const StrokeTime &end) -> int {
    StrokeTime s_end = end;
    auto ok = convert_orig2strokes(stroke_graph, s_end);
    // This original stroke is entirely deleted in the final graph
    assert(ok);
    assert(s_end.second == 0.0 || s_end.second == 1.0);
    const auto v = endpoint_to_vertex(
      stroke_graph, Endpoint(s_end.first, s_end.second == 0.0));
    assert(v.is_valid() && v.is_active());
    return v.index_;
  };
  return (stroke_graph.vertex(to_v_idx(junc.points[0])).valence() > 1);
}

void expand_corner_candidates(
  const StrokeGraph &stroke_graph,
  const std::vector<Junction> &corner_predictions,
  std::vector<Junction> &expanded_corner_predictions) {
  expanded_corner_predictions.clear();
  for (const auto &junc : corner_predictions) {
    size_t before_count = expanded_corner_predictions.size();
    std::vector<Junction> expanded_juncs;
    expand_adjacent_junction(stroke_graph, junc, expanded_juncs);
    // Assume all inputs are indeed corner junctions
    assert(
      (((expanded_juncs.empty()) ? junc.type : expanded_juncs.front().type) ==
         JunctionType::R &&
       (expanded_juncs.empty() || expanded_juncs.size() > 1)) ||
      (((expanded_juncs.empty()) ? junc.type : expanded_juncs.front().type) ==
         JunctionType::T &&
       (expanded_juncs.empty() || expanded_juncs.size() == 1)));
    expanded_corner_predictions.insert(expanded_corner_predictions.end(),
                                       expanded_juncs.begin(),
                                       expanded_juncs.end());
  }
}

void connect_graph(StrokeGraph &stroke_graph,
                   std::vector<Junction> &candidates) {
  // Connect. Note that we are hard-coded to use straight line connection here.
  // All connections shouldn't introduce new intersections since they already
  // passed the visibility check.
  for (auto &f_junc : candidates) {
    Junction junc = f_junc;
    bool to_swap = junc.points[0].first < junc.points[1].first;
    if (!original_stroke_to_stroke_indexing(stroke_graph, junc)) {
      goto failure;
    }

    {
      StrokeGraph::VertexView v1, v2;
      StrokeGraph varying_graph = stroke_graph.clone();
      varying_graph.strokes_[junc.points[0].first].ensure_arclengths();
      varying_graph.strokes_[junc.points[1].first].ensure_arclengths();
      add_vertex(varying_graph, junc.points[0].first,
                 junc.points[0].second *
                   varying_graph.strokes_[junc.points[0].first].length(),
                 v1);
      add_vertex(varying_graph, junc.points[1].first,
                 junc.points[1].second *
                   varying_graph.strokes_[junc.points[1].first].length(),
                 v2);
      if (to_swap)
        std::swap(v1, v2);

      // Snap the edge.
      if (v1 != v2 &&
          // !is_snap_spurious(varying_graph, v1.index_, v2.index_) &&
          graph_visible(varying_graph, f_junc) &&
          varying_graph.snap_vertices(v1.index_, v2.index_)) {
        // Find the next valid junction index
        int new_id = -1;
        for (size_t j = 0; j < varying_graph.vertices_.size(); ++j) {
          for (auto const &vid : varying_graph.vertices_[j].ids_) {
            if (vid.connection_type_ == StrokeGraph::VertexID::Junction)
              new_id = std::max(new_id, (int)vid.connection_index_);
          }
        }
        const auto &vid = varying_graph.vertices_[v2.index_].ids_.emplace_back(
          StrokeGraph::VertexID{StrokeGraph::VertexID::Junction,
                                (size_t)++new_id});
        f_junc.repr = vid.repr();

        stroke_graph = std::move(varying_graph);
        v1 = stroke_graph.vertex(v1.index_);
        v2 = stroke_graph.vertex(v2.index_);
      } else
        goto failure;
    }

    continue;

failure:
    std::stringstream ss;
    ss << std::setprecision(2) << "(" << f_junc.points[0].first << ", "
       << f_junc.points[0].second << " - " << f_junc.points[1].first << ", "
       << f_junc.points[1].second << ")";
    // SPDLOG_ERROR("Cannot snap: {}", ss.str());
  }
}

/// Visualization
void junc2pred(
  const StrokeGraph &graph, const std::vector<Junction> &junctions,
  std::vector<std::pair<ClassifierPrediction, VizProbabilities>> &predictions,
  const std::vector<Junction> &future_junctions) {
  predictions.reserve(junctions.size());
  for (size_t i = 0; i < junctions.size(); ++i) {
    const auto &junc = junctions[i];
    if (junc.type == JunctionType::R) {
      continue;
    }
    predictions.emplace_back();
    ClassifierPrediction &pred = predictions.back().first;
    pred.orig_a = junc.points[1];
    pred.orig_b = junc.points[0];
    pred.key.type = junc.type;
    pred.key.corner_type = junc.corner_type;
    pred.key.cand1 = (int)orig2endpoint(graph, junc.points[1]).index_;
    auto strokes_junc = junc;
    const auto ok = original_stroke_to_stroke_indexing(graph, strokes_junc);
    if (!ok)
      continue;
    pred.key.cand2 = strokes_junc.points[0].first;

    pred.p_a = graph.strokes_[strokes_junc.points[1].first].pos_norm(
      strokes_junc.points[1].second);
    pred.p_b = graph.strokes_[strokes_junc.points[0].first].pos_norm(
      strokes_junc.points[0].second);
    pred.prob = junc.probability;
    pred.fea = junc.fea;
    pred.alt_prob = junc.orig_dist;
    predictions.back().second.prob_ = junc.probability;
    predictions.back().second.repr = junc.repr;
    pred.junction_repr = junc.repr;
    if (!future_junctions.empty()) {
      predictions.back().second.future_prob_ =
        get_final_prediction(junc, future_junctions);
    }

    Float env_dist = junction_distance(graph, junc);
    predictions.back().second.env_dist_ = env_dist;
  }
  for (size_t i = 0; i < junctions.size(); ++i) {
    const auto &junc = junctions[i];
    if (junc.type != JunctionType::R) {
      continue;
    }
    predictions.emplace_back();
    ClassifierPrediction &pred = predictions.back().first;
    pred.orig_a = junc.points[1];
    pred.orig_b = junc.points[0];
    pred.key.type = junc.type;
    pred.key.corner_type = junc.corner_type;
    pred.key.cand1 = (int)orig2endpoint(graph, junc.points[1]).index_;
    pred.key.cand2 = (int)orig2endpoint(graph, junc.points[0]).index_;
    auto strokes_junc = junc;
    const auto ok = original_stroke_to_stroke_indexing(graph, strokes_junc);
    if (!ok)
      continue;
    pred.p_a = graph.strokes_[strokes_junc.points[1].first].pos_norm(
      strokes_junc.points[1].second);
    pred.p_b = graph.strokes_[strokes_junc.points[0].first].pos_norm(
      strokes_junc.points[0].second);
    pred.prob = junc.probability;
    pred.fea = junc.fea;
    pred.alt_prob = junc.orig_dist;
    predictions.back().second.prob_ = junc.probability;
    predictions.back().second.repr = junc.repr;
    pred.junction_repr = junc.repr;
    if (!future_junctions.empty()) {
      predictions.back().second.future_prob_ =
        get_final_prediction(junc, future_junctions);
    }

    Float env_dist = junction_distance(graph, junc);
    predictions.back().second.env_dist_ = env_dist;
  }
}

Junction pred2junc(const ClassifierPrediction &pred) {
  auto out = Junction(pred.key.type);
  out.points.emplace_back(pred.orig_a);
  out.points.emplace_back(pred.orig_b);
  out.probability = pred.prob;
  out.repr = pred.junction_repr;
  return out;
}

void viz_connections(const GraphState &state,
                     const std::vector<Junction> &future_junctions,
                     std::vector<std::string> &connected_junc_strs) {
  connected_junc_strs.reserve(state.candidates_.size());
  for (size_t i = 0; i < state.candidates_.size(); ++i) {
    if (!state.connectivity_[i])
      continue;
    assert(!state.candidates_[i].repr.empty() &&
           state.junc_distance_map_.count(state.candidates_[i].repr));
    std::stringstream junc_ss;
    junc_ss << std::setprecision(2) << "("
            << state.candidates_[i].points[0].first << ", "
            << state.candidates_[i].points[0].second << " - "
            << state.candidates_[i].points[1].first << ", "
            << state.candidates_[i].points[1].second << ")";

    std::stringstream ss;
    ss << state.candidates_[i].repr << ": " << junc_ss.str() << "\t"
       << ((state.candidates_[i].type == JunctionType::R) ? "R" : "T") << "\t"
       << state.candidates_[i].probability << "\t"
       << get_final_prediction(state.candidates_[i], future_junctions) << "\t"
       << state.junc_distance_map_.at(state.candidates_[i].repr) << std::endl;
    connected_junc_strs.emplace_back(ss.str());
  }
}

void viz_connected_junctions(
  const StrokeGraph &graph,
  std::vector<std::pair<std::string, Vec2>> &junctions) {
  // Dissolved junctions on edge
  for (size_t i = 0; i < graph.strokes_.size(); ++i)
    for (const auto &[t, vid] : graph.strokes2vid_[i]) {
      if (vid.repr().empty() || vid.repr()[0] == 'i')
        continue;
      if (graph.strokes_[i].size() > 0) {
        Vec2 vec(graph.strokes_[i].pos_norm(t).x(),
                 graph.strokes_[i].pos_norm(t).y());
        junctions.emplace_back(vid.repr(), vec);
      } else {
        SPDLOG_WARN("viz_connected_junctions: Empty stroke {}", i);
      }
    }

  // Vertices
  for (size_t i = 0; i < graph.vertices_.size(); ++i) {
    bool init_only = true;
    std::string vid_str;
    for (auto const &vid : graph.vertex(i).vertex_ids()) {
      vid_str += vid.repr() + ";";
      if (vid.connection_type_ == StrokeGraph::VertexID::Type::Junction)
        init_only = false;
    }

    if (init_only) {
      continue;
    }

    Vec2 vec(graph.vertex(i).pos().x(), graph.vertex(i).pos().y());
    junctions.emplace_back(vid_str, vec);
  }
}

} // namespace sketching
