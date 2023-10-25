// Include this header file for the functions in this file.
#include "stroke_graph_extra.h"

#include "detail/util.h"
#include "force_assert.h"
#include "incremental.h" // should_dissolve_vertex
#include "resample.h" // remove_duplicate_vertices

namespace sketching {

using HedgeView = StrokeGraph::HedgeView;
using VertexView = StrokeGraph::VertexView;

static constexpr auto invalid = StrokeGraph::invalid;
static constexpr auto NaN = std::numeric_limits<Float>::quiet_NaN();

/** Directly connect hedge hi's prev and next. */
static void skip_over(StrokeGraph& graph, const size_t hi) {
  const auto he = graph.hedge(hi);
  const auto prev = he.prev();
  const auto next = he.next();
  auto& prev_r = graph.hedges_[prev.index_];
  auto& next_r = graph.hedges_[next.index_];
  prev_r.next_ = next.index_;
  next_r.prev_ = prev.index_;
}

size_t dissolve_vertex(StrokeGraph& graph, const size_t vi) {
  Float join_arclen;
  auto dissolved = dissolved_stroke(graph, vi, join_arclen);
  if (dissolved.size() == 0) {
    return invalid;
  }

  const auto v = graph.vertex(vi);
  const auto save_vids = v.vertex_ids();
  const auto next = v.hedge();
  const auto he = next.prev();
  // We'll keep he.stroke() and remove next.stroke().
  const auto& edge_s = graph.strokes_[he.stroke_idx()];
  const auto& next_s = graph.strokes_[next.stroke_idx()];
  assert((edge_s.has_time() && next_s.has_time()) ||
         (!edge_s.has_time() && !next_s.has_time()));

  size_t s1 = he.stroke_idx();
  size_t s2 = next.stroke_idx();
  Float s1_len = edge_s.length();
  Float s2_len = next_s.length();
  std::pair<Float, Float> s1_to_s, s2_to_s;
  if (he.forward()) {
    s1_to_s = std::pair<Float, Float>(0., join_arclen / (s1_len + s2_len));
    s2_to_s = std::pair<Float, Float>(join_arclen / (s1_len + s2_len), 1.);
    if (!next.forward()) {
      std::swap(s2_to_s.first, s2_to_s.second);
    }
  } else { // !he.forward() case
    s1_to_s = std::pair<Float, Float>(join_arclen / (s1_len + s2_len), 1.);
    s2_to_s = std::pair<Float, Float>(0., join_arclen / (s1_len + s2_len));
    if (next.forward()) {
      std::swap(s2_to_s.first, s2_to_s.second);
    }
  }

  graph.strokes_[he.stroke_idx()] = std::move(dissolved);
  graph.bvh_.full_update(he.stroke_idx());

  if (!graph.faces_.empty()) {
    // Fix remaining face references.
    for (auto& cycle : graph.faces_[next.face_idx()].cycles_) {
      if (cycle == next.index_) {
        cycle = he.index_;
        break;
      }
    }
    for (auto& cycle : graph.faces_[next.twin().face_idx()].cycles_) {
      if (cycle == next.twin().index_) {
        cycle = he.twin().index_;
        break;
      }
    }
  }

  if (next.continuity() != invalid) {
    graph.hedges_[next.continuity()].continuity_ = invalid;
  }
  if (next.twin().continuity() != invalid) {
    graph.hedges_[next.twin().continuity()].continuity_ = he.twin().index_;
    graph.hedges_[he.twin().index_].continuity_ = next.twin().continuity();
  }

  const auto new_dest_i = next.dest().index_;
  graph.vertices_[new_dest_i].hedge_ = he.twin().index_;
  graph.hedges_[he.twin().index_].origin_ = new_dest_i;
  skip_over(graph, next.index_);
  skip_over(graph, next.twin().index_);
  graph.vertices_[v.index_].deactivate();

  // Merge mappings
  merge_stroke_mapping(graph, s1, s1_to_s, s2, s2_to_s, s1, save_vids);

  // Can't delete the vertex; incremental needs stable vertex indices.
  // delete_vertex(graph, v.index_);
  // Can't delete the edge; orienting vertex stars needs stable hedge indices.
  // delete_stroke(graph, next.index_ / 2);
  graph.bvh_.clear(next.index_ / 2);
  graph.hedges_[next.index_].deactivate();
  graph.hedges_[next.twin().index_].deactivate();

  return he.stroke_idx();
}

Stroke dissolved_stroke(const StrokeGraph& graph, const size_t vi, Float& join_arclen) {
  const auto v = graph.vertex(vi);
  assert(v.is_active());
  if (v.valence() != 2) {
    return Stroke();
  }
  auto next = v.hedge();
  auto he = next.prev();
  if (he.stroke_idx() == next.stroke_idx()) {
    return Stroke(); // Closed loop.
  }
  return concatenate_edges(graph, he.index_, next.index_, join_arclen);
}

Stroke concatenate_edges(const StrokeGraph& graph, const size_t hi1, const size_t hi2,
                         Float& join_arclen) {
  const auto he = graph.hedge(hi1);
  const auto next = graph.hedge(hi2);
  assert(he.is_active());
  assert(next.is_active());
  const auto& edge_s = graph.strokes_[he.stroke_idx()];
  const auto& next_s = graph.strokes_[next.stroke_idx()];
  assert((edge_s.has_time() && next_s.has_time()) ||
         (!edge_s.has_time() && !next_s.has_time()));

  auto out = Stroke(edge_s.size() + next_s.size() - 1, edge_s.has_time());
  if (he.forward()) {
    if (next.forward()) {
      copy_run(out, edge_s.size() - 1, out.size() - 1, next_s, 0, next_s.size() - 1);
    } else {
      copy_run(out, edge_s.size() - 1, out.size() - 1, next_s, next_s.size() - 1, 0);
    }
    assert(out.xy(edge_s.size() - 1).isApprox(edge_s.xy(Back)) &&
           "tried to concatenate geometrically disconnected edges");
    out.width(edge_s.size() - 1) =
      std::max(out.width(edge_s.size() - 1), edge_s.width(Back));
    copy_run(out, 0, edge_s.size() - 2, edge_s, 0, edge_s.size() - 2);
    join_arclen = edge_s.length();
  } else {
    if (next.forward()) {
      copy_run(out, 0, next_s.size() - 1, next_s, next_s.size() - 1, 0);
    } else {
      copy_run(out, 0, next_s.size() - 1, next_s, 0, next_s.size() - 1);
    }
    join_arclen = next_s.length();
    assert(out.xy(next_s.size() - 1).isApprox(edge_s.xy(0)) &&
           "tried to concatenate geometrically disconnected edges");
    out.width(next_s.size() - 1) =
      std::max(out.width(next_s.size() - 1), edge_s.width(0));
    copy_run(out, next_s.size(), out.size() - 1, edge_s, 1, edge_s.size() - 1);
  }
  return out;
}

bool get_forward_chain_original_indices(
  const StrokeGraph& graph, const size_t he_idx,
  std::vector<std::pair<size_t, bool>>& orig_indices, bool to_dissolve) {

  // Initialize frontiers
  auto frontier = graph.hedge(he_idx);
  assert(frontier);

  auto iterations = 0;
  do {
    auto he = frontier;
    frontier = HedgeView();

    const auto si = he.stroke_idx();
    force_assert(graph.strokes2orig_[si].size() == 1);
    const auto orig_si = graph.strokes2orig_[si][0];
    if (orig_indices.empty() ||
        (!orig_indices.empty() && orig_si != orig_indices.back().first &&
         orig_si != orig_indices[0].first)) {
      auto found = false;
      auto forward = false;
      for (const auto& [si2, map] : graph.orig2strokes_[orig_si]) {
        if (si == si2) {
          found = true;
          force_assert(map.range_arclens_.size() == 2);
          forward = map.range_arclens_[0] < map.range_arclens_[1];
          break;
        }
      }
      force_assert(found && "mapping references inconsistent");
      if (!he.forward()) {
        forward = !forward; // Match the direction of the half edge.
      }
      orig_indices.emplace_back(orig_si, forward);
    }

    if (he.twin().continuity_edge().is_valid()) {
      frontier = he.twin().continuity_edge();
    } else if (to_dissolve && should_dissolve_vertex(he.dest())) {
      frontier = he.next();
    }

    iterations++;
    force_assert(iterations < 1024 && "likely infinite loop found");
  } while (frontier.is_valid() && frontier.index_ != he_idx);

  force_assert(!orig_indices.empty());
  return frontier.index_ == he_idx;
}

Stroke get_chain(const HedgeView he, std::vector<std::pair<size_t, bool>>& orig_indices) {
  orig_indices.clear();
  const auto& graph = *he.graph_;
  const auto twin = he.twin();
  const auto cyclic =
    get_forward_chain_original_indices(graph, twin.index_, orig_indices, true);
  std::reverse(orig_indices.begin(), orig_indices.end());
  for (auto& pair : orig_indices) {
    pair.second = !pair.second; // Reverse directions.
  }
  if (!cyclic) {
    orig_indices.pop_back();
    const auto cyclic2 =
      get_forward_chain_original_indices(graph, he.index_, orig_indices, true);
    force_assert(!cyclic2 && "cyclic status inconsistent");
  }

#ifndef NDEBUG
  // Check for duplicates.
  auto sorted_copy = orig_indices;
  std::sort(sorted_copy.begin(), sorted_copy.end());
  for (size_t i = 1; i < sorted_copy.size(); ++i) {
    force_assert(sorted_copy[i - 1].first != sorted_copy[i].first);
  }
#endif

  const auto& strokes = graph.orig_strokes_;
  auto new_size = Index(0);
  for (size_t i = 0; i < orig_indices.size(); ++i) {
    const auto [orig_si, forward] = orig_indices[i];
    new_size += strokes[orig_si].size();
  }

  constexpr auto eps = 1e-6;
  auto out_stroke = Stroke(new_size, strokes[orig_indices[0].first].has_time());
  auto out_idx = Index(0);
  const auto n_orig_indices = (Index)orig_indices.size();
  for (Index i = 0; i < n_orig_indices; ++i) {
    const auto [orig_si, forward] = orig_indices[i];
    const auto& orig_stroke = strokes[orig_si];
    auto src_start = Index(0);
    auto src_start_fractional = Float(0);
    auto src_stop = orig_stroke.size() - 1;
    auto src_stop_fractional = Float(0);
    if (!forward) {
      std::swap(src_start, src_stop);
      std::swap(src_start_fractional, src_stop_fractional);
    }
    if (i > 0) {
      const auto [prev_orig_si, prev_forward] = orig_indices[i - 1];
      const auto prev_endpoint =
        (prev_forward ? strokes[prev_orig_si].xy(Back) : strokes[prev_orig_si].xy(0));
      force_assert(&strokes[prev_orig_si] != &orig_stroke);
      orig_stroke.ensure_arclengths();
      Vec2 proj;
      Float s = NaN;
      closest_point(orig_stroke, prev_endpoint, proj, s);
      force_assert(0.0 <= s && s <= orig_stroke.length());
      std::tie(src_start, src_start_fractional) = orig_stroke.fractional_index(s);
    }
    if (i + 1 < n_orig_indices) {
      const auto [next_orig_si, next_forward] = orig_indices[i + 1];
      const auto next_endpoint =
        (next_forward ? strokes[next_orig_si].xy(0) : strokes[next_orig_si].xy(Back));
      force_assert(&strokes[next_orig_si] != &orig_stroke);
      orig_stroke.ensure_arclengths();
      Vec2 proj;
      Float s = NaN;
      closest_point(orig_stroke, next_endpoint, proj, s);
      force_assert(0.0 <= s && s <= orig_stroke.length());
      std::tie(src_stop, src_stop_fractional) = orig_stroke.fractional_index(s);
    }

    if (src_start_fractional > 1.0 - eps) {
      src_start++;
      src_start_fractional = 0.0;
    } else if (src_start_fractional > 0.0) {
      out_stroke.x(out_idx) = //
        fast_lerp(orig_stroke.x(src_start), orig_stroke.x(src_start + 1),
                  src_start_fractional);
      out_stroke.y(out_idx) = //
        fast_lerp(orig_stroke.y(src_start), orig_stroke.y(src_start + 1),
                  src_start_fractional);
      out_stroke.width(out_idx) =
        fast_lerp(orig_stroke.width(src_start), orig_stroke.width(src_start + 1),
                  src_start_fractional);
      if (out_stroke.has_time()) {
        out_stroke.time(out_idx) =
          fast_lerp(orig_stroke.time(src_start), orig_stroke.time(src_start + 1),
                    src_start_fractional);
      }
      out_idx++;
      if (src_start < src_stop)
        src_start++;
      src_start_fractional = 0.0;
    }
    if (src_stop_fractional > 1.0 - eps) {
      src_stop++;
      src_stop_fractional = 0.0;
    }
    if (src_stop_fractional > eps && src_stop < src_start) {
      src_stop++; // Temporarily modify to make copy_run correct.
    }

    assert(out_idx + std::abs(src_stop - src_start) < out_stroke.size());
    copy_run(out_stroke, out_idx, out_idx + std::abs(src_stop - src_start), //
             orig_stroke, src_start, src_stop);
    out_idx += std::abs(src_stop - src_start) + 1;

    if (src_stop_fractional > eps && src_stop < src_start) {
      src_stop--; // Revert.
    }
    if (src_stop_fractional > eps) {
      out_stroke.x(out_idx) = //
        fast_lerp(orig_stroke.x(src_stop), orig_stroke.x(src_stop + 1),
                  src_stop_fractional);
      out_stroke.y(out_idx) = //
        fast_lerp(orig_stroke.y(src_stop), orig_stroke.y(src_stop + 1),
                  src_stop_fractional);
      out_stroke.width(out_idx) =
        fast_lerp(orig_stroke.width(src_stop), orig_stroke.width(src_stop + 1),
                  src_stop_fractional);
      if (out_stroke.has_time()) {
        out_stroke.time(out_idx) =
          fast_lerp(orig_stroke.time(src_stop), orig_stroke.time(src_stop + 1),
                    src_start_fractional);
      }
      out_idx++;
    }
  }
  assert(out_idx <= out_stroke.size());
  out_stroke.resize(out_idx);
  remove_duplicate_vertices(out_stroke);
  return out_stroke;
}

std::vector<Stroke>
chained_drawing(const StrokeGraph& graph,
                std::vector<std::vector<std::pair<size_t, bool>>>& mapping) {
  auto already_added_original_indices =
    std::vector<bool>(graph.orig_strokes_.size(), false);
  auto out_strokes = std::vector<Stroke>();

  for (size_t i = 0; i < graph.orig_strokes_.size(); ++i) {
    if (!already_added_original_indices[i]) {
      auto& o2s = graph.orig2strokes_[i];
      if (!o2s.empty()) {
        const auto hi = 2 * o2s[0].first;
        auto& orig_indices = mapping.emplace_back();
        out_strokes.emplace_back(get_chain(graph.hedge(hi), orig_indices));
        assert(!orig_indices.empty());

        for (const auto& [orig_si, forward] : orig_indices) {
          force_assert(!already_added_original_indices[orig_si]);
          already_added_original_indices[orig_si] = true;
        }
      }
    }
  }

  assert(mapping.size() == out_strokes.size());
  return out_strokes;
}

} // namespace sketching
