#include "stroke_graph_extra.h"

#include "deform.h"
#include "detail/alloca.h"
#include "detail/util.h"
#include "fairing.h"
#include "force_assert.h"
#include "incremental.h"
#include "intersect.h"
#include "mapping.h"
#include "stroke_graph.h"

#include <iomanip>

namespace sketching {

using HedgeView = StrokeGraph::HedgeView;
using VertexView = StrokeGraph::VertexView;

namespace {

constexpr auto invalid = StrokeGraph::invalid;
constexpr auto NaN = std::numeric_limits<Float>::quiet_NaN();

template <typename T>
void ensure_size(std::vector<T>& v, const size_t min_size) {
  if (v.size() < min_size) {
    v.resize(min_size);
  }
}

void insert_without_duplicates(std::vector<size_t>& v, const size_t x) {
  if (std::find(v.begin(), v.end(), x) == v.end()) {
    v.emplace_back(x);
  }
}

/**
 * Update the strokes2orig_ mapping to reflect an index change or removal in the
 * orig2strokes_ mapping.
 *
 * Performing an index update from old_edge_idx to new_edge_idx correctly will typically
 * look like this:
 *
 *     auto& map_pair = graph.orig2strokes_[orig_stroke_idx][...];
 *     const auto old_edge_idx = map_pair->first;
 *     map_pair->first = new_edge_idx;
 *     insert_without_duplicates(graph.strokes2orig_[new_edge_idx], orig_stroke_idx);
 *     remove_invalid_reverse_reference(graph, old_edge_idx, orig_stroke_idx);
 *
 * Note that we cannot indiscriminately remove orig_stroke_idx from
 * strokes2orig_[old_edge_idx] because old_edge_idx might equal new_edge_idx, or, more
 * importantly, there may have been more than one map_pair in
 * orig2strokes_[orig_stroke_idx] that refers to old_edge_idx.
 */
void remove_invalid_reverse_reference(StrokeGraph& graph, const size_t old_edge_idx,
                                      const size_t orig_stroke_idx) {
  assert(old_edge_idx < graph.strokes_.size() &&
         "did you get the order of old_edge_idx and orig_stroke_idx swapped?");
  assert(orig_stroke_idx < graph.orig_strokes_.size() &&
         "did you get the order of old_edge_idx and orig_stroke_idx swapped?");
  auto found = false;
  for (const auto& [i, _] : graph.orig2strokes_[orig_stroke_idx]) {
    if (i == old_edge_idx) {
      found = true;
      break;
    }
  }
  if (!found) {
    auto& original_indices = graph.strokes2orig_[old_edge_idx];
    original_indices.erase(
      std::remove(original_indices.begin(), original_indices.end(), orig_stroke_idx),
      original_indices.end());
  }
}

/// Is he and its twin adjacent to invalid face?
bool is_he_update_valid(const HedgeView he) {
  auto twin = he.twin();
  return (he.face_idx() == StrokeGraph::invalid ||
          twin.face_idx() == StrokeGraph::invalid);
}

/// Mark all strokes reachable from `he` with label `current_id`.
void connected_components_edge_he(const HedgeView he, const int current_id,
                                  const span<int> component_ids, bool to_update = false) {
  // TODO: Should rewrite to be non-recursive to avoid blowing the stack.
  const auto si = he.stroke_idx();
  if (component_ids[si] != -1) {
    return; // Already visited.
  }
  component_ids[si] = current_id;
  auto it = he.twin().next();
  // Assign umbrella around he's origin.
  auto valence = 0;
  while (it != he && (!to_update || is_he_update_valid(it))) {
    connected_components_edge_he(it, current_id, component_ids);
    it = it.twin().next();
    valence++;
    assert(valence < 1024 && "likely infinite loop found");
  }
  const auto next = he.next();
  it = next;
  // Assign umbrella around he's dest.
  valence = 0;
  do {
    if (!to_update || is_he_update_valid(it))
      connected_components_edge_he(it, current_id, component_ids);
    it = it.twin().next();
    valence++;
    assert(valence < 1024 && "likely infinite loop found");
  } while (it != next);
}

} // namespace

Endpoint vertex_to_endpoint(const VertexView v) {
  // We have a choice of which edge to use.
  // For consistency, we always pick the endpoint with the smallest value of `as_int`.
  const auto he = v.hedge();
  auto it = he;
  auto endp = Endpoint(std::numeric_limits<Endpoint::IdType>::max());
  auto valence = 0;
  do {
    const auto new_endp = Endpoint(it.stroke_idx(), it.forward());
    if (new_endp.as_int() < endp.as_int()) {
      endp = new_endp;
    }
    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
  assert(endp.stroke_idx() < v.graph_->strokes_.size());
  return endp;
}

static bool would_form_sliver_he_vert(const StrokeGraph& graph, const size_t hi,
                                      const size_t vi) {
  (void)graph;
  (void)hi;
  (void)vi;

#if 0
  const auto he = graph.hedge(hi);
  // If v1 and v2 share a common third vertex, avoid forming a sliver face.
  // See if any of the half-edges around v are similar to `he`.
  const auto v = graph.vertex(vi);
  const auto v_he = v.hedge();
  auto it = v_he;
  do {
    if (it.dest() == he.dest()) {
      const auto& s1 = he.stroke();
      const auto& s2 = it.stroke();
      auto s_start = 0.0, s_end = s1.length();
      if (!he.forward())
        std::swap(s_start, s_end);
      auto t_start = 0.0, t_end = s2.length();
      if (!it.forward())
        std::swap(t_start, t_end);

      if (similar_gaps_antisliver(s1, s_start, s_end, s2, t_start, t_end)) {
        return true;
      }
    }
    it = it.twin().next();
  } while (it != v_he);
#endif

  return false;
}

bool would_form_sliver_vert_interior(const StrokeGraph& graph, const size_t vi,
                                     const StrokeTime stroke_time) {
  assert(stroke_time.second <= 1.0 && "expected normalized arc length");

  (void)graph;
  (void)vi;

#if 0
  auto he = graph.hedge(2 * stroke_time.first);
  if (!he.forward())
    he = he.twin();
  // See if any of the half-edges around v are similar to `he`.
  const auto v = graph.vertex(vi);
  const auto v_he = v.hedge();
  auto it = v_he;
  auto valence = 0;
  do {
    const auto& s1 = he.stroke();
    const auto& s2 = it.stroke();
    // TODO: Eventually support anti-sliver logic when doing self-snaps.
    if (&s1 != &s2 && (it.dest() == he.dest() || it.dest() == he.origin())) {
      auto s_start = stroke_time.second * s1.length(), s_end = s1.length();
      if (it.dest() == he.origin()) {
        s_end = 0;
      }
      auto t_start = 0.0, t_end = s2.length();
      if (!it.forward())
        std::swap(t_start, t_end);

      if (similar_gaps_antisliver(s1, s_start, s_end, s2, t_start, t_end)) {
        return true;
      }
    }
    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != v_he);
#endif

  return false;
}

/**
 * Return true iff snapping vertex vi1 to vi2 would result in a spurious self-connection
 * for some edge.
 */
bool is_snap_spurious(const StrokeGraph& graph, const size_t vi1, const size_t vi2) {
  if (vi1 == vi2) {
    return true;
  }
  const auto dist = (graph.vertex(vi1).pos() - graph.vertex(vi2).pos()).norm();
  const auto he = graph.vertex(vi1).hedge();
  const auto he2 = graph.vertex(vi2).hedge();
  auto it = he;
  auto iterations = 0;
  do {
    // Check if snap is spurious from the perspective of the current edges.
    if (it.dest().index_ == vi2 &&
        // Keep this synced with others (grep for [spurious-condition]).
        !(std::max(3 * dist, M_PI * it.stroke().pen_width()) < it.stroke().length())) {
      return true;
    }

    // Check if snap is spurious from the perspective of the original strokes.

    if (!(it.flags() & StrokeGraph::HedgeRecord::Bridge)) {
      auto [orig_idx1, orig_arclen1] =
        as_orig_position(graph, Endpoint(it.stroke_idx(), it.forward()).as_pair());
      orig_arclen1 *= graph.orig_strokes_[orig_idx1].length();

      auto it2 = he2;
      do {
        if (!(it2.flags() & StrokeGraph::HedgeRecord::Bridge)) {
          auto [orig_idx2, orig_arclen2] =
            as_orig_position(graph, Endpoint(it2.stroke_idx(), it2.forward()).as_pair());
          if (orig_idx1 == orig_idx2) {
            const auto& orig_stroke = graph.orig_strokes_[orig_idx2];
            orig_arclen2 *= orig_stroke.length();
            const auto orig_p1 = orig_stroke.pos(orig_arclen1);
            const auto orig_p2 = orig_stroke.pos(orig_arclen2);
            const auto orig_dist = (orig_p1 - orig_p2).norm();
            // Keep this synced with others (grep for [spurious-condition]).
            if (!(std::max(3 * orig_dist, M_PI * orig_stroke.pen_width()) <
                  std::abs(orig_arclen1 - orig_arclen2))) {
              return true;
            }
          }
        }

        it2 = it2.twin().next();
        iterations++;
        force_assert(iterations < 2048 && "likely infinite loop found");
      } while (it2 != he2);
    }

    it = it.twin().next();
    iterations++;
    force_assert(iterations < 2048 && "likely infinite loop found");
  } while (it != he);
  return false;
}

bool is_bridge_vertex(const StrokeGraph& stroke_graph, const StrokeGraph::VertexView v) {
  if (v.is_valid()) {
    auto he = v.hedge();
    auto hit = he;
    do {
      if (hit.flags() & StrokeGraph::HedgeRecord::Bridge)
        return true;
      hit = hit.twin().next();
    } while (hit != he);
  }

  return false;
}

/**
 * Return true iff snapping vertex vi1 to the specified stroke interior would result in a
 * spurious self-connection for some edge.
 *
 * `stroke_time` is specified as (stroke index, normalized arc length).
 */
[[nodiscard]] static bool is_snap_spurious(const StrokeGraph& graph, const size_t vi1,
                                           const StrokeTime& stroke_time) {
  assert(stroke_time.second <= 1.0 && "expected normalized arc length");

  const auto si = stroke_time.first;
  const auto arclen2 = stroke_time.second * graph.strokes_[si].length();
  const auto dist = (graph.vertex(vi1).pos() - graph.strokes_[si].pos(arclen2)).norm();
  auto [orig_idx2, orig_arclen2] = as_orig_position(graph, stroke_time);
  const auto& orig_stroke = graph.orig_strokes_[orig_idx2];
  orig_arclen2 *= orig_stroke.length();
  const auto he = graph.vertex(vi1).hedge();
  auto it = he;
  auto valence = 0;
  do {
    // Check if snap is spurious from the perspective of the current edges.
    const auto arclen1 = (it.forward() ? 0.0 : it.stroke().length());
    if (it.stroke_idx() == si &&
        // Keep this synced with others (grep for [spurious-condition]).
        !(std::max(3 * dist, M_PI * it.stroke().pen_width()) <
          std::abs(arclen1 - arclen2))) {
      return true;
    }

    // Check if snap is spurious from the perspective of the original strokes.
    auto [orig_idx1, orig_arclen1] =
      as_orig_position(graph, Endpoint(it.stroke_idx(), it.forward()).as_pair());
    if (orig_idx1 == orig_idx2) {
      orig_arclen1 *= orig_stroke.length();
      const auto orig_dist =
        (orig_stroke.pos(orig_arclen1) - orig_stroke.pos(orig_arclen2)).norm();
      // Note this is a little stricter than the classifier needs---the classifier
      // reprojects onto the original find the best position.
      // Keep this synced with others (grep for [spurious-condition]).
      if (!(std::max(3 * orig_dist, M_PI * orig_stroke.pen_width()) <
            std::abs(orig_arclen1 - orig_arclen2))) {
        return true;
      }
    }

    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
  return false;
}

bool is_snap_valid(const StrokeGraph& graph, const size_t vi1, const size_t vi2) {
  if (is_snap_spurious(graph, vi1, vi2)) {
    return false;
  }

  const auto he = graph.vertex(vi1).hedge();
  auto it = he;
  do {
    if (would_form_sliver_he_vert(graph, it.index_, vi2)) {
      return false;
    }
    it = it.twin().next();
  } while (it != he);

  return true;
}

bool is_snap_valid(const StrokeGraph& graph, const size_t vi1,
                   const StrokeTime stroke_time) {
  if (is_snap_spurious(graph, vi1, stroke_time)) {
    return false;
  }
  if (would_form_sliver_vert_interior(graph, vi1, stroke_time)) {
    return false;
  }
  return true;
}

bool has_continuous_edge(const StrokeGraph::VertexView v) {
  const auto he = v.hedge();
  auto it = he;
  auto valence = 0;
  do {
    if (it.continuity() != invalid) {
      return true;
    }
    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
  return false;
}

int connected_components_edge(const StrokeGraph& graph, const span<int> component_ids) {
  assert(graph.strokes_.size() == component_ids.size());
  std::fill(component_ids.begin(), component_ids.end(), -1); // Mark unvisited.
  auto current_id = 0;
  for (size_t hi = 0; hi < graph.hedges_.size(); hi += 2) {
    const auto he = graph.hedge(hi);
    if (component_ids[he.stroke_idx()] == -1 && he.is_active() &&
        is_he_update_valid(he)) {
      connected_components_edge_he(he, current_id, component_ids);
      current_id++;
    }
  }
  return current_id;
}

void cycle_positions(const HedgeView he, std::vector<Vec2>& out_vertices) {
  auto n_points = Index(0);
  auto it = he;
  do {
    n_points += it.stroke().size() - 1;
    it = it.next();
  } while (it != he);
  out_vertices.reserve(n_points);
  auto iter = 0;
  do {
    const auto& s = it.stroke();
    if (it.forward()) {
      const auto stop = s.size() - 1;
      for (Index i = 0; i < stop; ++i) {
        out_vertices.emplace_back(s.xy(i));
      }
    } else {
      for (Index i = s.size() - 1; i > 0; --i) {
        out_vertices.emplace_back(s.xy(i));
      }
    }
    it = it.next();
    iter++;
    force_assert(iter < 1048576 && "likely infinite loop found");
  } while (it != he);
}

void dehook_dangling_edges(StrokeGraph& graph) {
  const auto num_vertices = graph.vertices_.size();
  for (size_t vi = 0; vi < num_vertices; ++vi) {
    const auto v = graph.vertex(vi);
    if (v.is_active() && v.is_dangling()) {
      dehook_dangling_vertex(graph, v);
    }
  }
}

void dehook_dangling_vertex(StrokeGraph& graph, const StrokeGraph::VertexView vert) {
  assert(vert.is_valid() && vert.is_active() && vert.is_dangling());

  StrokeMapping dehooked_strokes;
  StrokeMapping dehooked_orig_strokes;
  const auto he = vert.hedge();
  auto& stroke = graph.strokes_[he.stroke_idx()];
  auto endp = vertex_to_endpoint(vert).as_pair();
  const auto ok = convert_strokes2orig(graph, endp);
  force_assert(ok && "couldn't map from strokes to orig");
  const auto orig_si = endp.first;
  auto& orig_stroke = graph.orig_strokes_[orig_si];
  const auto [start, end] = dehooked_range(stroke); // @optimize We only use one end.
  Index trim_amount;
  Float before_s_len = stroke.length();

  // Currently assume we only trim one end
  dehooked_strokes.domain_arclens_.push_back(0.0);
  dehooked_strokes.domain_arclens_.push_back(1.0);

  if (he.forward()) {
    Float arc = stroke.arclength(start);
    graph.vertices_[vert.index_].p_ = stroke.xy(start);
    trim_amount = start;
    stroke.trim(start, stroke.size());
    stroke.compute_arclengths();

    dehooked_strokes.range_arclens_.push_back(arc / before_s_len);
    dehooked_strokes.range_arclens_.push_back(1.0);
  } else {
    Float arc = stroke.arclength(end - 1);
    graph.vertices_[vert.index_].p_ = stroke.xy(end - 1);
    trim_amount = stroke.size() - end;
    stroke.trim(0, end);
    stroke.compute_arclengths();

    dehooked_strokes.range_arclens_.push_back(0.0);
    dehooked_strokes.range_arclens_.push_back(arc / before_s_len);
  }
  graph.bvh_.full_update(he.stroke_idx());

  Float before_orig_s_len = orig_stroke.length();
  // Currently assume we only trim one end
  dehooked_orig_strokes.domain_arclens_.push_back(0.0);
  dehooked_orig_strokes.domain_arclens_.push_back(1.0);

  if (endp.second == 0) {
    Float arc = orig_stroke.arclength(trim_amount);
    orig_stroke.trim(trim_amount, orig_stroke.size());
    orig_stroke.compute_arclengths();

    dehooked_orig_strokes.range_arclens_.push_back(arc / before_orig_s_len);
    dehooked_orig_strokes.range_arclens_.push_back(1.0);
  } else {
    Float arc = orig_stroke.arclength(orig_stroke.size() - trim_amount - 1);
    orig_stroke.trim(0, orig_stroke.size() - trim_amount);
    orig_stroke.compute_arclengths();

    dehooked_orig_strokes.range_arclens_.push_back(0.0);
    dehooked_orig_strokes.range_arclens_.push_back(arc / before_orig_s_len);
  }

  // Update mapping
  for (size_t i = 0; i < graph.orig2strokes_.size(); ++i) {
    auto& mapping_arr = graph.orig2strokes_[i];
    for (auto& mapping : mapping_arr) {
      // Update the original
      if (i == orig_si) {
        mapping.second.domain_arclens_[0] =
          dehooked_orig_strokes.inv_map_clamp_from_start(
            mapping.second.domain_arclens_[0]);
        mapping.second.domain_arclens_[1] =
          dehooked_orig_strokes.inv_map_clamp_from_start(
            mapping.second.domain_arclens_[1]);
      }
      // Update the stroke
      size_t si = mapping.first;
      if (si == he.stroke_idx()) {
        mapping.second.range_arclens_[0] =
          dehooked_strokes.inv_map_clamp_from_start(mapping.second.range_arclens_[0]);
        mapping.second.range_arclens_[1] =
          dehooked_strokes.inv_map_clamp_from_start(mapping.second.range_arclens_[1]);
      }
    }
  }
}

void delete_stroke(StrokeGraph& graph, const size_t si) {
  assert(si < graph.strokes_.size() && "out of bounds");
  // We can delete from a std::vector in constant time by swapping with the last element
  // and popping the last element.
  const auto si_end = graph.strokes_.size() - 1;
  if (si != si_end) {
    graph.bvh_.swap_strokes(si, si_end);
    std::swap(graph.hedges_[2 * si], graph.hedges_[2 * si_end]);
    std::swap(graph.hedges_[2 * si + 1], graph.hedges_[2 * si_end + 1]);
    std::swap(graph.strokes2orig_[si], graph.strokes2orig_[si_end]);
    std::swap(graph.strokes2vid_[si], graph.strokes2vid_[si_end]);
  }
  graph.bvh_.pop_back();
  graph.hedges_.pop_back();
  graph.hedges_.pop_back();
  graph.strokes2orig_.pop_back();
  graph.strokes2vid_.pop_back();
  if (si != si_end && graph.strokes_[si].size() > 0) {
    // However, we need to update any references to the old last element.
    for (auto orig_si : graph.strokes2orig_[si]) {
      for (auto& sj_mapping : graph.orig2strokes_[orig_si]) {
        if (sj_mapping.first == si_end) {
          sj_mapping.first = si;
        }
      }
    }
    auto& edge_record = graph.hedges_[2 * si];
    auto& twin_record = graph.hedges_[2 * si + 1];
    if (edge_record.next_ == 2 * si_end) {
      edge_record.next_ = 2 * si;
    } else if (edge_record.next_ == 2 * si_end + 1) {
      edge_record.next_ = 2 * si + 1;
    }
    if (edge_record.prev_ == 2 * si_end) {
      edge_record.prev_ = 2 * si;
    } else if (edge_record.prev_ == 2 * si_end + 1) {
      edge_record.prev_ = 2 * si + 1;
    }
    if (edge_record.continuity_ == 2 * si_end) {
      edge_record.continuity_ = 2 * si;
    } else if (edge_record.continuity_ == 2 * si_end + 1) {
      edge_record.continuity_ = 2 * si + 1;
    }

    if (twin_record.next_ == 2 * si_end) {
      twin_record.next_ = 2 * si;
    } else if (twin_record.next_ == 2 * si_end + 1) {
      twin_record.next_ = 2 * si + 1;
    }
    if (twin_record.prev_ == 2 * si_end) {
      twin_record.prev_ = 2 * si;
    } else if (twin_record.prev_ == 2 * si_end + 1) {
      twin_record.prev_ = 2 * si + 1;
    }
    if (twin_record.continuity_ == 2 * si_end) {
      twin_record.continuity_ = 2 * si;
    } else if (twin_record.continuity_ == 2 * si_end + 1) {
      twin_record.continuity_ = 2 * si + 1;
    }

    graph.hedges_[edge_record.next_].prev_ = 2 * si;
    graph.hedges_[edge_record.prev_].next_ = 2 * si;
    graph.hedges_[twin_record.next_].prev_ = 2 * si + 1;
    graph.hedges_[twin_record.prev_].next_ = 2 * si + 1;
    graph.vertices_[edge_record.origin_].hedge_ = 2 * si;
    graph.vertices_[twin_record.origin_].hedge_ = 2 * si + 1;
    if (edge_record.continuity_ != StrokeGraph::invalid)
      graph.hedges_[edge_record.continuity_].continuity_ = 2 * si;
    if (twin_record.continuity_ != StrokeGraph::invalid)
      graph.hedges_[twin_record.continuity_].continuity_ = 2 * si + 1;
    if (edge_record.face_ != StrokeGraph::invalid) {
      for (auto& hi : graph.faces_[edge_record.face_].cycles_) {
        if (hi == 2 * si_end) {
          hi = 2 * si;
        }
      }
    }
    if (twin_record.face_ != StrokeGraph::invalid) {
      for (auto& hi : graph.faces_[twin_record.face_].cycles_) {
        if (hi == 2 * si_end + 1) {
          hi = 2 * si + 1;
        }
      }
    }
  }
}

void delete_vertex(StrokeGraph& graph, const size_t vi) {
  assert(vi < graph.vertices_.size() && "out of bounds");
  // We can delete from a std::vector in constant time by swapping with the last element
  // and popping the last element.
  const auto vi_end = graph.vertices_.size() - 1;
  if (vi != vi_end) {
    std::swap(graph.vertices_[vi], graph.vertices_[vi_end]);
  }
  graph.vertices_.pop_back();
  {
    // Delete mappings to and from vi.
    const auto it = graph.vertex2endpoint_.find(vi);
    if (it != graph.vertex2endpoint_.end()) {
      const auto endp = it->second;
      graph.vertex2endpoint_.erase(vi);
      graph.endpoint2vertex_.erase(endp);
    }
  }
  if (vi != vi_end) {
    // Modify mappings to and from vi_end.
    const auto it = graph.vertex2endpoint_.find(vi_end);
    if (it != graph.vertex2endpoint_.end()) {
      const auto endp = it->second;
      graph.vertex2endpoint_.erase(vi_end);
      graph.endpoint2vertex_.erase(endp);
      graph.vertex2endpoint_.insert({vi, endp});
      graph.endpoint2vertex_.insert({endp, vi});
    }
  }
  if (vi != vi_end && graph.vertex(vi).is_active()) {
    // However, we need to update any references to the old last element.
    auto valence = 0;
    const auto he = graph.vertex(vi).hedge();
    auto it = he;
    do {
      auto& origin = graph.hedges_[it.index_].origin_;
      assert(origin == vi_end);
      origin = vi;
      valence++;
      it = it.twin().next();
      force_assert(valence < 1024 && "likely infinite loop found");
    } while (it != he);
  }
}

bool move_vertex(StrokeGraph& graph, const size_t vi, const Vec2 new_pos,
                 const StrokeGraph::SnappingType type,
                 const StrokeGraph::VertexView move_to_v) {
  assert(&graph == move_to_v.graph_);
  using ST = StrokeGraph::SnappingType;
  const auto v = graph.vertex(vi);
  if (v.pos() != new_pos) {
    auto valence = 0;
    const auto he = v.hedge();
    auto it = he;
    do {
      size_t si = it.stroke_idx();
      auto& stroke = graph.strokes_[si];
      if (it.forward()) {
        if (type == ST::Deformation) {
          non_intersecting_smooth_deform(move_to_v, stroke, new_pos, stroke.xy(Back));
        } else if (type == ST::LinearSolve) {
          non_intersecting_smart_deform(move_to_v, stroke, new_pos, stroke.xy(Back));
        } else if (type == ST::Connection) {
          stroke.insert(0, new_pos.x_, new_pos.y_, stroke.width(0));
        } else {
          std::abort(); // Unknown snap method.
        }
        graph.bvh_.shallow_update(si, StrokeBVH::UpdateHint::AddToHead);
      } else {
        if (type == ST::Deformation) {
          non_intersecting_smooth_deform(move_to_v, stroke, stroke.xy(0), new_pos);
        } else if (type == ST::LinearSolve) {
          non_intersecting_smart_deform(move_to_v, stroke, stroke.xy(0), new_pos);
        } else if (type == ST::Connection) {
          stroke.push_back(new_pos.x_, new_pos.y_, stroke.width(Back));
        } else {
          std::abort(); // Unknown snap method.
        }
        graph.bvh_.shallow_update(si, StrokeBVH::UpdateHint::AddToTail);
      }
      it = it.twin().next();
      valence++;
      force_assert(valence < 1024 && "likely infinite loop found");
    } while (it != he);

    graph.vertices_[vi].p_ = new_pos;
  }

  return true;
}

size_t bridge_vertices(StrokeGraph& graph, const size_t vi1, const size_t vi2) {
  force_assert(vi1 != vi2);
  const auto v1 = graph.vertex(vi1);
  const auto v2 = graph.vertex(vi2);
  force_assert(v1.is_active());
  force_assert(v2.is_active());

  auto bridge_edge = Stroke(2, false);
  bridge_edge.x(0) = v1.pos().x_;
  bridge_edge.y(0) = v1.pos().y_;
  bridge_edge.x(1) = v2.pos().x_;
  bridge_edge.y(1) = v2.pos().y_;
  bridge_edge.width(0) = 0;
  bridge_edge.width(1) = 0;
  graph.bvh_.add(std::move(bridge_edge));
  if (!graph.strokes2vid_.empty()) {
    graph.strokes2vid_.emplace_back();
  }
  if (!graph.strokes2orig_.empty()) {
    graph.strokes2orig_.emplace_back();
  }

  graph.hedges_.emplace_back();
  graph.hedges_.emplace_back();
  const auto new_edge_idx = graph.hedges_.size() - 2;
  const auto new_twin_idx = new_edge_idx + 1;
  auto& edge_rec = graph.hedges_[new_edge_idx];
  auto& twin_rec = graph.hedges_[new_twin_idx];
  edge_rec.origin_ = vi1;
  twin_rec.origin_ = vi2;
  edge_rec.next_ = edge_rec.prev_ = new_twin_idx;
  twin_rec.next_ = twin_rec.prev_ = new_edge_idx;
  edge_rec.flags_ |= StrokeGraph::HedgeRecord::Bridge;
  twin_rec.flags_ |= StrokeGraph::HedgeRecord::Bridge;

  insert_into_star(graph, new_edge_idx, vi1);
  insert_into_star(graph, new_twin_idx, vi2);

  // Update faces.
  if (!graph.faces_.empty()) {
    auto edge_walks_back_to_self = false;
    auto twin_walks_back_to_self = false;
    auto it = graph.hedge(new_edge_idx);
    auto iters = 0;
    while (true) {
      it = it.next();
      if (it.index_ == new_edge_idx) {
        edge_walks_back_to_self = true;
        break;
      }
      if (it.index_ == new_twin_idx) {
        break;
      }
      iters++;
      force_assert(iters < 1024 && "likely infinite loop found");
    }
    it = graph.hedge(new_twin_idx);
    while (true) {
      it = it.next();
      if (it.index_ == new_twin_idx) {
        twin_walks_back_to_self = true;
        break;
      }
      if (it.index_ == new_edge_idx) {
        break;
      }
      iters++;
      force_assert(iters < 1024 && "likely infinite loop found");
    }

    if (edge_walks_back_to_self || twin_walks_back_to_self) {
      size_t edge_to_get_a_new_face = new_edge_idx;
      size_t edge_to_not_get_a_new_face = new_twin_idx;
      if (twin_walks_back_to_self) {
        std::swap(edge_to_get_a_new_face, edge_to_not_get_a_new_face);
      }

      const auto he = graph.hedge(edge_to_get_a_new_face);
      auto& new_face = graph.faces_.emplace_back();
      const auto new_face_idx = graph.faces_.size() - 1;
      new_face.cycles_.emplace_back(he.index_);
      it = he;
      do {
        auto& face = graph.hedges_[it.index_].face_;
        if (face != invalid) {
          auto& cycles = graph.faces_[face].cycles_;
          auto reverse_ref = std::find(cycles.begin(), cycles.end(), it.index_);
          if (reverse_ref != cycles.end()) {
            *reverse_ref = edge_to_not_get_a_new_face;
          }
        }
        face = new_face_idx;
        it = it.next();
      } while (it != he);
    }
    graph.hedges_[new_edge_idx].face_ = graph.hedge(new_edge_idx).next().face_idx();
    graph.hedges_[new_twin_idx].face_ = graph.hedge(new_twin_idx).next().face_idx();
  }

  return new_edge_idx / 2;
}

void add_stroke_mapping(StrokeGraph& graph, size_t orig_si,
                        std::pair<Float, Float> orig_domain, size_t si,
                        std::pair<Float, Float> s_range, int si_remain,
                        std::pair<Float, Float> s_remain_range) {
  orig_domain.first = std::max(orig_domain.first, 0.);
  orig_domain.second = std::min(orig_domain.second, 1.);
  s_range.first = std::clamp(s_range.first, 0.0, 1.0);
  s_range.second = std::clamp(s_range.second, 0.0, 1.0);

  assert(orig_domain.first < orig_domain.second);

  // Find the mapping to split
  int to_chop = -1;
  for (size_t i = 0; i < graph.orig2strokes_[orig_si].size(); ++i) {
    const auto& mapping = graph.orig2strokes_[orig_si][i].second;

    assert(mapping.domain_arclens_.size() == 2);

    // Don't allow chopping a domain from the middle
    if ((orig_domain.first > mapping.domain_arclens_[0] &&
         orig_domain.second < mapping.domain_arclens_[1]) ||
        (orig_domain.first < mapping.domain_arclens_[0] &&
         orig_domain.second > mapping.domain_arclens_[1])) {
      std::abort();
    }

    if (std::abs(mapping.domain_arclens_[0] - orig_domain.first) <
          std::numeric_limits<Float>::epsilon() ||
        std::abs(mapping.domain_arclens_[1] - orig_domain.second) <
          std::numeric_limits<Float>::epsilon()) {
      to_chop = (int)i;
      break;
    }
  }

  ensure_size(graph.strokes2orig_, (size_t)std::max((int)si, si_remain) + 1);

  // Adding a new mapping without chopping
  if (to_chop < 0) {
    StrokeMapping mapping;
    mapping.domain_arclens_.push_back(orig_domain.first);
    mapping.domain_arclens_.push_back(orig_domain.second);
    mapping.range_arclens_.push_back(s_range.first);
    mapping.range_arclens_.push_back(s_range.second);
    graph.orig2strokes_[orig_si].emplace_back(si, mapping);

    insert_without_duplicates(graph.strokes2orig_[si], orig_si);
    if (si_remain >= 0) {
      insert_without_duplicates(graph.strokes2orig_[si_remain], orig_si);
    }
    return;
  }

  size_t si_replace = graph.orig2strokes_[orig_si][to_chop].first;
  const auto& mapping = graph.orig2strokes_[orig_si][to_chop].second;

  // Match exactly
  // When si_remain < 0 and the match is not exact, we are updating the existing record.
  if ((mapping.domain_arclens_[0] == orig_domain.first &&
       mapping.domain_arclens_[1] == orig_domain.second) ||
      si_remain < 0) {
    graph.orig2strokes_[orig_si][to_chop].first = si;
    graph.orig2strokes_[orig_si][to_chop].second.domain_arclens_[0] = orig_domain.first;
    graph.orig2strokes_[orig_si][to_chop].second.domain_arclens_[1] = orig_domain.second;
    graph.orig2strokes_[orig_si][to_chop].second.range_arclens_[0] = s_range.first;
    graph.orig2strokes_[orig_si][to_chop].second.range_arclens_[1] = s_range.second;
  }
  // Split from the front
  else if (mapping.domain_arclens_[0] == orig_domain.first) {
    if (s_remain_range.first < s_remain_range.second) {
      s_remain_range.first = std::max(s_remain_range.first, 0.);
      s_remain_range.second = std::min(s_remain_range.second, 1.);
    } else {
      s_remain_range.second = std::max(s_remain_range.second, 0.);
      s_remain_range.first = std::min(s_remain_range.first, 1.);
    }
    assert(si_remain >= 0 && std::min(s_remain_range.first, s_remain_range.second) >= 0 &&
           std::max(s_remain_range.first, s_remain_range.second) > 0);
    graph.orig2strokes_[orig_si][to_chop].first = si_remain;
    graph.orig2strokes_[orig_si][to_chop].second.domain_arclens_[0] = orig_domain.second;
    graph.orig2strokes_[orig_si][to_chop].second.range_arclens_[0] = s_remain_range.first;
    graph.orig2strokes_[orig_si][to_chop].second.range_arclens_[1] =
      s_remain_range.second;
    auto insert_itr = graph.orig2strokes_[orig_si].begin() + to_chop;
    StrokeMapping new_mapping;
    new_mapping.domain_arclens_.push_back(orig_domain.first);
    new_mapping.domain_arclens_.push_back(orig_domain.second);
    new_mapping.range_arclens_.push_back(s_range.first);
    new_mapping.range_arclens_.push_back(s_range.second);
    graph.orig2strokes_[orig_si].insert(insert_itr, std::make_pair(si, new_mapping));
  }
  // Split from the back
  else if (mapping.domain_arclens_[1] == orig_domain.second) {
    if (s_remain_range.first < s_remain_range.second) {
      s_remain_range.first = std::max(s_remain_range.first, 0.);
      s_remain_range.second = std::min(s_remain_range.second, 1.);
    } else {
      s_remain_range.second = std::max(s_remain_range.second, 0.);
      s_remain_range.first = std::min(s_remain_range.first, 1.);
    }
    assert(si_remain >= 0 && std::min(s_remain_range.first, s_remain_range.second) >= 0 &&
           std::max(s_remain_range.first, s_remain_range.second) > 0);
    graph.orig2strokes_[orig_si][to_chop].first = si_remain;
    graph.orig2strokes_[orig_si][to_chop].second.domain_arclens_[1] = orig_domain.first;
    graph.orig2strokes_[orig_si][to_chop].second.range_arclens_[0] = s_remain_range.first;
    graph.orig2strokes_[orig_si][to_chop].second.range_arclens_[1] =
      s_remain_range.second;
    auto insert_itr = graph.orig2strokes_[orig_si].begin() + (to_chop + 1);
    StrokeMapping new_mapping;
    new_mapping.domain_arclens_.push_back(orig_domain.first);
    new_mapping.domain_arclens_.push_back(orig_domain.second);
    new_mapping.range_arclens_.push_back(s_range.first);
    new_mapping.range_arclens_.push_back(s_range.second);
    graph.orig2strokes_[orig_si].insert(insert_itr, std::make_pair(si, new_mapping));
  }

  graph.strokes2orig_[si_replace].erase(
    std::remove(graph.strokes2orig_[si_replace].begin(),
                graph.strokes2orig_[si_replace].end(), orig_si),
    graph.strokes2orig_[si_replace].end());
  insert_without_duplicates(graph.strokes2orig_[si], orig_si);
  if (si_remain >= 0)
    insert_without_duplicates(graph.strokes2orig_[si_remain], orig_si);
}

void split_stroke_mapping(StrokeGraph& graph, size_t si, Float split, size_t si1,
                          int si2) {
  assert(split > 0 && split < 1);

  struct VertexIDPos {
    StrokeGraph::VertexID vid;
    std::pair<int, Float> orig_pos;
  };
  // Transfer the stored vertex ids
  std::vector<VertexIDPos> vid_mapping;
  vid_mapping.reserve(graph.strokes2vid_[si].size());
  for (const auto& vid : graph.strokes2vid_[si]) {
    VertexIDPos vid_pos;
    vid_pos.vid = vid.second;
    vid_pos.orig_pos = std::make_pair((int)si, vid.first);
    const auto ok = convert_strokes2orig(graph, vid_pos.orig_pos);
    force_assert(ok && "couldn't map from strokes to orig");
    vid_mapping.emplace_back(vid_pos);
  }

  ensure_size(graph.strokes2orig_, (size_t)std::max((int)si1, si2) + 1);

  struct MappingSection {
    size_t orig_idx_ = (size_t)-1;
    std::pair<size_t, StrokeMapping>* map_pair_ = nullptr;
  };
  MappingSection to_split;
  std::vector<MappingSection> before_split, after_split;

  // Prevent numerical instability: if `split` is very close to an existing range value,
  // then snap to that value.
  for (auto orig_si : graph.strokes2orig_[si]) {
    auto& orig_mappings = graph.orig2strokes_[orig_si];
    for (auto& pair : orig_mappings) {
      if (pair.first == si) {
        auto& mapping = pair.second;
        for (int i = 0; i < 2; ++i) {
          if (std::abs(split - mapping.range_arclens_[i]) < 1e-10) {
            split = mapping.range_arclens_[i];
            goto done;
          }
        }
      }
    }
  }
done:

  // 1. Categorize the original strokes based on the splitting position
  for (auto orig_si : graph.strokes2orig_[si]) {
    auto& orig_mappings = graph.orig2strokes_[orig_si];
    // Note that an edge can appear multiple times in the range.
    for (auto& pair : orig_mappings) {
      if (pair.first == si) {
        auto& mapping = pair.second;
        auto range_s = std::min(mapping.range_arclens_[0], mapping.range_arclens_[1]);
        auto range_l = std::max(mapping.range_arclens_[0], mapping.range_arclens_[1]);
        if (split > range_s && split < range_l) {
          assert(!to_split.map_pair_);
          to_split = {orig_si, &pair};
        } else if (split >= range_l) {
          before_split.emplace_back(MappingSection{orig_si, &pair});
        } else if (split <= range_s && si2 >= 0) {
          after_split.emplace_back(MappingSection{orig_si, &pair});
        }
      }
    }
  }

  // 2. Update before and after
  for (auto [orig_idx, map_pair] : before_split) {
    const auto old_edge_idx = map_pair->first;
    map_pair->first = si1;
    insert_without_duplicates(graph.strokes2orig_[map_pair->first], orig_idx);
    remove_invalid_reverse_reference(graph, old_edge_idx, orig_idx);
    // Recompute the range after splitting
    auto& range = map_pair->second.range_arclens_;
    range[0] /= split;
    range[1] /= split;
    // Domain stays the same.
  }
  for (auto [orig_idx, map_pair] : after_split) {
    const auto old_edge_idx = map_pair->first;
    map_pair->first = (si2 >= 0 ? si2 : si);
    insert_without_duplicates(graph.strokes2orig_[map_pair->first], orig_idx);
    remove_invalid_reverse_reference(graph, old_edge_idx, orig_idx);
    // Recompute the range after splitting
    auto& range = map_pair->second.range_arclens_;
    range[0] = (range[0] - split) / (1 - split);
    range[1] = (range[1] - split) / (1 - split);
    // Domain stays the same.
  }

  // 3. Split a mapping if necessary
  if (to_split.map_pair_) {
    const auto& mapping = to_split.map_pair_->second;
    // Split position, from stroke to original stroke
    Float orig_split = mapping.inv_map_clamp_from_start(split);

    // Orig -> si -> si1, si2
    auto range_s = std::min(mapping.range_arclens_[0], mapping.range_arclens_[1]);
    auto range_l = std::max(mapping.range_arclens_[0], mapping.range_arclens_[1]);
    bool flipped = (mapping.range_arclens_[0] > mapping.range_arclens_[1]);

    Float si1_l = range_s / split;

    if (si2 >= 0) {
      Float si2_r = (range_l - split) / (1 - split);
      add_stroke_mapping(
        graph, to_split.orig_idx_,
        (!flipped) ? std::make_pair(orig_split, mapping.domain_arclens_[1])
                   : std::make_pair(mapping.domain_arclens_[0], orig_split),
        si2, (!flipped) ? std::make_pair(0., si2_r) : std::make_pair(si2_r, 0.), (int)si1,
        (!flipped) ? std::make_pair(si1_l, 1.) : std::make_pair(1., si1_l));
    } else {
      add_stroke_mapping(
        graph, to_split.orig_idx_,
        (!flipped) ? std::make_pair(mapping.domain_arclens_[0], orig_split)
                   : std::make_pair(orig_split, mapping.domain_arclens_[1]),
        si1, (!flipped) ? std::make_pair(si1_l, 1.) : std::make_pair(1., si1_l));
    }
  }

  // 4. Store the vid and compute their positions based on the new mapping
  graph.strokes2vid_.resize(std::max((int)si1, si2) + 1);
  graph.strokes2vid_[si1].clear();
  if (si2 >= 0)
    graph.strokes2vid_[si2].clear();
  for (const auto& vid_pos : vid_mapping) {
    auto s_pos = vid_pos.orig_pos;
    const auto ok = convert_orig2strokes(graph, s_pos);
    force_assert(ok && "couldn't map from orig to strokes");
    assert(s_pos.first == si1 || s_pos.first == si2);
    graph.strokes2vid_[s_pos.first].emplace(s_pos.second, vid_pos.vid);
  }
}

void merge_stroke_mapping(StrokeGraph& graph, size_t si1, std::pair<Float, Float> s1_to_s,
                          size_t si2, std::pair<Float, Float> s2_to_s, size_t si,
                          const std::vector<StrokeGraph::VertexID>& save_vids) {
  struct VertexIDPos {
    StrokeGraph::VertexID vid;
    std::pair<int, Float> orig_pos;
  };
  // Transfer the stored vertex ids
  Float split_t =
    (s1_to_s.first == 0 || s1_to_s.first == 1) ? s1_to_s.second : s1_to_s.first;

  std::vector<VertexIDPos> vid_mapping;
  vid_mapping.reserve(graph.strokes2vid_[si1].size() + graph.strokes2vid_[si2].size());
  for (const auto& vid : graph.strokes2vid_[si1]) {
    VertexIDPos vid_pos;
    vid_pos.vid = vid.second;
    vid_pos.orig_pos = std::make_pair((int)si1, vid.first);
    const auto ok = convert_strokes2orig(graph, vid_pos.orig_pos);
    force_assert(ok && "couldn't map from strokes to orig");
    vid_mapping.emplace_back(vid_pos);
  }
  for (const auto& vid : graph.strokes2vid_[si2]) {
    VertexIDPos vid_pos;
    vid_pos.vid = vid.second;
    vid_pos.orig_pos = std::make_pair((int)si2, vid.first);
    const auto ok = convert_strokes2orig(graph, vid_pos.orig_pos);
    force_assert(ok && "couldn't map from strokes to orig");
    vid_mapping.emplace_back(vid_pos);
  }

  // 0. Build the mapping from s1, s2 to merged s
  StrokeMapping mapping1;
  mapping1.domain_arclens_.push_back(0);
  mapping1.domain_arclens_.push_back(1);
  mapping1.range_arclens_.push_back(s1_to_s.first);
  mapping1.range_arclens_.push_back(s1_to_s.second);
  StrokeMapping mapping2;
  mapping2.domain_arclens_.push_back(0);
  mapping2.domain_arclens_.push_back(1);
  mapping2.range_arclens_.push_back(s2_to_s.first);
  mapping2.range_arclens_.push_back(s2_to_s.second);

  // 1. Determine the merged original stroke set
  std::unordered_set<size_t> origs;
  for (auto orig_si : graph.strokes2orig_[si1]) {
    origs.emplace(orig_si);
  }
  for (auto orig_si : graph.strokes2orig_[si2]) {
    origs.emplace(orig_si);
  }
  // We only overwrite. Don't delete since other functions will do that.
  if (si == si1)
    graph.strokes2orig_[si1].clear();
  if (si == si2)
    graph.strokes2orig_[si2].clear();

  ensure_size(graph.strokes2orig_, si + 1);
  graph.strokes2orig_[si].insert(graph.strokes2orig_[si].end(), origs.begin(),
                                 origs.end());

  // 2. Update the stored mappings (from orig to strokes)
  auto map2s = [](const StrokeMapping& mapping, Float v) -> Float {
    return mapping.map_clamp(v);
  };
  for (auto orig_si : graph.strokes2orig_[si]) {
    auto& orig_mappings = graph.orig2strokes_[orig_si];
    bool modified = false;
    for (auto& si_mapping : orig_mappings) {
      if ((si_mapping.first == si1) || (si_mapping.first == si2)) {
        const auto old_edge_idx = si_mapping.first;
        insert_without_duplicates(graph.strokes2orig_[si], orig_si);
        remove_invalid_reverse_reference(graph, old_edge_idx, orig_si);
        auto& mapping = si_mapping.second;
        mapping.range_arclens_[0] = map2s((si_mapping.first == si1) ? mapping1 : mapping2,
                                          mapping.range_arclens_[0]);
        mapping.range_arclens_[1] = map2s((si_mapping.first == si1) ? mapping1 : mapping2,
                                          mapping.range_arclens_[1]);
        si_mapping.first = si;
        modified = true;
      }
    }
    assert(modified);
  }

  // Merge mapping records if two adjacent entries are continuous and point to
  // the same stroke.
  for (auto orig_si : graph.strokes2orig_[si]) {
    auto& orig_mappings = graph.orig2strokes_[orig_si];
    for (size_t i = 1; i < orig_mappings.size(); ++i) {
      if (orig_mappings[i - 1].first == orig_mappings[i].first) {
        auto& domain_before = orig_mappings[i - 1].second.domain_arclens_;
        auto& range_before = orig_mappings[i - 1].second.range_arclens_;
        auto& domain_current = orig_mappings[i].second.domain_arclens_;
        auto& range_current = orig_mappings[i].second.range_arclens_;
        // Domains are continuous
        assert(domain_before[1] == domain_current[0]);
        if (range_before[1] == range_current[0]) { // Ranges are also continuous
          domain_current[0] = domain_before[0];
          // Either both ranges are increasing or both ranges are decreasing.
          assert((range_before[0] < range_before[1]) ==
                 (range_current[0] < range_current[1]));
          range_current[0] = range_before[0];

          orig_mappings[i - 1].first = (size_t)-1; // Mark for deletion.
        }
        // Note that if the mapping is like
        //     [0, 0.5] -> [0.3, 1],   [0.5, 1] -> [0, 0.3],
        // we shouldn't merge them because the ranges aren't continuous.  This can happen
        // if we dissolve the original endpoint vertex on a closed stroke with one other
        // vertex in the interior.
      }
    }

    orig_mappings.erase(
      std::remove_if(orig_mappings.begin(), orig_mappings.end(),
                     [](const auto& pair) { return pair.first == (size_t)-1; }),
      orig_mappings.end());
  }

  // 3. Update the vertex ID record
  graph.strokes2vid_[si].clear();
  for (const auto& vid_pos : vid_mapping) {
    auto s_pos = vid_pos.orig_pos;
    const auto ok = convert_orig2strokes(graph, s_pos);
    force_assert(ok && "couldn't map from orig to strokes");
    assert(s_pos.first == si);
    graph.strokes2vid_[s_pos.first].emplace(s_pos.second, vid_pos.vid);
  }
  for (auto const& vid : save_vids)
    graph.strokes2vid_[si].emplace(split_t, vid);
}

void clamp_edge_mapping(StrokeGraph& graph, const size_t edge_idx,
                        const std::pair<Float, Float> new_domain) {
  assert(new_domain.first < new_domain.second);

  const auto [start, stop] = new_domain;
  for (const auto orig_si : graph.strokes2orig_[edge_idx]) {
    auto& o2s = graph.orig2strokes_[orig_si];
    for (auto& [si, mapping] : o2s) {
      if (si == edge_idx) {
        assert(mapping.range_arclens_.size() == 2);
        mapping.range_arclens_[0] =
          std::clamp((mapping.range_arclens_[0] - start) / (stop - start), 0.0, 1.0);
        mapping.range_arclens_[1] =
          std::clamp((mapping.range_arclens_[1] - start) / (stop - start), 0.0, 1.0);
      }
    }
    // Remove now-empty mappings.
    auto remove_it = std::remove_if(o2s.begin(), o2s.end(), [](const auto& m) {
      return m.second.range_arclens_[0] == m.second.range_arclens_[1];
    });
    const auto n_removed = std::distance(remove_it, o2s.end());
    if (n_removed > 0) {
      auto removed_stroke_indices = ALLOCA_SPAN(size_t, n_removed);
      auto removed_stroke_indices_inserter = removed_stroke_indices.begin();
      auto it = remove_it;
      for (; it != o2s.end(); ++it, ++removed_stroke_indices_inserter) {
        *removed_stroke_indices_inserter = it->first;
      }
      o2s.erase(remove_it, o2s.end());
      for (const auto si : removed_stroke_indices) {
        remove_invalid_reverse_reference(graph, si, orig_si);
      }
    }
  }
}

bool check_mapping(const StrokeGraph& graph) {
  bool no_dangling = true;

  // Check the continuity on original strokes
  std::unordered_set<size_t> seen_strokes;
  for (size_t i = 0; i < graph.orig2strokes_.size(); ++i) {
    const auto& orig_mappings = graph.orig2strokes_[i];
    if (orig_mappings.empty()) {
      continue;
    }

    int prev_stroke = -1;
    std::unordered_map<Float, size_t> interval_count;
    for (size_t j = 0; j < orig_mappings.size(); ++j) {
      const auto& mapping = orig_mappings[j];
      seen_strokes.emplace(mapping.first);

      // TODO: This is no longer required.
      if (prev_stroke == (int)mapping.first) {
        force_assert(orig_mappings[j - 1].second.range_arclens_[1] !=
                     mapping.second.range_arclens_[0]);
      }

      interval_count[mapping.second.domain_arclens_[0]]++;
      interval_count[mapping.second.domain_arclens_[1]]++;
      prev_stroke = (int)mapping.first;
    }

    force_assert(interval_count.count(0) && interval_count.count(1));

    for (auto [v, c] : interval_count) {
      if (v == 0 || v == 1)
        continue;
      force_assert(c == 2);
    }
  }

  // Check the continuity on strokes
  for (size_t si = 0; si < graph.strokes2orig_.size(); ++si) {
    if (!seen_strokes.count(si)) {
      no_dangling = false;
      continue;
    }
    const auto& origs = graph.strokes2orig_[si];
    std::unordered_map<Float, size_t> interval_count;

    for (auto orig_si : origs) {
      const auto& orig_mappings = graph.orig2strokes_[orig_si];
      for (const auto& pair : orig_mappings) {
        // Note that an edge can appear multiple times in the range.
        if (pair.first == si) {
          const auto& mapping = pair.second;
          if (!interval_count.count(mapping.range_arclens_[0]))
            interval_count[mapping.range_arclens_[0]] = 0;
          interval_count[mapping.range_arclens_[0]]++;
          if (!interval_count.count(mapping.range_arclens_[1]))
            interval_count[mapping.range_arclens_[1]] = 0;
          interval_count[mapping.range_arclens_[1]]++;
        }
      }
    }

    force_assert(interval_count.count(0) && interval_count.count(1));

    for (auto [v, c] : interval_count) {
      if (v == 0 || v == 1)
        continue;
      force_assert(c == 2);
    }
  }

  for (size_t hi = 0; hi < graph.hedges_.size(); ++hi) {
    const auto he = graph.hedge(hi);
    if (he.is_active() && !(he.flags() & StrokeGraph::HedgeRecord::Bridge) &&
        he.continuity() == invalid) {
      auto end = StrokeTime((int)he.stroke_idx(), he.forward() ? 0.0 : 1.0);
      const auto ok = convert_strokes2orig(graph, end);
      force_assert(ok && "couldn't map from strokes to orig");
      force_assert(end.second == 0.0 || end.second == 1.0);
    }
  }

  // Check if positions match (this only works when there's no deformation)
  // for (size_t i = 0; i < graph.strokes_.size(); ++i) {
  //  auto orig_p = std::make_pair((int)i, 0.5);
  //  Vec2 s_pos = graph.strokes_[orig_p.first].pos_norm(orig_p.second);
  //  convert_strokes2orig(graph, orig_p);
  //  Vec2 orig_s_pos = graph.orig_strokes_[orig_p.first].pos_norm(orig_p.second);
  //  assert((s_pos - orig_s_pos).norm() < 1e-6);
  //}

  return no_dangling;
}

void check_interior_intersection(const StrokeGraph& graph) {
  // 1. Build bboxes
  std::vector<PolylineBVHLeaf> stroke_nodes;
  stroke_nodes.reserve(graph.strokes_.size());
  for (const auto& s : graph.strokes_) {
    s.ensure_arclengths();
    stroke_nodes.emplace_back(s, bounds(s));
  }

  // 2. Check the interior intersection between all pairs
  for (size_t i = 0; i + 1 < graph.strokes_.size(); ++i) {
    for (size_t j = i + 1; j < graph.strokes_.size(); ++j) {
      std::vector<Vec2> out;
      intersect_different(stroke_nodes[i], stroke_nodes[j], out);
      for (const auto& out_v : out) {
        assert(out_v.y() <= 1e-10 ||
               std::abs(graph.strokes_[j].length() - out_v.y()) <= 1e-10);
      }
    }
  }
}

bool convert_strokes2orig(const StrokeGraph& graph, StrokeTime& end) {
  size_t si = end.first;
  Float s_time = end.second;
  if (si >= graph.strokes2orig_.size())
    return false;
  for (auto orig_si : graph.strokes2orig_[si]) {
    for (const auto& [sii, mapping] : graph.orig2strokes_[orig_si]) {
      if (sii == si &&
          s_time >= std::min(mapping.range_arclens_[0], mapping.range_arclens_[1]) &&
          s_time <= std::max(mapping.range_arclens_[0], mapping.range_arclens_[1])) {
        end.first = (int)orig_si;
        end.second = mapping.inv_map_clamp_from_end(s_time);
        return true;
      }
    }
  }
  return false;
}

bool convert_orig2strokes(const StrokeGraph& graph, StrokeTime& end) {
  size_t orig_si = end.first;
  Float orig_s_time = end.second;
  if (orig_si >= graph.orig2strokes_.size())
    return false;
  for (const auto& si_mapping : graph.orig2strokes_[orig_si]) {
    if (orig_s_time >= std::min(si_mapping.second.domain_arclens_[0],
                                si_mapping.second.domain_arclens_[1]) &&
        orig_s_time <= std::max(si_mapping.second.domain_arclens_[0],
                                si_mapping.second.domain_arclens_[1])) {
      end.first = (int)si_mapping.first;
      end.second = si_mapping.second.map_clamp(orig_s_time);
      return true;
    }
  }
  return false;
}

StrokeTime as_orig_position(const StrokeGraph& graph, const StrokeTime& edge_pos) {
  auto pos = edge_pos;
  const auto ok = convert_strokes2orig(graph, pos);
  force_assert(ok && "couldn't map from strokes to orig");
  return pos;
}

StrokeTime as_edge_position(const StrokeGraph& graph, const StrokeTime& orig_pos) {
  auto pos = orig_pos;
  const auto ok = convert_orig2strokes(graph, pos);
  force_assert(ok && "couldn't map from orig to strokes");
  return pos;
}

bool stroke_to_original_stroke_indexing(const StrokeGraph& stroke_graph, Junction& junc) {
  const auto p0 = junc.points[0]; // Store the old point so we can revert on failure.
  if (!convert_strokes2orig(stroke_graph, junc.points[0])) {
    return false;
  }
  if (!convert_strokes2orig(stroke_graph, junc.points[1])) {
    junc.points[0] = p0;
    return false;
  }
  return true;
}

bool original_stroke_to_stroke_indexing(const StrokeGraph& stroke_graph, Junction& junc) {
  const auto p0 = junc.points[0]; // Store the old point so we can revert on failure.
  if (!convert_orig2strokes(stroke_graph, junc.points[0])) {
    return false;
  }
  if (!convert_orig2strokes(stroke_graph, junc.points[1])) {
    junc.points[0] = p0;
    return false;
  }
  return true;
}

StrokeGraph::VertexView orig2endpoint(const StrokeGraph& graph, const StrokeTime& end) {
  auto s_end = end;
  if (!convert_orig2strokes(graph, s_end)) {
    return VertexView();
  }

  // Dissolved vertex
  if (!(s_end.second == 0.0 || s_end.second == 1.0))
    return StrokeGraph::VertexView();

  auto he = graph.hedge(2 * s_end.first);
  return (s_end.second == 0.0) ? he.origin() : he.dest();
}

static void drop_hedge(StrokeGraph& graph, const size_t hi) {
  const auto he = graph.hedge(hi);
  const auto prev_i = he.prev().index_;
  const auto next_i = he.twin().next().index_;
  const auto origin_i = he.origin().index_;

  if (he.continuity() != invalid) {
    graph.hedges_[he.continuity()].continuity_ = invalid;
  }
  graph.hedges_[prev_i].next_ = next_i;
  graph.hedges_[next_i].prev_ = prev_i;
  if (graph.hedge(next_i).is_active()) {
    graph.vertices_[origin_i].hedge_ = next_i;
  } else {
    graph.vertices_[origin_i].deactivate();
  }

  if (!graph.faces_.empty()) {
    for (auto& cycle_hi : graph.faces_[he.face_idx()].cycles_) {
      if (cycle_hi == he.index_) {
        cycle_hi = next_i;
        break;
      }
    }
  }

  graph.hedges_[he.index_].deactivate();
}

StrokeTime reproject_cand2_on_original(const StrokeGraph& graph,
                                       const StrokeTime& orig_pos1,
                                       const StrokeTime& cand2) {
  const auto query_arclen =
    orig_pos1.second * graph.orig_strokes_[orig_pos1.first].length();
  const auto xy = graph.orig_strokes_[orig_pos1.first].pos(query_arclen);

  auto best_dist = Float(INFINITY);
  auto best_st = StrokeTime();
  for (const auto orig_si : graph.strokes2orig_[cand2.first]) {
    const auto& orig_stroke = graph.orig_strokes_[orig_si];
    for (const auto& [si, map] : graph.orig2strokes_[orig_si]) {
      if ((int)si == cand2.first) {
        const auto start_arclen = map.domain_arclens_[0] * orig_stroke.length();
        const auto end_arclen = map.domain_arclens_[1] * orig_stroke.length();
        auto proj = Vec2::Empty();
        Float s;
        auto dist = Float(INFINITY);
        if ((int)orig_si == orig_pos1.first) {
          assert(!(start_arclen < query_arclen && query_arclen < end_arclen));
          dist = closest_point_to_own_substroke(orig_stroke, start_arclen, end_arclen,
                                                query_arclen, proj, s);
        } else {
          dist =
            closest_point_substroke(orig_stroke, start_arclen, end_arclen, xy, proj, s);
        }
        if (dist < best_dist) {
          best_dist = dist;
          best_st = StrokeTime((int)orig_si, s / orig_stroke.length());
        }
      }
    }
  }

  force_assert(best_dist < INFINITY && "could not perform a non-spurious reprojection");
  return best_st;
}

void drop_edge(StrokeGraph& graph, const size_t si) {
  const auto hedge = graph.hedge(2 * si);
  const auto twin = graph.hedge(2 * si + 1);
  drop_hedge(graph, hedge.index_);
  drop_hedge(graph, twin.index_);

  // Remove this segment from the mapping structures.
  {
    auto& orig_stroke_indices = graph.strokes2orig_[si];
    // TODO: Technically, this can happen, but it seems unrealistic if we are only calling
    //       this function from the anti-sliver logic...
    force_assert(orig_stroke_indices.size() == 1 && "not implemented yet");
    const auto orig_si = orig_stroke_indices[0];
    auto& mapping = graph.orig2strokes_[orig_si];
    if (mapping.size() == 1) {
      // Drop the only edge for this original stroke; it's like we skipped adding the
      // original stroke in the first place.
      mapping.clear();
    } else if (mapping[0].first == si) {
      // Restart the original stroke at the cut point.
      mapping[1].second.domain_arclens_[0] = 0.0;
      mapping.erase(mapping.begin());
    } else if (mapping.back().first == si) {
      // End the original stroke at the cut point.
      mapping[mapping.size() - 2].second.domain_arclens_[1] = 1.0;
      mapping.pop_back();
    } else {
      // A dangling edge has to be the first or last section of an original stroke.
      std::abort();
    }
    orig_stroke_indices.clear();
    // We could cut the original stroke too, but it's not clear if that is desirable...
    // The current mapping behaviour is to act like we shrank the corresponding end
    // segment when we constructed the edge.
  }

  graph.strokes_[si].clear(); // Mark for deletion.
}

void drop_vertex_and_neighbors(StrokeGraph& graph, const size_t vi,
                               const span<Junction> junctions_to_update) {
  const auto v = graph.vertex(vi);
  force_assert(v.is_active());
  force_assert(v.is_dangling() && "non-dangling case not implemented yet");
  const auto he = v.hedge();
  const auto dest = he.dest();
  force_assert(!dest.is_dangling() && "floating case not implemented yet");

  auto& twin_rec = graph.hedges_[he.twin().index_];
  if (twin_rec.continuity_ != invalid) {
    graph.hedges_[twin_rec.continuity_].continuity_ = invalid;
  }

  const auto prev = he.prev().prev().index_;
  const auto next = he.next().index_;
  graph.hedges_[prev].next_ = next;
  graph.hedges_[next].prev_ = prev;
  graph.vertices_[dest.index_].hedge_ = next;

  if (!graph.faces_.empty()) {
    for (auto& cycle_hi : graph.faces_[he.face_idx()].cycles_) {
      if (cycle_hi == he.index_) {
        cycle_hi = next;
        break;
      }
    }
  }

  // Remove this segment from the mapping structures.
  {
    const auto si = he.stroke_idx();
    for (const auto orig_si : graph.strokes2orig_[si]) {
      auto& mapping = graph.orig2strokes_[orig_si];
      auto valid_domain_begin = Float(0.0);
      auto valid_domain_end = Float(1.0);
      auto remap_domain_begin = NaN;
      auto remap_domain_end = NaN;
      auto remap_range_begin = NaN;
      auto remap_range_end = NaN;
      if (mapping.size() == 1) {
        // Drop the only edge for this original stroke; it's like we skipped adding the
        // original stroke in the first place.
        mapping.clear();
        // Setting an empty interval leads to no valid positions corresponding to this
        // original stroke, which is what we want.
        valid_domain_begin = 1.0;
        valid_domain_end = 0.0;
      }
      if (!mapping.empty() && mapping[0].first == si) {
        // Restart the original stroke at the cut point.
        for (size_t i = 1; i < mapping.size(); ++i) {
          if (mapping[i].first != si) {
            auto& domain = mapping[i].second.domain_arclens_;
            remap_domain_begin = domain[0];
            remap_domain_end = domain[1];

            valid_domain_begin = domain[0];
            domain[0] = 0.0;

            remap_range_begin = domain[0];
            remap_range_end = domain[1];

            mapping.erase(mapping.begin(), mapping.begin() + i);
            break;
          }
        }
      }
      if (!mapping.empty() && mapping.back().first == si) {
        // End the original stroke at the cut point.
        for (auto i = (Index)mapping.size() - 2; i >= 0; --i) {
          if (mapping[i].first != si) {
            auto& domain = mapping[i].second.domain_arclens_;

            remap_domain_begin = domain[0];
            remap_domain_end = domain[1];

            valid_domain_end = domain[1];
            domain[1] = 1.0;

            remap_range_begin = domain[0];
            remap_range_end = domain[1];

            mapping.erase(mapping.begin() + i + 1, mapping.end());
            break;
          }
        }
      }

      force_assert(
        (valid_domain_begin > 0 || valid_domain_end < 1) &&
        "a dangling edge has to be the first or last section of an original stroke");

      for (auto& junc : junctions_to_update) {
        for (auto& stroke_time : junc.points) {
          if (stroke_time.first == (int)orig_si) {
            if (stroke_time.second < valid_domain_begin ||
                stroke_time.second > valid_domain_end) {
              // Signal that this point no longer maps to a valid point.
              stroke_time.second = -1;
            } else if (remap_domain_begin <= stroke_time.second &&
                       stroke_time.second <= remap_domain_end) {
              stroke_time.second =
                std::clamp(linear_remap(stroke_time.second, //
                                        remap_domain_begin, remap_domain_end,
                                        remap_range_begin, remap_range_end),
                           0.0, 1.0);
            }
          }
        }
      }

      // We could cut the original stroke too, but it's not clear if that is desirable...
      // The current mapping behaviour is to act like we shrank the corresponding end
      // segment when we constructed the edge.
    }
  }

  graph.hedges_[he.index_].deactivate();
  twin_rec.deactivate();
  graph.vertices_[vi].deactivate();
  graph.strokes_[he.stroke_idx()].clear(); // Mark for deletion.
}

Float vertex_radius(const VertexView v) {
  const auto he = v.hedge();
  if (!he) { // Normally this would be an assert, but given the way this function is used,
             // it's more convenient if we return 0 here.
    return 0.0;
  }
  auto it = he;
  auto valence = size_t(0);
  auto v_width = 0.0;
  do {
    valence++;
    const auto width = (it.forward() ? it.stroke().width(0) : it.stroke().width(Back));
    v_width = std::max(v_width, width);
    it = it.twin().next();
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
  return 0.5 * v_width;
}

std::string orig2strokes_mapping_repr(const StrokeGraph& graph) {
  auto ss = std::stringstream();
  ss << std::setprecision(4);
  for (size_t orig_si = 0; orig_si < graph.orig2strokes_.size(); ++orig_si) {
    ss << 'o' << orig_si << ": ";
    const auto& piecewise_mappings = graph.orig2strokes_[orig_si];
    for (size_t i = 0; i < piecewise_mappings.size(); ++i) {
      const auto& linear_mapping = piecewise_mappings[i];
      const auto& domain = linear_mapping.second.domain_arclens_;
      const auto& range = linear_mapping.second.range_arclens_;
      ss << "([" << domain[0] << ", " << domain[1] << "] -> " //
         << 'e' << linear_mapping.first << '[' << range[0] << ", " << range[1] << "]) ";
    }
    ss << '\n';
  }
  return ss.str();
}

std::string strokes2orig_mapping_repr(const StrokeGraph& graph) {
  auto ss = std::stringstream();
  ss << std::setprecision(4);
  for (size_t si = 0; si < graph.strokes2orig_.size(); ++si) {
    ss << 'e' << si << ": ";
    const auto& orig_stroke_indices = graph.strokes2orig_[si];
    for (size_t i = 0; i < orig_stroke_indices.size(); ++i) {
      const auto orig_stroke_idx = orig_stroke_indices[i];
      const auto& piecewise_mappings = graph.orig2strokes_[orig_stroke_idx];
      for (size_t j = 0; j < piecewise_mappings.size(); ++j) {
        const auto& linear_mapping = piecewise_mappings[j];
        if (linear_mapping.first == si) {
          const auto& domain = linear_mapping.second.domain_arclens_;
          const auto& range = linear_mapping.second.range_arclens_;
          if (range[0] < range[1]) {
            ss << "([" << range[0] << ", " << range[1] << "] -> " //
               << 'o' << orig_stroke_idx << '[' << domain[0] << ", " << domain[1]
               << "]) ";
          } else {
            ss << "([" << range[1] << ", " << range[0] << "] -> " //
               << 'o' << orig_stroke_idx << '[' << domain[1] << ", " << domain[0]
               << "]) ";
          }
        }
      }
    }
    ss << '\n';
  }
  return ss.str();
}

bool similar_gaps_antisliver( //
  const Stroke& stroke1, const Float s_start, const Float s_end, //
  const Stroke& stroke2, const Float t_start, const Float t_end) {

  (void)stroke1;
  (void)stroke2;
  (void)s_start;
  (void)t_start;
  (void)s_end;
  (void)t_end;

  return false;

#if 0
  constexpr auto n_samples = 4;
  auto max_width = 0.0;
  for (int i = 1; i < n_samples; ++i) {
    const auto u = i / Float(n_samples);
    const auto s = clamped_lerp(s_start, s_end, u);
    const auto t = clamped_lerp(t_start, t_end, u);
    max_width = std::max(max_width, stroke1.width_at(s));
    max_width = std::max(max_width, stroke2.width_at(t));
  }
  const auto sq_tol = 1.5 * max_width * max_width;
  for (int i = 1; i < n_samples; ++i) {
    const auto u = i / Float(n_samples);
    const auto s = clamped_lerp(s_start, s_end, u);
    const auto t = clamped_lerp(t_start, t_end, u);
    if ((stroke1.pos(s) - stroke2.pos(t)).squaredNorm() >= sq_tol) {
      return false;
    }
  }
  return true;
#endif
}

bool is_corner(const StrokeGraph& graph, size_t v_idx, bool include_t) {
  auto v = graph.vertex(v_idx);
  std::set<std::pair<size_t, size_t>> continuous_pairs;
  size_t discontinous_valence = 0;

  auto valence = 0;
  const auto he = v.hedge();
  auto it = he;
  do {
    if (it.continuity_edge().is_valid())
      continuous_pairs.insert(std::make_pair(std::min(it.index_, it.continuity()),
                                             std::max(it.index_, it.continuity())));
    else
      discontinous_valence++;
    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);

  // Check if is dissolveable
  if (discontinous_valence == 2 && continuous_pairs.empty())
    return !should_dissolve_vertex(v);

  // TODO: Do we consider a T-junction a corner?
  return (discontinous_valence > 1) ||
         (include_t && discontinous_valence > 0 && !continuous_pairs.empty());
}

void get_corner_original_positions(const StrokeGraph& graph, size_t v_idx,
                                   std::vector<Vec2>& positions, bool include_t) {
  if (!is_corner(graph, v_idx, include_t))
    return;

  auto v = graph.vertex(v_idx);
  auto valence = 0;
  const auto he = v.hedge();
  auto it = he;
  do {
    if (!it.continuity_edge().is_valid()) {
      auto endp = StrokeTime((int)it.stroke_idx(), it.forward() ? 0.0 : 1.0);
      const auto ok = convert_strokes2orig(graph, endp);
      force_assert(ok && "couldn't map from strokes to orig");
      positions.emplace_back(graph.orig_strokes_[endp.first].pos_norm(endp.second));
    }
    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
}

/**
 * Get the original stroke positions on both sides of a mapping discontinuity.
 *
 * si is a stroke index in the graph (as opposed to an original stroke index).
 * norm_stroke_arclen is a normalized (between [0, 1]) arc length position on stroke si.
 * (si, norm_stroke_arclen) needs to correspond exactly to one of the (stroke index,
 * range) values in the orig2strokes_ mapping for some original stroke.
 */
static std::pair<StrokeTime, StrokeTime>
original_stroke_positions_on_both_sides(const StrokeGraph& graph, const size_t si,
                                        const Float norm_stroke_arclen) {
  assert(si < graph.strokes_.size() && "expected a stroke index");

  auto before = StrokeTime(-1, 0.0);
  auto after = StrokeTime(-1, 0.0);
  for (const auto orig_si : graph.strokes2orig_[si]) {
    for (const auto& map : graph.orig2strokes_[orig_si]) {
      if (map.first == si) {
        const auto& range = map.second.range_arclens_;
        const auto& domain = map.second.domain_arclens_;
        if (range[0] == norm_stroke_arclen) {
          assert(after.first == -1);
          after = StrokeTime((int)orig_si, domain[0]);
        }
        if (range[1] == norm_stroke_arclen) {
          assert(before.first == -1);
          before = StrokeTime((int)orig_si, domain[1]);
        }
      }
    }
  }
  assert(before.first != -1 || after.first != -1);
  return {before, after};
}

Connections::Connections(const StrokeGraph& graph)
  : max_orig_stroke_idx_((int)graph.orig2strokes_.size() - 1) {

  // Find connections at vertices.
  auto buffer = std::vector<StrokeTime>();
  for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
    const auto v = graph.vertex(vi);
    if (v.is_active()) {
      // Collect connections at this vertex.
      const auto he = v.hedge();
      auto it = he;
      do {
        if (it.continuity() >= it.index_) { // Don't add continuous points twice.
          auto endp = StrokeTime((int)it.stroke_idx(), it.forward() ? 0.0 : 1.0);
          const auto ok = convert_strokes2orig(graph, endp);
          force_assert(ok && "couldn't map from strokes to orig");
          buffer.push_back(endp);
        }
        it = it.twin().next();
      } while (it != he);

      // Insert all connected pairs.
      for (size_t i = 0; i < buffer.size(); ++i) {
        for (size_t j = i + 1; j < buffer.size(); ++j) {
          auto location1 = buffer[i];
          auto location2 = buffer[j];
          if (location1.second != 0 && location1.second != 1) {
            // Try to put the stroke interior second.
            std::swap(location1, location2);
          }

          // We only care if one of them corresponds to an original endpoint.
          if (location1.second == 0 || location1.second == 1) {
            if (location2.second != 0 && location2.second != 1) {
              location2.second = 0.5; // For assessing with a stroke interior, we do not
                                      // care where on the stroke the connection happens.
            } else {
              // Allow bidirectional lookup for end-end.
              connections_.emplace_back(location2, location1);
            }
            connections_.emplace_back(location1, location2);
          }
        }
      }
      buffer.clear();
    }
  }

  // Find connections at dissolved ex-vertices.
  for (size_t orig_si = 0; orig_si < graph.orig2strokes_.size(); ++orig_si) {
    const auto& o2s = graph.orig2strokes_[orig_si];
    if (!o2s.empty()) {
      const auto first_range_start = o2s[0].second.range_arclens_[0];
      if (first_range_start != 0 && first_range_start != 1) {
        const auto si = o2s[0].first;
        const auto [stroke_time_before, stroke_time_after] =
          original_stroke_positions_on_both_sides(graph, si, first_range_start);
        assert(stroke_time_after.first == (int)orig_si);
        assert(stroke_time_after.second == 0);
        connections_.emplace_back(stroke_time_before, stroke_time_after);
      }

      const auto last_range_end = o2s.back().second.range_arclens_.back();
      if (last_range_end != 0 && last_range_end != 1) {
        const auto si = o2s.back().first;
        const auto [stroke_time_before, stroke_time_after] =
          original_stroke_positions_on_both_sides(graph, si, last_range_end);
        assert(stroke_time_before.first == (int)orig_si);
        assert(stroke_time_before.second == 1);
        connections_.emplace_back(stroke_time_before, stroke_time_after);
      }
    }
  }

  std::sort(connections_.begin(), connections_.end());
  // Remove duplicates.  This can happen if we have a connection at a connection at a
  // self-intersection point.
  connections_.erase(std::unique(connections_.begin(), connections_.end()),
                     connections_.end());
}

Connections::Connections(const EnvelopeBVH& strokes)
  : max_orig_stroke_idx_((int)strokes.nodes.size() - 1) {

  for (auto& node : strokes.nodes) {
    node.geometry->ensure_arclengths();
  }

  const auto n_strokes = strokes.nodes.size();
  for (size_t i = 0; i < n_strokes; ++i) {
    const auto& stroke = *strokes.nodes[i].geometry;
    if (stroke.size() > 1) {
      for (size_t j = 0; j < n_strokes; ++j) {
        for (const auto arclen : {0.0, 1.0}) {
          Float closest_env_dist = 0.0;
          auto [dist, other_arclen] = find_snap_point(
            stroke, /*head=*/arclen == 0.0, strokes.nodes[j], &closest_env_dist);
          if (std::isfinite(dist)) {
            other_arclen /= strokes.nodes[j].geometry->length();
            if (other_arclen != 0.0 && other_arclen != 1.0) {
              // For assessing with a stroke interior, we do not care where on the stroke
              // the connection happens.
              other_arclen = 0.5;
            }
            connections_.emplace_back(StrokeTime((int)i, arclen),
                                      StrokeTime((int)j, other_arclen));
          }
        }
      }
    }
  }

  std::sort(connections_.begin(), connections_.end());
  // Remove duplicates.
  connections_.erase(std::unique(connections_.begin(), connections_.end()),
                     connections_.end());
}

span<const Connections::Connection>
Connections::associated_connections(const StrokeTime endp) const {
  force_assert(endp.second == 0 || endp.second == 1);
  assert(endp.first <= max_orig_stroke_idx_ && "are you using original stroke indexing?");

  using index_t = std::numeric_limits<decltype(StrokeTime::first)>;
  const auto min_si = index_t::min();
  const auto max_si = index_t::max();

  const auto key_begin = Connection(endp, StrokeTime(min_si, -infinity));
  const auto begin =
    std::lower_bound(connections_.begin(), connections_.end(), key_begin);
  if (begin == connections_.end() || begin->first != endp) {
    return {nullptr, 0};
  }

  const auto key_end = Connection(endp, StrokeTime(max_si, infinity));
  const auto end = std::upper_bound(connections_.begin(), connections_.end(), key_end);
  return {begin, end};
}

int Connections::largest_connected_stroke_idx(const StrokeTime endp) const {
  force_assert(endp.second == 0 || endp.second == 1);
  assert(endp.first <= max_orig_stroke_idx_ && "are you using original stroke indexing?");

  auto max_si = endp.first;
  for (const auto& [_, other] : associated_connections(endp)) {
    max_si = std::max(max_si, other.first);
  }
  return max_si;
}

bool Connections::contains(StrokeTime a, StrokeTime b) const {
  assert(a.first <= max_orig_stroke_idx_ && b.first <= max_orig_stroke_idx_ &&
         "are you using original stroke indexing?");
  if (a.second != 0 && a.second != 1) {
    std::swap(a, b);
    force_assert(
      (a.second == 0 || a.second == 1) &&
      "one of the locations of a connection lookup must be an original endpoint");
  }

  if (b.second == 0 || b.second == 1) {
    // Endpoint-endpoint lookup.
    return std::binary_search(connections_.begin(), connections_.end(),
                              a < b ? std::make_pair(a, b) : std::make_pair(b, a));
  } else {
    // Endpoint and original stroke interior lookup.
    b.second = 0.5;
    return std::binary_search(connections_.begin(), connections_.end(),
                              std::make_pair(a, b));
  }
}

void Connections::insert(const span<const StrokeTime> a, const span<const StrokeTime> b) {
  force_assert(a.size() == b.size());
  for (size_t i = 0; i < a.size(); ++i) {
    const auto p = a[i];
    force_assert((p.second == 0 || p.second == 1) &&
                 "array a can only contain endpoints");
    auto q = b[i];
    if (q.second != 0 && q.second != 1) {
      q.second = 0.5;
    } else {
      connections_.emplace_back(q, p);
    }
    connections_.emplace_back(p, q);

    max_orig_stroke_idx_ = std::max({max_orig_stroke_idx_, p.first, q.first});
  }
  std::sort(connections_.begin(), connections_.end());
  // Remove duplicates.
  connections_.erase(std::unique(connections_.begin(), connections_.end()),
                     connections_.end());
}

Connections Connections::subset_involving(size_t si) const {
  Connections new_intersection_set;
  new_intersection_set.connections_.reserve(new_intersection_set.connections().size());

  // Filter out allowed/enforced set given the newest original stroke index
  for (const auto& conn : connections()) {
    if (conn.first.first == si || conn.second.first == si)
      new_intersection_set.connections_.emplace_back(conn);
  }

  return new_intersection_set;
}

} // namespace sketching
