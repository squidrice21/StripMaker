#include "stroke_graph.h"

#include "deform.h"
#include "detail/alloca.h"
#include "detail/suppress_warning.h"
#include "detail/util.h"
#include "eigen_compat.h"
#include "force_assert.h"
#include "intersect.h"
#include "junction.h"
#include "junction_type.h"
#include "resample.h"
#include "sketching.h"
#include "stroke_graph_extra.h"

#include <unordered_set>

namespace sketching {

namespace {

using VertexView = StrokeGraph::VertexView;
using HedgeView = StrokeGraph::HedgeView;
using FaceView = StrokeGraph::FaceView;
constexpr auto invalid = StrokeGraph::invalid;

std::pair<size_t, size_t> make_twins(std::vector<StrokeGraph::HedgeRecord>& hedges) {
  hedges.emplace_back(); // Do not use this reference, it will become invalid if the
                         // vector resizes on the next line.
  hedges.emplace_back();
  const auto new_edge_idx = hedges.size() - 2;
  const auto new_twin_idx = new_edge_idx + 1;
  return {new_edge_idx, new_twin_idx};
}

/**
 * Determine a scale for which it is safe to compute outgoing tangents for the purpose of
 * finding a consistent, sensible sort angle.
 */
double squared_scale(const StrokeGraph& graph, const span<size_t> hedge_indices) {
  auto closest_furthest = infinity;
  const auto p = graph.hedge(hedge_indices[0]).origin().pos();
  for (const auto hi : hedge_indices) {
    const auto& he = graph.hedge(hi);
    const auto& stroke = he.stroke();
    auto furthest = -infinity;
    for (Index i = 0; i < stroke.size(); ++i) {
      furthest = std::max((stroke.xy(i) - p).squaredNorm(), furthest);
    }
    closest_furthest = std::min(furthest, closest_furthest);
  }
  // Two geometrically distinct hedges might end at the same point, so we want to stop
  // just before that point.
  return 0.9 * closest_furthest;
}

double squared_scale(const VertexView v) {
  auto valence = 0;
  auto closest_furthest = infinity;
  const auto p = v.pos();
  const auto he = v.hedge();
  auto it = he;
  do {
    const auto& stroke = it.stroke();
    assert(stroke.size() > 1);
    auto furthest = -infinity;
    for (Index i = 0; i < stroke.size(); ++i) {
      furthest = std::max((stroke.xy(i) - p).squaredNorm(), furthest);
    }
    closest_furthest = std::min(furthest, closest_furthest);
    valence++;
    it = it.twin().next();
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
  // Two geometrically distinct hedges might end at the same point, so we want to stop
  // just before that point.
  return 0.9 * closest_furthest;
}

double tangent_angle(const HedgeView he, const Float squared_scale) {
  const auto p = he.origin().pos();
  const auto& stroke = he.stroke();
  if (he.forward()) {
    for (Index j = 1; j < stroke.size(); ++j) {
      const auto last = j - 1;
      const auto dist = (stroke.xy(j) - p).squaredNorm();
      if (dist >= squared_scale) {
        const auto last_dist = (stroke.xy(last) - p).norm();
        const auto u =
          (std::sqrt(squared_scale) - last_dist) / (std::sqrt(dist) - last_dist);
        assert(0.0 <= u && u <= 1.0);
        const auto q = lerp(stroke.xy(last), stroke.xy(j), u);
        const auto tangent = q - p;
        assert(tangent != Vec2(0, 0));
        return std::atan2(tangent.y(), tangent.x());
      }
    }
  } else {
    for (Index j = stroke.size() - 2; j >= 0; --j) {
      const auto last = j + 1;
      const auto dist = (stroke.xy(j) - p).squaredNorm();
      if (dist >= squared_scale) {
        const auto last_dist = (stroke.xy(last) - p).norm();
        const auto u =
          (std::sqrt(squared_scale) - last_dist) / (std::sqrt(dist) - last_dist);
        assert(0.0 <= u && u <= 1.0);
        const auto q = lerp(stroke.xy(last), stroke.xy(j), u);
        const auto tangent = q - p;
        assert(tangent != Vec2(0, 0));
        return std::atan2(tangent.y(), tangent.x());
      }
    }
  }
  std::abort(); // Shouldn't happen.
}

/// Use this for consistently sorting by tangent angle if some tangent angles are equal.
/// This avoids Mobius-like cycles that can occur otherwise.
std::pair<double, double> augmented_tangent_angle(const HedgeView he,
                                                  const double squared_scale) {
  return {tangent_angle(he, squared_scale),
          (he.forward() ? double(he.index_) : -double(he.index_))};
}

/**
 * Sort a vertex star by angle.  Robust to tiny intersections near the beginning of edges.
 *
 * @param graph
 * @param hedge_indices Indices of outbound hedges which share an origin vertex.  Will be
 *                      sorted by angle when this function returns.
 * @param angle_buffer Buffer of same size as `hedge_indices`. Buffer contents are
 *                     overwritten.
 * @param index_buffer Buffer of same size as `hedge_indices`. Buffer contents are
 *                     overwritten.
 */
void sort_outbound_hedges(const StrokeGraph& graph, const span<size_t> hedge_indices,
                          const span<std::pair<double, double>> angle_buffer,
                          const span<size_t> index_buffer) {
  assert(hedge_indices.size() == angle_buffer.size());
  assert(hedge_indices.size() == index_buffer.size());
  if (hedge_indices.empty())
    return;

  const auto n = hedge_indices.size();
  for (size_t i = 0; i < n; ++i) {
    index_buffer[i] = i;
  }

  // There may be intersections between the edges near the beginning, so determine a safe
  // scale.
  const auto sqscale = squared_scale(graph, hedge_indices);

  // Compute angles.
  for (size_t i = 0; i < n; ++i) {
    const auto& he = graph.hedge(hedge_indices[i]);
    angle_buffer[i] = augmented_tangent_angle(he, sqscale);
  }

  const auto* const angles = angle_buffer.data();
  std::sort(index_buffer.begin(), index_buffer.end(),
            [=](size_t a, size_t b) { return angles[a] < angles[b]; });
  for (size_t i = 0; i < n; ++i) {
    index_buffer[i] = hedge_indices[index_buffer[i]];
  }
  for (size_t i = 0; i < n; ++i) {
    hedge_indices[i] = index_buffer[i];
  }
}

double turns_det(const StrokeGraph::HedgeView he) {
  const auto& s = he.stroke();
  const auto prev = he.prev();
  const auto& pr_s = prev.stroke();
  const Vec2 a = (prev.forward() ? pr_s.xy(pr_s.size() - 2) : pr_s.xy(1));
  const Vec2 b = (prev.forward() ? pr_s.xy(pr_s.size() - 1) : pr_s.xy(0));
  const Vec2 c = (he.forward() ? s.xy(1) : s.xy(s.size() - 2));
  return counter_clockwise(a, b, c);
}

struct ContainsResult {
  std::int8_t res_;
  explicit constexpr operator bool() const { return res_ == 1; }
  constexpr bool operator==(ContainsResult other) const { return res_ == other.res_; }
  constexpr bool operator!=(ContainsResult other) const { return !(*this == other); }
};
constexpr auto Inside = ContainsResult{1};
constexpr auto OnBoundary = ContainsResult{0};
constexpr auto Outside = ContainsResult{-1};
static_assert(Inside && !OnBoundary && !Outside);

/**
 * Check if point p is inside the closed polygon formed by the cycle obtained by
 * traversing `he = he.next()` repeatedly.  Orientation does not matter.
 */
ContainsResult point_in_cycle(const Vec2 p, const HedgeView he) {
  // Based on https://github.com/rowanwins/point-in-polygon-hao (MIT).

  auto it = he;
  auto* current_stroke = &he.stroke();
  auto current_stroke_size = current_stroke->size();
  auto currentP = (he.forward() ? current_stroke->xy(0) : current_stroke->xy(Back));
  auto u1 = currentP.x_ - p.x_;
  auto v1 = currentP.y_ - p.y_;
  auto k = 0;
  do {
    for (Index i = 0; i < current_stroke_size - 1; ++i) {
      // Find the next point on the cycle.
      auto nextP = Vec2::Empty();
      const auto next_i = (it.forward() ? i + 1 : current_stroke_size - (i + 1) - 1);
      nextP = current_stroke->xy(next_i);

      const auto v2 = nextP.y_ - p.y_;

      if ((v1 < 0 && v2 < 0) || (v1 > 0 && v2 > 0)) {
        currentP = nextP;
        v1 = v2;
        u1 = currentP.x_ - p.x_;
        continue;
      }

      const auto u2 = nextP.x_ - p.x_;

      if (v2 > 0 && v1 <= 0) {
        const auto f = (u1 * v2) - (u2 * v1);
        if (f > 0.0) {
          k++;
        } else if (f == 0) {
          return OnBoundary;
        }
      } else if (v1 > 0 && v2 <= 0) {
        const auto f = (u1 * v2) - (u2 * v1);
        if (f < 0.0) {
          k++;
        } else if (f == 0) {
          return OnBoundary;
        }
      } else if (v2 == 0 && v1 < 0) {
        const auto f = (u1 * v2) - (u2 * v1);
        if (f == 0.0) {
          return OnBoundary;
        }
      } else if (v1 == 0 && v2 < 0) {
        const auto f = u1 * v2 - u2 * v1;
        if (f == 0.0) {
          return OnBoundary;
        }
      } else if (v1 == 0 && v2 == 0) {
        if ((u2 <= 0 && u1 >= 0) || (u1 <= 0 && u2 >= 0)) {
          return OnBoundary;
        }
      }
      currentP = nextP;
      v1 = v2;
      u1 = u2;
    }

    it = it.next();
    current_stroke = &it.stroke();
    current_stroke_size = current_stroke->size();
  } while (it != he);

  return (k % 2 == 0 ? Outside : Inside);
}

ContainsResult point_in_face(const Vec2 p, const FaceView f) {
  if (f.cycles().empty()) {
    return Outside;
  }
  const auto& graph = *f.graph_;
  // Outer polygon.
  const auto boundary_condition = point_in_cycle(p, graph.hedge(f.cycles()[0]));
  if (boundary_condition != Inside) {
    return boundary_condition; // Outside or on boundary.
  }
  for (size_t i = 1; i < f.cycles().size(); ++i) {
    const auto condition = point_in_cycle(p, graph.hedge(f.cycles()[i]));
    if (condition == Inside) {
      return Outside; // Inside a hole -> outside.
    } else if (condition == OnBoundary) {
      return OnBoundary;
    }
  }
  return Inside;
}

} // namespace

bool point_in_face(const FaceView f, const Vec2 p) {
  return point_in_face(p, f) != Outside;
}

void check_continuity(const StrokeGraph& graph) {
  for (size_t hi = 0; hi < graph.hedges_.size(); hi++) {
    const auto& he = graph.hedges_[hi];
    if (he.is_active() && he.continuity_ != invalid) {
      auto& cont_he = graph.hedges_[he.continuity_];
      force_assert(cont_he.is_active());
      force_assert(cont_he.continuity_ == hi);
      force_assert(cont_he.origin_ == he.origin_);
    }
  }
}

void check_consistency(const StrokeGraph& graph) {
  force_assert(graph.hedges_.size() == 2 * graph.strokes_.size());
  force_assert(graph.strokes2vid_.size() == graph.strokes_.size());
  for (size_t vi = 0; vi < graph.vertices_.size(); ++vi) {
    const auto v = graph.vertex(vi);
    if (v.is_active()) {
      const auto hedge = v.hedge();
      force_assert(hedge.index_ < graph.hedges_.size());
      force_assert(hedge.is_valid());
      force_assert(hedge.origin() == v);
    }
  }
  for (size_t hi = 0; hi < graph.hedges_.size(); ++hi) {
    const auto hedge = graph.hedge(hi);
    if (hedge.is_active()) {
      force_assert(hedge.is_valid());
      force_assert(hedge.next().prev() == hedge);
      force_assert(hedge.prev().next() == hedge);
      force_assert(hedge.origin().index_ < graph.vertices_.size());
      force_assert(hedge.origin().is_valid());
      const auto v = hedge.origin();
      if (hedge.forward()) {
        force_assert((v.pos() - hedge.stroke().xy(0)).norm() < 1e-6);
      } else {
        force_assert((v.pos() - hedge.stroke().xy(Back)).norm() < 1e-6);
      }
      force_assert(hedge.face_idx() != StrokeGraph::invalid);
      force_assert(hedge.face_idx() == hedge.next().face_idx());
      force_assert(hedge.face_idx() == hedge.prev().face_idx());
    }
  }
  for (size_t fi = 0; fi < graph.faces_.size(); ++fi) {
    for (const auto hi : graph.faces_[fi].cycles_) {
      const auto he = graph.hedge(hi);
      auto it = he;
      do {
        force_assert(it.face().index_ == fi);
        it = it.next();
      } while (it != he);
    }
  }
}

void construct_faces(StrokeGraph& graph) {
  const auto bounding_boxes = graph.bvh_.centerline_bbs();
  assert(bounding_boxes.size() == graph.strokes_.size());
  auto component_ids = std::make_unique<int[]>(graph.strokes_.size());
  const auto num_components =
    connected_components_edge(graph, {component_ids.get(), graph.strokes_.size()});
  if (num_components == 0) {
    return;
  }

  struct CComponent {
    int id; // Use this to index into `component_ids`.
    Vec2 leftmost = {infinity, infinity};
    size_t leftmost_hi = invalid;
    Index leftmost_pi = -1;
    BoundingBox bb;
  };
  auto ccomponents = std::make_unique<CComponent[]>(num_components);
  for (int i = 0; i < num_components; ++i) {
    ccomponents[i].id = i;
  }
  for (size_t si = 0; si < graph.strokes_.size(); ++si) {
    const auto hi = 2 * si;
    const auto he = graph.hedge(hi);
    if (graph.hedges_[hi].is_active()) {
      auto& ccomp = ccomponents[component_ids[si]];
      const auto& s = he.stroke();
      ccomp.bb.unite(bounding_boxes[si]);
      for (Index pi = 0; pi < s.size(); ++pi) {
        const auto p = Vec2{s.xy(pi)};
        if (p < ccomp.leftmost) {
          ccomp.leftmost = p;
          ccomp.leftmost_hi = hi;
          ccomp.leftmost_pi = pi;
        }
      }
    }
  }

  // To ensure outer components are processed before inner components, we sort them by
  // leftmost position.
  std::sort(
    &ccomponents[0], &ccomponents[num_components],
    [](const CComponent& c1, const CComponent& c2) { return c1.leftmost < c2.leftmost; });

  graph.faces_.emplace_back(); // Create boundary face.

  for (int i = 0; i < num_components; ++i) {
    // Ensure each half-edge lies on the boundary's outside.
    auto& ccomp = ccomponents[i];
    auto he = graph.hedge(ccomp.leftmost_hi); // Our desired (not-yet) boundary half-edge.
    assert(he.forward()); // Remember, by construction `he` faces forward.
    const auto& stroke = he.stroke();
    if (ccomp.leftmost_pi == 0 || ccomp.leftmost_pi == stroke.size() - 1) {
      if (ccomp.leftmost_pi == stroke.size() - 1) {
        he = he.twin();
      }
      // Go around in a circle until we either get back where we started or find a
      // clockwise cycle.  If we get back where we started, it's probably a colinear
      // stroke, which can be a boundary too.
      auto it = he;
      do {
        if (turns_det(it) < 0.0) {
          he = it;
          break;
        }
        it = it.twin().next();
      } while (it != he);
    } else {
      const auto det =
        counter_clockwise(stroke.xy(ccomp.leftmost_pi - 1), stroke.xy(ccomp.leftmost_pi),
                          stroke.xy(ccomp.leftmost_pi + 1));
      if (det > 0.0) {
        // Half-edge turns counter-clockwise, but we want the boundary to turn clockwise.
        he = he.twin();
      }
    }

    // Determine if this component is inside any of the previous components.
    // If not, then it is part of the boundary.
    auto face_idx = graph.boundary_face_;
    for (int j = 0; j < i; ++j) {
      if (ccomponents[j].bb.contains(ccomp.bb)) {
        const auto p = ccomp.leftmost; // Pick an arbitrary point in ccomp.
        for (size_t fi = 1; fi < graph.faces_.size(); ++fi) {
          const auto face = graph.face(fi);
          auto relevant = false;
          for (const auto hi : face.cycles()) {
            if (component_ids[hi / 2] == ccomponents[j].id) {
              relevant = true;
              break;
            }
          }
          if (relevant && point_in_face(p, graph.face(fi))) {
            face_idx = fi;
            break;
          }
        }
      }
    }
    graph.faces_[face_idx].cycles_.push_back(he.index_);
    {
      auto it = he;
      do {
        assert(graph.hedges_[it.index_].face_ == StrokeGraph::invalid);
        graph.hedges_[it.index_].face_ = face_idx;
        it = it.next();
      } while (it != he);
    }

    // Assign all the other faces in this component.
    for (size_t hi = 0; hi < graph.hedges_.size(); ++hi) {
      const auto& he_r = graph.hedges_[hi];
      if (he_r.is_active() && component_ids[hi / 2] == ccomp.id &&
          he_r.face_ == invalid) {
        const auto it_start = graph.hedge(hi);
        const auto new_face_idx = graph.faces_.size();
        graph.faces_.emplace_back().cycles_.push_back(hi);
        auto it = it_start;
        do {
          assert(graph.hedges_[it.index_].face_ == StrokeGraph::invalid);
          graph.hedges_[it.index_].face_ = new_face_idx;
          it = it.next();
        } while (it != it_start);
      }
    }
  }
}

void insert_into_star(StrokeGraph& graph, const size_t hi, const size_t vi) {
  const auto he = graph.hedge(hi);
  const auto v = graph.vertex(vi);
  auto sqscale = squared_scale(v);
  {
    const auto& stroke = he.stroke();
    auto furthest = -infinity;
    for (Index i = 0; i < stroke.size(); ++i) {
      furthest = std::max((stroke.xy(i) - v.pos()).squaredNorm(), furthest);
    }
    sqscale = std::min(0.9 * furthest, sqscale);
  }
  assert(sqscale > 0);

  const auto ta = augmented_tangent_angle(he, sqscale);
  auto valence = 0;
  const auto begin = v.hedge();
  auto it = begin;
  do {
    auto next = it.twin().next();
    const auto ta_it = augmented_tangent_angle(it, sqscale);
    const auto ta_next = augmented_tangent_angle(next, sqscale);
    if ((ta_it > ta && ta > ta_next) ||
        (ta_it <= ta_next && (ta > ta_next || ta < ta_it))) {
      // Place `he` in-between `it` and `next` in the umbrella around `v`.
      auto& next_rec = graph.hedges_[next.index_];
      auto& it_twin_rec = graph.hedges_[it.twin().index_];
      auto& edge_rec = graph.hedges_[he.index_];
      auto& twin_rec = graph.hedges_[he.twin().index_];
      it_twin_rec.next_ = he.index_;
      edge_rec.prev_ = it.twin().index_;
      twin_rec.next_ = next.index_;
      next_rec.prev_ = he.twin().index_;
      return;
    }

    valence++;
    it = it.twin().next();
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != begin);
}

static VertexView insert_vertex_in_middle_of_edge( //
  StrokeGraph& graph, const size_t stroke_idx, const Float arclen) {

  auto& stroke = graph.strokes_[stroke_idx];
  // We rely on stroke i corresponding to hedges 2i and 2i + 1.
  // This invariant must be upheld if you wish to use this function.
  const auto edge = HedgeView(&graph, 2 * stroke_idx);
  const auto twin = HedgeView(&graph, 2 * stroke_idx + 1);
  force_assert(&edge.stroke() == &stroke);
  force_assert(&twin.stroke() == &stroke);
  force_assert(edge.forward() && !twin.forward());

  // Save the face indices for setting the faces of the new edges.
  auto edge_face = edge.face_idx();
  auto twin_face = twin.face_idx();

  auto slice = stroke.fractional_slice(arclen, stroke.length());
  if (slice.size() <= 1) {
    return VertexView(); // Would create an invalid edge.
  }
  Float split = arclen / stroke.length();
  auto tmp_stroke = stroke.clone();
  tmp_stroke.trim(0.0, arclen);
  if (tmp_stroke.size() <= 1) {
    return VertexView(); // Would create an invalid edge.
  }
  slice.compute_arclengths();
  tmp_stroke.compute_arclengths();
  force_assert(slice.length() > 0.0);
  force_assert(tmp_stroke.length() > 0.0);
  stroke = std::move(tmp_stroke);

  graph.bvh_.full_update(stroke_idx);
  const auto new_stroke_idx = (size_t)graph.bvh_.add(std::move(slice));
  auto& new_stroke = graph.strokes_[new_stroke_idx];
  split_stroke_mapping(graph, stroke_idx, split, stroke_idx, new_stroke_idx);
  // `stroke` can now be a dangling reference. Don't use it any more.

  if (!(edge.flags() & StrokeGraph::HedgeRecord::Bridge)) {
    StrokeTime test_end1((int)stroke_idx, 0.0), test_end2((int)stroke_idx, 1.0),
      test_end3((int)new_stroke_idx, 0.0), test_end4((int)new_stroke_idx, 1.0);
    const auto ok1 = convert_strokes2orig(graph, test_end1);
    const auto ok2 = convert_strokes2orig(graph, test_end2);
    const auto ok3 = convert_strokes2orig(graph, test_end3);
    const auto ok4 = convert_strokes2orig(graph, test_end4);
    force_assert(ok1 && ok2 && ok3 && ok4 && "graph corrupted");
  }

  size_t twin_cont = twin.continuity();

  const auto [new_edge_idx, new_twin_idx] = make_twins(graph.hedges_);
  auto& new_edge = graph.hedges_[new_edge_idx];
  auto& new_twin = graph.hedges_[new_twin_idx];
  graph.hedges_[edge.next().index_].prev_ = new_edge_idx;
  new_edge.next_ = edge.next().index_;
  new_edge.prev_ = edge.index_;
  new_twin.next_ = twin.index_;
  new_twin.prev_ = twin.prev().index_;
  new_edge.flags_ = edge.flags();
  new_twin.flags_ = twin.flags();
  graph.hedges_[twin.prev().index_].next_ = new_twin_idx;

  new_edge.face_ = edge_face;
  new_twin.face_ = twin_face;

  const auto old_dest = twin.origin();
  const auto new_vertex_idx = graph.vertices_.size();
  auto& new_vertex = graph.vertices_.emplace_back(new_stroke.xy(0));
  graph.vertices_[old_dest.index_].hedge_ = new_twin_idx;
  new_vertex.hedge_ = new_edge_idx;
  new_edge.origin_ = new_vertex_idx;
  new_twin.origin_ = old_dest.index_;

  if (twin_cont != StrokeGraph::invalid) {
    graph.hedges_[twin_cont].continuity_ = new_twin_idx;
    new_twin.continuity_ = twin_cont;
  }

  auto& edge_rec = graph.hedges_[edge.index_];
  auto& twin_rec = graph.hedges_[twin.index_];
  edge_rec.next_ = new_edge_idx;
  twin_rec.prev_ = new_twin_idx;
  twin_rec.origin_ = new_vertex_idx;

  twin_rec.continuity_ = new_edge_idx;
  new_edge.continuity_ = twin.index_;

  return graph.vertex(new_vertex_idx);
}

bool add_vertex_precise(StrokeGraph& graph, const size_t stroke_idx, const Float arclen,
                        StrokeGraph::VertexView& new_vertexview) {
  auto& stroke = graph.strokes_[stroke_idx];
  force_assert(0.0 <= arclen && arclen <= stroke.length());
  const auto edge = HedgeView(&graph, 2 * stroke_idx);
  force_assert(&edge.stroke() == &stroke);
  force_assert(edge.forward());

  if (arclen == 0.0) {
    new_vertexview = edge.origin();
    return false;
  } else if (arclen == stroke.length()) {
    new_vertexview = edge.dest();
    return false;
  }

  new_vertexview = insert_vertex_in_middle_of_edge(graph, stroke_idx, arclen);
  return new_vertexview.is_valid();
}

bool add_vertex(StrokeGraph& graph, const size_t stroke_idx, Float arclen,
                VertexView& new_vertexview) {
  auto& stroke = graph.strokes_[stroke_idx];
  arclen = std::clamp(arclen, 0.0, stroke.length());
  const auto edge = HedgeView(&graph, 2 * stroke_idx);
  force_assert(&edge.stroke() == &stroke);
  force_assert(edge.forward());

  // Too close, re-use existing vertex.
  const auto fudge_amount = std::max(0.5 * stroke.width_at(arclen), 1e-5);
  if (arclen <= std::min(fudge_amount, 0.5 * stroke.length())) {
    new_vertexview = edge.origin();
    return false;
  } else if (arclen >= stroke.length() - fudge_amount) {
    new_vertexview = edge.dest();
    return false;
  }

  new_vertexview = insert_vertex_in_middle_of_edge(graph, stroke_idx, arclen);
  return new_vertexview.is_valid();
}

size_t StrokeGraph::VertexView::valence() const {
  const auto he = hedge();
  auto it = he;
  auto valence = size_t(0);
  do {
    valence++;
    it = it.twin().next();
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);
  return valence;
}

bool StrokeGraph::VertexView::is_dangling() const {
  const auto he = hedge();
  return he == he.twin().next();
}

size_t HedgeView::orig_stroke_idx() const {
  VertexView v = origin();
  const Stroke& s = stroke();
  size_t si = stroke_idx();

  // Warning: this function can only retrieve a single original stroke if we have an
  // edge with the same vertex as the origin and destination
  // TODO: let the caller could specify where on the edge they want to look up the
  // original index.
  Float s_time =
    ((s.pos_norm(0) - v.pos()).norm() < (s.pos_norm(1) - v.pos()).norm()) ? 0 : 1;
  for (auto orig_si : graph_->strokes2orig_[stroke_idx()]) {
    const auto& orig_mappings = graph_->orig2strokes_[orig_si];
    for (auto orig_itr = orig_mappings.begin(); orig_itr != orig_mappings.end();
         ++orig_itr) {
      if (orig_itr->first != si)
        continue;
      const auto& mapping = orig_itr->second;
      auto range_s = std::min(mapping.range_arclens_[0], mapping.range_arclens_[1]);
      auto range_l = std::max(mapping.range_arclens_[0], mapping.range_arclens_[1]);
      if (s_time >= range_s && s_time <= range_l) {
        return orig_si;
      }
    }
  }

  // Must find the original stroke
  std::abort();
}

size_t StrokeGraph::FaceView::n_neighbors() const {
  auto seen = std::unordered_set<size_t>();
  auto iterations = 0;
  for (const auto hi : cycles()) {
    const auto he = graph_->hedge(hi);
    auto it = he;
    do {
      seen.insert(it.twin().face_idx());
      it = it.next();
      iterations++;
      force_assert(iterations < 1024 && "likely infinite loop found");
    } while (it != he);
  }
  seen.erase(index_);
  return seen.size();
}

size_t StrokeGraph::FaceView::n_edges() const {
  auto count = size_t(0);
  for (const auto hi : cycles()) {
    const auto he = graph_->hedge(hi);
    auto it = he;
    do {
      count++;
      it = it.next();
    } while (it != he);
  }
  return count;
}

StrokeGraph::StrokeGraph(span<const Stroke> strokes, SnappingType type)
  : strokes_(this)
  , snapping_method_type_(type) {
  orig_strokes_.reserve(strokes.size());
  for (auto& s : strokes) {
    orig_strokes_.emplace_back(s.clone());
    orig_strokes_.back().ensure_arclengths();
  }
  orig_bvh_ = std::make_unique<PolylineBVH>(orig_strokes_);
  strokes = orig_strokes_;

  auto polyline_bvh = PolylineBVH(strokes);
  auto junctions = std::vector<ExtendedJunction>();
  intersection_junctions(polyline_bvh, junctions);
  // 1. Centerline intersecting junctions
  init(strokes, StrokeCoverage(junctions));

  // Intersecting if consider width but not intersecting along centerline

  // Perform snapping for loose endpoints whose envelopes overlap without a centerline
  // intersection.
  const auto n = vertices_.size();
  for (size_t vi = 0; vi < n; ++vi) {
    const auto v = vertex(vi);
    if (v.is_dangling()) {
      auto he = v.hedge();
      auto& he_s = strokes_[he.stroke_idx()];
      auto best_dist = infinity;
      auto stroke_idx = invalid;
      auto stroke_arclen = -infinity;
      for (size_t si = 0; si < strokes_.size(); ++si) {
        const auto& stroke = strokes_[si];
        if (stroke.size() > 1) {
          stroke.ensure_arclengths();
          // FIXME: Switch to find_snap_point once it is well-tested.
          const auto [dist, arclen] =
            find_snap_point_legacy(he_s, he.forward(), bvh_.envelope_bvh_leaf(si));
          if (dist < best_dist) {
            best_dist = dist;
            stroke_idx = si;
            stroke_arclen = arclen;
          }
        }
      }

      if (best_dist < infinity) {
        snap_endpoint_to_edge_legacy(v, stroke_idx, stroke_arclen, snapping_method_type_);
      }
    }
  }

  for (size_t vi = 0; vi < vertices_.size();) {
    if (vertices_[vi].is_active()) {
      vi++; // Move on to next.
    } else {
      delete_vertex(*this, vi);
      // At this point vertex vi has changed, so we should try again.
    }
  }
  for (size_t si = 0; si < strokes_.size();) {
    if (strokes_[si].size() == 0) {
      delete_stroke(*this, si);
    } else {
      si++;
    }
  }

  construct_faces(*this);

  // Renew the junction types. Changes may happen due to deletion.
  {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      auto& v = vertices_[i];
      auto he = HedgeView(this, v.hedge_);
      auto hit = he;
      size_t t_bar = invalid;
      size_t discont_count = 0;
      do {
        if (hit.continuity() != invalid) {
          t_bar = hit.continuity();
        } else {
          discont_count++;
        }
        hit = hit.twin().next();
      } while (hit != he);

      // T-junction
      if (discont_count > 0 && t_bar != invalid) {
        v.junc_type_ = JunctionType::T;
        v.hedge_ = t_bar;
      }
      // End-end
      else if (discont_count > 1 && t_bar == invalid) {
        v.junc_type_ = JunctionType::R;
      }
      // X-junction
      else {
        v.junc_type_ = JunctionType::X;
        auto valence = vertex(i).valence();
        if (t_bar != invalid && valence <= 2) {
          size_t t_bar_pair = hedges_[t_bar].continuity_;
          hedges_[t_bar].continuity_ = invalid;
          hedges_[t_bar_pair].continuity_ = invalid;
        }
      }
    }
  }

#ifndef NDEBUG
  check_consistency(*this);
  check_continuity(*this);
#endif

  // Determine and save junction types; assign unique vertex IDs
  {
    for (size_t i = 0; i < vertices_.size(); ++i) {
      auto& v = vertices_[i];
      auto he = HedgeView(this, v.hedge_);
      auto hit = he;
      size_t t_bar = invalid;
      size_t discont_count = 0;
      do {
        if (hit.continuity() != invalid) {
          t_bar = hit.continuity();
        } else {
          discont_count++;
        }
        hit = hit.twin().next();
      } while (hit != he);

      // T-junction
      if (discont_count > 0 && t_bar != invalid) {
        v.junc_type_ = JunctionType::T;
      }
      // End-end
      else if (discont_count > 1 && t_bar == invalid) {
        v.junc_type_ = JunctionType::R;
      }
      // X-junction
      else {
        v.junc_type_ = JunctionType::X;
      }

      v.ids_.emplace_back(
        StrokeGraph::VertexID({StrokeGraph::VertexID::Initialization, i}));
    }
  }
}

StrokeGraph StrokeGraph::StrokeGraph::clone() const {
  auto dst = StrokeGraph();
  dst.boundary_face_ = boundary_face_;
  dst.bvh_ = bvh_;
  dst.vertices_ = vertices_;
  dst.hedges_ = hedges_;
  dst.faces_ = faces_;
  dst.orig_strokes_.reserve(orig_strokes_.size());
  for (const auto& os : orig_strokes_) {
    dst.orig_strokes_.emplace_back(os.clone());
  }
  if (orig_bvh_) {
    dst.orig_bvh_ = std::make_unique<PolylineBVH>(dst.orig_strokes_);
    dst.orig_bvh_->masked_nodes = orig_bvh_->masked_nodes;
    dst.orig_bvh_->combine_he_indices = orig_bvh_->combine_he_indices;
  }
  dst.endpoint2vertex_ = endpoint2vertex_;
  dst.vertex2endpoint_ = vertex2endpoint_;
  dst.orig2strokes_ = orig2strokes_;
  dst.strokes2orig_ = strokes2orig_;
  dst.snapping_method_type_ = snapping_method_type_;
  dst.strokes2vid_ = strokes2vid_;
  dst.snap_history_ = snap_history_;
  dst.parallel_endpoints_cache_ = parallel_endpoints_cache_;

#ifndef NDEBUG
  check_continuity(*this);
  check_continuity(dst);
#endif

  return dst;
}

StrokeGraph::StrokeGraph(StrokeGraph&& other) noexcept
  : strokes_(this) {
  *this = std::move(other);
}

StrokeGraph& StrokeGraph::operator=(StrokeGraph&& other) noexcept {
  // Move everything except `strokes_` which must still reference `this`.
  boundary_face_ = other.boundary_face_;
  bvh_ = std::move(other.bvh_);
  vertices_ = std::move(other.vertices_);
  hedges_ = std::move(other.hedges_);
  faces_ = std::move(other.faces_);
  orig_strokes_ = std::move(other.orig_strokes_);
  orig_bvh_ = std::move(other.orig_bvh_);
  endpoint2vertex_ = std::move(other.endpoint2vertex_);
  vertex2endpoint_ = std::move(other.vertex2endpoint_);
  orig2strokes_ = std::move(other.orig2strokes_);
  strokes2orig_ = std::move(other.strokes2orig_);
  snapping_method_type_ = other.snapping_method_type_;
  strokes2vid_ = std::move(other.strokes2vid_);
  snap_history_ = std::move(other.snap_history_);
  parallel_endpoints_cache_ = std::move(other.parallel_endpoints_cache_);
  return *this;
}

void StrokeGraph::init(const span<const Stroke> strokes, const StrokeCoverage& coverage) {
  orig2strokes_.resize(strokes.size());
  // TODO: Assess how good these guesses are.
  vertices_.reserve(2 * strokes.size());

  static constexpr auto grazing_vid = (size_t)-2;

  auto cut_points = // Mapping[stroke, (arc length, vertex index)]
    std::unordered_map<const Stroke*, std::vector<std::pair<Float, size_t>>>();
  auto head2vertex = std::unordered_map<const Stroke*, size_t>();
  auto tail2vertex = std::unordered_map<const Stroke*, size_t>();
  auto junc2vertex = std::unordered_map<const ExtendedJunction*, size_t>();

  // Create all junction vertices.
  for (const auto* junc : coverage.junctions()) {
    // auto skip = true;
    // for (const auto& [si, _] : junc->ranges_) {
    //   if (!skip_stroke(strokes[si])) {
    //     skip = false;
    //     break;
    //   }
    // }
    if (junc->is_grazing()) {
      junc2vertex.insert({junc, grazing_vid});
    } else {
      const auto vertex_pos = junc->position(strokes);
      const auto new_vert_idx = vertices_.size();
      vertices_.emplace_back(vertex_pos); // Must set `new_vert.hedge_` later.
      junc2vertex.insert({junc, new_vert_idx});
    }
  }
  for (auto si = 0; si < (int)strokes.size(); ++si) {
    const auto& stroke = strokes[si];
    if (skip_stroke(stroke))
      continue;
    const auto it = coverage.coverage_.find(si);
    auto& cuts = cut_points[&stroke];
    if (it != coverage.coverage_.end()) {
      for (const auto& terr : it->second) {
        const auto vid = junc2vertex[terr.junction_];
        const auto grazing = terr.junction_->is_grazing();
        assert(grazing == (vid == grazing_vid));
        auto pushed = false;
        if (terr.range_.start_ < 1e-5) {
          assert(!grazing);
          head2vertex[&stroke] = vid;
          cuts.emplace_back(terr.range_.mid() * stroke.length(), vid);
          pushed = true;
        }
        if (terr.range_.end_ > 1.0 - 1e-5) {
          assert(!grazing);
          tail2vertex[&stroke] = vid;
          cuts.emplace_back(terr.range_.mid() * stroke.length(), vid);
          pushed = true;
        }
        if (!pushed) {
          cuts.emplace_back(terr.range_.mid() * stroke.length(), vid);
        }
      }
    }
    if (head2vertex.find(&stroke) == head2vertex.end()) {
      const auto new_vertex_idx = vertices_.size();
      vertices_.emplace_back(stroke.xy(0));
      cuts.emplace_back(0.0, new_vertex_idx);
      head2vertex[&stroke] = new_vertex_idx;
    }
    if (tail2vertex.find(&stroke) == tail2vertex.end()) {
      const auto new_vertex_idx = vertices_.size();
      vertices_.emplace_back(stroke.xy(Back));
      cuts.emplace_back(stroke.length(), new_vertex_idx);
      tail2vertex[&stroke] = new_vertex_idx;
    }
    std::sort(cuts.begin(), cuts.end());
  }

  for (auto& [p_stroke, cuts] : cut_points) {
    for (auto& [arclen, vid] : cuts) {
      if (vid == grazing_vid) {
        vid = vertices_.size();
        auto& v = vertices_.emplace_back(p_stroke->pos(arclen));
        v.flags_ |= VertexRecord::Grazing;
      }
    }
  }

  // Cut strokes at junctions.
  auto vertex2outbound = std::unordered_multimap<size_t, size_t>();
  for (size_t si = 0; si < strokes.size(); ++si) {
    const auto& stroke = strokes[si];
    if (skip_stroke(stroke))
      continue;
    const auto cut_it = cut_points.find(&stroke);
    assert(cut_it != cut_points.end());
    const auto& cuts = cut_it->second;
    force_assert(cuts.size() >= 2);
    if (cuts.size() == 2) { // Stroke does not need splitting.
      // Ignore strokes that collapse to a single point.
      if (cuts[1].first - cuts[0].first > 1e-6) {
        const auto new_stroke_idx =
          bvh_.add(stroke.fractional_slice(cuts[0].first, cuts[1].first));
        const auto [new_edge_idx, new_twin_idx] = make_twins(hedges_);
        auto& new_edge = hedges_[new_edge_idx];
        auto& new_twin = hedges_[new_twin_idx];
        const auto it = head2vertex.find(&stroke);
        assert(it != head2vertex.end());
        const auto head_v = it->second;
        const auto it2 = tail2vertex.find(&stroke);
        assert(it2 != tail2vertex.end());
        const auto tail_v = it2->second;
        new_edge.origin_ = head_v;
        new_twin.origin_ = tail_v;
        vertices_[head_v].hedge_ = new_edge_idx;
        vertices_[tail_v].hedge_ = new_twin_idx;
        // Will set next_ and prev_ later (we need to know all the outbound edges first).
        // For now, make a record.
        vertex2outbound.insert({head_v, new_edge_idx});
        vertex2outbound.insert({tail_v, new_twin_idx});

        endpoint2vertex_.insert({Endpoint(si, true), head_v});
        vertex2endpoint_.insert({head_v, Endpoint(si, true)});
        endpoint2vertex_.insert({Endpoint(si, false), tail_v});
        vertex2endpoint_.insert({tail_v, Endpoint(si, false)});
        add_stroke_mapping(*this, si,
                           std::make_pair(cuts[0].first / stroke.length(),
                                          cuts[1].first / stroke.length()),
                           new_stroke_idx, std::make_pair(0., 1.));
      }
    } else { // Stroke needs cutting.
      auto split_values = std::vector<Float>();
      split_values.reserve(cuts.size());
      for (const auto& cut : cuts)
        split_values.push_back(cut.first);
      Float old_stroke_length = stroke.length();
      std::vector<std::pair<Float, Float>> domain_intervals;
      domain_intervals.reserve(split_values.size() - 1);
      for (size_t i = 0; i + 1 < split_values.size(); i++) {
        domain_intervals.emplace_back(split_values[i] / old_stroke_length,
                                      split_values[i + 1] / old_stroke_length);
      }

      auto split_strokes = std::vector<Stroke>();
      split_strokes.reserve(split_values.size() - 1);
      stroke.split(split_values, split_strokes);
      const auto first_idx = strokes_.size(); // First stroke index of splits.
      for (auto& split_stroke : split_strokes) {
        bvh_.add(std::move(split_stroke));
      }
      split_strokes.clear();

      // The outbound hedge from the previous vertex
      size_t prev_twin_idx = invalid;

      for (auto stroke_idx = first_idx; stroke_idx < strokes_.size(); ++stroke_idx) {
        size_t head_v, tail_v;
        if (stroke_idx == first_idx) { // First slice.
          const auto it = head2vertex.find(&stroke);
          force_assert(it != head2vertex.end());
          head_v = it->second;
        } else {
          head_v = cuts[stroke_idx - first_idx].second;
        }
        if (stroke_idx == strokes_.size() - 1) { // Last slice.
          const auto it = tail2vertex.find(&stroke);
          force_assert(it != tail2vertex.end());
          tail_v = it->second;
        } else {
          tail_v = cuts[stroke_idx - first_idx + 1].second;
        }
        const auto [new_edge_idx, new_twin_idx] = make_twins(hedges_);
        auto& new_edge = hedges_[new_edge_idx];
        auto& new_twin = hedges_[new_twin_idx];
        new_edge.origin_ = head_v;
        new_twin.origin_ = tail_v;

        // Save the continuity info at the cut vertex
        if (prev_twin_idx != invalid) {
          new_edge.continuity_ = prev_twin_idx;
          hedges_[prev_twin_idx].continuity_ = new_edge_idx;
        }
        prev_twin_idx = new_twin_idx;

        if (stroke_idx != strokes_.size() - 1) {
          new_edge.next_ = hedges_.size(); // Not created yet.
          new_twin.next_ = new_edge_idx - 1;
        }
        if (stroke_idx != first_idx) {
          new_edge.prev_ = new_edge_idx - 2;
          new_twin.prev_ = hedges_.size() + 1;
        }
        vertices_[head_v].hedge_ = new_edge_idx;
        vertices_[tail_v].hedge_ = new_twin_idx;
        vertex2outbound.insert({head_v, new_edge_idx});
        vertex2outbound.insert({tail_v, new_twin_idx});

        if (stroke_idx == first_idx) {
          endpoint2vertex_.insert({Endpoint(si, true), head_v});
          vertex2endpoint_.insert({head_v, Endpoint(si, true)});
        } else if (stroke_idx == strokes_.size() - 1) {
          endpoint2vertex_.insert({Endpoint(si, false), tail_v});
          vertex2endpoint_.insert({tail_v, Endpoint(si, false)});
        }
        auto domain = domain_intervals[stroke_idx - first_idx];
        add_stroke_mapping(*this, si, std::make_pair(domain.first, domain.second),
                           stroke_idx, std::make_pair(0., 1.));
      }
    }
  }

  // This deforms the strokes.
  snap_endpoints_legacy();

  // Create vertex stars.
  auto outbound_hedges = std::vector<size_t>();
  auto angle_buffer = std::vector<std::pair<double, double>>();
  auto index_buffer = std::vector<size_t>();
  for (size_t vi = 0; vi < vertices_.size(); ++vi) {
    auto [it, end] = vertex2outbound.equal_range(vi);
    for (; it != end; ++it) {
      outbound_hedges.push_back(it->second);
      force_assert(hedge(it->second).origin().index_ == vi);
    }
    angle_buffer.resize(outbound_hedges.size());
    index_buffer.resize(outbound_hedges.size());
    sort_outbound_hedges(*this, outbound_hedges, angle_buffer, index_buffer);

    const auto n = outbound_hedges.size();
    for (size_t i = 0; i < n; ++i) {
      const auto current = outbound_hedges[i];
      const auto next = outbound_hedges[(i + 1) % n]; // Next outbound.
      auto& curr_he = hedges_[current];
      const auto prev_he_idx = hedge(next).twin().index_;
      auto& prev_he = hedges_[prev_he_idx];
      curr_he.prev_ = prev_he_idx;
      prev_he.next_ = current;
    }
    outbound_hedges.clear();
  }

  for (size_t i = 0; i < orig2strokes_.size(); ++i) {
    auto& orig_mappings = orig2strokes_[i];
    if (orig_mappings.empty())
      continue;
    orig_mappings[0].second.domain_arclens_[0] = 0;
    orig_mappings.back().second.domain_arclens_[1] = 1;
  }

#ifndef NDEBUG
  check_continuity(*this);
  check_mapping(*this);
#endif

  // Determine and save junction types
  for (auto& v : vertices_) {
    auto he = HedgeView(this, v.hedge_);
    auto hit = he;
    size_t t_bar = invalid;
    size_t discont_count = 0;
    do {
      if (hit.continuity() != invalid) {
        t_bar = hit.continuity();
      } else {
        discont_count++;
      }
      hit = hit.twin().next();
    } while (hit != he);

    // T-junction
    if (discont_count > 0 && t_bar != invalid) {
      v.junc_type_ = JunctionType::T;
      v.hedge_ = t_bar;
    }
    // End-end
    else if (discont_count > 1 && t_bar == invalid) {
      v.junc_type_ = JunctionType::R;
    }
    // X-junction
    else {
      v.junc_type_ = JunctionType::X;
    }
  }

  strokes2vid_.resize(strokes_.size());
}

std::vector<Vec2> StrokeGraph::face_positions(size_t face_index) const {
  std::vector<Vec2> poly;
  const auto& face = faces_[face_index];
  // TODO: Handle holes.
  const auto hi = face.cycles_[0];
  const auto he = hedge(hi);
  auto out = std::vector<Vec2>();
  cycle_positions(he, out);
  return out;
}

void StrokeGraph::snap_endpoints() {
  for (size_t hi = 0; hi < hedges_.size(); ++hi) {
    const auto he = hedge(hi);
    if (he.is_active() && he.forward()) {
      const auto si = he.stroke_idx();
      Vec2 start = he.origin().pos();
      Vec2 end = he.dest().pos();
      auto& stroke = strokes_[si];
      force_assert(stroke.xy(0).isApprox(start) &&
                   "snap_endpoints cannot do large deformations");
      force_assert(stroke.xy(Back).isApprox(end) &&
                   "snap_endpoints cannot do large deformations");
      stroke.x(0) = start.x_;
      stroke.y(0) = start.y_;
      stroke.x(Back) = end.x_;
      stroke.y(Back) = end.y_;
    }
  }
}

void StrokeGraph::snap_endpoints_legacy() {
  for (size_t hi = 0; hi < hedges_.size(); ++hi) {
    const auto he = hedge(hi);
    if (he.is_active() && he.forward()) {
      const auto si = he.stroke_idx();
      Vec2 start = he.origin().pos();
      Vec2 end = he.dest().pos();
      if (strokes_[si].size() == 1) {
        strokes_[si].resize(2);
        strokes_[si].x(1) = strokes_[si].x(0);
        strokes_[si].y(1) = strokes_[si].y(0);
        strokes_[si].width(1) = strokes_[si].width(0);
        strokes_[si].compute_arclengths();
      }
      smooth_deform(strokes_[si], start, end);
      bvh_.shallow_update(si, StrokeBVH::UpdateHint::AddToHead);
    }
  }
}

VertexView StrokeGraph::snap_endpoint_to_edge_legacy( //
  const VertexView v, const size_t stroke_idx, const Float stroke_arclen,
  const SnappingType type) {

  auto he = v.hedge();
  const auto forward = he.forward();

  VertexView new_vertex;
  bool vertex_added = add_vertex(*this, stroke_idx, stroke_arclen, new_vertex);

  // If the function finds a close vertex to the input vertex and the found one is the
  // same as the input one. We do nothing.
  if (new_vertex.index_ == v.index_)
    return v;

  // Modify the type of the vertex based on the snap position
  JunctionType::Type v_type =
    (new_vertex.is_dangling()) ? JunctionType::R : JunctionType::T;

  // Only move the he if any stroke gets cut in add_vertex.
  if (vertex_added && he.stroke_idx() == stroke_idx && !forward) {
    // We cut our own stroke, and our he_s is now at the end.
    he = hedge(hedges_.size() - 1);
  }

  // he_s could now be a dangling reference, so get a new one.
  MSVC_WARNING_SUPPRESS(4456)
  auto& he_s = strokes_[he.stroke_idx()];

  // Modify the type of the vertex based on the snap position
  vertices_[v.index_].junc_type_ = v_type;

  if (new_vertex.index_ == he.dest().index_ && he_s.length() < he_s.pen_width()) {
    // Drop this edge.
    auto& edge_rec = hedges_[he.index_];
    auto& twin_rec = hedges_[he.twin().index_];
    vertices_[edge_rec.origin_].deactivate();
    if (he.next().stroke_idx() != he.stroke_idx()) {
      vertices_[he.dest().index_].hedge_ = he.next().index_;
    } else if (he.twin().prev().stroke_idx() != he.stroke_idx()) {
      vertices_[he.dest().index_].hedge_ = he.twin().prev().index_;
    } else {
      std::abort(); // Shouldn't happen.
    }

    // Remove the continuity reference
    if (edge_rec.continuity_ != invalid)
      hedges_[edge_rec.continuity_].continuity_ = invalid;
    edge_rec.continuity_ = invalid;
    if (twin_rec.continuity_ != invalid)
      hedges_[twin_rec.continuity_].continuity_ = invalid;
    twin_rec.continuity_ = invalid;

    hedges_[edge_rec.prev_].next_ = edge_rec.next_;
    hedges_[edge_rec.next_].prev_ = edge_rec.prev_;
    hedges_[twin_rec.prev_].next_ = twin_rec.next_;
    hedges_[twin_rec.next_].prev_ = twin_rec.prev_;
    bvh_.clear(he.stroke_idx()); // Mark for later deletion.
    edge_rec.deactivate();
    twin_rec.deactivate();
  } else {
    if (he.forward()) {
      if (he_s.xy(0) != new_vertex.pos()) {
        if (type == Deformation) {
          non_intersecting_smooth_deform(new_vertex, he_s, new_vertex.pos(),
                                         he_s.xy(he_s.size() - 1));
        } else if (type == LinearSolve) {
          non_intersecting_smart_deform(new_vertex, he_s, new_vertex.pos(),
                                        he_s.xy(he_s.size() - 1));
        } else if (type == Connection) {
          he_s.insert(0, new_vertex.pos().x_, new_vertex.pos().y_, he_s.width(0));
        } else {
          std::abort(); // Unknown snap method.
        }
        bvh_.shallow_update(he.stroke_idx(), StrokeBVH::UpdateHint::AddToHead);
      }
    } else {
      if (he_s.xy(Back) != new_vertex.pos()) {
        if (type == Deformation) {
          non_intersecting_smooth_deform(new_vertex, he_s, he_s.xy(0), new_vertex.pos());
        } else if (type == LinearSolve) {
          non_intersecting_smart_deform(new_vertex, he_s, he_s.xy(0), new_vertex.pos());
        } else if (type == Connection) {
          he_s.push_back(new_vertex.pos().x_, new_vertex.pos().y_, he_s.width(Back));
        } else {
          std::abort(); // Unknown snap method.
        }
        bvh_.shallow_update(he.stroke_idx(), StrokeBVH::UpdateHint::AddToTail);
      }
    }
    if (he.origin() != new_vertex) {
      vertices_[he.origin().index_].deactivate();
      vertices_[new_vertex.index_].ids_.insert(vertices_[new_vertex.index_].ids_.end(),
                                               vertices_[he.origin().index_].ids_.begin(),
                                               vertices_[he.origin().index_].ids_.end());
      vertices_[he.origin().index_].ids_.clear();
    }
    hedges_[he.index_].origin_ = new_vertex.index_;
    insert_into_star(*this, he.index_, new_vertex.index_);
  }

  // Update the affected faces if faces are already built.
  if (!faces_.empty()) {
    faces_.clear();
    for (size_t i = 0; i < hedges_.size(); ++i) {
      hedges_[i].face_ = StrokeGraph::invalid;
    }
    construct_faces(*this);
  }

#ifndef NDEBUG
  if (!faces_.empty())
    check_consistency(*this);
  check_continuity(*this);
#endif

  return new_vertex;
}

bool StrokeGraph::snap_vertices(const size_t vi1, const size_t vi2,
                                const SnappingType type, const Float prob,
                                const Vec2* snap_pos) {
  assert(vi1 != vi2);

  const auto v1 = vertex(vi1);
  const auto v2 = vertex(vi2);
  auto new_pos = Vec2::Empty();
  if (!snap_pos) {
    if (v1.is_dangling()) {
      new_pos = v2.pos(); // We assume there is a bias towards v2 if both are dangling.
    } else if (v2.is_dangling()) {
      new_pos = v1.pos();
    } else {
      new_pos = 0.5 * (v1.pos() + v2.pos());
    }
  } else {
    new_pos = *snap_pos;
  }

  // Connect the geometry.
  StrokeGraph graph_tmp = this->clone();
  move_vertex(graph_tmp, vi1, new_pos, type, graph_tmp.vertex(vi2));
  move_vertex(graph_tmp, vi2, new_pos, type, graph_tmp.vertex(vi1));

  if (is_bridge_vertex(graph_tmp, graph_tmp.vertex(vi1)) ||
      is_bridge_vertex(graph_tmp, graph_tmp.vertex(vi2)))
    return false;

  // Check if there's any new intersection once we deformed all strokes within the two
  // 1-rings.
  const std::vector<size_t> vv{vi1, vi2};
  std::vector<Vec2> out;
  for (const auto vi : vv) {
    const auto v = graph_tmp.vertex(vi);
    auto valence = 0;
    const auto he = v.hedge();
    auto it = he;
    do {
      size_t si = it.stroke_idx();
      auto& stroke = graph_tmp.strokes_[si];

      auto node1 = PolylineBVHLeaf(stroke, bounds(stroke));
      node1.geometry->ensure_arclengths();

      // Check for self-intersection
      {
        out.clear();
        intersect_self(node1, out);
        Float check_length = node1.geometry->length();
        for (const auto& out_v : out) {
          if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
            return false;
          }
        }
      }

      for (size_t i = 0; i < v.graph_->strokes_.size(); ++i) {
        if (node1.geometry == &v.graph_->strokes_[i] || v.graph_->strokes_[i].size() == 0)
          continue;
        out.clear();
        const Stroke& s = v.graph_->strokes_[i];
        s.ensure_arclengths();
        Float check_length = s.length();
        const auto node2 = PolylineBVHLeaf(s, bounds(s));
        intersect_different(node1, node2, out);
        for (const auto& out_v : out) {
          if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
            return false;
          }
        }
      }

      it = it.twin().next();
      valence++;
      force_assert(valence < 1024 && "likely infinite loop found");
    } while (it != he);
  }

  // Accept the snap.

  {
    const auto he1 = v1.hedge();
    auto it1 = he1;
    const auto he2 = v2.hedge();
    do {
      auto it2 = he2;
      do {
        // TODO: Make this more detailed, and more accurate.
        auto& pred = graph_tmp.snap_history_.emplace_back();
        pred.key.type = JunctionType::R;
        pred.key.cand1 = (int)vi1;
        pred.key.cand2 = (int)vi2;
        pred.p_a = v1.pos();
        pred.p_b = v2.pos();
        pred.orig_a = Endpoint(it1.stroke_idx(), it1.forward()).as_pair();
        pred.orig_b = Endpoint(it2.stroke_idx(), it2.forward()).as_pair();
        auto ok = convert_strokes2orig(*this, pred.orig_a);
        force_assert(ok && "could not convert from strokes to orig");
        ok = convert_strokes2orig(*this, pred.orig_b);
        force_assert(ok && "could not convert from strokes to orig");
        pred.prob = prob;

        it2 = it2.twin().next();
      } while (it2 != he2);

      it1 = it1.twin().next();
    } while (it1 != he1);
  }

  *this = std::move(graph_tmp);

  // Collect indices of half-edges to snap.
  const auto valence = (int)v1.valence();
  auto hi_arr = ALLOCA_SPAN(size_t, valence);
  {
    const auto he1 = v1.hedge();
    auto it = he1;
    auto e_out = 0;
    do {
      hi_arr[e_out++] = it.index_;
      it = it.twin().next();
    } while (it != he1);
    assert(e_out == valence);
  }

  // Connect the topology.
  if (vi1 != vi2) {
    vertices_[vi1].deactivate();
    if (vertices_[vi1].is_originally_dangling()) {
      vertices_[vi2].flags_ |= VertexRecord::OriginallyDangling;
    }
    vertices_[vi2].ids_.insert(vertices_[vi2].ids_.end(), vertices_[vi1].ids_.begin(),
                               vertices_[vi1].ids_.end());
    vertices_[vi1].ids_.clear();
  }
  for (const auto hi : hi_arr) {
    hedges_[hi].origin_ = v2.index_;
    insert_into_star(*this, hi, v2.index_);
  }

  // Update the affected faces if faces are already built.
  if (!faces_.empty()) {
    faces_.clear();
    for (size_t i = 0; i < hedges_.size(); ++i) {
      hedges_[i].face_ = StrokeGraph::invalid;
    }
    construct_faces(*this);
  }

#ifndef NDEBUG
  if (!faces_.empty())
    check_consistency(*this);
  check_continuity(*this);
#endif

  return true;
}

VertexView StrokeGraph::snap_endpoint_to_edge( //
  const VertexView v, const size_t stroke_idx, Float stroke_arclen,
  const SnappingType type) {

  VertexView new_vertex;
  add_vertex(*this, stroke_idx, stroke_arclen, new_vertex);
  if (v == new_vertex)
    return new_vertex;
  if (!is_snap_valid(*this, v.index_, new_vertex.index_)) {
    return VertexView();
  }

  // Modify the type of the vertex based on the snap position
  // TODO: For high-valence junctions this may be wrong, but I'm not sure what downstream
  //       code expects here.  Revisit if we actually need this information.
  vertices_[v.index_].junc_type_ =
    (new_vertex.is_dangling() ? JunctionType::R : JunctionType::T);

  if (snap_vertices(v.index_, new_vertex.index_, type)) {
    return new_vertex;
  }

  // Snapping failed, probably due to line of sight.
  // Retry with more precise vertex insertion.
  // Careful: add_vertex_precise returns false if either no vertex was inserted due to
  // there already being one there, or if we failed to create two valid edges.
  stroke_arclen = std::clamp(stroke_arclen, 0.0, strokes_[stroke_idx].length());
  if (add_vertex_precise(*this, stroke_idx, stroke_arclen, new_vertex)) {
    force_assert(is_snap_valid(*this, v.index_, new_vertex.index_));
    if (snap_vertices(v.index_, new_vertex.index_, type)) {
      return new_vertex;
    }
  }

  return VertexView(); // Failure.
}

std::string StrokeGraph::repr() const {
  auto ss = std::stringstream();
  ss << "StrokeGraph(\n";
  for (size_t vi = 0; vi < vertices_.size(); ++vi) {
    const auto& v = vertices_[vi];
    if (v.is_active()) {
      ss << "  Vertex " << vi << ": {"
         << " p: (" << v.p_.x_ << ", " << v.p_.y_ << "), "
         << " he: " << v.hedge_ << " }\n";
    }
  }
  for (size_t hi = 0; hi < hedges_.size(); ++hi) {
    const auto& he = hedges_[hi];
    if (he.is_active()) {
      ss << "  Hedge " << hi << ": {"
         << " next: " << he.next_ << ", "
         << " prev: " << he.prev_ << ", "
         << " origin: " << he.origin_ << ", "
         << " face: " << he.face_ << " } \n";
    }
  }
  for (size_t fi = 0; fi < faces_.size(); ++fi) {
    const auto& f = faces_[fi];
    ss << "  Face " << fi << ": {"
       << " cycles: [";
    for (const auto hi : f.cycles_) {
      ss << hi << ", ";
    }
    ss << "] }\n";
  }
  ss << ')';
  return ss.str();
}

bool skip_stroke(const Stroke& s) {
  return s.size() < 2 || s.pen_width() == 0.0 || s.length() < 1.5 * s.pen_width();
}

} // namespace sketching
