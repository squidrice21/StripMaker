#include "bvh.h"

#include "detail/alloca.h"
#include "detail/util.h"
#include "force_assert.h"
#include "resample.h"

namespace sketching {

/**
 * Construct oriented bounding boxes for the stroke's envelope and centerline from vertex
 * `l2.start_index_` to vertex `l2.end_index_`.
 */
static void construct_bbs(Depth2BVH::L2Node& l2, const Stroke& stroke) {
  const auto start = l2.start_index_;
  const auto end = l2.end_index_;
  assert(start < end);
  auto axis = Vec2();
  auto axis_norm = Float(0.0);
  for (Index i = start + 1; i <= end; ++i) {
    const auto try_axis = stroke.xy(i) - stroke.xy(start);
    const auto try_axis_norm = try_axis.norm();
    if (try_axis_norm > axis_norm) {
      axis = try_axis;
      axis_norm = try_axis_norm;
    }
  }
  if (axis_norm < 1e-7) {
    // Degenerate segment.
    l2.centerline_orbb_.axis_ = Vec2(1, 0); // Arbitrary.
    l2.centerline_orbb_.extents_ = Vec2(1e-7, 1e-7);
    l2.centerline_orbb_.center_ = stroke.xy(start);
    l2.envelope_orbb_.axis_ = Vec2(1, 0);
    auto width = Float(0);
    for (int j = start; j <= end; ++j)
      width = std::max(width, stroke.width(j));
    l2.envelope_orbb_.extents_ = Vec2(0.5 * width, 0.5 * width);
    l2.envelope_orbb_.center_ = stroke.xy(start);
    return;
  }
  axis /= axis_norm;
  const auto perp = Vec2(-axis.y_, axis.x_);
  auto centerline_axis_min = infinity;
  auto centerline_axis_max = -infinity;
  auto centerline_perp_min = infinity;
  auto centerline_perp_max = -infinity;
  auto envelope_axis_min = infinity;
  auto envelope_axis_max = -infinity;
  auto envelope_perp_min = infinity;
  auto envelope_perp_max = -infinity;

  for (int j = start; j <= end; ++j) {
    const auto r = 0.5 * stroke.width(j);
    const auto p = stroke.xy(j) - stroke.xy(start);
    // Re-express p in terms of local axes of the oriented bounding box.
    const auto u = p.dot(axis);
    centerline_axis_min = std::min(centerline_axis_min, u);
    centerline_axis_max = std::max(centerline_axis_max, u);
    envelope_axis_min = std::min(envelope_axis_min, u - r);
    envelope_axis_max = std::max(envelope_axis_max, u + r);
    const auto v = p.dot(perp);
    centerline_perp_min = std::min(centerline_perp_min, v);
    centerline_perp_max = std::max(centerline_perp_max, v);
    envelope_perp_min = std::min(envelope_perp_min, v - r);
    envelope_perp_max = std::max(envelope_perp_max, v + r);
  }

  // Adjust for round-off error (try to get the bounding box to actually bound).
  // (No guarantee.)
  centerline_axis_max += 1e-14 * std::abs(centerline_axis_max);
  centerline_axis_min -= 1e-14 * std::abs(centerline_axis_min);
  envelope_axis_max += 1e-14 * std::abs(envelope_axis_max);
  envelope_axis_min -= 1e-14 * std::abs(envelope_axis_min);

  l2.centerline_orbb_.axis_ = axis;
  l2.centerline_orbb_.extents_ = 0.5 * Vec2(centerline_axis_max - centerline_axis_min,
                                            centerline_perp_max - centerline_perp_min);
  l2.centerline_orbb_.center_ =
    (stroke.xy(start) + //
     0.5 * (centerline_axis_min + centerline_axis_max) * axis +
     0.5 * (centerline_perp_min + centerline_perp_max) * perp);
  assert(!l2.centerline_orbb_.hasNaN());
  l2.envelope_orbb_.axis_ = axis;
  l2.envelope_orbb_.extents_ = 0.5 * Vec2(envelope_axis_max - envelope_axis_min,
                                          envelope_perp_max - envelope_perp_min);
  l2.envelope_orbb_.center_ = //
    (stroke.xy(start) + //
     0.5 * (envelope_axis_min + envelope_axis_max) * axis +
     0.5 * (envelope_perp_min + envelope_perp_max) * perp);
  assert(!l2.envelope_orbb_.hasNaN());
}

Index Depth2BVH::add(Stroke in_stroke) {
  const auto new_index = strokes_.size();
  strokes_.emplace_back(std::move(in_stroke));
  full_update(new_index);
  return new_index;
}

void Depth2BVH::shallow_update(const Index si, const UpdateHint hint) {
  auto& l1 = l1nodes_[si];
  const auto& stroke = strokes_[si];
  centerline_bbs_[si] = bounds(stroke);
  envelope_bbs_[si] = visual_bounds(stroke);

  // Update ranges if necessary.
  auto& last_s2 = l2nodes_[l1.start_index_ + l1.size_ - 1];
  if (last_s2.end_index_ != stroke.size()) {
    // Can't do a shallow update when trimming (yet).
    force_assert(last_s2.end_index_ < stroke.size() && "not implemented yet");
    if (hint == UpdateHint::AddToTail) {
      last_s2.end_index_ = int(stroke.size() - 1);
    } else if (hint == UpdateHint::AddToHead) {
      const auto shift = int(stroke.size() - last_s2.end_index_ - 1);
      for (int i = l1.start_index_; i < l1.start_index_ + l1.size_; ++i) {
        l2nodes_[i].start_index_ += shift;
        l2nodes_[i].end_index_ += shift;
      }
      assert(l2nodes_[l1.start_index_].end_index_ != 0);
      l2nodes_[l1.start_index_].start_index_ = 0;
    } else {
      std::abort(); // Unknown hint.
    }
  }

  for (int i = l1.start_index_; i < l1.start_index_ + l1.size_; ++i) {
    construct_bbs(l2nodes_[i], stroke);

#ifdef BVH_PARANOID
    for (int j = l2nodes_[i].start_index_; j <= l2nodes_[i].end_index_; ++j) {
      const auto p = stroke.xy(j);
      const auto dist = l2nodes_[i].centerline_orbb_.squared_distance_to(p);
      force_assert(dist <= 1e-20);
    }
#endif
  }
}

void Depth2BVH::full_update(const Index si) {
  const auto& stroke = strokes_[si];

  assert(centerline_bbs_.size() == l1nodes_.size());
  assert(si <= (Index)l1nodes_.size());
  if (si < (Index)l1nodes_.size()) {
    centerline_bbs_[si] = bounds(stroke);
    envelope_bbs_[si] = visual_bounds(stroke);
  } else {
    // This case is needed for `add`.
    centerline_bbs_.emplace_back(bounds(stroke));
    envelope_bbs_.emplace_back(visual_bounds(stroke));
    l1nodes_.emplace_back();
  }

  auto& l1 = l1nodes_[si];
  if (stroke.size() == 0) {
    // No level 2 nodes associated with this stroke.
    l1.size_ = 0;
  } else if (stroke.size() == 1) {
    l1.size_ = 1;
    if (l1.size_ > l1.capacity_) {
      l1.start_index_ = (int)l2nodes_.size();
      l2nodes_.resize(l2nodes_.size() + l1.size_);
      l1.capacity_ = l1.size_;
    }
    auto& l2 = l2nodes_[l1.start_index_];
    l2.start_index_ = 0;
    l2.end_index_ = 0;
    l2.centerline_orbb_.axis_ = Vec2(1, 0); // Arbitrary.
    l2.centerline_orbb_.extents_ = Vec2(0, 0);
    l2.centerline_orbb_.center_ = stroke.xy(0);
    l2.envelope_orbb_.axis_ = Vec2(1, 0);
    l2.envelope_orbb_.extents_ = Vec2(0.5 * stroke.width(0), 0.5 * stroke.width(0));
    l2.envelope_orbb_.center_ = stroke.xy(0);
  } else {
    const auto tolerance = 2 * stroke.pen_width(); // TODO: Tune me for best performance.
    auto include = ALLOCA_SPAN(bool, stroke.size());
    ramer_douglas_peucker(stroke, tolerance, &include[0]);

    l1.size_ = 0;
    for (size_t i = 1; i < include.size(); ++i) {
      if (include[i]) {
        l1.size_++;
      }
    }
    if (l1.size_ > l1.capacity_) {
      l1.start_index_ = (int)l2nodes_.size();
      l2nodes_.resize(l2nodes_.size() + l1.size_);
      l1.capacity_ = l1.size_;
    }

    auto start_i = size_t(0);
    auto l2i = l1.start_index_;
    for (size_t i = 1; i < include.size(); ++i) {
      if (include[i]) {
        auto& l2 = l2nodes_[l2i++];
        l2.start_index_ = (int)start_i;
        l2.end_index_ = (int)i;
        construct_bbs(l2, stroke);
        start_i = i;
      }
    }
  }
}

void Depth2BVH::clear(const Index si) {
  strokes_[si].clear();
  centerline_bbs_[si] = BoundingBox();
  envelope_bbs_[si] = BoundingBox();
  l1nodes_[si].size_ = 0;
}

void Depth2BVH::clear_all() {
  strokes_.clear();
  centerline_bbs_.clear();
  envelope_bbs_.clear();
  l1nodes_.clear();
  l2nodes_.clear();
}

void Depth2BVH::swap_strokes(const Index si, const Index sj) {
  std::swap(strokes_[si], strokes_[sj]);
  std::swap(centerline_bbs_[si], centerline_bbs_[sj]);
  std::swap(envelope_bbs_[si], envelope_bbs_[sj]);
  std::swap(l1nodes_[si], l1nodes_[sj]);
}

void Depth2BVH::pop_back() {
  strokes_.pop_back();
  centerline_bbs_.pop_back();
  envelope_bbs_.pop_back();
  l1nodes_.pop_back();
}

void Depth2BVH::swap(Depth2BVH& other) {
  centerline_bbs_.swap(other.centerline_bbs_);
  envelope_bbs_.swap(other.envelope_bbs_);
  l1nodes_.swap(other.l1nodes_);
  l2nodes_.swap(other.l2nodes_);
  strokes_.swap(other.strokes_);
}

PolylineBVH Depth2BVH::polyline_bvh() const {
  auto pbvh = PolylineBVH();
  auto pbvh_bb = BoundingBox();
  pbvh.nodes.reserve(strokes_.size());
  for (size_t i = 0; i < strokes_.size(); ++i) {
    pbvh.nodes.emplace_back(strokes_[i], centerline_bbs_[i]);
    pbvh_bb.unite(centerline_bbs_[i]);
  }
  pbvh.bb = pbvh_bb;

  return pbvh;
}

void Depth2BVH::check_consistency() const {
  force_assert(strokes_.size() == centerline_bbs_.size());
  force_assert(strokes_.size() == envelope_bbs_.size());
  force_assert(strokes_.size() == l1nodes_.size());
  const auto n_strokes = strokes_.size();
  for (size_t si = 0; si < n_strokes; ++si) {
    const auto& l1 = l1nodes_[si];
    const auto n_vertices = strokes_[si].size();
    if (n_vertices == 0) {
      force_assert(l1.size_ == 0);
    } else {
      force_assert(l1.size_ > 0);
      force_assert(l1.capacity_ >= l1.size_);
      force_assert(l1.start_index_ >= 0);
      const auto bb = bounds(strokes_[si]);
      force_assert(bb == centerline_bbs_[si]);
      const auto visual_bb = visual_bounds(strokes_[si]);
      force_assert(visual_bb == envelope_bbs_[si]);

      const auto end = l1.start_index_ + l1.size_;
      force_assert((size_t)end <= l2nodes_.size());
      auto boundary = 0;
      for (int i = l1.start_index_; i < end; ++i) {
        const auto& l2 = l2nodes_[i];
        force_assert(!l2.centerline_orbb_.hasNaN());
        force_assert(!l2.envelope_orbb_.hasNaN());
        force_assert(l2.start_index_ == boundary);
        force_assert(l2.end_index_ < n_vertices);
        boundary = l2.end_index_;
      }
    }
  }
}

Depth2BVH::Depth2BVH(const Depth2BVH& other)
  : centerline_bbs_(other.centerline_bbs_)
  , envelope_bbs_(other.envelope_bbs_)
  , l1nodes_(other.l1nodes_)
  , l2nodes_(other.l2nodes_) { // TODO: We can defragment l2nodes here.

  strokes_.reserve(other.strokes_.size());
  for (const auto& s : other.strokes_) {
    strokes_.emplace_back(s.clone());
  }
}

Depth2BVH& Depth2BVH::operator=(const Depth2BVH& other) {
  Depth2BVH tmp(other);
  swap(tmp);
  return *this;
}

} // namespace sketching
