#pragma once

#include "bounding_box.h"
#include "sketching.h"

#include <vector>

namespace sketching {

struct PolylineBVHLeaf {
  const Stroke* geometry = nullptr;
  BoundingBox bb;

  PolylineBVHLeaf(const Stroke& g, const BoundingBox& b)
    : geometry(&g)
    , bb(b) {}
};

struct EnvelopeBVHLeaf {
  const Stroke* geometry = nullptr;
  BoundingBox bb;

  explicit EnvelopeBVHLeaf(const Stroke& g)
    : geometry(&g)
    , bb(visual_bounds(g)) {}

  EnvelopeBVHLeaf(const Stroke& g, BoundingBox new_bb)
    : geometry(&g)
    , bb(new_bb) {}
};

/**
 * Flat bounding volume "hierarchy" representing the centerlines of strokes.
 */
struct PolylineBVH {
  std::vector<PolylineBVHLeaf> nodes;
  BoundingBox bb;
  std::vector<size_t> masked_nodes;
  std::vector<size_t> combine_he_indices;

  PolylineBVH() = default;

  /**
   * The lifetime of the BVH must not exceed the lifetime of the strokes.
   */
  explicit PolylineBVH(const span<const Stroke> strokes) { init(strokes); }

  void init(span<const Stroke> strokes) {
    constexpr auto infinity_ = std::numeric_limits<Float>::infinity();
    for (const auto& stroke : strokes) {
      auto xmin = infinity_, ymin = infinity_, xmax = -infinity_, ymax = -infinity_;
      for (auto i = 0; i < stroke.size(); ++i) {
        const auto p = stroke.xy(i);
        xmin = std::min(xmin, p.x());
        ymin = std::min(ymin, p.y());
        xmax = std::max(xmax, p.x());
        ymax = std::max(ymax, p.y());
      }
      nodes.emplace_back(stroke, BoundingBox(xmin, xmax, ymin, ymax));
      bb.unite(nodes.back().bb);
    }
  }
};

/**
 * Flat bounding volume "hierarchy" representing the envelopes of strokes, which take the
 * thickness of strokes into account.
 */
struct EnvelopeBVH {
  std::vector<EnvelopeBVHLeaf> nodes;
  BoundingBox bb;

  /**
   * The lifetime of the BVH must not exceed the lifetime of the strokes.
   */
  explicit EnvelopeBVH(const span<const Stroke> strokes) {
    for (const auto& stroke : strokes) {
      nodes.emplace_back(stroke);
      bb.unite(nodes.back().bb);
    }
  }
};

/**
 * Bounding volume hierarchy with two levels.
 *
 *  - In the first (highest) level, we store an axis-aligned bounding box per stroke.
 *  - In the second level, we store a collection of oriented bounding boxes for each
 *    disjoint section of the stroke.
 *  - In the last level we store the raw geometry.
 */
struct Depth2BVH {
  Depth2BVH() = default;

  span<const Stroke> strokes() const { return strokes_; }

  /**
   * After modifying a stroke, you must call either shallow_update or full_update in order
   * to update the associated bounding boxes.
   */
  Stroke& stroke(Index si) { return strokes_[si]; }
  const Stroke& stroke(Index si) const { return strokes_[si]; }

  [[nodiscard]] Stroke coarse_stroke(Index si) const;

  /**
   * Add a stroke to the BVH.
   *
   * Indices assigned to each stroke will be the same as the insertion order, i.e. the
   * first stroke inserted will receive index 0, the second will receive index 1, etc.
   *
   * Returns the index associated with the newly inserted stroke.
   */
  Index add(Stroke stroke);

  enum class UpdateHint {
    AddToHead,
    AddToTail,
  };

  /**
   * Perform a "shallow" update of the bounding boxes associated with stroke si.
   *
   * A shallow update will preserve the number of L2 nodes associated with this stroke,
   * which makes it suitable for small deformations of the stroke, where corners are more
   * or less preserved.
   *
   * If the number of stroke vertices has changed, make sure to specify the hint depending
   * on whether new vertices were added to the front (`UpdateHint::AddToHead`) or back
   * (`UpdateHint::AddToTail`).  If the number of vertices hasn't changed, then the hint
   * is ignored.
   */
  void shallow_update(Index si, UpdateHint hint);

  /**
   * Update of the bounding boxes associated with stroke si, possibly creating or removing
   * L2 nodes as necessary.  This function should be called if the stroke is trimmed,
   * extended, or otherwise has its geometry substantially changed.
   *
   * This function usually increases the fragmentation of the L2 nodes array.
   */
  void full_update(Index si);

  /** Defragment the L2 nodes array so that no space is wasted storing deleted nodes. */
  // void defragment(); // TODO: Implement this.

  /**
   * Clear stroke si.
   */
  void clear(Index si);

  /**
   * Remove all strokes from this BVH.
   */
  void clear_all();

  /** Swap strokes si and sj and their corresponding bounding boxes. */
  void swap_strokes(Index si, Index sj);

  /**
   * Remove the stroke associated with the highest index from this BVH.
   *
   * Note that you may delete any stroke by swapping it with the last one, then calling
   * this function.
   */
  void pop_back();

  void swap(Depth2BVH& other);

  span<const BoundingBox> centerline_bbs() const { return centerline_bbs_; }

  span<const BoundingBox> envelope_bbs() const { return envelope_bbs_; }

  /** Compatibility function. */
  EnvelopeBVHLeaf envelope_bvh_leaf(Index si) const {
    return {strokes_[si], envelope_bbs_[si]};
  }

  /** Compatibility function. */
  PolylineBVHLeaf polyline_bvh_leaf(Index si) const {
    return {strokes_[si], centerline_bbs_[si]};
  }

  /** Compatibility function. */
  PolylineBVH polyline_bvh() const;

  void check_consistency() const;

  // Copy constructors.
  Depth2BVH(const Depth2BVH&);
  Depth2BVH& operator=(const Depth2BVH&);

  // Move constructors.
  Depth2BVH(Depth2BVH&& other) noexcept = default;
  Depth2BVH& operator=(Depth2BVH&& other) noexcept = default;

  struct L1Node;
  struct L2Node;

private:
  std::vector<Stroke> strokes_;

  // Level 1.

  /// Axis-aligned bounding boxes, one for each stroke centerline.
  std::vector<BoundingBox> centerline_bbs_;
  /// Axis-aligned bounding boxes, one for each stroke envelope.
  std::vector<BoundingBox> envelope_bbs_;
  /// Level 1 nodes, one for each stroke in `strokes_`.
  std::vector<L1Node> l1nodes_;

  // Level 2.

  /// Level 2 nodes.
  std::vector<L2Node> l2nodes_;

public:
  struct L1Node {
    /// Index into l2nodes_ array representing first L2Node for this stroke.
    int start_index_ = 0;
    /// Number of L2Nodes associated with this stroke.
    int size_ = 0;
    /// We have space for this many L2Nodes in this section of the l2nodes_ array.
    int capacity_ = 0;
  };

  struct L2Node {
    /// Bounding box of the envelope of the stroke segment covered by this node.
    OrientedBoundingBox envelope_orbb_;
    /// Bounding box of the centerline of the stroke segment covered by this node.
    OrientedBoundingBox centerline_orbb_;
    /// Index of first vertex of stroke segment covered by this node.
    int start_index_ = 0;
    /// Index of last vertex of stroke segment covered by this node (inclusive).
    int end_index_ = 0;
  };
};

using StrokeBVH = Depth2BVH;

} // namespace sketching
