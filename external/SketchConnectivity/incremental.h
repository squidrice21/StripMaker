#pragma once

#include "classifier.h"
#include "stroke_graph_extra.h"

#include <optional>

namespace sketching {

struct ClassifierPrediction;

struct IncrementalCache {
  ////////////
  // Caches //
  ////////////

  // All original strokes including future ones.
  // Must be set if allowed_connections_ is set.
  std::shared_ptr<Stroke[]> orig_strokes_;
  // Bounding boxes of strokes in `orig_strokes_`.
  std::shared_ptr<BoundingBox[]> orig_strokes_bb_;
  // Size of `orig_strokes_` and `orig_strokes_bb_` arrays.
  size_t n_orig_strokes_ = 0;

  // Temporary buffers.
  std::vector<Vec2> intersections_;
  std::vector<ClassifierPrediction> predictions;

  ////////////////
  // Parameters //
  ////////////////

  // Enforce these connections.
  Connections enforced_connections_;

  // Limit connections to these plus those in enforced_connections_.
  Connections allowed_connections_;

  // Immediately accept connections with probability greater than the threshold.
  Float accept_threshold_ = 0.0;

  bool trim_overshoots_ = true;

  ///////////
  // State //
  ///////////

  // Candidates from the region solve.
  std::vector<Junction> candidates_;
};

struct SnapInfo {
  enum Status {
    Dangling,
    Overlap,
    PredictionPos,
    PredictionDelay,
    PredictionNeg,
    ReservedPos,
    ReservedDelay,
    Skip,
  };

  Status status = Dangling;
  /// Only valid if status is PredictionPos, PredictionDelay, or PredictionNeg.
  Float max_prob = -1.0;
  JunctionType::Type prediction_type = JunctionType::X;
  size_t other_vertex = StrokeGraph::invalid;
  StrokeTime stroke_pos1, stroke_pos2;

  [[nodiscard]] bool is_connected() const {
    return status == Overlap || status == PredictionPos || status == ReservedPos;
  }
};

struct StrokeSnapInfo {
  SnapInfo head;
  SnapInfo tail;

  std::vector<ClassifierPrediction> predictions;
};

std::unique_ptr<IncrementalCache> make_incremental_cache(Float accept_threshold = -1.0);

void set_enforced_connections(IncrementalCache* cache, const EnvelopeBVH& strokes);

/**
 * Set enforced and allowed connections based on the input strokes.
 */
void set_future_constraints(IncrementalCache* cache, span<const Stroke> strokes,
                            const size_t num_candidates, bool train_time);

/**
 * @param cache Persistent information needed across iterations.  Keep this value around.
 */
StrokeSnapInfo add_stroke_incremental(StrokeGraph& graph, IncrementalCache* cache,
                                      const Stroke& new_stroke, bool to_dissolve = true);

struct CutPoint {
  Float s; // Arc length value on new_stroke.
  size_t vertex; // Index of pre-existing vertex.
  bool is_self_intersection;
};

StrokeSnapInfo
add_stroke_incremental_topological(StrokeGraph& graph, IncrementalCache* cache,
                                   const Stroke& new_stroke, bool to_dissolve = true);

bool end_end_pre_filter(StrokeGraph::VertexView vertex,
                        StrokeGraph::VertexView other_vertex, bool train_time);

void snap_candidates(const StrokeGraph& graph, size_t vi,
                     const IncrementalCache* const cache,
                     const Connections& intersection_set, FeatureType feature_type,
                     std::vector<SnapInfo>& candidates, size_t cand_count, bool to_expand,
                     bool train_time);
void snap_candidates(const StrokeGraph& graph, size_t vi,
                     const IncrementalCache* const cache, FeatureType feature_type,
                     std::vector<SnapInfo>& candidates, size_t cand_count, bool to_expand,
                     bool train_time);
void snap_corner_candidates(const StrokeGraph& graph, size_t vi,
                            const IncrementalCache* const cache, FeatureType feature_type,
                            std::vector<SnapInfo>& candidates, size_t cand_count = 3);

void all_snap_candidates(const StrokeGraph& graph, size_t orig_si,
                         const IncrementalCache* const cache,
                         const Connections& intersection_set, FeatureType feature_type,
                         std::vector<std::vector<Junction>>& junction_sets,
                         size_t cand_count, bool train_time);

std::vector<ClassifierPrediction> finalize_incremental(StrokeGraph& graph,
                                                       IncrementalCache* cache,
                                                       bool to_dissolve = true);

/** Find the projection of the vertex onto stroke si. */
Float find_projection(StrokeGraph::VertexView vertex, size_t si, Vec2& proj, Float& s);

/**
 * Return the best cand2 value for the endpoint at `v` for a connection with any stroke
 * corresponding to original stroke `orig_si`.
 *
 * `last_orig_si` is the index of the most recently added original stroke.
 *
 * Returns std::nullopt if there are no good candidates.
 */
[[nodiscard]] std::optional<StrokeTime>
reprojected_vertex_stroke_candidate(StrokeGraph::VertexView v, int orig_si,
                                    int last_orig_si, const IncrementalCache& cache);

bool should_dissolve_vertex(const StrokeGraph::VertexView v);

/** Determine if two positions are too far. */
bool endpoint_too_far(const Stroke& s1, const StrokeTime& end1, const Stroke& s2,
                      const StrokeTime& end2, Float* out_cen_dist = nullptr);

} // namespace sketching
