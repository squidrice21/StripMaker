#pragma once

#include "endpoint.h"
#include "types.h"

#include <vector>

namespace sketching {

constexpr auto env_distance_over_pen_width_hard_threshold = 10;

using StrokeTime = std::pair<int, double>;

struct JunctionFeature;
struct PolylineBVH;
struct StrokeGraph;

////////////////////
// High-level API //
////////////////////

struct JunctionType {
  enum Type : std::uint8_t {
    R = 2, // Regular junction
    T = 1, // T-junction, where the _ of the T is a single stroke
    X = 0, // intersection
  };
};

/// A single feature vector for the end-end classifier.
struct EndEndFeatures {
  static constexpr size_t n_features_ = 14;
  Float data_[n_features_];

  explicit operator span<Float>() { return {data_, n_features_}; }
  explicit operator span<const Float>() const { return {data_, n_features_}; }
};

/// A single feature vector for the end-stroke classifier.
struct EndStrokeFeatures {
  static constexpr size_t n_features_ = 8;
  Float data_[n_features_];

  explicit operator span<Float>() { return {data_, n_features_}; }
  explicit operator span<const Float>() const { return {data_, n_features_}; }
};

// Polymorphic type (tagged union).
struct FeatureVector {
  enum Type {
    Uninitialized = 0,
    EndEnd = JunctionType::R,
    EndStroke = JunctionType::T,
  } type_;
  union {
    EndEndFeatures ee_fea_;
    EndStrokeFeatures es_fea_;
  };

  explicit operator span<const Float>() const {
    if (type_ == EndEnd) {
      return span<const Float>(ee_fea_);
    } else if (type_ == EndStroke) {
      return span<const Float>(es_fea_);
    } else {
      std::abort();
    }
  }

  explicit operator span<Float>() {
    auto data = span<const Float>(*this);
    return {(Float *)data.data(), data.size()};
  }
};

struct ClassifierPrediction {
  /// The tuple (type, cand1, cand2) is used to identify a connection.
  struct {
    JunctionType::Type type = JunctionType::X;
    JunctionType::Type corner_type = JunctionType::X;
    /// Index of a vertex.
    int cand1 = 0;
    /// Either index of a vertex (for end-end) or a stroke index (for
    /// end-stroke).
    int cand2 = 0;
  } key;
  Vec2 p_a; ///< Position of cand1 at prediction time.
  Vec2 p_b; ///< Position of cand2 at prediction time.
  /// Location of cand1 on the original drawing as (index, normalized arc
  /// length).
  StrokeTime orig_a{0, 0};
  /// Location of cand2 on the original drawing as (index, normalized arc
  /// length).
  StrokeTime orig_b{0, 0};
  Float prob = 0; ///< Prediction probability.
  Float alt_prob = -1; ///< Alternative prediction probability.
  FeatureVector fea;
  /// If this candidate is already connected in the corresponding graph.
  bool connected = false;
  /// `repr` field of the Junction this was created from (if any).
  std::string junction_repr;

  constexpr bool operator<(const ClassifierPrediction &other) const {
    if (key.type < other.key.type) {
      return true;
    } else if (key.type > other.key.type) {
      return false;
    }
    if (key.cand1 < other.key.cand1) {
      return true;
    } else if (key.cand1 > other.key.cand1) {
      return false;
    }
    return key.cand2 < other.key.cand2;
  }
};

/// Choose feature computation method
enum FeatureType : std::uint32_t { //
  Graph = 0,
  OrigStroke = 1,
};

/**
 * Batch-compute the probabilities of several potential connections.
 *
 * This function supports high-valence junctions, and works for both
 * vertex-vertex and vertex-stroke connections.
 *
 * @param graph The drawing to operate on.  cand1 and cand2 are expected to
 * index into the graph's strokes_ array (rather than orig_strokes_).
 * @param bvh BVH for the stroke graph.
 * @param cand1 Connecting vertex, expressed as an endpoint.  High-valence
 * endpoints will be handled for you, meaning that a single connection
 * probability may be the aggregation of multiple classifier evaluations.
 * @param cand2 Connected-to point, as (stroke index, normalized arc length).
 * For connecting to a vertex, the normalized arc length value must be exactly 0
 *              or 1.  High-valence endpoints will be handled for you.
 * @param out_prob Array to output probabilities to.
 * @param out_predictions Optional extra information about each classifier
 * evaluation.
 * @param expand Expand the candidates.
 */
void compute_junction_probabilities(
  const StrokeGraph &graph, span<const Endpoint> cand1,
  span<const StrokeTime> cand2, FeatureType feature_type, span<Float> out_prob,
  std::vector<ClassifierPrediction> *out_predictions = nullptr,
  bool expand = true,
  const std::vector<JunctionType::Type> &junction_types =
    std::vector<JunctionType::Type>());

///////////////////
// Low-level API //
///////////////////

Float clf_endpoint(const EndEndFeatures &feature_vec);

Float clf_tjunction(const EndStrokeFeatures &feature_vec);

std::vector<std::unique_ptr<JunctionFeature>> get_end_end_features();

std::vector<std::unique_ptr<JunctionFeature>> get_end_stroke_features();

std::vector<std::string> get_end_end_feature_descriptions();

std::vector<std::string> get_end_stroke_feature_descriptions();

void human_readable_end_end_features(span<Float> feature_vec);

void human_readable_end_stroke_features(span<Float> feature_vec);

/**
 * Return the number of classifier evaluations needed to assess the probability
 * of a potentially high-valence connection between two vertices.
 */
Index n_classifier_evaluations_vv(const StrokeGraph &graph, const Endpoint e1,
                                  const Endpoint e2);
/**
 * Return the number of classifier evaluations needed to assess the probability
 * of a potentially high-valence connection between a vertex and the interior of
 * an edge.
 */
Index n_classifier_evaluations_vs(const StrokeGraph &graph, const Endpoint e1);
void expand_candidates_vv(const StrokeGraph &graph, const Endpoint e1,
                          const Endpoint e2, const span<Endpoint> out_cand1,
                          const span<StrokeTime> out_cand2);
void expand_candidates_vs(const StrokeGraph &graph, const Endpoint e1,
                          const span<Endpoint> out_cand);

/** Not for high-valence junctions. */
// Run on graph
void compute_features(const StrokeGraph &graph, span<const Endpoint> cand1,
                      span<const StrokeTime> cand2,
                      span<FeatureVector> out_features,
                      const std::vector<JunctionType::Type> &junction_types =
                        std::vector<JunctionType::Type>());
// Run on stroke
void compute_features(const StrokeGraph &graph, const PolylineBVH &bvh,
                      const span<const Endpoint> cand1,
                      const span<const StrokeTime> cand2,
                      const span<FeatureVector> out_features);

} // namespace sketching
