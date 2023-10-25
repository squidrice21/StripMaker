#pragma once

#include "../external/SketchConnectivity/classifier.h"
#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Fitting.h"
#include "SketchConnectivity/bvh.h"
#include "SketchConnectivity/sketching.h"
#include "glm/detail/type_vec.hpp"

#include <vector>

// Interface for the junction prediction code in external/SketchConnectivity.
// (https://github.com/enjmiah/SketchConnectivity/tree/main).

/**
 * @brief A struct representing a junction between two strokes.
 *
 */
struct Junction {
  // Stroke index and normalized (i.e. in [0, 1]) arc length pairs.
  // Non-endpoints come first, followed by endpoints.  Within these regions,
  // they are sorted using `operator<` of `StrokeTime`.
  sketching::StrokeTime from;
  sketching::StrokeTime to;

  glm::dvec2 from_pos;
  glm::dvec2 to_pos;

  sketching::JunctionType::Type type;

  // The probability that this junction is connected. Since we only save
  // positive junctions for now. This value would be > 0.5 (unless we change the
  // threshold).
  double probability = -1.0;
};

/**
 * @brief Convert a fitted curve to a sketching stroke.
 *
 * @param fit The fitted curve to convert.
 * @param junc_stroke The converted sketching stroke.
 */
void convert_stroke(const FittedCurve &fit, sketching::Stroke &junc_stroke);

/**
 * @brief Construct a bounding volume hierarchy for a set of fitted curves.
 *
 * @param fits The set of fitted curves to construct the BVH for.
 * @param bvh The constructed BVH.
 */
void construct_bvh(const std::vector<FittedCurve> &fits,
                   sketching::StrokeBVH &bvh);

/**
 * @brief Prepare a set of fitted curves for junction prediction.
 *
 * @param clusters The set of stroke clusters to prepare the fits for.
 * @param fits The prepared set of fitted curves.
 */
void prepare_fits(const std::map<size_t, Cluster> &clusters,
                  std::vector<FittedCurve> &fits);

/**
 * @brief Predict a candidate junction between a stroke and a fitted curve.
 *
 * @param stroke The stroke to predict the junction for.
 * @param is_head Whether the junction is at the head of the stroke.
 * @param fit The fitted curve to predict the junction for.
 * @param bvh The BVH for the set of fitted curves.
 * @param cand The position of the other end of the predicted junction.
 * @param junc The predicted junction.
 */
void predict_candidate(const sketching::Stroke &stroke, bool is_head,
                       const FittedCurve &fit, sketching::StrokeBVH &bvh,
                       sketching::StrokeTime cand, Junction &junc);

/**
 * @brief Predict all junctions between a fitted curve and a set of strokes.
 *
 * @param fit The fitted curve to predict the junctions for.
 * @param is_head Whether the junction is at the head of the stroke.
 * @param knn The number of nearest neighbors to consider.
 * @param bvh The BVH for the set of fitted curves.
 * @param junctions The predicted junctions.
 */
void predict_junctions(const FittedCurve &fit, bool is_head, size_t knn,
                       sketching::StrokeBVH &bvh,
                       std::vector<Junction> &junctions);
