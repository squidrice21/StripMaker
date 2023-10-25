#pragma once

#include "SketchConnectivity/bvh.h"
#include "SketchConnectivity/classifier.h"
#include "SketchConnectivity/sketching.h"

#include "glm/detail/type_vec.hpp"

#include <vector>

extern double stroke_stroke_distance_ratio;
extern double endpoint_distance_ratio;

// Check if there is a line of sight between two endpoints
bool line_of_sight_junc(
  sketching::Vec2 p, sketching::Vec2 q,
  nonstd::span<const sketching::Stroke> strokes,
  nonstd::span<const sketching::BoundingBox> centerline_bbs);

// Check if a stroke endpoint is roughly an endpoint
bool is_roughly_end(const sketching::Stroke &s1, double t1);

// Propose end-end intersections for a given stroke cluster
void propose_end_end_intersections(
  size_t cur_c_idx, const sketching::StrokeBVH &bvh,
  std::vector<std::pair<sketching::StrokeTime, sketching::StrokeTime>>
    &out_intersections);

// Propose end-end envelope intersections for a given stroke cluster
void propose_end_end_envelope_intersections(
  size_t cur_c_idx, const sketching::StrokeBVH &bvh,
  std::vector<std::pair<sketching::StrokeTime, sketching::StrokeTime>>
    &out_intersections);

// Propose end-end candidates for a given stroke cluster and query point
void propose_end_end_candidates(const glm::dvec2 &v, size_t cur_c_idx,
                                size_t knn, sketching::StrokeBVH &bvh,
                                std::vector<sketching::StrokeTime> &candidates,
                                int to_c_idx = -1);

// Propose end-stroke candidates for a given stroke cluster and query point
void propose_end_stroke_candidates(
  const glm::dvec2 &v, size_t cur_c_idx, size_t knn,
  const sketching::StrokeBVH &bvh,
  std::vector<sketching::StrokeTime> &candidates, int to_c_idx = -1);
