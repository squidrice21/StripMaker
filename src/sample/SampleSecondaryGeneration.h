#pragma once

#include <map>
#include <vector>

#include "SampleGeneration.h"

/**
 * Generates positive global stroke-strip training example candidates
 * based on the given input and matching pairs.
 *
 * @param input The input data to generate samples from.
 * @param pos_cluster_samples The vector to store the generated positive
 * samples.
 * @param matching_pair_ptr A pointer to the vector of matching pairs to use for
 * generating samples.
 */
void generate_secondary_stroke_positive_samples(
  const Input &input, std::vector<Sample> &pos_cluster_samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr);

/**
 * @brief Generates positive global strip-strip training example candidates
 * based on the input and matching pairs.
 *
 * @param input The input sketch.
 * @param pos_cluster_samples The vector of positive cluster samples to be
 * generated.
 * @param matching_pair_ptr The pointer to the vector of matching pairs.
 */
void generate_secondary_cluster_positive_samples(
  const Input &input, std::vector<Sample> &pos_cluster_samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr);

/**
 * @brief Generates negative global strip-strip training example candidates
 * based on the input and matching pairs.
 *
 * @param input The input sketch.
 * @param neg_cluster_samples The vector of negative cluster samples to be
 * generated.
 * @param matching_pair_ptr The pointer to the vector of matching pairs.
 */
void generate_secondary_cluster_negative_samples(
  const Input &input, std::vector<Sample> &neg_cluster_samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr);
