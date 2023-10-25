#pragma once

#include <map>
#include <vector>

#include "../../stroke_strip_src/Cluster.h"

struct Input;
typedef std::pair<std::vector<int>, std::vector<int>> Sample;
typedef std::vector<int> TypicalSample;

template <typename T>
bool next_combination(const T first, const T last, int k);
bool check_overlap_two_cluster(const Input &input, Cluster cluster,
                               std::vector<int> cluster_stroke_ind,
                               std::vector<int> sub_cluster_stroke_ind);

/**
 * Generates positive stroke-strip training example candidates from the given
 * input and stores them in the provided vector.
 *
 * @param input The input (with GT labels) to generate samples from.
 * @param pos_stroke_samples The vector to store the generated samples in.
 */
void generate_stroke_positive_samples(const Input &input,
                                      std::vector<Sample> &pos_stroke_samples);
/**
 * Generates negative stroke-strip training example candidates based on the
 * given input.
 *
 * @param input The input (with GT labels) to generate negative training example
 * candidates from.
 * @param neg_stroke_samples The vector to store the generated negative stroke
 * samples in.
 * @param matching_pair_ptr A pointer to a vector of matching pairs to use in
 * the generation process.
 */
void generate_stroke_negative_samples(
  const Input &input, std::vector<Sample> &neg_stroke_samples,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr);
/**
 * Generates positive strip-strip training example candidates for clustering
 * from the given input and stores them in the provided vector.
 *
 * @param input The input data (with GT labels) to generate samples from.
 * @param pos_cluster_samples The vector to store the generated positive
 * training example candidates in.
 */
void generate_cluster_positive_samples(
  const Input &input, std::vector<Sample> &pos_cluster_samples);
/**
 * Generates negative strip-strip training example candidates for clustering by
 * selecting nearby pairs of strips.
 *
 * @param input The input data (with GT labels) to generate samples from.
 * @param neg_cluster_samples The vector to store the generated negative
 * training example candidates.
 * @param matching_pair_ptr A pointer to a vector of matching pairs of sketches.
 *                          This is used to ensure that negative samples are not
 * generated from matching pairs.
 */
void generate_cluster_negative_samples(
  const Input &input, std::vector<Sample> &neg_cluster_samples,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr);
