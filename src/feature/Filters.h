/**
 * @brief Filters for sample generation and inference.
 *
 */

#pragma once

#include <set>
#include <string>
#include <vector>

struct FeatureVector;
struct Cluster;

/**
 * @brief Filter to check if the two strips are too far to be considered.
 *
 * @param cluster1 First cluster to compare.
 * @param selected_strokes1 Set of selected strokes for the first cluster.
 * @param cluster2 Second cluster to compare.
 * @param selected_strokes2 Set of selected strokes for the second cluster.
 * @return true if the two clusters pass the filter, false otherwise.
 */
bool filter_sample(const Cluster &cluster1,
                   const std::set<int> &selected_strokes1,
                   const Cluster &cluster2,
                   const std::set<int> &selected_strokes2);

/**
 * @brief Filter to check if the two strips are too far to be considered.
 *
 * @param cluster1 First cluster to compare.
 * @param selected_strokes1 Vector of selected strokes for the first cluster.
 * @param cluster2 Second cluster to compare.
 * @param selected_strokes2 Vector of selected strokes for the second cluster.
 * @return true if the two clusters pass the filter, false otherwise.
 */
bool filter_sample(const Cluster &cluster1,
                   const std::vector<int> &selected_strokes1,
                   const Cluster &cluster2,
                   const std::vector<int> &selected_strokes2);

/**
 * @brief Filter to check if the two strips have enough overlapping to be
 * considered for sample generation.
 *
 * @param cluster1 First cluster to compare.
 * @param selected_strokes1 Vector of selected strokes for the first cluster.
 * @param cluster2 Second cluster to compare.
 * @param selected_strokes2 Vector of selected strokes for the second cluster.
 * @param matching_pair_ptr Pointer to a vector of matching stroke pairs.
 * @param cache_folder Path to the cache folder.
 * @param vis_folder Path to the visualization folder.
 * @return true if the two clusters pass the filter, false otherwise.
 */
bool filter_sample_overlapping(
  const Cluster &cluster1, const std::vector<int> &selected_strokes1,
  const Cluster &cluster2, const std::vector<int> &selected_strokes2,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  const std::string &cache_folder, const std::string &vis_folder);

/**
 * @brief Filters to check if the two strips are parallel enough to be
 * considered.
 *
 * @param cluster1 First cluster to compare.
 * @param selected_strokes1 Vector of selected strokes for the first cluster.
 * @param cluster2 Second cluster to compare.
 * @param selected_strokes2 Vector of selected strokes for the second cluster.
 * @param matching_pair_ptr Pointer to a vector of matching stroke pairs.
 * @param cache_folder Path to the cache folder.
 * @return true if the two clusters pass the filter, false otherwise.
 */
bool filter_sample_angle(
  const Cluster &cluster1, const std::vector<int> &selected_strokes1,
  const Cluster &cluster2, const std::vector<int> &selected_strokes2,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  const std::string &cache_folder);

/**
 * @brief Filter to check if the two strips are overlapping at all.
 *
 * @param fea Feature vector to filter.
 * @return true if the feature vector passes the filter, false otherwise.
 */
bool filter_feature(const FeatureVector &fea);

/**
 * @brief Filter to check if the two strips have enough overlapping to be
 * considered for inference.
 *
 * @param cluster1 First cluster to compare.
 * @param cluster2 Second cluster to compare.
 * @param merged_cluster Merged cluster of the first and second clusters.
 * @return true if the two clusters pass the filter, false otherwise.
 */
bool filter_test_overlapping(const Cluster &cluster1, const Cluster &cluster2,
                             const Cluster &merged_cluster);

/**
 * @brief Filter for strip merging.
 *
 * @param cluster1 First cluster to compare.
 * @param cluster2 Second cluster to compare.
 * @param vis_folder Path to the visualization folder.
 * @return true if the two clusters pass the filter, false otherwise.
 */
bool filter_merging(const Cluster &cluster1, const Cluster &cluster2,
                    const std::string &vis_folder);
