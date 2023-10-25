/**
 * This file contains functions used to filter training sample candidates that
 * visually doesn't align with the labels. This may happen because we generate
 * multiple candidates from a single final example by taking snapshots at
 * different time stamps. A positive or negative example at the final time stamp
 * may not always be a visually correct posiitve or negative example.
 */
#pragma once

#include "SampleGeneration.h"

#include <set>
#include <vector>

struct Cluster;

struct FilteredSample {
  Sample sample;
  bool is_pos;
  enum FilterCause : std::uint32_t {
    Parameterization = 1,
    Overlapping = 2,
    Orientation = 3,
    FillIn = 4,
    Branch = 5,
    Continuity = 6,
    OverlappingInconsistence = 7,
    ShortOverlapping = 8,
    LargeAngle = 9,
  } cause;
};

/**
 * Filters training sample candidates whose orientations differ from the
 * parameterization of the strip at the final time stamp.
 *
 * @param merged_cluster The final strip with parameterization.
 * @param pos_cluster The positive cluster to filter against.
 * @return True if the sample should be filtered, false otherwise.
 */
bool filter_orientation_sample(const Cluster &merged_cluster,
                               const Cluster &pos_cluster);
/**
 * Filters positive training sample candidates which have a gap in the middle
 * (filled in at some later time stamp).
 *
 * @param merged_cluster The final strip with parameterization.
 * @param gt_cluster The ground truth cluster to compare against.
 * @return True if the sample should be filtered, false otherwise.
 */
bool filter_fillin_sample(const Cluster &merged_cluster,
                          const Cluster &gt_cluster);
/**
 * Filters negative training sample candidates which haven't formed a pronounced
 * branching.
 *
 * @param merged_cluster The final strip with parameterization.
 * @param neg_cluster The negative cluster.
 * @param selected_branches The selected branches to filter.
 * @param min_dist The minimum distance between the merged cluster and the
 * negative cluster.
 * @return True if the sample should be filtered, false otherwise.
 */
bool filter_branch_sample(
  const Cluster &merged_cluster, const Cluster &neg_cluster,
  std::vector<std::pair<std::set<size_t>, std::set<size_t>>> &selected_branches,
  const double min_dist);

/**
 * Filters the training sample candidates and returns filtered samples along
 * with positive and negative samples.
 *
 * @param input The input data to filter.
 * @param pos_samples The vector of positive samples.
 * @param neg_samples The vector of negative samples.
 * @param filtered_samples The vector of filtered samples.
 * @param matching_pair_ptr A pointer to a vector of matching pairs.
 * @param cache_folder The folder to cache the filtered samples.
 * @param output_vis_folder The folder to output visualization files.
 * @param disable_short_overlapping A flag to disable short overlapping.
 */
void filter_sample_extra(
  Input &input, std::vector<Sample> &pos_samples,
  std::vector<Sample> &neg_samples,
  std::vector<FilteredSample> &filtered_samples,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr = nullptr,
  const std::string &cache_folder = "",
  const std::string &output_vis_folder = "",
  bool disable_short_overlapping = false);
