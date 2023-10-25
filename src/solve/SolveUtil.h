/**
 * This file contains the context used for the solve steps (including the
 * feature set definitions) and helper functions used to debug solving.
 */
#pragma once

#include "../feature/Features.h"
#include "../stroke_strip_src/Cluster.h"

#include <set>

struct SolveContext {
  // Input info
  int width;
  int height;
  double norm_thickness;
  double input_thickness;

  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    matching_pair;

  double cutoff = 0.5;
  std::map<std::string, Cluster> merged_cluster_cache;
  double junc_cutoff = 0.7;
  std::vector<std::pair<size_t, size_t>> separation_sid;

  std::vector<double> split_merge_cutoffs{0.3, 0.4, 0.5};

  std::set<size_t> final_single_sid;
  double single_single_final_cutoff = 0.55;
  double single_final_cutoff = 0.45;

  // Features
  std::vector<std::string> feature_set = {
    "average_angle",
    "median_angle",
    "average_distance",
    "median_distance",
    "overlapping_length2combined_length",
    "maxC12_avg_narrowness",
    "maxC12_median_narrowness",
    "minC12_avg_narrowness",
    "minC12_median_narrowness",
    "avg_narrowness_combined",
    "median_narrowness_combined",
    "avg_distance2LocalMaxWidth",
    "median_distance2LocalMaxWidth",
    "avg_distance2LocalMinWidth",
    "median_distance2LocalMinWidth",
    "avg_distance2width_combined",
    "median_distance2width_combined",
    "maxC12_90th_LocalMedianGap2width",
    "minC12_90th_LocalMedianGap2width",
    "maxC12_10th_LocalMedianGap2width",
    "minC12_10th_LocalMedianGap2width",
    "maxC12_local_nonoverlapping2overlapping",
    "minC12_local_nonoverlapping2overlapping",
    "maxC12_AvgNonoverlapping2AvgOverlapping",
    "minC12_AvgNonoverlapping2AvgOverlapping",
    "velocity",
    "alignment"};
  std::vector<std::string> sec_feature_set = {
    "average_angle",
    "median_angle",
    "average_distance",
    "median_distance",
    "overlapping_length2combined_length",
    "maxC12_avg_narrowness",
    "maxC12_median_narrowness",
    "minC12_avg_narrowness",
    "minC12_median_narrowness",
    "avg_narrowness_combined",
    "median_narrowness_combined",
    "avg_distance2LocalMaxWidth",
    "median_distance2LocalMaxWidth",
    "avg_distance2LocalMinWidth",
    "median_distance2LocalMinWidth",
    "avg_distance2width_combined",
    "median_distance2width_combined",
    "maxC12_90th_LocalMedianGap2width",
    "minC12_90th_LocalMedianGap2width",
    "maxC12_10th_LocalMedianGap2width",
    "minC12_10th_LocalMedianGap2width",
    "maxC12_local_nonoverlapping2overlapping",
    "minC12_local_nonoverlapping2overlapping",
    "maxC12_AvgNonoverlapping2AvgOverlapping",
    "minC12_AvgNonoverlapping2AvgOverlapping",
    "velocity",
    "alignment",
    "average_distance2global_average_average_distance",
    "median_distance2global_average_median_distance",
    "average_distance2global_median_average_distance",
    "median_distance2global_median_median_distance",
    "average_distance2global_p90th_average_distance",
    "median_distance2global_p90th_median_distance"};
  std::vector<std::string> typical_feature_set;
  std::map<size_t, FeatureVector> cluster_features;

  // Debug
  std::map<size_t, size_t>
    gt_sid2cid; // Note this is only used for the GT debugging mode!

  // Stages
  enum SolveStage : std::uint32_t {
    Initial = 1,
    Splitting = 2,
    Split_Merging = 3,
    Merging = 4,
    Endpoint = 5,
  } stage = SolveStage::Initial;
  bool to_junc_sep;

  // Use an epsilon value to determine similar predictions
  double tier_epsilon = 0.025;

  // Timing counters
  double feature_timer = 0;
  double classifier_timer = 0;

  bool count_skip = false;
  size_t ours_count = 0;
  size_t ours_call_count = 0;
  size_t skip_count = 0;
  size_t call_count = 0;

  size_t single_count = 0;
  size_t single_call_count = 0;
};

extern SolveContext solve_context;
extern size_t seed_xsec_count;
extern std::string tier_breaker_feature;
extern double seed_prob_threshold;

extern double width_split_threshold;

FeatureVector inf_feature(int width, int height);
std::vector<double> pick_features(const FeatureVector &fea,
                                  const std::vector<std::string> &feature_set,
                                  std::vector<std::string> &feature_set_out);

FeatureVector compute_feature(Cluster &cluster1, Cluster &cluster2,
                              Cluster &merged_cluster, int width, int height,
                              double norm_thickness, std::string vis_folder);

void print_csv_header(const std::vector<std::string> &feature_set, bool typical,
                      std::string &csv_str);
void print_csv_row(const std::vector<std::string> &feature_set,
                   const std::vector<double> &feature_vec, double score,
                   double typical, int base_number, int target_number,
                   bool is_merging, bool decision, std::string fit_output,
                   std::string isoline, std::string &csv_str);

bool to_merge(const Cluster &cluster1, const Cluster &cluster2, double score);
