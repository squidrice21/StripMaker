#pragma once

#include "../feature/Features.h"
#include "../stroke_strip_src/Cluster.h"

struct PredictionInfo {
  std::vector<std::string>
    feature_set; /**< Vector of feature names used in the prediction. */
  std::vector<double>
    feature_vec; /**< Vector of feature values used in the prediction. */
  double score; /**< Prediction score. */
  int base_number; /**< Number of strokes in the base strip. */
  int target_number; /**< Number of strokes in the target strip. */
  bool is_merging; /**< Flag indicating whether the prediction is for merging
                      existing strips. */
  bool decision; /**< Prediction decision. */
  std::string hash_str_comb; /**< Hash string of the combined strip. */
};

/**
 * @brief Compares a strip with other strips and updates the scores vector
 * with the prediction results.
 *
 * @param curr_c The current strip to compare.
 * @param in_cid The ID of the current strip.
 * @param clusters Map of all strips.
 * @param later_cluster_only Flag indicating whether to compare with later
 * clusters only.
 * @param split_within Flag indicating whether to allow splitting within
 * clusters.
 * @param merged_cluster_ptr Pointer to the merged strip.
 * @param width Width of the input image.
 * @param height Height of the input image.
 * @param norm_thickness Normalized thickness of the strokes.
 * @param scores Vector of pairs of prediction scores and prediction
 * information.
 * @param vis_folder Path to the visualization folder.
 */
void compare_clusters(Cluster &curr_c, int in_cid,
                      std::map<size_t, Cluster> &clusters,
                      bool later_cluster_only, bool split_within,
                      Cluster *merged_cluster_ptr, int width, int height,
                      double norm_thickness,
                      std::vector<std::pair<double, PredictionInfo>> &scores,
                      std::string vis_folder);

/**
 * @brief Compares a single strip with another strip and updates the
 * prediction information and score.
 *
 * @param cmpr_c The strip to compare with.
 * @param curr_c The current strip.
 * @param merged_cluster The merged strip.
 * @param width Width of the input image.
 * @param height Height of the input image.
 * @param norm_thickness Normalized thickness of the strokes.
 * @param prediction_info Prediction information to update.
 * @param prediction Prediction score to update.
 * @param vis_folder Path to the visualization folder.
 * @param seed_picking Flag indicating whether to use seed picking.
 * @param split_within Flag indicating whether to allow splitting within
 * strips.
 */
void compare_cluster_single(Cluster &cmpr_c, Cluster &curr_c,
                            Cluster &merged_cluster, int width, int height,
                            double norm_thickness,
                            PredictionInfo &prediction_info, double &prediction,
                            std::string vis_folder, bool seed_picking = false,
                            bool split_within = false);

/**
 * @brief Recursively compares a strip with other strips and merges them until
 * the probability is lower than threshold or the max recursion level is
 * reached.
 *
 * @param input Input sketch.
 * @param cid ID of the strip to compare.
 * @param clusters Map of all strips.
 * @param merge_history Vector of strip IDs that have been merged.
 * @param vis_folder Path to the visualization folder.
 * @param csv_str String to store the CSV data.
 * @param in_place Flag indicating whether to perform the comparison in place.
 * @param merged_cluster_ptr Pointer to the merged strip.
 * @return int The ID of the merged strip.
 */
int compare_cluster_cluster_recursive(const Input &input, size_t cid,
                                      std::map<size_t, Cluster> &clusters,
                                      std::vector<size_t> &merge_history,
                                      std::string vis_folder,
                                      std::string &csv_str,
                                      bool in_place = false,
                                      Cluster *merged_cluster_ptr = nullptr);

/**
 * @brief The initial local temporal consolidation solve. This step follows the
 * input drawing order. Note that it delays decisions on strip pairs with
 * insufficient overlapping.
 *
 * @param input Input sketch.
 * @param clusters Map of all strips.
 * @param lazy_cluster_comparison Flag indicating whether to use lazy strip
 * comparison.
 * @param vis_folder Path to the visualization folder.
 */
void solve(const Input &input, std::map<size_t, Cluster> &clusters,
           bool lazy_cluster_comparison, std::string vis_folder);
