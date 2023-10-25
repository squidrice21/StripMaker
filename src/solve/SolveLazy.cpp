#include "SolveLazy.h"

#include "../Logger.h"
#include "../Util.h"
#include "../classifier/forest.h"
#include "../classifier/forest_sec.h"
#include "../endpoint/JunctionSeparation.h"
#include "../feature/FeatureIO.h"
#include "../feature/Filters.h"
#include "../stroke_strip_src/Serialization.h"
#include "SolveUtil.h"

#include <map>
#include <oneapi/tbb/parallel_for_each.h>

#include <limits>
#include <regex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

SolveContext solve_context;

size_t stroke_step = 0;

bool sum_comparisons = false;
bool single_comparisons = false;

void compare_cluster_single(Cluster &cmpr_c, Cluster &curr_c,
                            Cluster &merged_cluster, int width, int height,
                            double norm_thickness,
                            PredictionInfo &prediction_info, double &prediction,
                            std::string vis_folder, bool seed_picking,
                            bool split_within) {
  if (sum_comparisons)
    solve_context.call_count++;
  else if (solve_context.stage == SolveContext::SolveStage::Initial)
    solve_context.ours_call_count++;

  if (single_comparisons &&
      solve_context.stage == SolveContext::SolveStage::Initial)
    solve_context.single_call_count++;

  // Only apply this filter when fitting exists (aka not in the spliting stage)
  if (!seed_picking &&
      filter_sample(cmpr_c, std::vector<int>(), curr_c, std::vector<int>())) {
    return;
  }

  // Assume the comparison cluster is parameterized successfully
  assert(!cmpr_c.strokes.empty());

  // Only apply this filter when fitting exists (aka not in the spliting stage)
  if (!split_within && !curr_c.fit.centerline.empty() &&
      filter_sample_overlapping(cmpr_c, std::vector<int>(), curr_c,
                                std::vector<int>(),
                                &solve_context.matching_pair, "", vis_folder)) {
    return;
  }

  // Apply this filter in the final merging
  if (solve_context.stage == SolveContext::SolveStage::Merging) {
    if (filter_merging(cmpr_c, curr_c, vis_folder)) {
      logger().info("Skip: {} vs {}", curr_c.strokes.front().stroke_ind,
                    cmpr_c.strokes.front().stroke_ind);
      return;
    }
  }

  // Compute feature
  FeatureVector fea = compute_feature(cmpr_c, curr_c, merged_cluster, width,
                                      height, norm_thickness, vis_folder);
  if (sum_comparisons)
    solve_context.skip_count++;
  else if (solve_context.stage == SolveContext::SolveStage::Initial)
    solve_context.ours_count++;

  if (single_comparisons &&
      solve_context.stage == SolveContext::SolveStage::Initial)
    solve_context.single_count++;

  // Directly determine prob = 0 if the side-by-side section does not
  // exist or is too short

  double alignment = -1;
  if (merged_cluster.obj_term_values.size() >= 2 &&
      merged_cluster.obj_term_values.count("alignment"))
    alignment = merged_cluster.obj_term_values.at("alignment");
  else
    return;

  std::vector<double> fea_vec =
    (solve_context.stage == SolveContext::SolveStage::Initial)
      ? pick_features(fea, solve_context.feature_set,
                      prediction_info.feature_set)
      : pick_features(fea, solve_context.sec_feature_set,
                      prediction_info.feature_set);
  prediction_info.feature_set =
    (solve_context.stage == SolveContext::SolveStage::Initial)
      ? solve_context.feature_set
      : solve_context.sec_feature_set;
  prediction_info.feature_vec = fea_vec;

  // Skip overlapping filter when splitting within a cluster
  if (fea_vec.front() < std::numeric_limits<double>::infinity() &&
      (split_within ||
       !filter_test_overlapping(cmpr_c, curr_c, merged_cluster))) {

    std::vector<std::string> alignment_str, avg_angle_str;
    auto alignments =
      pick_features(fea, std::vector<std::string>{"alignment"}, alignment_str);
    double alignment = alignments[0];
    auto avg_angles = pick_features(
      fea, std::vector<std::string>{"average_angle"}, avg_angle_str);
    double avg_angle = avg_angles[0];
    double prob = 0;

    // We set the probability to 0 and pick the furthest pair when splitting
    if (seed_picking &&
        filter_sample(cmpr_c, std::vector<int>(), curr_c, std::vector<int>())) {
      prediction_info.score = 0;
      return;
    }

    double cur_alignment_threshold = alignment_threshold;
    if (solve_context.stage == SolveContext::SolveStage::Merging)
      cur_alignment_threshold = alignment_threshold / 2;
    if (alignment <= cur_alignment_threshold &&
        avg_angle <= large_angle_threshold) {
      if (solve_context.gt_sid2cid.empty()) {
        auto begin_clf = std::chrono::high_resolution_clock::now();
        if (solve_context.stage == SolveContext::SolveStage::Initial) {
          normalize(fea_vec);
          prob = clf(fea_vec);
        } else if (solve_context.stage >= SolveContext::SolveStage::Splitting) {
          normalize_sec(fea_vec);
          prob = clf_sec(fea_vec);
        }
        auto end_clf = std::chrono::high_resolution_clock::now();
        solve_context.classifier_timer +=
          std::chrono::duration_cast<std::chrono::milliseconds>(end_clf -
                                                                begin_clf)
            .count() /
          1000.0;
      } else {
        assert(
          solve_context.gt_sid2cid.count(curr_c.strokes.front().stroke_ind) &&
          solve_context.gt_sid2cid.count(cmpr_c.strokes.front().stroke_ind));
        // It's sufficient to compare a single arbitrary stroke from each
        // cluster since we'll never cluster strokes from different
        // clusters in the first place, given GT.
        if (solve_context.gt_sid2cid.at(curr_c.strokes.front().stroke_ind) ==
            solve_context.gt_sid2cid.at(cmpr_c.strokes.front().stroke_ind))
          prob = 1;
        else
          prob = 0;
      }
    }
    prediction = prob;
  } else if (filter_test_overlapping(cmpr_c, curr_c, merged_cluster)) {
    // Debug printout
    prediction_info.score = -0.1;
    prediction_info.decision = to_merge(cmpr_c, curr_c, prediction);
    std::string hash_str_comb = get_hash_str_comb(cmpr_c, curr_c);
    prediction_info.hash_str_comb = hash_str_comb;
  } else {
    prediction = 0;
  }

  if (prediction > 0 || (seed_picking && prediction >= 0)) {
    prediction_info.score = prediction;
    prediction_info.decision = to_merge(cmpr_c, curr_c, prediction);
    std::string hash_str_comb = get_hash_str_comb(cmpr_c, curr_c);
    prediction_info.hash_str_comb = hash_str_comb;
  }
}

void compare_clusters(Cluster &curr_c, int in_cid,
                      std::map<size_t, Cluster> &clusters,
                      bool later_cluster_only, bool split_within,
                      Cluster *merged_cluster_ptr, int width, int height,
                      double norm_thickness,
                      std::vector<std::pair<double, PredictionInfo>> &scores,
                      std::string vis_folder) {
  std::vector<size_t> keys;
  std::map<size_t, double> predictions;
  std::map<size_t, PredictionInfo> prediction_info;
  std::set<std::pair<double, int>> dedup_scores;
  keys.reserve(clusters.size());
  for (auto &cid2c : clusters) {
    keys.emplace_back(cid2c.first);
    predictions[cid2c.first] = -1;
    prediction_info[cid2c.first] = PredictionInfo();
  }

  assert(in_cid < 0 || clusters.count(in_cid));

  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<size_t>(0u, keys.size()),
    [&](const oneapi::tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        int cid = keys[i];
        // Skip if it's the same cluster
        if (in_cid >= 0 && in_cid == cid)
          continue;
        if (later_cluster_only && cid < in_cid)
          continue;
        assert(clusters.count(cid));
        Cluster &cmpr_c = clusters[cid];

        std::string hash_str_comb = get_hash_str_comb(cmpr_c, curr_c);
        Cluster &merged_cluster =
          (!merged_cluster_ptr)
            ? solve_context.merged_cluster_cache[hash_str_comb]
            : *merged_cluster_ptr;

        compare_cluster_single(cmpr_c, curr_c, merged_cluster, width, height,
                               norm_thickness, prediction_info[cid],
                               predictions[cid], vis_folder, false,
                               split_within);

        if (predictions[cid] > 0) {
          prediction_info[cid].base_number =
            (in_cid >= 0) ? in_cid : curr_c.strokes.front().stroke_ind;
          prediction_info[cid].target_number = cid;
          prediction_info[cid].is_merging = (in_cid >= 0);
          logger().info("compare_cluster_single: {} - {} => {}",
                        prediction_info[cid].base_number,
                        prediction_info[cid].target_number, predictions[cid]);
        } else if (prediction_info[cid].score < 0) {
          // Debug printout for short overlapping filter
          prediction_info[cid].base_number =
            (in_cid >= 0) ? in_cid : curr_c.strokes.front().stroke_ind;
          prediction_info[cid].target_number = cid;
          prediction_info[cid].is_merging = (in_cid >= 0);
        }
      }
    });

  int key_offset = 10000;
  for (auto const &p : prediction_info) {
    if (p.second.score != 0) {
      int cid_key = -1;
      const auto &c = clusters[p.first];
      for (auto const &s : c.strokes)
        cid_key = std::max(cid_key, (int)s.stroke_ind);
      cid_key = cid_key * key_offset + p.first;
      cid_key = -cid_key;
      dedup_scores.emplace(std::make_pair(-p.second.score, cid_key));
    }
  }

  for (auto const &p_c : dedup_scores) {
    int cid = p_c.second;
    cid = -cid;
    cid = cid % key_offset;
    assert(prediction_info.count(cid));
    assert(!prediction_info[cid].feature_set.empty());
    scores.emplace_back(std::make_pair(p_c.first, prediction_info[cid]));
  }

  size_t scores_i = 0;
  for (auto &s : scores) {
    if (scores_i > 0) {
      s.second.decision = false;
    }
    scores_i++;
  }
}

int compare_cluster_cluster(const Input &input, size_t cid,
                            std::map<size_t, Cluster> &clusters, bool in_place,
                            bool later_cluster_only,
                            Cluster *merged_cluster_ptr, std::string vis_folder,
                            std::string &csv_str) {
  double norm_thickness = 1;
  std::vector<std::pair<double, PredictionInfo>> all_scores;
  bool split_within =
    solve_context.stage == SolveContext::SolveStage::Splitting ||
    solve_context.stage == SolveContext::SolveStage::Split_Merging;
  compare_clusters(clusters[cid], cid, clusters, later_cluster_only,
                   split_within, merged_cluster_ptr, input.width, input.height,
                   norm_thickness, all_scores, vis_folder);

  if (sum_comparisons)
    return -1;

  std::map<size_t, size_t> s2c;
  std::vector<std::set<size_t>> separation_cid;
  for (auto const &cc : clusters) {
    for (auto const &s : cc.second.strokes) {
      s2c[s.stroke_ind] = cc.first;
    }
  }
  for (auto const &ss : solve_context.separation_sid) {
    assert((s2c.count(ss.first) && s2c.count(ss.second)) ||
           (!s2c.count(ss.first) && !s2c.count(ss.second)));
    if (s2c.count(ss.first)) {
      separation_cid.emplace_back();
      separation_cid.back().emplace(s2c[ss.first]);
      separation_cid.back().emplace(s2c[ss.second]);
    }
  }
  std::vector<std::pair<double, PredictionInfo>> scores;
  scores.reserve(all_scores.size());
  for (auto const &p_score : all_scores) {
    bool seen_sep = false;
    for (auto const &cc : separation_cid) {
      if (cc.count(cid) && cc.count(p_score.second.target_number)) {
        seen_sep = true;
        break;
      }
    }
    // Junction separation
    if (!seen_sep && solve_context.stage == SolveContext::SolveStage::Merging &&
        to_merge(clusters[p_score.second.target_number], clusters[cid],
                 std::abs(p_score.first))) {
      std::map<size_t, Cluster> sep_clusters;
      sep_clusters[p_score.second.target_number] =
        clusters[p_score.second.target_number];
      sep_clusters[cid] = clusters[cid];
      std::string hash_str_comb = get_hash_str_comb(
        clusters[p_score.second.target_number], clusters[cid]);
      Cluster merged_cluster =
        solve_context.merged_cluster_cache[hash_str_comb];
      merged_cluster.fit.cluster_idx = 0;
      bool are_separable = are_separable_clusters(input, clusters, sep_clusters,
                                                  merged_cluster, vis_folder);
      if (are_separable) {
        seen_sep = true;
      }
    }
    if (!seen_sep)
      scores.emplace_back(p_score);
  }

  size_t merged_c2 =
    (!scores.empty()) ? scores.begin()->second.target_number : 0;
  if (!scores.empty()) {
    logger().info("{} vs {}: {}", clusters[cid].strokes.front().stroke_ind,
                  clusters[merged_c2].strokes.front().stroke_ind,
                  std::abs(scores.begin()->first));
  }
  if (!scores.empty() && to_merge(clusters[merged_c2], clusters[cid],
                                  std::abs(scores.begin()->first))) {
    std::string hash_str_comb =
      get_hash_str_comb(clusters[merged_c2], clusters[cid]);
    size_t updated_c;

    // Also update the fitting width in the final merging. This is needed for
    // filtering
    if (solve_context.stage == SolveContext::SolveStage::Merging) {
      Input input_fit;
      input_fit.clusters[0] = solve_context.merged_cluster_cache[hash_str_comb];
      context.widths = true;
      input_fit.clusters[0].fit = fea_fitting.fit(&input_fit).begin()->second;
      solve_context.merged_cluster_cache[hash_str_comb].fit =
        input_fit.clusters[0].fit;
      context.widths = false;
    }

    if (!merged_cluster_ptr) {
      assert(solve_context.merged_cluster_cache.count(hash_str_comb));
      if (!in_place) {
        clusters[merged_c2] = solve_context.merged_cluster_cache[hash_str_comb];
        clusters.erase(cid);
        updated_c = merged_c2;
      } else {
        clusters[cid] = solve_context.merged_cluster_cache[hash_str_comb];
        clusters.erase(merged_c2);
        updated_c = cid;
      }
    } else { // Use a given common parameterization
      if (!in_place) {
        for (size_t i = 0; i < clusters[cid].strokes.size(); ++i) {
          clusters[merged_c2].strokes.emplace_back(clusters[cid].strokes[i]);
          clusters[merged_c2].original_input_strokes.emplace_back(
            clusters[cid].original_input_strokes[i]);
        }
        clusters.erase(cid);
        updated_c = merged_c2;

        std::sort(clusters[merged_c2].strokes.begin(),
                  clusters[merged_c2].strokes.end(),
                  [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                    return a.stroke_ind < b.stroke_ind;
                  });
        std::sort(clusters[merged_c2].original_input_strokes.begin(),
                  clusters[merged_c2].original_input_strokes.end(),
                  [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                    return a.stroke_ind < b.stroke_ind;
                  });
      } else {
        for (size_t i = 0; i < clusters[merged_c2].strokes.size(); ++i) {
          clusters[cid].strokes.emplace_back(clusters[merged_c2].strokes[i]);
          clusters[cid].original_input_strokes.emplace_back(
            clusters[merged_c2].original_input_strokes[i]);
        }
        clusters.erase(merged_c2);
        updated_c = cid;

        std::sort(clusters[cid].strokes.begin(), clusters[cid].strokes.end(),
                  [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                    return a.stroke_ind < b.stroke_ind;
                  });
        std::sort(clusters[cid].original_input_strokes.begin(),
                  clusters[cid].original_input_strokes.end(),
                  [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                    return a.stroke_ind < b.stroke_ind;
                  });
      }
    }

    if (!vis_folder.empty()) {
      double typical_score = -1;
      std::string fit_output = vis_folder + "/fit/" +
                               scores.begin()->second.hash_str_comb +
                               "_fit.svg";
      std::string isoline = vis_folder + "/isoline/" +
                            scores.begin()->second.hash_str_comb +
                            "_isolines.svg";
      if (!in_place)
        print_csv_row(
          scores.begin()->second.feature_set,
          scores.begin()->second.feature_vec, scores.begin()->second.score,
          typical_score, scores.begin()->second.base_number,
          scores.begin()->second.target_number,
          scores.begin()->second.is_merging, scores.begin()->second.decision,
          fit_output, isoline, csv_str);
      else
        print_csv_row(
          scores.begin()->second.feature_set,
          scores.begin()->second.feature_vec, scores.begin()->second.score,
          typical_score, scores.begin()->second.target_number,
          scores.begin()->second.base_number, scores.begin()->second.is_merging,
          scores.begin()->second.decision, fit_output, isoline, csv_str);
    }

    return merged_c2;
  } else if (!scores.empty()) {
    if (!vis_folder.empty()) {
      for (auto const &score : scores) {
        std::string fit_output =
          vis_folder + "/fit/" + score.second.hash_str_comb + "_fit.svg";
        std::string isoline = vis_folder + "/isoline/" +
                              score.second.hash_str_comb + "_isolines.svg";
        print_csv_row(score.second.feature_set, score.second.feature_vec,
                      score.second.score, -1, score.second.base_number,
                      score.second.target_number, score.second.is_merging,
                      score.second.decision, fit_output, isoline, csv_str);
      }
    }
  }

  return -1;
}

int compare_cluster_cluster_recursive(const Input &input, size_t cid,
                                      std::map<size_t, Cluster> &clusters,
                                      std::vector<size_t> &merge_history,
                                      std::string vis_folder,
                                      std::string &csv_str, bool in_place,
                                      Cluster *merged_cluster_ptr) {

  size_t max_recursion = 3;
  if (solve_context.stage == SolveContext::SolveStage::Merging)
    max_recursion = 1;
  if (solve_context.stage == SolveContext::SolveStage::Split_Merging)
    max_recursion = 5;
  int recur_cid = cid;
  int ret_cid = -1;
  // For the separate merging step, we first only compare to later clusters to
  // avoid duplicate comparisons
  bool later_cluster_only = in_place;
  for (size_t i = 0; i < max_recursion && clusters.size() > 1; ++i) {
    ret_cid = compare_cluster_cluster(input, recur_cid, clusters, in_place,
                                      later_cluster_only, merged_cluster_ptr,
                                      vis_folder, csv_str);
    logger().info("compare_cluster_cluster_recursive: {} - return: {}",
                  recur_cid, ret_cid);
    if (ret_cid < 0)
      break;

    // Compare to all clusters since the current one got updated
    later_cluster_only = false;
    merge_history.emplace_back(ret_cid);
    recur_cid = (in_place) ? recur_cid : ret_cid;
  }
  return ret_cid;
}

void solve(const Input &input, std::map<size_t, Cluster> &clusters,
           bool lazy_cluster_comparison, std::string vis_folder) {
  solve_context.stage = SolveContext::SolveStage::Initial;

  std::string csv_str;
  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : input.clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      sid2s[s.stroke_ind] = s;
    }
  }

  int last_updated_cluster = -1;
  double norm_thickness = 1;
  std::vector<std::pair<double, PredictionInfo>> scores;

  for (auto const &i2s : sid2s) {
    size_t i = i2s.first;
    stroke_step = i;
    scores.clear();

    Cluster curr_c;
    curr_c.strokes.emplace_back(i2s.second);
    curr_c.original_input_strokes.emplace_back(i2s.second);
    bool to_orient = false;
    bool test_time = true;
    fit_cluster(norm_thickness, input.width, input.height, curr_c, to_orient,
                test_time, &solve_context.matching_pair);

    // Compare to all existing clusters
    Cluster *merged_cluster_ptr = nullptr; // No common cached parameterization
    bool later_cluster_only = false;
    bool split_within = false;
    if (solve_context.count_skip)
      single_comparisons = true;
    compare_clusters(curr_c, -1, clusters, later_cluster_only, split_within,
                     merged_cluster_ptr, input.width, input.height,
                     norm_thickness, scores, vis_folder);
    if (solve_context.count_skip)
      single_comparisons = false;

    size_t merged_c2 =
      (!scores.empty()) ? scores.begin()->second.target_number : 0;
    int updated_cid = -1;
    // Make decision: no mering, or merging to an existing cluster.
    if (scores.empty() || !to_merge(clusters[merged_c2], curr_c,
                                    std::abs(scores.begin()->first))) {
      size_t new_c = (clusters.empty()) ? 0 : (clusters.rbegin()->first + 1);
      clusters[new_c] = curr_c;

      updated_cid = new_c;

      if (!scores.empty() && !vis_folder.empty()) {
        for (auto const &score : scores) {
          std::string fit_output =
            vis_folder + "/fit/" + score.second.hash_str_comb + "_fit.svg";
          std::string isoline = vis_folder + "/isoline/" +
                                score.second.hash_str_comb + "_isolines.svg";

          print_csv_row(score.second.feature_set, score.second.feature_vec,
                        score.second.score,
                        (!solve_context.typical_feature_set.empty()) ? 1 : -1,
                        score.second.base_number, score.second.target_number,
                        score.second.is_merging, score.second.decision,
                        fit_output, isoline, csv_str);
        }
      }
      if (!vis_folder.empty()) {
        std::string fit_output = vis_folder + "/fit/";
        std::string isoline = vis_folder + "/isoline/";
        FeatureVector fea = inf_feature(input.width, input.height);
        std::vector<std::string> feature_set_out;
        std::vector<double> fea_vec =
          (solve_context.stage == SolveContext::SolveStage::Initial)
            ? pick_features(fea, solve_context.feature_set, feature_set_out)
            : pick_features(fea, solve_context.sec_feature_set,
                            feature_set_out);
        print_csv_row(feature_set_out, fea_vec, -1,
                      (!solve_context.typical_feature_set.empty()) ? 1 : -1,
                      curr_c.strokes.front().stroke_ind, new_c, true, true,
                      fit_output, isoline, csv_str);
      }
    } else {
      size_t merged_c = scores.begin()->second.target_number;
      std::string hash_str_comb = get_hash_str_comb(clusters[merged_c], curr_c);
      assert(solve_context.merged_cluster_cache.count(hash_str_comb));
      clusters[merged_c] = solve_context.merged_cluster_cache[hash_str_comb];
      updated_cid = merged_c;

      if (!vis_folder.empty()) {
        double typical_score = -1;
        std::string fit_output = vis_folder + "/fit/" +
                                 scores.begin()->second.hash_str_comb +
                                 "_fit.svg";
        std::string isoline = vis_folder + "/isoline/" +
                              scores.begin()->second.hash_str_comb +
                              "_isolines.svg";
        print_csv_row(
          scores.begin()->second.feature_set,
          scores.begin()->second.feature_vec, scores.begin()->second.score,
          typical_score, scores.begin()->second.base_number,
          scores.begin()->second.target_number,
          scores.begin()->second.is_merging, scores.begin()->second.decision,
          fit_output, isoline, csv_str);
      }
    }

    {
      std::vector<size_t> merge_history;

      // If a cluster gets merged into, try to merge this updated cluster.
      // Run cluster-cluster comparisons recursively.
      if (!lazy_cluster_comparison || i == sid2s.rbegin()->first) {
        bool extra_compare =
          (lazy_cluster_comparison && i == sid2s.rbegin()->first &&
           last_updated_cluster >= 0 && last_updated_cluster != updated_cid);

        if (extra_compare) {
          compare_cluster_cluster_recursive(input, last_updated_cluster,
                                            clusters, merge_history, vis_folder,
                                            csv_str, false);
        }

        if (clusters.count(updated_cid)) {
          merge_history.clear();
          compare_cluster_cluster_recursive(input, updated_cid, clusters,
                                            merge_history, vis_folder, csv_str,
                                            false);
          if (!merge_history.empty())
            updated_cid = merge_history.back();
        }
      } else { // Lazy strategy: Only update the last updated cluster if the new
        // stroke is added into a new cluster
        if (solve_context.count_skip) {
          sum_comparisons = true;
          compare_cluster_cluster_recursive(input, updated_cid, clusters,
                                            merge_history, vis_folder, csv_str,
                                            false);
          merge_history.clear();
          sum_comparisons = false;
        }
        if (last_updated_cluster >= 0 && last_updated_cluster != updated_cid) {
          compare_cluster_cluster_recursive(input, last_updated_cluster,
                                            clusters, merge_history, vis_folder,
                                            csv_str, false);

          // This is for the case that the latest updated cluster is merged by
          // the current lazy comparison iteration (this would leave the latest
          // updated cluster empty in the record)
          if (!clusters.count(updated_cid)) {
            assert(!merge_history.empty());
            last_updated_cluster = merge_history.back();
          }
        }
      }
    }

    for (const auto &c : clusters)
      assert(!c.second.strokes.empty());

    if (clusters.count(updated_cid))
      last_updated_cluster = updated_cid;
  }

  if (!vis_folder.empty()) {
    std::ofstream csv_ofs(vis_folder + "/features.csv");
    csv_ofs.write(csv_str.c_str(), csv_str.size());
    csv_ofs.close();
  }
}
