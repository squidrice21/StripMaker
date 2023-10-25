#include "SolveEndSplit.h"

#include "../Logger.h"
#include "../Util.h"
#include "../endpoint/JunctionSeparation.h"
#include "../feature/Features.h"
#include "SolveEndMerge.h"
#include "SolveLazy.h"
#include "SolveUtil.h"
#include "glm/detail/type_vec.hpp"

#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

size_t seed_xsec_count = 3;
std::string tier_breaker_feature = "average_distance";
double seed_prob_threshold = 0.6;

double width_split_threshold = 5;

int larger_split_xsec_ind(const Cluster &cluster, size_t sid1, size_t sid2) {
  assert(sid1 != sid2);

  std::map<double, size_t> overlapping_xsec_dist;
  for (size_t i = 0; i < cluster.xsecs.size(); ++i) {
    const auto &xsec = cluster.xsecs[i];
    if (xsec.points.size() < 2)
      continue;
    std::map<size_t, glm::dvec2> sid2points;
    double gap = -1;
    for (size_t j = 0; j < xsec.points.size(); ++j) {
      const auto &p = xsec.points[j];
      if (p.stroke_ind == sid1 || p.stroke_ind == sid2)
        sid2points.emplace(p.stroke_ind, p.point);
      if (j > 0)
        gap = std::max(
          gap, glm::distance(xsec.points[j].point, xsec.points[j - 1].point));
    }

    if (sid2points.size() == 2) {
      double dist = gap;
      overlapping_xsec_dist[-dist] = i;
    }
  }

  if (overlapping_xsec_dist.empty())
    return -1;
  return overlapping_xsec_dist.begin()->second;
}

bool is_cut_stroke(const Cluster::Stroke &s) {
  for (auto const &p : solve_context.matching_pair) {
    if (p.first.first == s.stroke_ind || p.second.first == s.stroke_ind)
      return true;
  }
  return false;
}

void find_seed_xsecs_probability(
  Cluster &cluster, std::vector<std::pair<size_t, size_t>> &seeds) {
  if (cluster.strokes.empty())
    return;
  size_t in_cid = cluster.strokes.front().cluster_ind;
  std::map<size_t, Cluster::Stroke> sid2s, sid2orig_s;
  for (const auto &s : cluster.strokes) {
    sid2s[s.stroke_ind] = s;
  }
  for (const auto &s : cluster.original_input_strokes) {
    sid2orig_s[s.stroke_ind] = s;
  }

  std::vector<std::pair<std::pair<double, double>, std::pair<size_t, size_t>>>
    sort_seeds;
  int pick_i = -1;
  for (size_t i = 0; i + 1 < cluster.strokes.size(); ++i) {
    for (size_t j = i + 1; j < cluster.strokes.size(); ++j) {
      size_t sid1 = cluster.strokes[i].stroke_ind;
      size_t sid2 = cluster.strokes[j].stroke_ind;
      // Don't use cut stroke as seed stroke since if it were cut due to high
      // curvature, it would be hard to merge it back.
      if (is_cut_stroke(cluster.strokes[i]) ||
          is_cut_stroke(cluster.strokes[j]))
        continue;
      if (!is_overlapping_stroke(cluster, sid1, sid2))
        continue;
      Cluster curr_c;
      curr_c.strokes.emplace_back(sid2s[sid1]);
      curr_c.original_input_strokes.emplace_back(sid2orig_s[sid1]);

      Cluster cmp_c;
      cmp_c.strokes.emplace_back(sid2s[sid2]);
      cmp_c.original_input_strokes.emplace_back(sid2orig_s[sid2]);

      double prediction = -1;
      PredictionInfo prediction_info;
      prediction_info.score = -1;
      bool seed_picking = true;
      bool split_within = false;
      compare_cluster_single(cmp_c, curr_c, cluster, solve_context.width,
                             solve_context.height, solve_context.norm_thickness,
                             prediction_info, prediction, "", seed_picking,
                             split_within);

      // Pick this as a potential seed if it suggests not merging
      if (prediction_info.score >= 0 && !prediction_info.feature_vec.empty()) {
        if (pick_i < 0)
          for (pick_i = 0; pick_i < (int)prediction_info.feature_set.size();
               ++pick_i) {
            if (prediction_info.feature_set[pick_i] == tier_breaker_feature)
              break;
          }

        sort_seeds.emplace_back(std::make_pair(
          std::make_pair(prediction_info.score,
                         prediction_info.feature_vec[pick_i]),
          std::make_pair(std::min(sid1, sid2), std::max(sid1, sid2))));
      }
    }
  }

  std::sort(
    sort_seeds.begin(), sort_seeds.end(),
    [](const std::pair<std::pair<double, double>, std::pair<size_t, size_t>> &a,
       const std::pair<std::pair<double, double>, std::pair<size_t, size_t>>
         &b) {
      double prob_a = a.first.first;
      double prob_b = b.first.first;
      double dist_a = a.first.second;
      double dist_b = b.first.second;
      return prob_a < prob_b || (prob_a == prob_b && dist_a > dist_b) ||
             (prob_a == prob_b && dist_a == dist_b &&
              std::abs(int(a.second.first) - int(a.second.second)) >
                std::abs(int(b.second.first) - int(b.second.second)));
    });

  std::vector<std::pair<std::pair<double, double>, std::pair<size_t, size_t>>>
    sort_epsilon_seeds = sort_seeds;

  for (size_t i = 0; i < seed_xsec_count && i < sort_epsilon_seeds.size();
       ++i) {
    if (sort_epsilon_seeds[i].first.first < seed_prob_threshold)
      seeds.emplace_back(sort_epsilon_seeds[i].second);

    if (sort_epsilon_seeds[i].first.first < seed_prob_threshold && i == 0)
      logger().info("Try splitting cluster {} by {},{} with p={}", in_cid,
                    sort_epsilon_seeds[i].second.first,
                    sort_epsilon_seeds[i].second.second,
                    sort_epsilon_seeds[i].first.first);
    else if (i == 0)
      logger().info("Cannot split cluster {} by {},{} with p={}", in_cid,
                    sort_epsilon_seeds[i].second.first,
                    sort_epsilon_seeds[i].second.second,
                    sort_epsilon_seeds[i].first.first);
  }
}

void split_from_seed(Cluster &cluster, const std::pair<size_t, size_t> &seed,
                     std::map<size_t, Cluster> &clusters,
                     std::string vis_folder, std::string &csv_str,
                     size_t csv_cluster_offset) {
  std::map<size_t, Cluster::Stroke> sid2s, sid2orig_s;
  for (const auto &s : cluster.strokes) {
    sid2s[s.stroke_ind] = s;
  }
  for (const auto &s : cluster.original_input_strokes) {
    sid2orig_s[s.stroke_ind] = s;
  }
  std::set<size_t> processed_sids;

  // Find the two seeds
  std::vector<size_t> seed_sid{std::min(seed.first, seed.second),
                               std::max(seed.first, seed.second)};
  for (auto seed : seed_sid) {
    Cluster curr_c1;
    curr_c1.strokes.emplace_back(sid2s[seed]);
    curr_c1.original_input_strokes.emplace_back(sid2orig_s[seed]);
    size_t new_c = (clusters.empty()) ? 0 : (clusters.rbegin()->first + 1);
    clusters[new_c] = curr_c1;
    processed_sids.emplace(seed);

    if (!vis_folder.empty()) {
      std::string fit_output = vis_folder + "/fit/";
      std::string isoline = vis_folder + "/isoline/";
      FeatureVector fea =
        inf_feature(solve_context.width, solve_context.height);
      std::vector<std::string> feature_set_out;
      std::vector<double> fea_vec =
        (solve_context.stage == SolveContext::SolveStage::Initial)
          ? pick_features(fea, solve_context.feature_set, feature_set_out)
          : pick_features(fea, solve_context.sec_feature_set, feature_set_out);
      print_csv_row(feature_set_out, fea_vec, -1,
                    (!solve_context.typical_feature_set.empty()) ? 1 : -1,
                    curr_c1.strokes.front().stroke_ind,
                    new_c + csv_cluster_offset, true, true, fit_output, isoline,
                    csv_str);
    }
  }

  // Find the larger split location
  int key_xsec_ind = larger_split_xsec_ind(cluster, seed_sid[0], seed_sid[1]);
  assert(key_xsec_ind >= 0);

  // Iteratively reassign other strokes
  int last_updated_cid = -1;
  while (processed_sids.size() < cluster.strokes.size()) {
    std::vector<size_t> unprocessed_sids;
    unprocessed_sids.reserve(sid2s.size());

    // Get unprocessed strokes
    for (auto const &ss : sid2s) {
      if (processed_sids.count(ss.first))
        continue;

      // Check if this stroke overlaps with both seed clusters
      if (is_overlapping_cluster(cluster, ss.first, clusters[0],
                                 key_xsec_ind) &&
          is_overlapping_cluster(cluster, ss.first, clusters[1], key_xsec_ind))
        unprocessed_sids.emplace_back(ss.first);
    }

    // No remaining stroke overlaps with both seeds
    if (unprocessed_sids.empty()) {
      for (auto const &ss : sid2s) {
        if (processed_sids.count(ss.first))
          continue;

        // Check if this stroke overlaps with one seed cluster
        if (is_overlapping_cluster(cluster, ss.first, clusters[0]) ||
            is_overlapping_cluster(cluster, ss.first, clusters[1]))
          unprocessed_sids.emplace_back(ss.first);
      }
    }

    if (unprocessed_sids.empty()) {
      for (auto const &ss : sid2s) {
        if (processed_sids.count(ss.first))
          continue;
        bool overlaps_any = false;
        for (auto const &cc : clusters) {
          if (is_overlapping_cluster(cluster, ss.first, cc.second)) {
            overlaps_any = true;
            break;
          }
        }
        if (overlaps_any)
          unprocessed_sids.emplace_back(ss.first);
      }
    }

    // For each one, predict wrt each existing cluster
    auto merge_to_cluster = [&](const PredictionInfo &merge_prediction_info) {
      Cluster &target_cluster = clusters[merge_prediction_info.target_number];
      target_cluster.strokes.emplace_back(
        sid2s[merge_prediction_info.base_number]);
      target_cluster.original_input_strokes.emplace_back(
        sid2orig_s[merge_prediction_info.base_number]);
      std::sort(target_cluster.strokes.begin(), target_cluster.strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });
      std::sort(target_cluster.original_input_strokes.begin(),
                target_cluster.original_input_strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });

      processed_sids.emplace(merge_prediction_info.base_number);

      if (merge_prediction_info.target_number <= 1)
        last_updated_cid = merge_prediction_info.target_number;

      if (!vis_folder.empty()) {
        double typical_score = -1;
        std::string fit_output = vis_folder + "/fit/" +
                                 merge_prediction_info.hash_str_comb +
                                 "_fit.svg";
        std::string isoline = vis_folder + "/isoline/" +
                              merge_prediction_info.hash_str_comb +
                              "_isolines.svg";
        print_csv_row(
          merge_prediction_info.feature_set, merge_prediction_info.feature_vec,
          merge_prediction_info.score, typical_score,
          merge_prediction_info.base_number,
          merge_prediction_info.target_number + csv_cluster_offset,
          merge_prediction_info.is_merging, merge_prediction_info.decision,
          fit_output, isoline, csv_str);
      }
    };

    std::vector<PredictionInfo> prediction_infos, prediction_infos_smaller;
    std::vector<PredictionInfo> prediction_epsilon_infos;
    for (auto sid : unprocessed_sids) {
      Cluster curr_c;
      curr_c.strokes.emplace_back(sid2s[sid]);
      curr_c.original_input_strokes.emplace_back(sid2orig_s[sid]);

      for (auto cc : clusters) {
        double prediction = -1;
        PredictionInfo prediction_info;
        bool seed_picking = false;
        bool split_within = true;
        compare_cluster_single(
          cc.second, curr_c, cluster, solve_context.width, solve_context.height,
          solve_context.norm_thickness, prediction_info, prediction, vis_folder,
          seed_picking, split_within);
        if (prediction >= 0 && !prediction_info.feature_vec.empty()) {
          prediction_infos.emplace_back(prediction_info);
          prediction_infos.back().base_number =
            curr_c.strokes.front().stroke_ind;
          prediction_infos.back().target_number = cc.first;
          prediction_infos.back().is_merging = false;
        }
      }
    }

    // Find the pair with highest probability
    prediction_epsilon_infos = prediction_infos;
    if (!prediction_infos.empty()) {
      size_t pick_i = 0;
      for (; pick_i < prediction_infos.begin()->feature_set.size(); ++pick_i) {
        if (prediction_infos.begin()->feature_set[pick_i] ==
            tier_breaker_feature)
          break;
      }
      std::sort(prediction_infos.begin(), prediction_infos.end(),
                [&pick_i, &last_updated_cid](const PredictionInfo &a,
                                             const PredictionInfo &b) {
                  double score1 = a.score;
                  double score2 = b.score;

                  double fea1 = a.feature_vec[pick_i];
                  double fea2 = b.feature_vec[pick_i];

                  return score1 > score2 ||
                         (score1 == score2 &&
                          a.feature_vec[pick_i] < b.feature_vec[pick_i]);
                });

      if (solve_context.tier_epsilon > 0) {
        prediction_epsilon_infos.clear();
        for (auto const &info : prediction_infos) {
          if (std::abs(prediction_infos.front().score - info.score) <
              solve_context.tier_epsilon) {
            prediction_epsilon_infos.emplace_back(info);
          } else {
            prediction_infos_smaller.emplace_back(info);
          }
        }
        std::sort(prediction_epsilon_infos.begin(),
                  prediction_epsilon_infos.end(),
                  [&pick_i, &last_updated_cid](const PredictionInfo &a,
                                               const PredictionInfo &b) {
                    double score1 = a.score;
                    double score2 = b.score;

                    double fea1 = a.feature_vec[pick_i];
                    double fea2 = b.feature_vec[pick_i];

                    return a.feature_vec[pick_i] < b.feature_vec[pick_i] ||
                           (score1 > score2 &&
                            a.feature_vec[pick_i] == b.feature_vec[pick_i]);
                  });
        prediction_epsilon_infos.insert(prediction_epsilon_infos.end(),
                                        prediction_infos_smaller.begin(),
                                        prediction_infos_smaller.end());
      }
    }

    if (!prediction_epsilon_infos.empty()) {
      PredictionInfo merge_prediction_info = prediction_epsilon_infos.front();
      for (const auto &info : prediction_epsilon_infos) {
        if (info.score >= solve_context.cutoff) {
          merge_prediction_info = info;
          break;
        }
      }
      if (merge_prediction_info.score >= solve_context.cutoff) {
        logger().info("Merging {} into c{} with prob={}:",
                      merge_prediction_info.base_number,
                      merge_prediction_info.target_number,
                      merge_prediction_info.score);
        for (auto const &s :
             clusters[merge_prediction_info.target_number].strokes) {
          logger().info("\ts{}", s.stroke_ind);
        }
        // Merge
        merge_to_cluster(merge_prediction_info);
      } else {
        logger().info("Cannot merging {} into c{} with prob={}:",
                      merge_prediction_info.base_number,
                      merge_prediction_info.target_number,
                      merge_prediction_info.score);
        for (auto const &s :
             clusters[merge_prediction_info.target_number].strokes) {
          logger().info("\ts{}", s.stroke_ind);
        }
        // New cluster
        size_t new_cluster_sid = merge_prediction_info.base_number;
        Cluster curr_c;
        curr_c.strokes.emplace_back(sid2s[new_cluster_sid]);
        curr_c.original_input_strokes.emplace_back(sid2orig_s[new_cluster_sid]);

        size_t new_c = (clusters.empty()) ? 0 : (clusters.rbegin()->first + 1);
        clusters[new_c] = curr_c;

        processed_sids.emplace(new_cluster_sid);

        if (!vis_folder.empty()) {
          for (auto const &info : prediction_epsilon_infos) {
            std::string fit_output =
              vis_folder + "/fit/" + info.hash_str_comb + "_fit.svg";
            std::string isoline =
              vis_folder + "/isoline/" + info.hash_str_comb + "_isolines.svg";
            print_csv_row(
              info.feature_set, info.feature_vec, info.score, -1,
              info.base_number, info.target_number + csv_cluster_offset,
              info.is_merging, info.decision, fit_output, isoline, csv_str);
          }

          std::string fit_output = vis_folder + "/fit/";
          std::string isoline = vis_folder + "/isoline/";
          FeatureVector fea =
            inf_feature(solve_context.width, solve_context.height);
          std::vector<std::string> feature_set_out;
          std::vector<double> fea_vec =
            (solve_context.stage == SolveContext::SolveStage::Initial)
              ? pick_features(fea, solve_context.feature_set, feature_set_out)
              : pick_features(fea, solve_context.sec_feature_set,
                              feature_set_out);
          print_csv_row(feature_set_out, fea_vec, -1,
                        (!solve_context.typical_feature_set.empty()) ? 1 : -1,
                        curr_c.strokes.front().stroke_ind,
                        new_c + csv_cluster_offset, true, true, fit_output,
                        isoline, csv_str);
        }
      }
    } else {
      logger().warn("Cannot find any overlapping stroke");

      // Pick one stroke closest to any existing one and add it to that cluster
      std::map<double, std::pair<size_t, size_t>> sort_u_dist;
      for (auto const &ss : sid2s) {
        if (processed_sids.count(ss.first))
          continue;

        size_t new_cluster_sid = ss.first;
        Cluster curr_c;
        curr_c.strokes.emplace_back(sid2s[new_cluster_sid]);
        curr_c.original_input_strokes.emplace_back(sid2orig_s[new_cluster_sid]);
        for (auto const &cc : clusters) {
          double u_dist = u_distance(curr_c, cc.second, cluster);
          sort_u_dist[u_dist] = std::make_pair(ss.first, cc.first);
        }
      }

      // Merge
      FeatureVector fea =
        inf_feature(solve_context.width, solve_context.height);
      std::vector<std::string> feature_set_out;
      std::vector<double> fea_vec =
        (solve_context.stage == SolveContext::SolveStage::Initial)
          ? pick_features(fea, solve_context.feature_set, feature_set_out)
          : pick_features(fea, solve_context.sec_feature_set, feature_set_out);
      PredictionInfo merge_prediction_info;
      merge_prediction_info.base_number = sort_u_dist.begin()->second.first;
      merge_prediction_info.target_number = sort_u_dist.begin()->second.second;
      merge_prediction_info.is_merging = false;
      merge_prediction_info.feature_set =
        (solve_context.stage == SolveContext::SolveStage::Initial)
          ? solve_context.feature_set
          : solve_context.sec_feature_set;
      merge_prediction_info.feature_vec = fea_vec;
      std::fill(merge_prediction_info.feature_vec.begin(),
                merge_prediction_info.feature_vec.end(), 0);
      merge_prediction_info.score = 1;
      merge_prediction_info.decision = true;
      merge_prediction_info.hash_str_comb = "u_dist";
      merge_to_cluster(merge_prediction_info);
    }
  }
}

void further_split(size_t in_cid, Cluster &cluster,
                   std::vector<Cluster> &clusters, std::string vis_folder,
                   std::string &csv_str, size_t csv_cluster_offset,
                   bool filtering) {
  if (cluster.strokes.size() < 2) {
    clusters.emplace_back(cluster);
    return;
  }

  bool reorient = false;
  FeatureVector fea(solve_context.width, solve_context.height);
  fea.compute_typical(cluster, solve_context.norm_thickness,
                      &solve_context.matching_pair, reorient, "");
  std::vector<std::string> max_width_str;
  auto max_widths = pick_features(
    fea, std::vector<std::string>{"average_distance"}, max_width_str);
  double max_width = max_widths[0];

  max_widths =
    pick_features(fea, std::vector<std::string>{"max_distance"}, max_width_str);
  double max_gap = max_widths[0];
  max_widths =
    pick_features(fea, std::vector<std::string>{"max_distance2width_combined"},
                  max_width_str);
  double max_gap_ratio = max_widths[0];

  double max_other_width = -1;
  for (const auto &c_fea : solve_context.cluster_features) {
    if (c_fea.first != in_cid) {
      std::vector<std::string> max_width_str;
      auto max_widths = pick_features(
        c_fea.second, std::vector<std::string>{"average_distance"},
        max_width_str);
      double max_width = max_widths[0];
      max_other_width = std::max(max_other_width, max_width);
    }
  }

  if (max_gap > 5)
    logger().info("Further: {} - {}", max_gap, max_gap_ratio);

  if (filtering && max_width < width_split_threshold * max_other_width &&
      !(max_gap > 5 && max_gap_ratio > 0.7)) {
    clusters.emplace_back(cluster);
    return;
  }

  // 1. Find xsec candidates to generate seeds
  std::vector<std::pair<size_t, size_t>> seed_sid;

  size_t c_tmp = cluster.strokes.front().cluster_ind;
  cluster.strokes.front().cluster_ind = in_cid;
  find_seed_xsecs_probability(cluster, seed_sid);
  cluster.strokes.front().cluster_ind = c_tmp;

  if (seed_sid.empty()) {
    clusters.emplace_back(cluster);
    return;
  }

  // 2. Start splitting from a seed. For now just split the first one (the
  // widest one).
  // 3. Iteratively reassign within group
  std::string new_csv_str = (!csv_str.empty()) ? "@" : "";
  std::map<size_t, Cluster> split_clusters;
  split_from_seed(cluster, seed_sid.front(), split_clusters, vis_folder,
                  new_csv_str, csv_cluster_offset);

  for (auto const &cc : split_clusters) {
    clusters.emplace_back(cc.second);
  }
}

bool disable_further_filter(Cluster &cluster) {
  bool reorient = false;
  FeatureVector fea(solve_context.width, solve_context.height);
  fea.compute_typical(cluster, solve_context.norm_thickness,
                      &solve_context.matching_pair, reorient, "");
  std::vector<std::string> max_width_str;
  auto max_widths = pick_features(
    fea, std::vector<std::string>{"average_distance"}, max_width_str);
  double max_width = max_widths[0];

  max_widths =
    pick_features(fea, std::vector<std::string>{"max_distance"}, max_width_str);
  double max_gap = max_widths[0];
  max_widths =
    pick_features(fea, std::vector<std::string>{"max_distance2width_combined"},
                  max_width_str);
  double max_gap_ratio = max_widths[0];
  logger().info("disable_further_filter: {} - {}", max_gap, max_gap_ratio);

  return (max_gap > 3.5 && max_gap_ratio > 0.6);
}

void solve_end_split(const Input &input, std::map<size_t, Cluster> &clusters,
                     std::map<size_t, Cluster> &intermediate_clusters,
                     std::string intermediate_folder, std::string vis_folder) {
  solve_context.stage = SolveContext::SolveStage::Splitting;

  std::map<size_t, std::map<size_t, Cluster>> to_merge_clusters;
  std::map<size_t, size_t> round_c_offsets;
  std::vector<Cluster> result_clusters;
  std::string csv_str = "";
  std::string merge_csv_str = "@";
  for (auto &cid_cluster : clusters) {
    auto &cluster = cid_cluster.second;

    if (cluster.strokes.size() <= 1) {
      size_t new_c = intermediate_clusters.size();
      intermediate_clusters[new_c] = cluster;
      intermediate_clusters[new_c].fit.cluster_idx = new_c;
      result_clusters.emplace_back(cluster);
      continue;
    }

    // 1. Find xsec candidates to generate seeds
    std::vector<std::pair<size_t, size_t>> seed_sid;

    size_t c_tmp = cluster.strokes.front().cluster_ind;
    cluster.strokes.front().cluster_ind = cid_cluster.first;
    find_seed_xsecs_probability(cluster, seed_sid);
    cluster.strokes.front().cluster_ind = c_tmp;

    if (seed_sid.empty()) {
      size_t new_c = intermediate_clusters.size();
      intermediate_clusters[new_c] = cluster;
      intermediate_clusters[new_c].fit.cluster_idx = new_c;
      result_clusters.emplace_back(cluster);
      continue;
    }

    // 2. Start splitting from a seed. For now just split the first one (the
    // widest one).
    // 3. Iteratively reassign within group
    std::string new_csv_str = (!csv_str.empty()) ? "@" : "";
    size_t round_c_offset = intermediate_clusters.size();
    std::map<size_t, Cluster> split_clusters_init;
    bool further_filter = !disable_further_filter(cluster);
    further_filter = true;
    split_from_seed(cluster, seed_sid.front(), split_clusters_init, vis_folder,
                    new_csv_str, round_c_offset);

    std::map<size_t, Cluster> split_clusters;
    for (auto &cc : split_clusters_init) {
      std::vector<Cluster> further_clusters;
      auto assign_cluster_index = [](const Cluster &ref_cluster, size_t cid,
                                     Cluster &ret_cluster) {
        std::map<size_t, size_t> s2c;
        for (const auto &s : ref_cluster.strokes) {
          s2c[s.stroke_ind] = cid;
        }

        for (auto &s : ret_cluster.strokes) {
          if (s2c.count(s.stroke_ind))
            s.cluster_ind = s2c[s.stroke_ind];
        }

        for (auto &xsec : ret_cluster.xsecs) {
          for (auto &p : xsec.points) {
            if (s2c.count(p.stroke_ind))
              p.cluster_ind = s2c[p.stroke_ind];
          }
        }
      };
      // Assign an invalid flag first
      cc.second.xsecs = cluster.xsecs;
      assign_cluster_index(cluster, std::numeric_limits<size_t>::max(),
                           cc.second);
      assign_cluster_index(cc.second, 0, cc.second);
      cc.second.xsecs = filter_xsecs(cc.second.xsecs);
      cc.second.obj_term_values = cluster.obj_term_values;
      further_split(cid_cluster.first, cc.second, further_clusters, vis_folder,
                    new_csv_str, round_c_offset, further_filter);
      if (!further_filter && further_clusters.size() > 1) {
        std::map<size_t, Cluster> further_merge_clusters;
        for (auto const &c : further_clusters)
          further_merge_clusters[further_merge_clusters.size()] = c;
        std::string merge_csv_str_tmp;
        double tmp_cutoff = solve_context.cutoff;
        solve_context.cutoff = 0.45;
        solve_end_merge(input, further_merge_clusters, vis_folder,
                        merge_csv_str_tmp, round_c_offset, &cc.second);
        further_clusters.clear();
        solve_context.cutoff = tmp_cutoff;
        for (auto const &fcc : further_merge_clusters)
          further_clusters.emplace_back(fcc.second);
      }
      for (auto const &c : further_clusters)
        split_clusters[split_clusters.size()] = c;
    }

    if (!csv_str.empty())
      new_csv_str = new_csv_str.substr(1);
    csv_str += new_csv_str;

    round_c_offsets[cid_cluster.first] = round_c_offset;

    for (auto &cc : split_clusters) {
      size_t new_c = intermediate_clusters.size();
      intermediate_clusters[new_c] = cc.second;
      intermediate_clusters[new_c].fit.cluster_idx = new_c;
    }
    to_merge_clusters[cid_cluster.first] = split_clusters;
  }

  // Optional: Junction separation
  if (solve_context.to_junc_sep) {
    std::vector<std::pair<size_t, size_t>> sep_cid;
    separable_clusters(input, clusters, intermediate_clusters, sep_cid,
                       intermediate_folder);
    for (auto const &cc : sep_cid) {
      solve_context.separation_sid.emplace_back(std::make_pair(
        std::min(intermediate_clusters[cc.first].strokes.front().stroke_ind,
                 intermediate_clusters[cc.second].strokes.front().stroke_ind),
        std::max(intermediate_clusters[cc.first].strokes.front().stroke_ind,
                 intermediate_clusters[cc.second].strokes.front().stroke_ind)));
    }
  }

  logger().info("Split-merge");
  for (auto &c_split_clusters : to_merge_clusters) {
    double tmp_cutoff = solve_context.cutoff;

    std::map<size_t, Cluster> split_clusters = c_split_clusters.second;
    auto &cluster = clusters[c_split_clusters.first];
    size_t round_c_offset = round_c_offsets[c_split_clusters.first];

    // 4. Merge resulting clusters again (including the two seeds)
    solve_end_merge(input, split_clusters, vis_folder, merge_csv_str,
                    round_c_offset, &cluster);

    // 5. Merge resulting clusters again without giving the init cluster
    // parameterization
    if (split_clusters.size() > 1) {
      for (auto &cc : split_clusters) {
        auto &cluster = cc.second;
        for (auto &s : cluster.strokes) {
          s.cluster_ind = cc.first;
        }

        // Recompute parameterization
        std::string hash_str = get_hash_str(cluster);
        if (solve_context.merged_cluster_cache.count(hash_str))
          cluster = solve_context.merged_cluster_cache[hash_str];
        else {
          bool test_time = true;
          bool to_orient = true; // Assume strokes are correctly oriented in
                                 // the initial cluster
          fit_cluster(solve_context.norm_thickness, solve_context.width,
                      solve_context.height, cluster, to_orient, test_time,
                      &solve_context.matching_pair);
          solve_context.merged_cluster_cache[hash_str] = cluster;
        }
      }
      solve_end_merge(input, split_clusters, vis_folder, merge_csv_str, 0,
                      nullptr);
    }
    solve_context.cutoff = tmp_cutoff;

    c_split_clusters.second = split_clusters;
  }

  /////////////////
  solve_context.stage = SolveContext::SolveStage::Split_Merging;
  for (auto const &c_split_clusters : to_merge_clusters) {
    double tmp_cutoff = solve_context.cutoff;

    std::map<size_t, Cluster> split_clusters = c_split_clusters.second;
    auto &cluster = clusters[c_split_clusters.first];
    size_t round_c_offset = round_c_offsets[c_split_clusters.first];

    // 4. Merge resulting clusters again (including the two seeds)
    if (split_clusters.size() > 1) {
      solve_end_merge(input, split_clusters, vis_folder, merge_csv_str,
                      round_c_offset, &cluster);

      // 5. Merge resulting clusters again without giving the init cluster
      // parameterization
      if (split_clusters.size() > 1) {
        for (auto &cc : split_clusters) {
          auto &cluster = cc.second;
          for (auto &s : cluster.strokes) {
            s.cluster_ind = cc.first;
          }

          // Recompute parameterization
          std::string hash_str = get_hash_str(cluster);
          if (solve_context.merged_cluster_cache.count(hash_str))
            cluster = solve_context.merged_cluster_cache[hash_str];
          else {
            bool test_time = true;
            bool to_orient = true; // Assume strokes are correctly oriented in
                                   // the initial cluster
            fit_cluster(solve_context.norm_thickness, solve_context.width,
                        solve_context.height, cluster, to_orient, test_time,
                        &solve_context.matching_pair);
            solve_context.merged_cluster_cache[hash_str] = cluster;
          }
        }
        solve_end_merge(input, split_clusters, vis_folder, merge_csv_str, 0,
                        nullptr);
      }
    }

    solve_context.cutoff = tmp_cutoff;

    // Save
    for (auto &cc : split_clusters) {
      result_clusters.emplace_back(cc.second);
    }
  }
  solve_context.stage = SolveContext::SolveStage::Splitting;
  /////////////////

  logger().info("Split-merge ends");

  merge_csv_str = merge_csv_str.substr(1);
  merge_csv_str = csv_str + merge_csv_str;

  // Read out and recompute parameterization
  clusters.clear();
  for (size_t i = 0; i < result_clusters.size(); ++i) {
    clusters[i] = result_clusters[i];
    for (auto &s : clusters[i].strokes) {
      s.cluster_ind = i;
    }

    // Recompute parameterization
    std::string hash_str = get_hash_str(clusters[i]);
    if (solve_context.merged_cluster_cache.count(hash_str))
      clusters[i] = solve_context.merged_cluster_cache[hash_str];
    else {
      bool test_time = true;
      bool to_orient =
        true; // Assume strokes are correctly oriented in the initial cluster
      fit_cluster(solve_context.norm_thickness, solve_context.width,
                  solve_context.height, clusters[i], to_orient, test_time,
                  &solve_context.matching_pair);
      solve_context.merged_cluster_cache[hash_str] = clusters[i];
    }
  }

  for (auto &cc : intermediate_clusters) {
    for (auto &s : cc.second.strokes) {
      s.cluster_ind = cc.first;
    }

    // Recompute parameterization
    std::string hash_str = get_hash_str(cc.second);
    if (solve_context.merged_cluster_cache.count(hash_str))
      cc.second = solve_context.merged_cluster_cache[hash_str];
    else {
      bool test_time = true;
      bool to_orient =
        true; // Assume strokes are correctly oriented in the initial cluster
      fit_cluster(solve_context.norm_thickness, solve_context.width,
                  solve_context.height, cc.second, to_orient, test_time,
                  &solve_context.matching_pair);
      solve_context.merged_cluster_cache[hash_str] = cc.second;
    }
  }

  if (!intermediate_folder.empty()) {
    std::ofstream csv_ofs(intermediate_folder + "/features_split.csv");
    csv_ofs.write(csv_str.c_str(), csv_str.size());
    csv_ofs.close();
    std::ofstream merge_csv_ofs(intermediate_folder + "/features_merge.csv");
    merge_csv_ofs.write(merge_csv_str.c_str(), merge_csv_str.size());
    merge_csv_ofs.close();
  }
}
