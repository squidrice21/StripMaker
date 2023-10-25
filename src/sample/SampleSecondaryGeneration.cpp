#include "SampleSecondaryGeneration.h"
#include "../../stroke_strip_src/Cluster.h"
#include "../../stroke_strip_src/SketchInfo.h"
#include "../Logger.h"
#include "../Util.h"
#include "../feature/Features.h"
#include "../feature/Filters.h"
#include "../solve/SolveLazy.h"
#include "../solve/SolveUtil.h"

#include <algorithm>
#include <limits>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

size_t sample_seed_xsec_count = 1;

bool is_cut_stroke(
  const Cluster::Stroke &s,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    &matching_pair) {
  for (auto const &p : matching_pair) {
    if (p.first.first == s.stroke_ind || p.second.first == s.stroke_ind)
      return true;
  }
  return false;
}

void find_seed_xsecs_distance(
  Cluster &cluster, std::vector<std::pair<size_t, size_t>> &seeds,
  size_t seed_xsec_count,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    &matching_pair) {
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

  std::vector<std::pair<double, std::pair<size_t, size_t>>> sort_seeds;
  int pick_i = -1;
  for (size_t i = 0; i + 1 < cluster.strokes.size(); ++i) {
    for (size_t j = i + 1; j < cluster.strokes.size(); ++j) {
      size_t sid1 = cluster.strokes[i].stroke_ind;
      size_t sid2 = cluster.strokes[j].stroke_ind;
      // Don't use cut stroke as seed stroke since if it were cut due to high
      // curvature, it would be hard to merge it back.
      if (is_cut_stroke(cluster.strokes[i], matching_pair) ||
          is_cut_stroke(cluster.strokes[j], matching_pair))
        continue;
      if (!is_overlapping_stroke(cluster, sid1, sid2))
        continue;
      Cluster curr_c;
      curr_c.strokes.emplace_back(sid2s[sid1]);
      curr_c.original_input_strokes.emplace_back(sid2orig_s[sid1]);

      Cluster cmp_c;
      cmp_c.strokes.emplace_back(sid2s[sid2]);
      cmp_c.original_input_strokes.emplace_back(sid2orig_s[sid2]);

      double norm_thickness = 1;
      FeatureVector fea(-1, -1);
      bool reorient = false;
      fea.compute_distance(cmp_c, curr_c, norm_thickness, cluster, nullptr,
                           reorient, "");
      std::vector<std::string> avg_dist_str;
      auto avg_dists = pick_features(
        fea, std::vector<std::string>{"average_distance"}, avg_dist_str);
      double avg_dist = avg_dists[0];

      if (avg_dist < std::numeric_limits<double>::infinity()) {
        sort_seeds.emplace_back(
          std::make_pair(avg_dist, std::make_pair(std::min(sid1, sid2),
                                                  std::max(sid1, sid2))));
      }
    }
  }

  std::sort(sort_seeds.begin(), sort_seeds.end(),
            [](const std::pair<double, std::pair<size_t, size_t>> &a,
               const std::pair<double, std::pair<size_t, size_t>> &b) {
              double dist_a = a.first;
              double dist_b = b.first;
              return dist_a > dist_b;
            });

  for (size_t i = 0; i < seed_xsec_count && i < sort_seeds.size(); ++i) {
    seeds.emplace_back(sort_seeds[i].second);
  }
}

void split_cluster(Cluster &cluster, const std::pair<size_t, size_t> &seed,
                   Sample &sample, size_t left_seed_count,
                   std::vector<Sample> *intermediate_samples_ptr = nullptr) {
  std::map<size_t, Cluster> clusters;
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
  }

  // Iteratively reassign other strokes
  while (processed_sids.size() < cluster.strokes.size()) {
    std::vector<size_t> unprocessed_sids;
    unprocessed_sids.reserve(sid2s.size());

    // Get unprocessed strokes
    for (auto const &ss : sid2s) {
      if (processed_sids.count(ss.first))
        continue;

      // Check if this stroke overlaps the left cluster if we haven't reached
      // desired size
      if (is_overlapping_cluster(cluster, ss.first, clusters[0]) &&
          clusters[0].strokes.size() < left_seed_count)
        unprocessed_sids.emplace_back(ss.first);
    }

    // Reached the desired size of left cluster
    if (unprocessed_sids.empty()) {
      for (auto const &ss : sid2s) {
        if (processed_sids.count(ss.first))
          continue;

        // Check if this stroke overlaps with the right cluster
        if (is_overlapping_cluster(cluster, ss.first, clusters[1]))
          unprocessed_sids.emplace_back(ss.first);
      }
    }

    // For each one, predict wrt each existing cluster
    auto merge_to_cluster = [&](size_t base_number, size_t target_number) {
      Cluster &target_cluster = clusters[target_number];

      if (intermediate_samples_ptr) {
        std::vector<int> pre_cluster_stroke_ind;
        for (auto const &s : target_cluster.strokes) {
          pre_cluster_stroke_ind.emplace_back(s.stroke_ind);
        }
        auto sample = std::make_pair(std::vector<int>{(int)base_number},
                                     pre_cluster_stroke_ind);
        intermediate_samples_ptr->emplace_back(sample);
      }

      target_cluster.strokes.emplace_back(sid2s[base_number]);
      target_cluster.original_input_strokes.emplace_back(
        sid2orig_s[base_number]);
      std::sort(target_cluster.strokes.begin(), target_cluster.strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });
      std::sort(target_cluster.original_input_strokes.begin(),
                target_cluster.original_input_strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });

      processed_sids.emplace(base_number);
    };

    std::vector<std::pair<double, std::pair<size_t, size_t>>> sort_strokes;
    for (auto sid : unprocessed_sids) {
      Cluster curr_c;
      curr_c.strokes.emplace_back(sid2s[sid]);
      curr_c.original_input_strokes.emplace_back(sid2orig_s[sid]);

      for (auto &cc : clusters) {
        if (clusters[0].strokes.size() < left_seed_count && cc.first != 0)
          continue;
        if (clusters[0].strokes.size() >= left_seed_count && cc.first == 0)
          continue;

        double norm_thickness = 1;
        Cluster merged_cluster;
        FeatureVector fea(-1, -1);
        bool reorient = false;
        fea.compute_distance(cc.second, curr_c, norm_thickness, cluster,
                             nullptr, reorient, "");
        std::vector<std::string> avg_dist_str;
        auto avg_dists = pick_features(
          fea, std::vector<std::string>{"average_distance"}, avg_dist_str);
        double avg_dist = avg_dists[0];

        if (avg_dist < std::numeric_limits<double>::infinity()) {
          sort_strokes.emplace_back(
            std::make_pair(avg_dist, std::make_pair(sid, cc.first)));
        }
      }
    }

    // Find the pair with highest probability
    if (!sort_strokes.empty()) {
      std::sort(sort_strokes.begin(), sort_strokes.end(),
                [](const std::pair<double, std::pair<size_t, size_t>> &a,
                   const std::pair<double, std::pair<size_t, size_t>> &b) {
                  return a.first < b.first;
                });
    }

    if (!sort_strokes.empty()) {
      // Merge
      merge_to_cluster(sort_strokes.front().second.first,
                       sort_strokes.front().second.second);
    } else {
      logger().warn("Cannot find any overlapping stroke");

      // Pick one stroke closest to any existing one and add it to that cluster
      bool updated;
      do {
        updated = false;
        for (auto const &ss : sid2s) {
          if (processed_sids.count(ss.first))
            continue;

          if (is_overlapping_cluster(cluster, ss.first, clusters[0])) {
            merge_to_cluster(ss.first, 0);
            updated = true;
          } else if (is_overlapping_cluster(cluster, ss.first, clusters[1])) {
            merge_to_cluster(ss.first, 1);
            updated = true;
          }
        }
      } while (processed_sids.size() < cluster.strokes.size() && updated);

      while (processed_sids.size() < cluster.strokes.size()) {
        std::map<double, std::pair<size_t, size_t>> sort_u_dist;
        for (auto const &ss : sid2s) {
          if (processed_sids.count(ss.first))
            continue;

          size_t new_cluster_sid = ss.first;
          Cluster curr_c;
          curr_c.strokes.emplace_back(sid2s[new_cluster_sid]);
          curr_c.original_input_strokes.emplace_back(
            sid2orig_s[new_cluster_sid]);
          for (auto const &cc : clusters) {
            double u_dist = u_distance(curr_c, cc.second, cluster);
            sort_u_dist[u_dist] = std::make_pair(ss.first, cc.first);
          }
        }

        merge_to_cluster(sort_u_dist.begin()->second.first,
                         sort_u_dist.begin()->second.second);
      }
    }
  }

  std::vector<int> pre_cluster_stroke_ind;
  std::vector<int> post_cluster_stroke_ind;
  for (auto const &s : clusters[0].strokes) {
    pre_cluster_stroke_ind.emplace_back(s.stroke_ind);
  }
  for (auto const &s : clusters[1].strokes) {
    post_cluster_stroke_ind.emplace_back(s.stroke_ind);
  }
  sample = std::make_pair(pre_cluster_stroke_ind, post_cluster_stroke_ind);
}

void generate_secondary_stroke_positive_samples(
  const Input &input, std::vector<Sample> &pos_cluster_samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  std::set<size_t> selected_strokes;
  std::set<std::string> seen_subclusters;
  auto clusters = input.clusters;
  for (auto &cluster : clusters) {
    // The parameterization fails
    if (cluster.second.fit.centerline.empty() ||
        cluster.second.strokes.size() < 2)
      continue;

    // 1. Find a seed
    std::vector<std::pair<size_t, size_t>> seeds;
    const std::vector<
      std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
      matching_pair =
        (matching_pair_ptr)
          ? *matching_pair_ptr
          : std::vector<
              std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>();
    find_seed_xsecs_distance(cluster.second, seeds, sample_seed_xsec_count,
                             matching_pair);

    // These may happen if we mostly have cut strokes in this cluster
    if (seeds.empty())
      continue;

    // 2. Split
    for (size_t i = 0; i < cluster.second.strokes.size(); ++i) {
      Sample sample;
      size_t left_seed_count = i + 1;
      std::vector<Sample> intermediate_samples;
      split_cluster(cluster.second, seeds.front(), sample, left_seed_count,
                    &intermediate_samples);

      std::string hash_str1 = example_hash(sample.first);
      std::string hash_str2 = example_hash(sample.second);
      seen_subclusters.emplace(hash_str1 + '-' + hash_str2);
      seen_subclusters.emplace(hash_str2 + '-' + hash_str1);

      for (auto const &s : intermediate_samples) {
        std::string hash_str1 = example_hash(s.first);
        std::string hash_str2 = example_hash(s.second);
        if (!seen_subclusters.count(hash_str1 + "-" + hash_str2) &&
            !seen_subclusters.count(hash_str2 + "-" + hash_str1)) {
          pos_cluster_samples.emplace_back(s);
          seen_subclusters.emplace(hash_str1 + '-' + hash_str2);
          seen_subclusters.emplace(hash_str2 + '-' + hash_str1);
        }
      }
    }
  }
}

void generate_secondary_cluster_positive_samples(
  const Input &input, std::vector<Sample> &pos_cluster_samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  std::set<size_t> selected_strokes;
  std::set<std::string> seen_subclusters;
  auto clusters = input.clusters;
  for (auto &cluster : clusters) {
    // The parameterization fails
    if (cluster.second.fit.centerline.empty() ||
        cluster.second.strokes.size() < 2)
      continue;

    // 1. Find a seed
    std::vector<std::pair<size_t, size_t>> seeds;
    const std::vector<
      std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
      matching_pair =
        (matching_pair_ptr)
          ? *matching_pair_ptr
          : std::vector<
              std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>();
    find_seed_xsecs_distance(cluster.second, seeds, sample_seed_xsec_count,
                             matching_pair);

    // These may happen if we mostly have cut strokes in this cluster
    if (seeds.empty())
      continue;

    // 2. Split
    for (size_t i = 0; i < cluster.second.strokes.size(); ++i) {
      Sample sample;
      size_t left_seed_count = i + 1;
      split_cluster(cluster.second, seeds.front(), sample, left_seed_count);

      std::string hash_str1 = example_hash(sample.first);
      std::string hash_str2 = example_hash(sample.second);
      if (!seen_subclusters.count(hash_str1) &&
          !seen_subclusters.count(hash_str2)) {
        pos_cluster_samples.emplace_back(sample);
        seen_subclusters.emplace(hash_str1);
        seen_subclusters.emplace(hash_str2);
      }
    }
  }
}

void generate_secondary_cluster_negative_samples(
  const Input &input, std::vector<Sample> &neg_cluster_samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  std::set<size_t> selected_strokes;
  for (const auto &cluster : input.clusters) {
    // The parameterization fails
    if (cluster.second.fit.centerline.empty())
      continue;

    for (const auto &cluster_sub : input.clusters) {
      // The parameterization fails
      if (cluster_sub.second.fit.centerline.empty())
        continue;

      if (cluster.first >= cluster_sub.first)
        continue;
      std::vector<int> cluster_stroke_ind;
      std::vector<int> sub_cluster_stroke_ind;
      for (const auto &stroke : cluster.second.strokes)
        cluster_stroke_ind.emplace_back(stroke.stroke_ind);
      for (const auto &stroke : cluster_sub.second.strokes)
        sub_cluster_stroke_ind.emplace_back(stroke.stroke_ind);

      // 0. Trivial distance filter
      if (filter_sample(cluster.second, cluster_stroke_ind, cluster_sub.second,
                        sub_cluster_stroke_ind))
        continue;

      // No need to check, the two clusters are continuous by construction
      // Only check if the two clusters overlap
      Cluster cluster1, cluster2, merged_cluster;

      selected_strokes.clear();
      selected_strokes.insert(cluster_stroke_ind.begin(),
                              cluster_stroke_ind.end());
      reuse_cluster(cluster.second, cluster1);

      selected_strokes.clear();
      selected_strokes.insert(sub_cluster_stroke_ind.begin(),
                              sub_cluster_stroke_ind.end());
      reuse_cluster(cluster_sub.second, cluster2);

      FeatureVector fea(input.width, input.height);
      fea.compute(cluster1, cluster2, input.thickness, merged_cluster, false,
                  matching_pair_ptr);
      if (fea.features_.front() >= std::numeric_limits<double>::infinity())
        continue;

      neg_cluster_samples.emplace_back(
        std::make_pair(cluster_stroke_ind, sub_cluster_stroke_ind));
    }
  }
}
