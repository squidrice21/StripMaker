#include "Filters.h"

#include "../Logger.h"
#include "../Util.h"
#include "../solve/SolveUtil.h"
#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Utils.h"
#include "Features.h"

#include <fstream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

double fit_angle_threshold = 60;
double filter_min_length = 20;
double large_angle_threshold = 35;

double filter_max_overlapping_ratio = 0.2;

size_t min_overlapping_xsec_count = 5;

double alignment_threshold = 2;

double final_merge_ratio_threshold = 1.5;

bool filter_sample(const Cluster &cluster1,
                   const std::set<int> &selected_strokes1,
                   const Cluster &cluster2,
                   const std::set<int> &selected_strokes2) {
  std::unordered_map<size_t, std::vector<Cluster::Stroke const *>> c2strokes;
  for (const auto &s : cluster1.original_input_strokes) {
    if (selected_strokes1.empty() || selected_strokes1.count(s.stroke_ind))
      c2strokes[0].emplace_back(&s);
  }
  for (const auto &s : cluster2.original_input_strokes) {
    if (selected_strokes2.empty() || selected_strokes2.count(s.stroke_ind))
      c2strokes[1].emplace_back(&s);
  }

  double min_dist = std::numeric_limits<double>::infinity();
  for (const auto &c_s1 : c2strokes) {
    for (const auto &c_s2 : c2strokes) {
      if (c_s1.first == c_s2.first)
        continue;

      for (const auto &s1 : c_s1.second) {
        for (const auto &s2 : c_s2.second) {
          double dist12 = stroke_sample_distance(*s1, *s2);
          double dist21 = stroke_sample_distance(*s2, *s1);
          double dist = std::min(dist12, dist21);
          min_dist = std::min(min_dist, dist);
          if (min_dist < pointwise_distance_threshold)
            return false;
        }
      }
    }
  }

  return min_dist >= pointwise_distance_threshold;
}

bool filter_sample(const Cluster &cluster1,
                   const std::vector<int> &selected_strokes1,
                   const Cluster &cluster2,
                   const std::vector<int> &selected_strokes2) {
  std::set<int> selected_strokes_set1, selected_strokes_set2;
  selected_strokes_set1.insert(selected_strokes1.begin(),
                               selected_strokes1.end());
  selected_strokes_set2.insert(selected_strokes2.begin(),
                               selected_strokes2.end());
  return filter_sample(cluster1, selected_strokes_set1, cluster2,
                       selected_strokes_set2);
}

bool filter_sample_overlapping(
  const Cluster &in_cluster1, const std::vector<int> &selected_strokes1,
  const Cluster &in_cluster2, const std::vector<int> &selected_strokes2,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  const std::string &cache_folder, const std::string &vis_folder) {
  Cluster cluster1, cluster2;
  std::string json_filename1 = "";
  std::string json_filename2 = "";
  if (!cache_folder.empty()) {
    std::vector<size_t> stroke_indices1, stroke_indices2;
    if (!selected_strokes1.empty()) {
      for (auto const &s : selected_strokes1) {
        stroke_indices1.emplace_back(s);
      }
    } else {
      for (auto const &s : in_cluster1.strokes) {
        stroke_indices1.emplace_back(s.stroke_ind);
      }
    }
    if (!selected_strokes2.empty()) {
      for (auto const &s : selected_strokes2) {
        stroke_indices2.emplace_back(s);
      }
    } else {
      for (auto const &s : in_cluster2.strokes) {
        stroke_indices2.emplace_back(s.stroke_ind);
      }
    }
    std::sort(stroke_indices1.begin(), stroke_indices1.end());
    std::sort(stroke_indices2.begin(), stroke_indices2.end());

    std::string hash_str1 = example_hash(stroke_indices1);
    std::string hash_str2 = example_hash(stroke_indices2);

    json_filename1 = cache_folder + "/" + hash_str1 + ".json";
    json_filename2 = cache_folder + "/" + hash_str2 + ".json";
    bool hit1 = read_cache(json_filename1, cluster1);
    bool hit2 = read_cache(json_filename2, cluster2);

    if (!hit1)
      cluster1 = in_cluster1;
    if (!hit2)
      cluster2 = in_cluster2;

    if (cluster1.strokes.size() == 1 && cluster1.fit.centerline.empty()) {
      for (const auto &p : cluster1.strokes.front().points)
        cluster1.fit.centerline.emplace_back(p);
    }
    if (cluster2.strokes.size() == 1 && cluster2.fit.centerline.empty()) {
      for (const auto &p : cluster2.strokes.front().points)
        cluster2.fit.centerline.emplace_back(p);
    }

    // Cache exists but the individual parameterization failed
    if ((hit1 && cluster1.fit.centerline.empty()) ||
        (hit2 && cluster2.fit.centerline.empty())) {
      return true;
    }
  } else {
    cluster1 = in_cluster1;
    cluster2 = in_cluster2;
  }

  if (cluster1.fit.centerline.size() < 2 || cluster2.fit.centerline.size() < 2)
    return true;

  // Determine overlapping using the two fitting curves
  Cluster fit_cluster1, fit_cluster2, fit_merged_cluster;
  reuse_cluster(cluster1, fit_cluster1, true);
  reuse_cluster(cluster2, fit_cluster2, true);

  FeatureVector fea(0, 0);
  bool reorient = true;
  fea.compute_angle(fit_cluster1, fit_cluster2, 1.0, fit_merged_cluster,
                    nullptr, reorient, vis_folder);

  // If parameterization fails, skip filtering and leave the work for the full
  // run
  if (fit_merged_cluster.strokes.empty())
    return false;

  if (fea.features_.front() >= std::numeric_limits<double>::infinity()) {
    // Skip if both are single-stroke clusters
    if (cluster1.strokes.size() == 1 && cluster2.strokes.size() == 1)
      return false;

    // Skip these two test if any stroke is too short (the sampling may miss the
    // overlapping)
    if (total_length(fit_cluster1.strokes.front().points) < filter_min_length ||
        total_length(fit_cluster2.strokes.front().points) < filter_min_length)
      return false;

    if (is_roughly_circular(fit_cluster1.strokes.front()) ||
        is_roughly_circular(fit_cluster2.strokes.front()))
      return false;

    return true;
  }

  // To avoid parameterization again, check angle here.
  for (size_t i = 0; i < fea.descriptions_.size(); ++i) {
    if (fea.descriptions_[i] == "median_angle") {
      if (fea.features_[i] > fit_angle_threshold) {
        return true;
      } else
        return false;
    }
  }

  return false;
}

bool filter_merging(const Cluster &cluster1, const Cluster &cluster2,
                    const std::string &vis_folder) {
  assert(!cluster1.fit.widths.empty() && !cluster2.fit.widths.empty());

  // Determine overlapping using the two fitting curves
  Cluster fit_cluster1, fit_cluster2, fit_merged_cluster;
  reuse_cluster(cluster1, fit_cluster1, true);
  reuse_cluster(cluster2, fit_cluster2, true);

  FeatureVector fea(0, 0);
  bool reorient = true;
  double norm_thickness = 1.0;
  fea.compute_distance(fit_cluster1, fit_cluster2, norm_thickness,
                       fit_merged_cluster, nullptr, reorient, vis_folder);

  // For the final merging, if the parameterization fails, we filter this pair.
  if (fit_merged_cluster.strokes.empty())
    return true;

  if (fea.features_.front() >= std::numeric_limits<double>::infinity())
    return true;

  double to_width1 =
    *std::max_element(cluster1.fit.widths.begin(), cluster1.fit.widths.end());
  double to_width2 =
    *std::max_element(cluster2.fit.widths.begin(), cluster2.fit.widths.end());
  std::vector<std::string> avg_dist_str;
  auto avg_dists = pick_features(
    fea, std::vector<std::string>{"average_distance"}, avg_dist_str);
  double avg_dist = avg_dists[0];
  avg_dist -= 0.5 * to_width1 + 0.5 * to_width2;
  avg_dist = std::max(avg_dist, 0.0);

  double max_gap = -1;
  for (const auto &c_fea : solve_context.cluster_features) {
    std::vector<std::string> max_width_str;
    auto max_gaps = pick_features(
      c_fea.second, std::vector<std::string>{"max_distance"}, max_width_str);
    double max_g = max_gaps[0];
    max_gap = std::max(max_gap, max_g);
  }
  max_gap += 1.0;

  return avg_dist > final_merge_ratio_threshold * max_gap;
}

bool filter_sample_angle(
  const Cluster &in_cluster1, const std::vector<int> &selected_strokes1,
  const Cluster &in_cluster2, const std::vector<int> &selected_strokes2,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  const std::string &cache_folder) {
  Cluster cluster1, cluster2;
  std::string json_filename1 = "";
  std::string json_filename2 = "";
  if (!cache_folder.empty()) {
    std::vector<size_t> stroke_indices1, stroke_indices2;
    for (auto const &s : selected_strokes1) {
      stroke_indices1.emplace_back(s);
    }
    for (auto const &s : selected_strokes2) {
      stroke_indices2.emplace_back(s);
    }
    std::sort(stroke_indices1.begin(), stroke_indices1.end());
    std::sort(stroke_indices2.begin(), stroke_indices2.end());

    std::string hash_str1 = example_hash(stroke_indices1);
    std::string hash_str2 = example_hash(stroke_indices2);

    json_filename1 = cache_folder + "/" + hash_str1 + ".json";
    json_filename2 = cache_folder + "/" + hash_str2 + ".json";
    bool hit1 = read_cache(json_filename1, cluster1);
    bool hit2 = read_cache(json_filename2, cluster2);

    if (!hit1)
      cluster1 = in_cluster1;
    if (!hit2)
      cluster2 = in_cluster2;

    if (cluster1.strokes.size() == 1 && cluster1.fit.centerline.empty()) {
      for (const auto &p : cluster1.strokes.front().points)
        cluster1.fit.centerline.emplace_back(p);
    }
    if (cluster2.strokes.size() == 1 && cluster2.fit.centerline.empty()) {
      for (const auto &p : cluster2.strokes.front().points)
        cluster2.fit.centerline.emplace_back(p);
    }

    // Cache exists but the individual parameterization failed
    if ((hit1 && cluster1.fit.centerline.empty()) ||
        (hit2 && cluster2.fit.centerline.empty())) {
      return true;
    }
  }

  // Determine overlapping using the two fitting curves
  Cluster fit_cluster1, fit_cluster2, fit_merged_cluster;
  reuse_cluster(cluster1, fit_cluster1);
  reuse_cluster(cluster2, fit_cluster2);
  FeatureVector fea(0, 0);
  fea.compute(fit_cluster1, fit_cluster2, 1.0, fit_merged_cluster, true,
              matching_pair_ptr);
  for (size_t i = 0; i < fea.descriptions_.size(); ++i) {
    if (fea.descriptions_[i] == "median_angle") {
      if (fea.features_[i] > fit_angle_threshold) {
        return true;
      } else
        return false;
    }
  }

  return false;
}

bool filter_sample_new(const Cluster &cluster1,
                       const std::set<int> &selected_strokes1,
                       const Cluster &cluster2,
                       const std::set<int> &selected_strokes2) {
  std::unordered_map<size_t, std::vector<Cluster::Stroke const *>> c2strokes;
  for (const auto &s : cluster1.strokes) {
    if (selected_strokes1.count(s.stroke_ind))
      c2strokes[0].emplace_back(&s);
  }
  for (const auto &s : cluster2.strokes) {
    if (selected_strokes2.count(s.stroke_ind))
      c2strokes[1].emplace_back(&s);
  }

  struct AABB {
    std::pair<double, double> x_range, y_range;
  };
  auto build_AABB = [](const Cluster::Stroke &s) -> AABB {
    double min_x = std::numeric_limits<double>::infinity();
    double max_x = -std::numeric_limits<double>::infinity();
    double min_y = std::numeric_limits<double>::infinity();
    double max_y = -std::numeric_limits<double>::infinity();

    for (const auto &point : s.points) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }

    AABB aabb;
    aabb.x_range = std::make_pair(min_x - pointwise_distance_threshold,
                                  max_x + pointwise_distance_threshold);
    aabb.y_range = std::make_pair(min_y - pointwise_distance_threshold,
                                  max_y + pointwise_distance_threshold);

    return aabb;
  };

  // Use the union of per-stroke AABBs as the bounding area for a cluster
  std::unordered_map<size_t, AABB> s2aabbs;
  for (const auto &c_s : c2strokes) {
    for (const auto &s : c_s.second) {
      s2aabbs[s->stroke_ind] = build_AABB(*s);
    }
  }

  bool hit = false;
  for (const auto &c_s1 : c2strokes) {
    for (const auto &c_s2 : c2strokes) {
      if (c_s1.first == c_s2.first)
        continue;

      for (const auto &s1 : c_s1.second) {
        for (const auto &s2 : c_s2.second) {
          for (const auto &p1 : s1->points) {
            if (p1.x > s2aabbs[s2->stroke_ind].x_range.first &&
                p1.x < s2aabbs[s2->stroke_ind].x_range.second &&
                p1.y > s2aabbs[s2->stroke_ind].y_range.first &&
                p1.y < s2aabbs[s2->stroke_ind].y_range.second) {
              hit = true;
              goto found_hit;
            }
          }
        }
      }

      for (const auto &s2 : c_s2.second) {
        for (const auto &s1 : c_s1.second) {
          for (const auto &p2 : s2->points) {
            if (p2.x > s2aabbs[s1->stroke_ind].x_range.first &&
                p2.x < s2aabbs[s1->stroke_ind].x_range.second &&
                p2.y > s2aabbs[s1->stroke_ind].y_range.first &&
                p2.y < s2aabbs[s1->stroke_ind].y_range.second) {
              hit = true;
              goto found_hit;
            }
          }
        }
      }
    }
  }

found_hit:
  if (!hit)
    return true;

  double min_dist = std::numeric_limits<double>::infinity();
  for (const auto &c_s1 : c2strokes) {
    for (const auto &c_s2 : c2strokes) {
      if (c_s1.first == c_s2.first)
        continue;

      for (const auto &s1 : c_s1.second) {
        for (const auto &s2 : c_s2.second) {
          double dist12 = stroke_sample_distance(*s1, *s2);
          double dist21 = stroke_sample_distance(*s2, *s1);
          double dist = std::min(dist12, dist21);
          min_dist = std::min(min_dist, dist);
          if (min_dist < pointwise_distance_threshold)
            return false;
        }
      }
    }
  }

  return min_dist >= pointwise_distance_threshold;
}

bool filter_feature(const FeatureVector &fea) {
  if (fea.features_.empty())
    return false;
  return fea.features_.front() >= std::numeric_limits<double>::infinity();
}

bool filter_test_overlapping(const Cluster &cluster1, const Cluster &cluster2,
                             const Cluster &merged_cluster) {
  size_t overlapping_xsec_count = 0;
  std::map<size_t, size_t> s2c;
  for (const auto &s : cluster1.strokes) {
    s2c[s.stroke_ind] = 0;
  }
  for (const auto &s : cluster2.strokes) {
    s2c[s.stroke_ind] = 1;
  }

  std::vector<double> continuous_overlapping_u;
  std::vector<bool> is_current_u;
  std::map<size_t, std::vector<double>> cluster_u;
  cluster_u[0] = std::vector<double>();
  cluster_u[1] = std::vector<double>();
  for (auto const &xsec : merged_cluster.xsecs) {
    if (is_overlapping(xsec, s2c)) {
      overlapping_xsec_count++;

      if (continuous_overlapping_u.empty() || !is_current_u.back()) {
        continuous_overlapping_u.emplace_back(xsec.u);
        is_current_u.emplace_back(true);
      }
    } else if (!continuous_overlapping_u.empty() && is_current_u.back()) {
      is_current_u.back() = false;
      continuous_overlapping_u.back() =
        std::abs(xsec.u - continuous_overlapping_u.back());
    }

    std::unordered_set<size_t> c_indices;
    for (const auto &p : xsec.points) {
      if (s2c.count(p.stroke_ind)) {
        c_indices.emplace(s2c.at(p.stroke_ind));
      }
    }
    for (auto c : c_indices)
      cluster_u[c].emplace_back(xsec.u);
  }
  if (!continuous_overlapping_u.empty() && is_current_u.back()) {
    is_current_u.back() = false;
    continuous_overlapping_u.back() =
      std::abs(merged_cluster.xsecs.back().u - continuous_overlapping_u.back());
  }

  if (cluster_u[0].empty() || cluster_u[1].empty())
    return true;

  auto min0 = *std::min_element(cluster_u[0].begin(), cluster_u[0].end());
  auto max0 = *std::max_element(cluster_u[0].begin(), cluster_u[0].end());
  auto min1 = *std::min_element(cluster_u[1].begin(), cluster_u[1].end());
  auto max1 = *std::max_element(cluster_u[1].begin(), cluster_u[1].end());

  // No overlapping
  if (continuous_overlapping_u.empty())
    return true;

  auto max_overlapping = *std::max_element(continuous_overlapping_u.begin(),
                                           continuous_overlapping_u.end());

  auto max_overlapping_ratio =
    max_overlapping / std::min(max0 - min0, max1 - min1);

  bool ret = (overlapping_xsec_count < min_overlapping_xsec_count &&
              max_overlapping_ratio < filter_max_overlapping_ratio);

  return ret;
}
