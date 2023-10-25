#include "SampleFilters.h"

#include "../../stroke_strip_src/Serialization.h"
#include "../Util.h"
#include "../feature/FeatureImplementations.h"
#include "../feature/Features.h"
#include "../feature/Filters.h"
#include "glm/detail/type_vec.hpp"

#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

double fillin_gap_change_threshold = 2;
double continuity_overlapping_threshold = 0.1;
double continuity_width_ratio = 1.5;

static double get_cluster_median_width(const Cluster &cluster) {
  PercentileMergedWidthStrokeClusterFeature fea(0.5, 1);
  auto w_vec = fea.operator()(cluster, FittedCurve(), FittedCurve());
  double median_w =
    (w_vec.empty() || w_vec.front() == std::numeric_limits<double>::infinity())
      ? 1.0
      : w_vec.front();

  return median_w;
}

static void find_center(const Cluster::Stroke &s, glm::dvec2 &center) {
  center = glm::dvec2(0, 0);
  for (const auto &p : s.points) {
    center += p;
  }
  center /= s.points.size();
}

// Assume merged_cluster contains cluster1 and cluster2 which are sub-clusters
// of pos_cluster (a GT cluster) and both are parameterized.
bool filter_orientation_sample(const Cluster &merged_cluster,
                               const Cluster &pos_cluster) {
  auto is_stroke_flipped = [](const Cluster &ref_cluster,
                              const Cluster::Stroke &s) -> bool {
    for (const auto &ref_s : ref_cluster.strokes) {
      if (ref_s.stroke_ind == s.stroke_ind) {
        glm::dvec2 ref_center, s_center;
        find_center(ref_s, ref_center);
        find_center(s, s_center);

        // assert(ref_s.points.size() == s.points.size());
        double ff_dist = glm::distance((ref_s.points.front() - ref_center),
                                       (s.points.front() - s_center));
        double fb_dist = glm::distance((ref_s.points.front() - ref_center),
                                       (s.points.back() - s_center));
        assert(ff_dist < 0.5 || fb_dist < 0.5);

        if (glm::distance(ref_s.points.front(), ref_s.points.back()) > 0.1)
          return ff_dist < fb_dist;
        else if (ref_s.points.size() > 1) {
          double ff_dist = glm::distance((ref_s.points[1] - ref_center),
                                         (s.points[1] - s_center));
          double fb_dist =
            glm::distance((ref_s.points[1] - ref_center),
                          (s.points[s.points.size() - 2] - s_center));
          assert(ff_dist < 0.1 || fb_dist < 0.1);
          return ff_dist < fb_dist;
        }
        return false;
      }
    }
    abort();
  };

  // Check if all strokes are flipped or not flipped.
  std::vector<bool> stroke_flipped;
  for (const auto &s : merged_cluster.strokes) {
    stroke_flipped.emplace_back(is_stroke_flipped(pos_cluster, s));
  }

  bool front_flipped = stroke_flipped.front();
  bool ori_consistent = true;
  for (auto f : stroke_flipped) {
    bool ff = (front_flipped) ? f : !f;
    ori_consistent &= ff;
    if (!ori_consistent)
      return true;
  }

  return false;
}

std::vector<double>
get_cluster_percentile_inner_gap(double percentile, const Cluster &gt_cluster) {
  std::vector<double> gt_dist_vec;
  {
    for (const auto &c_samples : gt_cluster.xsecs) {
      double max_gap = -std::numeric_limits<double>::infinity();
      for (size_t i = 0; i + 1 < c_samples.points.size(); ++i) {
        double dist = glm::distance(c_samples.points[i].point,
                                    c_samples.points[i + 1].point) -
                      1.0;
        dist = std::max(dist, 0.0);
        max_gap = std::max(dist, max_gap);
      }
      if (max_gap > -std::numeric_limits<double>::infinity()) {
        max_gap = std::max(max_gap, 0.0);
        gt_dist_vec.emplace_back(max_gap);
      }
    }
    std::sort(gt_dist_vec.begin(), gt_dist_vec.end());
  }

  if (gt_dist_vec.empty())
    return gt_dist_vec;

  int pos = (gt_dist_vec.size() - 1) * percentile;
  pos = std::max(pos, 0);
  auto p = gt_dist_vec.begin() + pos;
  std::nth_element(gt_dist_vec.begin(), p, gt_dist_vec.end());

  double gt_dist = gt_dist_vec[pos];
  gt_dist_vec.clear();
  gt_dist_vec.emplace_back(gt_dist);

  return gt_dist_vec;
}

bool filter_fillin_sample(const Cluster &merged_cluster,
                          const Cluster &gt_cluster) {
  // 1. Compute the current and the final largest gap size
  double percentile = 0.5;

  std::vector<double> cur_dist_vec =
    get_cluster_percentile_inner_gap(percentile, merged_cluster);
  std::vector<double> gt_dist_vec =
    get_cluster_percentile_inner_gap(percentile, gt_cluster);

  if (cur_dist_vec.empty() || gt_dist_vec.empty())
    return true;

  // 2. Compare to determine if we have fill-in
  double cur_dist = cur_dist_vec.front();
  double gt_dist = gt_dist_vec.front();

  cur_dist = std::max(1.0, cur_dist);
  gt_dist = std::max(1.0, gt_dist);

  return (cur_dist / gt_dist > fillin_gap_change_threshold);
}

struct LocalSplit {
  double distance;
  std::set<size_t> small, large;
  bool single = false;
  bool split = false;
};

bool split_branch(std::vector<LocalSplit> &xsec_splits,
                  const std::map<size_t, size_t> &s2c,
                  std::set<size_t> &selected_strokes1,
                  std::set<size_t> &selected_strokes2,
                  const double min_dist = -1) {
  // Split from the largest gap until the two subclusters are not divideable
  // Make sure to find the split position within parts that are not split yet.
  auto split_max = std::max_element(
    xsec_splits.begin(), xsec_splits.end(),
    [](LocalSplit const &x, LocalSplit const &y) {
      return (int(!x.split) * 1000 + x.distance - int(x.single) * 10000) <
             (int(!y.split) * 1000 + y.distance - int(y.single) * 10000);
    });

  // There's no overlapping or not divideable
  if (split_max->distance <= 0 || split_max->split || split_max->single)
    return false;

  // This split is not large enough
  if (split_max->distance < min_dist)
    return false;

  auto split_branch = [&xsec_splits,
                       &s2c](size_t start_i, std::set<size_t> &strokes1,
                             std::set<size_t> &strokes2, bool inc) {
    int prev = -1;
    int i;
    for (i = start_i; (inc) ? i < xsec_splits.size() : i >= 0;
         (inc) ? ++i : --i) {
      auto &split = xsec_splits[i];

      // Stop if it becomes undivideable or the order flipped
      if (split.distance <= 0 && !split.single)
        break;

      if (!split.single) {
        bool flipped = false;
        if (prev >= 0) {
          for (auto l : split.large) {
            if (strokes1.count(l)) {
              flipped = true;
              break;
            }
          }
          for (auto l : split.small) {
            if (strokes2.count(l)) {
              flipped = true;
              break;
            }
          }
        }
        if (flipped)
          break;
        strokes1.insert(split.small.begin(), split.small.end());
        strokes2.insert(split.large.begin(), split.large.end());
      }
      // else {
      //   assert(!strokes1.empty() && !strokes2.empty());
      //   if (!split.small.empty() &&
      //       s2c.at(*split.small.begin()) == s2c.at(*strokes1.begin())) {
      //     strokes1.insert(split.small.begin(), split.small.end());
      //   } else {
      //     strokes2.insert(split.small.begin(), split.small.end());
      //   }
      // }

      split.split = true;
      prev = i;
    }
  };

  size_t start_i = std::distance(xsec_splits.begin(), split_max);
  split_branch(start_i, selected_strokes1, selected_strokes2, true);
  split_branch(start_i, selected_strokes1, selected_strokes2, false);

  return true;
}

// selected_strokes1&2 could return arguments or inputs (from some past run)
bool filter_branch_sample(
  const Cluster &merged_cluster, const Cluster &neg_cluster,
  std::vector<std::pair<std::set<size_t>, std::set<size_t>>> &selected_branches,
  const double min_dist) {
  int c1 = -1, c2 = -1;
  std::map<size_t, size_t> s2c;
  for (const auto &s : neg_cluster.strokes) {
    s2c[s.stroke_ind] = s.cluster_ind;
  }

  if (selected_branches.empty()) {
    std::vector<LocalSplit> xsec_splits;
    for (const auto &xsec : neg_cluster.xsecs) {
      xsec_splits.emplace_back();

      std::vector<double> euclidean_clusterwise;
      int split_i = -1;
      for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
        if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
          split_i = i;
          double dis =
            glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
          dis = std::max(dis, 0.0);
          euclidean_clusterwise.emplace_back(dis);
          // The two clusters are not clearly divided
          if (euclidean_clusterwise.size() > 1) {
            euclidean_clusterwise.clear();
            euclidean_clusterwise.emplace_back(0);
            break;
          }
        }
      }

      if (!euclidean_clusterwise.empty()) {
        xsec_splits.back().distance = euclidean_clusterwise.front();
        for (size_t i = 0; i < xsec.points.size(); ++i) {
          if (i <= split_i)
            xsec_splits.back().small.emplace(xsec.points[i].stroke_ind);
          else
            xsec_splits.back().large.emplace(xsec.points[i].stroke_ind);
        }
      } else {
        xsec_splits.back().distance = -1;
        xsec_splits.back().single = true;
        for (size_t i = 0; i < xsec.points.size(); ++i) {
          xsec_splits.back().small.emplace(xsec.points[i].stroke_ind);
        }
      }
    }

    bool b_split = false;
    do {
      std::set<size_t> strokes1, strokes2;

      // If this is the first branch. Assuming the labeling is correct, we
      // always cut this first branch.
      double min_d = min_dist;
      if (selected_branches.empty())
        min_d = -1;
      b_split = split_branch(xsec_splits, s2c, strokes1, strokes2, min_d);

      // Add to the final results
      if (b_split) {
        selected_branches.emplace_back();
        auto &selected_strokes1 = selected_branches.back().first;
        auto &selected_strokes2 = selected_branches.back().second;

        if (c1 < 0 || c2 < 0) {
          selected_strokes1.insert(strokes1.begin(), strokes1.end());
          selected_strokes2.insert(strokes2.begin(), strokes2.end());
          if (!selected_strokes1.empty())
            c1 = s2c[*selected_strokes1.begin()];
          if (!selected_strokes2.empty())
            c2 = s2c[*selected_strokes2.begin()];
        } else if (!strokes1.empty() && !strokes2.empty()) {
          size_t cc1 = s2c[*strokes1.begin()];
          size_t cc2 = s2c[*strokes2.begin()];
          assert(cc1 != cc2);
          if (c1 == cc1) {
            selected_strokes1.insert(strokes1.begin(), strokes1.end());
            selected_strokes2.insert(strokes2.begin(), strokes2.end());
          } else if (c2 == cc1) {
            selected_strokes1.insert(strokes2.begin(), strokes2.end());
            selected_strokes2.insert(strokes1.begin(), strokes1.end());
          }
        }
        assert(!selected_branches.back().first.empty() &&
               !selected_branches.back().second.empty());
      }
    } while (b_split);
  }

  // Actually determine if we want to filter this case.
  // One of the two subclusters has to be in one of the selected stroke set
  std::map<size_t, std::vector<size_t>> cur_c2ss;
  for (const auto &s : merged_cluster.strokes) {
    cur_c2ss[s.cluster_ind].emplace_back(s.stroke_ind);
  }

  std::vector<std::map<size_t, bool>> inside_vec, covered_vec;
  for (const auto &branch : selected_branches) {
    inside_vec.emplace_back();
    covered_vec.emplace_back();

    for (const auto &c2ss : cur_c2ss) {
      const auto &selected_strokes1 = branch.first;
      const auto &selected_strokes2 = branch.second;

      auto &inside = inside_vec.back();
      auto &covered = covered_vec.back();

      inside[c2ss.first] = true;
      covered[c2ss.first] = false;
      for (auto s : c2ss.second) {
        if (!selected_strokes1.count(s) && !selected_strokes2.count(s)) {
          inside[c2ss.first] = false;
        }
        if (selected_strokes1.count(s) || selected_strokes2.count(s)) {
          covered[c2ss.first] = true;
        }
      }
    }
  }

  // bool is_valid = false;
  // for (const auto &c2in : inside) {
  //   is_valid |= c2in.second;
  // }
  bool is_valid_all = false;
  for (const auto &covered : covered_vec) {
    bool is_valid = true;
    for (const auto &c2in : covered) {
      is_valid &= c2in.second;
    }
    is_valid_all |= is_valid;
  }

  return !is_valid_all;
}

bool filter_overlapping_inconsist_sample(const Cluster &merged_cluster,
                                         const Cluster &final_cluster) {
  std::map<size_t, size_t> cur_s2c;
  for (const auto &s : merged_cluster.strokes) {
    cur_s2c[s.stroke_ind] = s.cluster_ind;
  }

  auto find_overlapping_stroke_pairs =
    [&cur_s2c](const std::vector<Cluster::XSec> &xsecs,
               std::set<std::pair<size_t, size_t>> &overlapping_set) {
      for (const auto &xsec : xsecs) {
        for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
          for (size_t j = i + 1; j < xsec.points.size(); ++j) {
            size_t i_s = xsec.points[i].stroke_ind;
            size_t j_s = xsec.points[j].stroke_ind;
            if (cur_s2c.count(i_s) && cur_s2c.count(j_s) &&
                cur_s2c[i_s] != cur_s2c[j_s]) {
              overlapping_set.emplace(
                std::make_pair(std::min(i_s, j_s), std::max(i_s, j_s)));
            }
          }
        }
      }
    };

  std::set<std::pair<size_t, size_t>> overlapping_final, overlapping_cur;
  find_overlapping_stroke_pairs(final_cluster.xsecs, overlapping_final);
  find_overlapping_stroke_pairs(merged_cluster.xsecs, overlapping_cur);
  // return overlapping_final != overlapping_cur; // Xsecs may differ slightly
  // for pairs that are barely overlapping
  return overlapping_cur.empty() != overlapping_final.empty();
}

bool filter_large_angle(
  const Cluster &cluster1, const Cluster &cluster2,
  const Cluster &merged_cluster,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  FeatureVector fea(0, 0);
  bool reorient = true;
  Cluster c1 = cluster1, c2 = cluster2, mc = merged_cluster;
  fea.compute_angle(c1, c2, 1.0, mc, matching_pair_ptr, reorient);

  assert(fea.features_.front() < std::numeric_limits<double>::infinity());

  for (size_t i = 0; i < fea.descriptions_.size(); ++i) {
    if (fea.descriptions_[i] == "average_angle") {
      if (fea.features_[i] > large_angle_threshold) {
        return true;
      } else
        return false;
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////
// Filter the generated samples. Assume clusters in input are fitted.
void filter_sample_extra(
  Input &input, std::vector<Sample> &pos_samples,
  std::vector<Sample> &neg_samples,
  std::vector<FilteredSample> &filtered_samples,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  const std::string &cache_folder, const std::string &output_vis_folder,
  bool disable_short_overlapping) {
  // Prepare data: conversion, indexing.
  std::map<size_t, size_t> s2c;
  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : input.clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      s2c[s.stroke_ind] = s.cluster_ind;
      sid2s[s.stroke_ind] = s;
    }
  }

  std::vector<std::pair<Cluster, Cluster>> clusters;
  std::vector<bool> is_pos;
  std::map<size_t, Sample const *> sample_indexing;
  std::map<std::pair<size_t, size_t>,
           std::vector<std::pair<std::set<size_t>, std::set<size_t>>>>
    neg_selected_cache;
  for (const auto &samples : pos_samples) {
    sample_indexing[clusters.size()] = &samples;

    clusters.emplace_back();
    for (const auto &s : samples.first) {
      clusters.back().first.strokes.emplace_back(sid2s[s]);
    }
    for (const auto &s : samples.second) {
      clusters.back().second.strokes.emplace_back(sid2s[s]);
    }
    clusters.back().first.original_input_strokes =
      clusters.back().first.strokes;
    clusters.back().second.original_input_strokes =
      clusters.back().second.strokes;
    is_pos.emplace_back(true);
  }
  for (const auto &samples : neg_samples) {
    sample_indexing[clusters.size()] = &samples;

    clusters.emplace_back();
    for (const auto &s : samples.first) {
      clusters.back().first.strokes.emplace_back(sid2s[s]);
    }
    for (const auto &s : samples.second) {
      clusters.back().second.strokes.emplace_back(sid2s[s]);
    }
    clusters.back().first.original_input_strokes =
      clusters.back().first.strokes;
    clusters.back().second.original_input_strokes =
      clusters.back().second.strokes;
    is_pos.emplace_back(false);
  }

  std::vector<Sample> out_pos_samples, out_neg_samples;

  auto record_filtered_sample = [&filtered_samples, &sample_indexing, &is_pos](
                                  size_t i, FilteredSample::FilterCause cause) {
    FilteredSample f_sample;
    f_sample.sample = *sample_indexing[i];
    f_sample.is_pos = is_pos[i];
    f_sample.cause = cause;
    filtered_samples.emplace_back(f_sample);
  };

  // Call filters (differ based on pos or neg)
  double thickness = 1;
  bool cache_exists = false;
  for (size_t i = 0; i < clusters.size(); ++i) {
    auto &cc = clusters[i];

    auto &gt_cluster1 =
      input.clusters.at(s2c[cc.first.strokes.front().stroke_ind]);
    if (gt_cluster1.fit.centerline.empty()) {
      record_filtered_sample(i, FilteredSample::Parameterization);
      continue;
    }
    auto &gt_cluster2 =
      input.clusters.at(s2c[cc.second.strokes.front().stroke_ind]);
    if (gt_cluster2.fit.centerline.empty()) {
      record_filtered_sample(i, FilteredSample::Parameterization);
      continue;
    }

    // Max median cluster width between the two corresponding GT clusters
    double max_gt_cluster_width =
      std::max(get_cluster_median_width(gt_cluster1),
               get_cluster_median_width(gt_cluster2));

    std::string stroke_str;
    for (auto const &s : cc.first.strokes) {
      if (!stroke_str.empty())
        stroke_str += ' ';
      stroke_str += std::to_string(s.stroke_ind);
    }
    stroke_str += ',';
    for (auto const &s : cc.second.strokes) {
      if (!stroke_str.empty())
        stroke_str += ' ';
      stroke_str += std::to_string(s.stroke_ind);
    }
    // logger().info("Compute on {}", stroke_str);
    std::cout << "Compute on " + stroke_str << std::endl;

    // Look up cached parameterization
    std::string hash_str = "";
    std::string svg_name = (is_pos[i]) ? "pos_" : "neg_";
    {
      std::vector<size_t> stroke_indices1, stroke_indices2;
      for (auto const &s : cc.first.strokes) {
        stroke_indices1.emplace_back(s.stroke_ind);
      }
      for (auto const &s : cc.second.strokes) {
        stroke_indices2.emplace_back(s.stroke_ind);
      }

      hash_str = example_hash(stroke_indices1, stroke_indices2);
    }
    std::string json_filename1, json_filename2, json_filename_all;

    Cluster merged_cluster, neg_gt_cluster;
    json_filename1 = cache_folder + "/" + svg_name + hash_str + "_c1.json";
    json_filename2 = cache_folder + "/" + svg_name + hash_str + "_c2.json";
    json_filename_all = cache_folder + "/" + svg_name + hash_str + "_comb.json";

    // This is for the negative examples
    std::string gt_hash_str = "";
    {
      std::vector<size_t> gt_stroke_indices1, gt_stroke_indices2;
      for (auto const &s : gt_cluster1.strokes) {
        gt_stroke_indices1.emplace_back(s.stroke_ind);
      }
      for (auto const &s : gt_cluster2.strokes) {
        gt_stroke_indices2.emplace_back(s.stroke_ind);
      }
      gt_hash_str = example_hash(gt_stroke_indices1, gt_stroke_indices2);
    }
    std::string gt_json_filename_all =
      cache_folder + "/" + svg_name + gt_hash_str + "_gt_comb.json";

    auto file_exists = [](const std::string &file_name) {
      std::ifstream f(file_name.c_str());
      return f.good();
    };

    if (!is_pos[i] && !cache_folder.empty() &&
        file_exists(gt_json_filename_all)) {
      read_json(gt_json_filename_all, neg_gt_cluster);

      for (const auto &s : gt_cluster1.original_input_strokes)
        neg_gt_cluster.original_input_strokes.emplace_back(s);
      for (const auto &s : gt_cluster2.original_input_strokes)
        neg_gt_cluster.original_input_strokes.emplace_back(s);

      // Combined parameterization failed
      if (neg_gt_cluster.strokes.empty() ||
          neg_gt_cluster.obj_term_values["alignment"] >=
            pos_parameterization_alignment_threshold) {
        record_filtered_sample(i, FilteredSample::Parameterization);
        continue;
      }
    } else if (!is_pos[i]) {
      bool reorient = false;
      FeatureVector fea(input.width, input.height);
      fea.compute(gt_cluster1, gt_cluster2, thickness, neg_gt_cluster, false,
                  matching_pair_ptr, reorient, output_vis_folder);
      if (neg_gt_cluster.fit.centerline.empty()) {
        record_filtered_sample(i, FilteredSample::Parameterization);
        continue;
      }
      if (!cache_folder.empty()) {
        write_json(gt_json_filename_all, neg_gt_cluster);
      }
    }
    //

    bool to_write = !cache_folder.empty() && !(file_exists(json_filename1) &&
                                               file_exists(json_filename2) &&
                                               file_exists(json_filename_all));
    if (!cache_folder.empty())
      cache_exists |= file_exists(json_filename1) &&
                      file_exists(json_filename2) &&
                      file_exists(json_filename_all);

    if (!cache_folder.empty() && file_exists(json_filename1) &&
        file_exists(json_filename2) && file_exists(json_filename_all)) {
      read_json(json_filename1, cc.first);
      read_json(json_filename2, cc.second);
      read_json(json_filename_all, merged_cluster);

      // Combined parameterization failed
      if (!(cc.first.fit.centerline.empty() ||
            cc.first.obj_term_values["alignment"] >=
              pos_parameterization_alignment_threshold) &&
          !(cc.second.fit.centerline.empty() ||
            cc.second.obj_term_values["alignment"] >=
              pos_parameterization_alignment_threshold) &&
          (merged_cluster.strokes.empty() ||
           merged_cluster.obj_term_values["alignment"] >=
             pos_parameterization_alignment_threshold)) {
        record_filtered_sample(i, FilteredSample::Parameterization);
        continue;
      }
    }
    else {
      // Fit the two individual clusters
      bool to_orient = true;
      fit_cluster(input.thickness, input.width, input.height, cc.first,
                  to_orient, matching_pair_ptr);
      fit_cluster(input.thickness, input.width, input.height, cc.second,
                  to_orient, matching_pair_ptr);
    }

    // Either parameterizations fail
    if ((cc.first.strokes.empty() ||
         cc.first.obj_term_values["alignment"] >=
           pos_parameterization_alignment_threshold) ||
        (cc.second.strokes.empty() ||
         cc.second.obj_term_values["alignment"] >=
           pos_parameterization_alignment_threshold)) {
      record_filtered_sample(i, FilteredSample::Parameterization);
      continue;
    }

    // To merge the two clusters if cached data doesn't exist.
    if (merged_cluster.fit.centerline.empty()) {
      bool reorient = false;
      FeatureVector fea(input.width, input.height);
      fea.compute(cc.first, cc.second, thickness, merged_cluster, false,
                  matching_pair_ptr, reorient, output_vis_folder);
    } else {
      for (const auto &s : cc.first.original_input_strokes)
        merged_cluster.original_input_strokes.emplace_back(s);
      for (const auto &s : cc.second.original_input_strokes)
        merged_cluster.original_input_strokes.emplace_back(s);
    }

    if (to_write) {
      write_json(json_filename1, cc.first);
      write_json(json_filename2, cc.second);
      write_json(json_filename_all, merged_cluster);
    }

    if (!(merged_cluster.fit.centerline.empty() ||
          merged_cluster.obj_term_values["alignment"] >=
            pos_parameterization_alignment_threshold)) {
      // Non overlapping
      bool overlapping_filter = true;
      for (const auto &xsec : merged_cluster.xsecs)
        if (is_overlapping(xsec)) {
          overlapping_filter = false;
          break;
        }
      if (overlapping_filter) {
        record_filtered_sample(i, FilteredSample::Overlapping);
        continue;
      }

      // Check if the orientation is consistent between the current example
      // and the final GT cluster
      {
        const auto &final_cluster = (is_pos[i]) ? gt_cluster1 : neg_gt_cluster;
        bool ori_filter =
          filter_orientation_sample(merged_cluster, final_cluster);
        if (is_pos[i] && ori_filter) {
          record_filtered_sample(i, FilteredSample::Orientation);
          continue;
        }
        // else if (!is_pos[i] && ori_filter) {
        //   out_pos_samples.emplace_back(*sample_indexing[i]);
        //   continue;
        // }
      }

      // Filter out the pairs that are not overlapping wrt the GT final
      // parameterization
      {
        const auto &final_cluster = (is_pos[i]) ? gt_cluster1 : neg_gt_cluster;
        bool overlapping_inconsist =
          filter_overlapping_inconsist_sample(merged_cluster, final_cluster);
        if (is_pos[i] && overlapping_inconsist) {
          record_filtered_sample(i, FilteredSample::OverlappingInconsistence);
          continue;
        }
        // else if (!is_pos[i] && overlapping_inconsist) {
        //   out_pos_samples.emplace_back(*sample_indexing[i]);
        //   continue;
        // }
      }

      // Short overlapping
      if (!disable_short_overlapping) {
        // This is the same filter as the test time one
        // The reasoning is that we don't trust the measurement if the
        // overlapping is too short both in the absolute sense and wrt the
        // individual cluster length
        bool short_overlapping_filter =
          filter_test_overlapping(cc.first, cc.second, merged_cluster);
        if (short_overlapping_filter) {
          record_filtered_sample(i, FilteredSample::ShortOverlapping);
          continue;
        }
      }

      // Large average angle (35 deg). This is only for positive training
      // examples.
      if (is_pos[i]) {
        bool large_angle_filter = filter_large_angle(
          cc.first, cc.second, merged_cluster, matching_pair_ptr);
        if (large_angle_filter) {
          record_filtered_sample(i, FilteredSample::LargeAngle);
          continue;
        }
      }

      // Filter
      if (is_pos[i]) {
        bool fillin_filter = filter_fillin_sample(merged_cluster, gt_cluster1);
        if (fillin_filter) {
          record_filtered_sample(i, FilteredSample::FillIn);
          continue;
        }

        out_pos_samples.emplace_back(*sample_indexing[i]);
      } else {
        std::vector<std::pair<std::set<size_t>, std::set<size_t>>>
          selected_branches;
        auto neg_key =
          std::make_pair(std::min(gt_cluster1.strokes.front().cluster_ind,
                                  gt_cluster2.strokes.front().cluster_ind),
                         std::max(gt_cluster1.strokes.front().cluster_ind,
                                  gt_cluster2.strokes.front().cluster_ind));
        if (neg_selected_cache.count(neg_key)) {
          selected_branches = neg_selected_cache[neg_key];
        }

        double min_distance = max_gt_cluster_width;
        assert(min_distance < std::numeric_limits<double>::infinity());

        bool branch_filter = filter_branch_sample(
          merged_cluster, neg_gt_cluster, selected_branches, min_distance);

        if (!neg_selected_cache.count(neg_key)) {
          neg_selected_cache[neg_key] = selected_branches;
        }

        if (branch_filter) {
          record_filtered_sample(i, FilteredSample::Branch);
          continue;
        }

        out_neg_samples.emplace_back(*sample_indexing[i]);
      }
    } else {
      record_filtered_sample(i, FilteredSample::Parameterization);
    }
  }

  pos_samples = std::move(out_pos_samples);
  neg_samples = std::move(out_neg_samples);
}
