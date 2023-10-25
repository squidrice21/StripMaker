#include "SampleGeneration.h"
#include "../../stroke_strip_src/Cluster.h"
#include "../../stroke_strip_src/SketchInfo.h"
#include "../Util.h"
#include "../feature/Features.h"
#include "../feature/Filters.h"

#include <algorithm>
#include <limits>
#include <set>
#include <unordered_set>
#include <vector>

template <typename T>
bool next_combination(const T first, const T last, int k) {
  const T subset = first + k;
  // empty container | k = 0 | k == n
  if (first == last || first == subset || last == subset) {
    return false;
  }
  T src = subset;
  while (first != src) {
    src--;
    if (*src < *(last - 1)) {
      T dest = subset;
      while (*src >= *dest) {
        dest++;
      }
      iter_swap(src, dest);
      rotate(src + 1, dest + 1, last);
      rotate(subset, subset + (last - dest) - 1, last);
      return true;
    }
  }
  // restore
  rotate(first, subset, last);
  return false;
}

static bool vector_finder(std::vector<int> vec, int number) {
  auto itr = std::find(vec.begin(), vec.end(), number);
  size_t index = std::distance(vec.begin(), itr);
  if (index != vec.size()) {
    return true;
  } else {
    return false;
  }
}

static bool is_overlapping(const std::vector<Cluster::XSec> &xsecs,
                           const std::vector<int> &pre_cluster_stroke_ind,
                           const std::vector<int> &post_cluster_stroke_ind) {
  bool pre_xsec_include_pre = false;
  bool pre_xsec_include_post = false;
  int number_change_pre = 0;
  int number_change_post = 0;
  bool is_first = true;
  for (const auto &xsec : xsecs) {
    bool is_include_pre = false;
    bool is_include_post = false;
    for (const auto &point : xsec.points) {
      bool vector_finder_result =
        vector_finder(pre_cluster_stroke_ind, point.stroke_ind);
      is_include_pre = is_include_pre || vector_finder_result;
      is_include_post = is_include_post || !(vector_finder_result);
    }
    if (is_first) {
      pre_xsec_include_pre = is_include_pre;
      pre_xsec_include_post = is_include_post;
      is_first = false;
      continue;
    }
    if (pre_xsec_include_pre != is_include_pre) {
      number_change_pre++;
      pre_xsec_include_pre = is_include_pre;
    }
    if (pre_xsec_include_post != is_include_post) {
      number_change_post++;
      pre_xsec_include_post = is_include_post;
    }
  }
  if (number_change_pre > 1 || number_change_post > 1) {
    return false;
  } else {
    return true;
  }
}

static void find_overlapping_strokes(const std::vector<Cluster::XSec> &xsecs,
                                     const std::set<size_t> &selected_strokes,
                                     std::set<size_t> &overlapping_set) {
  for (const auto &xsec : xsecs) {
    std::unordered_set<size_t> c_indices;
    for (const auto &p : xsec.points) {
      c_indices.emplace(p.stroke_ind);
    }
    for (auto s_idx : selected_strokes) {
      if (c_indices.count(s_idx)) {
        for (auto idx : c_indices)
          if (!selected_strokes.count(idx))
            overlapping_set.emplace(idx);
      }
    }
  }
}

static bool is_cluster_continuous(const Cluster &cluster,
                                  const std::set<size_t> &selected_strokes) {
  std::map<size_t, std::vector<std::pair<double, double>>> s2intervals;
  for (const auto &s : cluster.strokes) {
    if (!selected_strokes.count(s.stroke_ind))
      continue;

    assert(s.u.size() > 1);

    std::map<int, size_t> inc_dec_counter;
    inc_dec_counter[-1] = 0;
    inc_dec_counter[1] = 0;
    for (size_t i = 0; i < s.u.size(); ++i) {
      if (i > 0) {
        inc_dec_counter[(s.u[i] > s.u[i - 1]) ? 1 : -1]++;
      }
    }
    bool is_inc = inc_dec_counter[1] > inc_dec_counter[-1];
    for (size_t i = 0; i < s.u.size(); ++i) {
      // Generate two intervals if the cluster is periodic and the stroke across
      // the connected endpoint (this is determined by seeing a change that
      // doesn't agree with the general increasing/decreasing trend.
      if (i > 0 && is_inc != (s.u[i] > s.u[i - 1])) {
        double end_u =
          (s.u[i - 1] > cluster.max_u() * 0.5) ? cluster.max_u() : 0;
        std::pair<double, double> interval = s2intervals[s.stroke_ind].back();
        s2intervals[s.stroke_ind].back() = std::pair<double, double>(
          std::min(interval.first, end_u), std::max(interval.first, end_u));
        double start_u = (s.u[i] > cluster.max_u() * 0.5) ? cluster.max_u() : 0;
        std::pair<double, double> new_interval(start_u, start_u);
        s2intervals[s.stroke_ind].emplace_back(new_interval);
      } else if (s2intervals[s.stroke_ind].empty()) {
        std::pair<double, double> interval(s.u[i], s.u[i]);
        s2intervals[s.stroke_ind].emplace_back(interval);
      } else { // Update existing interval
        std::pair<double, double> interval = s2intervals[s.stroke_ind].back();
        s2intervals[s.stroke_ind].back() = std::pair<double, double>(
          std::min(interval.first, s.u[i]), std::max(interval.first, s.u[i]));
      }
    }
  }

  // Check if intervals are overlapping
  // Or if there is a single discontinuity, check if the cluster is periodic and
  // strokes cross the connected end.
  std::vector<std::pair<double, double>> sorted_interval;
  for (const auto &s_intervals : s2intervals) {
    sorted_interval.insert(sorted_interval.end(), s_intervals.second.begin(),
                           s_intervals.second.end());
  }

  std::sort(sorted_interval.begin(), sorted_interval.end(),
            [](const std::pair<double, double> &a,
               const std::pair<double, double> &b) -> bool {
              return a.first < b.first;
            });
  size_t discontinuity_count = 0;
  size_t end_count = 0;
  for (size_t i = 1; i < sorted_interval.size(); ++i) {
    double max_right = -1;
    for (size_t j = 0; j < i; ++j) {
      max_right = std::max(max_right, sorted_interval[j].second);
    }
    if (max_right < sorted_interval[i].first)
      discontinuity_count++;
    if (sorted_interval[i].first == 0)
      end_count++;
    if (sorted_interval[i].second == cluster.max_u())
      end_count++;
  }

  return (discontinuity_count == 0) ||
         (cluster.periodic && discontinuity_count == 1 && end_count == 2);
}

bool check_overlap_two_cluster(const Input &input, Cluster cluster,
                               Cluster cluster_sub) {
  Cluster cluster1, cluster2, merged_cluster;

  reuse_cluster(cluster, cluster1);
  reuse_cluster(cluster_sub, cluster2);

  FeatureVector fea(input.width, input.height);
  bool reorient = true;
  fea.compute_angle(cluster1, cluster2, input.thickness, merged_cluster,
                    nullptr, reorient);
  if (fea.features_.front() >= std::numeric_limits<double>::infinity())
    return false;
  return true;
}

bool check_angle_two_cluster(const Input &input, Cluster cluster,
                             Cluster cluster_sub) {
  Cluster cluster1, cluster2, merged_cluster;

  reuse_cluster(cluster, cluster1);
  reuse_cluster(cluster_sub, cluster2);

  FeatureVector fea(input.width, input.height);
  bool reorient = true;
  fea.compute_angle(cluster1, cluster2, input.thickness, merged_cluster,
                    nullptr, reorient);
  for (size_t i = 0; i < fea.descriptions_.size(); ++i) {
    if (fea.descriptions_[i] == "average_angle") {
      if (fea.features_[i] > fit_angle_threshold) {
        return false;
      }
    }
  }
  return true;
}

void generate_stroke_positive_samples(const Input &input,
                                      std::vector<Sample> &pos_stroke_samples) {
  std::set<size_t> selected_strokes;
  for (const auto &cluster : input.clusters) {
    // The parameterization fails
    if (cluster.second.fit.centerline.empty())
      continue;

    std::vector<int> cluster_stroke_ind;
    for (const auto &stroke : cluster.second.strokes) {
      cluster_stroke_ind.emplace_back(stroke.stroke_ind);
    }
    std::sort(cluster_stroke_ind.begin(), cluster_stroke_ind.end());
    for (int i = 1; i < cluster_stroke_ind.size(); i++) {
      std::vector<int> pre_cluster_stroke_ind;
      for (int j = 0; j < i; j++) {
        pre_cluster_stroke_ind.emplace_back(cluster_stroke_ind[j]);
      }
      std::vector<int> target_stroke_ind(1, cluster_stroke_ind[i]);

      // 0. Trivial distance filter
      if (filter_sample(cluster.second, pre_cluster_stroke_ind, cluster.second,
                        target_stroke_ind))
        continue;

      // 1. Check if the subcluster is continuous
      selected_strokes.clear();
      selected_strokes.insert(pre_cluster_stroke_ind.begin(),
                              pre_cluster_stroke_ind.end());
      bool is_subc_continuous =
        is_cluster_continuous(cluster.second, selected_strokes);
      if (!is_subc_continuous)
        continue;

      // 2. Check if the two parts overlap in the GT final cluster
      std::set<size_t> overlapping_set;
      find_overlapping_strokes(cluster.second.xsecs,
                               std::set<size_t>{(size_t)cluster_stroke_ind[i]},
                               overlapping_set);
      bool is_overlapping_eventually = false;
      for (auto s_idx : pre_cluster_stroke_ind) {
        if (overlapping_set.count(s_idx)) {
          is_overlapping_eventually = true;
          break;
        }
      }
      if (!is_overlapping_eventually)
        continue;

      pos_stroke_samples.emplace_back(
        std::make_pair(pre_cluster_stroke_ind, target_stroke_ind));
    }
  }
}

void generate_stroke_negative_samples(
  const Input &input, std::vector<Sample> &neg_stroke_samples,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  std::map<int, int> stroke_ind_to_cluster_ind;
  std::set<size_t> selected_strokes;
  for (const auto &cluster : input.clusters) {
    // The parameterization fails
    if (cluster.second.fit.centerline.empty())
      continue;
    for (const auto &stroke : cluster.second.strokes) {
      stroke_ind_to_cluster_ind[stroke.stroke_ind] = stroke.cluster_ind;
    }
  }

  for (const auto &ind_pair : stroke_ind_to_cluster_ind) {
    std::map<int, std::vector<int>> pre_dict_cluster_stroke_ind;
    for (const auto &ind_pair_sub : stroke_ind_to_cluster_ind) {
      if (ind_pair.first == ind_pair_sub.first)
        break;
      if (ind_pair.second == ind_pair_sub.second)
        continue;
      pre_dict_cluster_stroke_ind[ind_pair_sub.second].emplace_back(
        ind_pair_sub.first);
    }
    std::vector<int> target_stroke_ind(1, ind_pair.first);

    for (const auto &pre_cluster_stroke_ind : pre_dict_cluster_stroke_ind) {
      // 0. Trivial distance filter
      if (filter_sample(input.clusters.at(pre_cluster_stroke_ind.first),
                        pre_cluster_stroke_ind.second,
                        input.clusters.at(ind_pair.second),
                        std::vector<int>{ind_pair.first}))
        continue;

      // 1. Check if the cluster is continuous
      selected_strokes.clear();
      selected_strokes.insert(pre_cluster_stroke_ind.second.begin(),
                              pre_cluster_stroke_ind.second.end());
      bool is_subc_continuous = is_cluster_continuous(
        input.clusters.at(pre_cluster_stroke_ind.first), selected_strokes);
      if (!is_subc_continuous)
        continue;

      // 2. Check if two parts overlap
      Cluster cluster1, cluster2, merged_cluster;
      reuse_cluster(input.clusters.at(pre_cluster_stroke_ind.first), cluster1);
      reuse_cluster(input.clusters.at(ind_pair.second), cluster2);
      FeatureVector fea(input.width, input.height);
      fea.compute(cluster1, cluster2, input.thickness, merged_cluster, false,
                  matching_pair_ptr);
      if (fea.features_.front() >= std::numeric_limits<double>::infinity())
        continue;

      neg_stroke_samples.emplace_back(
        std::make_pair(pre_cluster_stroke_ind.second, target_stroke_ind));
    }
  }
}

void generate_cluster_positive_samples(
  const Input &input, std::vector<Sample> &pos_cluster_samples) {
  std::set<size_t> selected_strokes;
  for (const auto &cluster : input.clusters) {
    // The parameterization fails
    if (cluster.second.fit.centerline.empty())
      continue;

    std::vector<int> cluster_stroke_ind;
    for (const auto &stroke : cluster.second.strokes) {
      cluster_stroke_ind.emplace_back(stroke.stroke_ind);
    }
    std::sort(cluster_stroke_ind.begin(), cluster_stroke_ind.end());
    for (int i = 1; i < cluster_stroke_ind.size(); i++) {
      std::vector<int> pre_cluster_stroke_ind;
      std::vector<int> post_cluster_stroke_ind;
      for (int j = 0; j < i; j++) {
        pre_cluster_stroke_ind.emplace_back(cluster_stroke_ind[j]);
      }
      for (int j = i; j < cluster_stroke_ind.size(); j++) {
        post_cluster_stroke_ind.emplace_back(cluster_stroke_ind[j]);
      }

      if (post_cluster_stroke_ind.size() == 1)
        continue;

      // 0. Trivial distance filter
      if (filter_sample(cluster.second, pre_cluster_stroke_ind, cluster.second,
                        post_cluster_stroke_ind))
        continue;

      // 1. Check if the two subclusters are both continuous
      selected_strokes.clear();
      selected_strokes.insert(pre_cluster_stroke_ind.begin(),
                              pre_cluster_stroke_ind.end());
      bool is_subc_continuous =
        is_cluster_continuous(cluster.second, selected_strokes);
      if (!is_subc_continuous)
        continue;
      selected_strokes.clear();
      selected_strokes.insert(post_cluster_stroke_ind.begin(),
                              post_cluster_stroke_ind.end());
      is_subc_continuous =
        is_cluster_continuous(cluster.second, selected_strokes);
      if (!is_subc_continuous)
        continue;

      // 2. Check if the two parts overlap in the GT final cluster
      std::set<size_t> overlapping_set;
      find_overlapping_strokes(cluster.second.xsecs, selected_strokes,
                               overlapping_set);
      bool is_overlapping_eventually = false;
      for (auto s_idx : pre_cluster_stroke_ind) {
        if (overlapping_set.count(s_idx)) {
          is_overlapping_eventually = true;
          break;
        }
      }
      if (!is_overlapping_eventually)
        continue;

      pos_cluster_samples.emplace_back(
        std::make_pair(pre_cluster_stroke_ind, post_cluster_stroke_ind));
    }
  }
}

void generate_cluster_negative_samples(
  const Input &input, std::vector<Sample> &neg_cluster_samples,
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
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

      // If the single stroke cluster is after all strokes in the other cluster,
      // then we've seen this example
      auto is_duplicate = [](int s_idx,
                             const std::vector<int> &stroke_indices) -> bool {
        for (auto idx : stroke_indices)
          if (idx > s_idx)
            return false;
        return true;
      };

      if (cluster_stroke_ind.size() == 1 &&
          is_duplicate(cluster_stroke_ind[0], sub_cluster_stroke_ind))
        continue;
      if (sub_cluster_stroke_ind.size() == 1 &&
          is_duplicate(sub_cluster_stroke_ind[0], cluster_stroke_ind))
        continue;

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
