#include "OutlierRemoval.h"

#include "../Logger.h"
#include "../Util.h"
#include "../feature/Filters.h"
#include "../solve/SolveUtil.h"
#include "Uncutting.h"
#include "glm/detail/func_geometric.hpp"
#include "glm/detail/type_vec.hpp"

#include <algorithm>
#include <limits>
#include <set>
#include <vector>

double envelope_ratio = 1.5;
double cut_envelope_ratio = 2;
double close_envelope_ratio = 1;
double envelope_coverage = 0.9;
double cut_envelope_coverage = 0.8;
double close_envelope_coverage = 0.7;

double long_stroke_length = 30;
double long_envelope_ratio = 1.2;
double long_envelope_coverage = 0.8;

size_t max_stroke_num = 3;

bool is_point_close(const glm::dvec2 center, const glm::dvec2 p,
                    const Cluster &cluster, double ratio) {
  // Find the closest point on fitting
  double min_dist = std::numeric_limits<double>::infinity();
  size_t closest_i = 0;
  for (size_t i = 0; i < cluster.fit.centerline.size(); ++i) {
    double d = glm::distance(cluster.fit.centerline[i], p);
    if (d < min_dist) {
      min_dist = d;
      closest_i = i;
    }
  }
  min_dist /= solve_context.input_thickness;

  // Xsec
  glm::dvec2 transformed_p =
    cluster.fit.fit_sample_xsecs[closest_i].points.front().point;
  transformed_p *= solve_context.input_thickness;
  transformed_p += center;
  glm::dvec2 env_p =
    (glm::dot(p - cluster.fit.centerline[closest_i],
              transformed_p - cluster.fit.centerline[closest_i]) > 0)
      ? cluster.fit.fit_sample_xsecs[closest_i].points.front().point
      : cluster.fit.fit_sample_xsecs[closest_i].points.back().point;
  env_p *= solve_context.input_thickness;
  env_p += center;
  double env_d = glm::distance(env_p, cluster.fit.centerline[closest_i]) /
                 solve_context.input_thickness;
  env_d += 1;
  if (min_dist < ratio * env_d)
    return true;

  // Fit width
  assert(!cluster.fit.widths.empty());
  double w = cluster.fit.widths[closest_i];
  w /= solve_context.input_thickness;
  if (min_dist < ratio * 0.5 * w)
    return true;

  double max_w =
    *std::max_element(cluster.fit.widths.begin(), cluster.fit.widths.end());
  max_w /= solve_context.input_thickness;
  if (min_dist < ratio * 0.5 * max_w)
    return true;
  return false;
}

bool is_close(const glm::dvec2 center, const Cluster &single_cluster,
              const Cluster &cluster, bool is_cut,
              const std::vector<Junction> &junctions) {
  // Determine if it has a junction outside
  bool outside_junc = false;
  for (auto const &junc : junctions) {
    if ((junc.from.first == single_cluster.fit.cluster_idx &&
         junc.to.first != single_cluster.fit.cluster_idx) ||
        (junc.to.first == single_cluster.fit.cluster_idx &&
         junc.from.first != single_cluster.fit.cluster_idx)) {
      size_t other_cid = (junc.from.first == single_cluster.fit.cluster_idx)
                           ? junc.to.first
                           : junc.from.first;
      if (cluster.fit.cluster_idx == other_cid)
        continue;
      outside_junc = true;
      break;
    }
  }

  std::set<size_t> c_sid;
  for (auto const &s : cluster.strokes)
    c_sid.emplace(s.stroke_ind);
  for (auto const &s : single_cluster.strokes) {
    for (const auto &p : solve_context.matching_pair) {
      if ((p.first.first == s.stroke_ind && !c_sid.count(p.second.first)) ||
          (p.second.first == s.stroke_ind && !c_sid.count(p.first.first))) {
        outside_junc = true;
        break;
      }
    }
  }

  SketchUI::Polyline2D fit_poly;
  fit_poly.from_fits(single_cluster.fit, 0);

  // reparameterize
  double rate = 1;
  size_t min_num_samples = 50;
  fit_poly.reparameterize(
    std::min(rate, fit_poly.totalLen() / min_num_samples));

  double e_ratio = envelope_ratio;
  double coverage = envelope_coverage;
  if (is_cut) {
    e_ratio = cut_envelope_ratio;
    coverage = cut_envelope_coverage;
  }
  if (outside_junc) {
    e_ratio = close_envelope_ratio;
  }

  double length = total_length(single_cluster.fit.centerline);
  length /= solve_context.input_thickness;
  if (length > long_stroke_length) {
    e_ratio = long_envelope_ratio;
    coverage = long_envelope_coverage;
  }

  if (length > 5 && !is_cut) {
    Cluster fit_cluster1, fit_cluster2, fit_merged_cluster;
    reuse_cluster(single_cluster, fit_cluster1, true);
    reuse_cluster(cluster, fit_cluster2, true);
    FeatureVector fea(0, 0);
    bool reorient = true;
    double norm_thickness = 1.0;
    fea.compute_angle(fit_cluster1, fit_cluster2, norm_thickness,
                      fit_merged_cluster, nullptr, reorient, "");
    std::vector<std::string> avg_dist_str;
    auto max_angles =
      pick_features(fea, std::vector<std::string>{"max_angle"}, avg_dist_str);
    double max_angle = max_angles[0];
    if (max_angle < std::numeric_limits<double>::infinity() && max_angle > 80) {
      e_ratio = close_envelope_ratio;
      coverage = envelope_coverage;
    }
  }

  double close_p_num = 0;
  for (auto &p : fit_poly.points) {
    if (is_point_close(center, glm::dvec2(p.first.x, p.first.y), cluster,
                       e_ratio)) {
      close_p_num++;
    }
  }

  double ratio = close_p_num / fit_poly.points.size();
  bool loose_coverage = ratio > coverage;
  if (loose_coverage)
    return true;

  if (outside_junc)
    return false;

  close_p_num = 0;
  for (auto &p : fit_poly.points) {
    if (is_point_close(center, glm::dvec2(p.first.x, p.first.y), cluster,
                       close_envelope_ratio)) {
      close_p_num++;
    }
  }
  ratio = close_p_num / fit_poly.points.size();
  bool close_coverage = ratio > close_envelope_coverage;
  return close_coverage;
}

bool is_tight_close(const glm::dvec2 center, const Cluster &single_cluster,
                    const Cluster &cluster, bool is_cut,
                    const std::vector<Junction> &junctions) {
  // Determine if it has a junction outside
  bool outside_junc = false;
  for (auto const &junc : junctions) {
    if ((junc.from.first == single_cluster.fit.cluster_idx &&
         junc.to.first != single_cluster.fit.cluster_idx) ||
        (junc.to.first == single_cluster.fit.cluster_idx &&
         junc.from.first != single_cluster.fit.cluster_idx)) {
      size_t other_cid = (junc.from.first == single_cluster.fit.cluster_idx)
                           ? junc.to.first
                           : junc.from.first;
      if (cluster.fit.cluster_idx == other_cid)
        continue;
      outside_junc = true;
      break;
    }
  }

  std::set<size_t> c_sid;
  for (auto const &s : cluster.strokes)
    c_sid.emplace(s.stroke_ind);
  for (auto const &s : single_cluster.strokes) {
    for (const auto &p : solve_context.matching_pair) {
      if ((p.first.first == s.stroke_ind && !c_sid.count(p.second.first)) ||
          (p.second.first == s.stroke_ind && !c_sid.count(p.first.first))) {
        outside_junc = true;
        break;
      }
    }
  }

  SketchUI::Polyline2D fit_poly;
  fit_poly.from_fits(single_cluster.fit, 0);

  // reparameterize
  double rate = 1;
  size_t min_num_samples = 50;
  fit_poly.reparameterize(
    std::min(rate, fit_poly.totalLen() / min_num_samples));

  double e_ratio = envelope_ratio;
  double coverage = envelope_coverage;
  if (is_cut) {
    e_ratio = cut_envelope_ratio;
    coverage = cut_envelope_coverage;
  }
  if (outside_junc) {
    e_ratio = close_envelope_ratio;
  }

  double length = total_length(single_cluster.fit.centerline);
  length /= solve_context.input_thickness;
  if (length > long_stroke_length) {
    e_ratio = long_envelope_ratio;
    coverage = long_envelope_coverage;
  }

  if (length > 5 && !is_cut) {
    Cluster fit_cluster1, fit_cluster2, fit_merged_cluster;
    reuse_cluster(single_cluster, fit_cluster1, true);
    reuse_cluster(cluster, fit_cluster2, true);
    FeatureVector fea(0, 0);
    bool reorient = true;
    double norm_thickness = 1.0;
    fea.compute_angle(fit_cluster1, fit_cluster2, norm_thickness,
                      fit_merged_cluster, nullptr, reorient, "");
    std::vector<std::string> avg_dist_str;
    auto max_angles =
      pick_features(fea, std::vector<std::string>{"max_angle"}, avg_dist_str);
    double max_angle = max_angles[0];
    if (max_angle < std::numeric_limits<double>::infinity() && max_angle > 80) {
      e_ratio = close_envelope_ratio;
      coverage = envelope_coverage;
    }
  }

  double close_p_num = 0;
  for (auto &p : fit_poly.points) {
    if (is_point_close(center, glm::dvec2(p.first.x, p.first.y), cluster,
                       close_envelope_ratio)) {
      close_p_num++;
    }
  }
  double ratio = close_p_num / fit_poly.points.size();
  bool close_coverage = ratio > close_envelope_coverage;
  return close_coverage;
}

bool remove_outliers(const Input &input,
                     const std::map<size_t, Cluster> &in_clusters,
                     const std::map<int, FittedCurve> &in_fits,
                     std::map<size_t, Cluster> &clusters,
                     const std::vector<Junction> &junctions, bool skip_cut) {
  std::set<size_t> single_stroke_clusters;
  for (auto const &cc : in_clusters) {
    if (cc.second.strokes.size() <= max_stroke_num)
      single_stroke_clusters.emplace(cc.first);
  }
  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : in_clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      sid2s[s.stroke_ind] = s;
    }
  }

  // 1. Check if the single stroke cluster is part of a cut stroke
  auto find_cut_cluster = [&in_clusters, &sid2s](size_t in_cid,
                                                 size_t sid) -> int {
    for (auto const &cc : in_clusters) {
      if (cc.second.strokes.size() == 1 || cc.first == in_cid)
        continue;
      std::set<size_t> c_sid;
      for (auto const &s : cc.second.strokes)
        c_sid.emplace(s.stroke_ind);
      for (const auto &p : solve_context.matching_pair) {
        if ((p.first.first == sid && c_sid.count(p.second.first)) ||
            (p.second.first == sid && c_sid.count(p.first.first))) {
          size_t other_sid =
            (p.first.first == sid) ? p.second.first : p.first.first;
          if (total_length(sid2s[sid].points) <
              total_length(sid2s[other_sid].points))
            return cc.first;
          else
            return -1;
        }
      }
    }
    return -1;
  };
  auto is_cycle_cut = [&in_clusters](size_t in_cid, size_t sid) -> int {
    for (auto const &cc : in_clusters) {
      if (cc.second.strokes.size() == 1 || cc.first == in_cid)
        continue;
      std::set<size_t> c_sid;
      for (auto const &s : cc.second.strokes)
        c_sid.emplace(s.stroke_ind);
      for (const auto &p : solve_context.matching_pair) {
        if (p.first.first == sid && c_sid.count(p.second.first)) {
          for (const auto &p1 : solve_context.matching_pair) {
            if (p1.second.first == sid && p1.first.first == p.second.first) {
              return cc.first;
            }
          }
        }
        if (p.second.first == sid && c_sid.count(p.first.first)) {
          for (const auto &p1 : solve_context.matching_pair) {
            if (p1.first.first == sid && p1.second.first == p.first.first) {
              return cc.first;
            }
          }
        }
      }
    }
    return -1;
  };
  std::map<size_t, size_t> single_cut_pairs;
  std::set<size_t> skip_cid;
  for (auto single_c : single_stroke_clusters) {
    for (auto const &s : in_clusters.at(single_c).strokes) {
      int cid = find_cut_cluster(single_c, s.stroke_ind);
      int cyc_cid = is_cycle_cut(single_c, s.stroke_ind);
      if (cyc_cid >= 0) {
        skip_cid.emplace(single_c);
        break;
      }
      if (cid >= 0) {
        single_cut_pairs[single_c] = cid;
        break;
      }
    }
  }

  // 2. Make other close single-multi pairs
  std::map<size_t, size_t> single_multi_pairs;
  for (auto single_c : single_stroke_clusters) {
    if (skip_cid.count(single_c))
      continue;
    if (!skip_cut && single_cut_pairs.count(single_c)) {
      single_multi_pairs[single_c] = single_cut_pairs[single_c];
      continue;
    }

    // Min hausdorff
    std::map<double, size_t> s_c_dist;
    for (auto const &cc : in_clusters) {
      if (cc.second.strokes.size() == 1 || cc.first == single_c)
        continue;
      if (skip_cut && single_cut_pairs.count(single_c) &&
          single_cut_pairs[single_c] == cc.first)
        continue;
      if (skip_cut && !single_cut_pairs.count(single_c))
        continue;
      if (filter_sample(in_clusters.at(single_c), std::vector<int>(), cc.second,
                        std::vector<int>()))
        continue;
      double d = hausdorff_distance_one_sided(in_clusters.at(single_c).fit,
                                              cc.second.fit);
      Cluster fit_cluster1, fit_cluster2, fit_merged_cluster;
      reuse_cluster(in_clusters.at(single_c), fit_cluster1, true);
      reuse_cluster(cc.second, fit_cluster2, true);

      FeatureVector fea(0, 0);
      bool reorient = true;
      double norm_thickness = 1.0;
      fea.compute_distance(fit_cluster1, fit_cluster2, norm_thickness,
                           fit_merged_cluster, nullptr, reorient, "");
      std::vector<std::string> avg_dist_str;
      auto avg_dists = pick_features(
        fea, std::vector<std::string>{"average_distance"}, avg_dist_str);
      double avg_dist = avg_dists[0];

      double dist = std::min(avg_dist, d);
      dist /= solve_context.input_thickness;

      s_c_dist[dist] = cc.first;
    }
    double max_w = 0;
    if (!s_c_dist.empty()) {
      max_w = *std::max_element(
        in_clusters.at(s_c_dist.begin()->second).fit.widths.begin(),
        in_clusters.at(s_c_dist.begin()->second).fit.widths.end());
      max_w /= solve_context.input_thickness;
    }
    if (!s_c_dist.empty() &&
        (s_c_dist.begin()->first < pointwise_distance_threshold ||
         s_c_dist.begin()->first < 5 * max_w)) {
      single_multi_pairs[single_c] = s_c_dist.begin()->second;
    }
  }

  // 3. Check if the single stroke cluster is within a cluster envelope
  if (skip_cut)
    single_cut_pairs.clear();
  std::map<size_t, std::vector<size_t>> to_merge_singles;
  std::set<size_t> tight_cid;
  for (auto const &sc : single_multi_pairs) {
    bool is_cut = single_cut_pairs.count(sc.first) &&
                  (single_cut_pairs[sc.first] == sc.second);
    if (is_close(input.orig_center, in_clusters.at(sc.first),
                 in_clusters.at(sc.second), is_cut, junctions)) {
      to_merge_singles[sc.second].emplace_back(sc.first);
      if (is_tight_close(input.orig_center, in_clusters.at(sc.first),
                         in_clusters.at(sc.second), is_cut, junctions))
        tight_cid.emplace(sc.first);
    }
  }

  // 4. Merge
  auto try_merging = [&in_clusters, &input,
                      &tight_cid](const std::set<size_t> &c_indices,
                                  Cluster &cluster, size_t base_cid) -> bool {
    assert(!c_indices.empty());
    if (c_indices.size() == 1) {
      cluster = in_clusters.at(*c_indices.begin());
      return true;
    }
    for (auto cid : c_indices) {
      cluster.strokes.insert(cluster.strokes.end(),
                             in_clusters.at(cid).strokes.begin(),
                             in_clusters.at(cid).strokes.end());
      cluster.original_input_strokes.insert(
        cluster.original_input_strokes.end(),
        in_clusters.at(cid).original_input_strokes.begin(),
        in_clusters.at(cid).original_input_strokes.end());
      cluster.periodic |= in_clusters.at(cid).periodic;
    }

    bool to_orient = true;
    bool test_time = true;
    {
      // Turn cutting off
      context.ignore_endpoint = true;
      // fit_cluster(solve_context.norm_thickness, solve_context.width,
      //             solve_context.height, cluster, to_orient, test_time,
      //             &solve_context.matching_pair);
      fit_cluster(solve_context.norm_thickness, solve_context.width,
                  solve_context.height, cluster, to_orient, test_time, nullptr);
      context.ignore_endpoint = false;
      if (cluster.strokes.empty()) {
        fit_cluster(solve_context.norm_thickness, solve_context.width,
                    solve_context.height, cluster, to_orient, test_time,
                    &solve_context.matching_pair);
      }
    }
    if (cluster.strokes.empty())
      return false;

    // Fit with thickness
    bool in_width = context.widths;
    Input input_final;
    context.widths = true;
    input_final.clusters[0] = cluster;
    cluster.fit = fea_fitting.fit(&input_final).begin()->second;

    transform_fit_curve(cluster.fit, input.orig_center,
                        solve_context.input_thickness);

    context.widths = in_width;

    double hausdorff = -1;
    double max_w = -1;
    double input_length = 0;

    double new_length = total_length(cluster.fit.centerline);
    double new_w =
      *std::max_element(cluster.fit.widths.begin(), cluster.fit.widths.end());
    double add_length = 0;
    bool is_tight = false;
    for (auto cid : c_indices) {
      double w = *std::max_element(in_clusters.at(cid).fit.widths.begin(),
                                   in_clusters.at(cid).fit.widths.end());
      double d =
        hausdorff_distance_one_sided(in_clusters.at(cid).fit, cluster.fit);
      hausdorff = std::max(hausdorff, d);
      max_w = std::max(max_w, w);
      input_length += total_length(in_clusters.at(cid).fit.centerline);
      if (cid != base_cid) {
        add_length = std::max(add_length,
                              total_length(in_clusters.at(cid).fit.centerline));
        is_tight |= tight_cid.count(cid);
      }
      double width_change = std::max(new_w, max_w) / std::min(new_w, max_w);
      logger().info("{} => {}, {}, {} => {}", cid, hausdorff,
                    max_w / solve_context.input_thickness, new_w, width_change);
    }

    double width_change = std::max(new_w, max_w) / std::min(new_w, max_w);
    add_length /= solve_context.input_thickness;
    if (!is_tight && add_length > 8 && hausdorff > add_length)
      return false;

    return (hausdorff < 2 && width_change < 1.3) || (width_change < 1.1);
  };
  std::set<size_t> merged_cids;
  for (auto const &c_singles : to_merge_singles) {
    Cluster last_cluster = in_clusters.at(c_singles.first);
    std::vector<size_t> c_indices;
    for (size_t i = 0; i < c_singles.second.size(); ++i) {
      std::set<size_t> c_indices_all;
      for (auto cid : c_indices) {
        c_indices_all.emplace(cid);
      }
      c_indices_all.emplace(c_singles.second[i]);
      c_indices_all.emplace(c_singles.first);

      Cluster updated_cluster;
      bool is_merge_ok =
        try_merging(c_indices_all, updated_cluster, c_singles.first);

      if (is_merge_ok &&
          last_cluster.strokes.size() != updated_cluster.strokes.size()) {
        updated_cluster.fit.cluster_idx = c_singles.first;
        last_cluster = updated_cluster;
        c_indices.emplace_back(c_singles.second[i]);
        for (auto cid : c_indices) {
          merged_cids.emplace(cid);
        }
      }
    }

    // Save to results
    clusters[c_singles.first] = last_cluster;
  }

  bool has_cut_left = false;

  for (auto const &cc : in_clusters) {
    if (!clusters.count(cc.first) && !merged_cids.count(cc.first)) {
      has_cut_left |= single_cut_pairs.count(cc.first);
      clusters[cc.first] = cc.second;
    }
  }

  return has_cut_left;
}
