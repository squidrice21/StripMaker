#include "JunctionEndMerge.h"

#include "../Util.h"
#include "../measure/Measurement.h"
#include "../solve/SolveUtil.h"
#include "../uncutting/Uncutting.h"
#include "Junction.h"
#include "JunctionProposal.h"
#include "JunctionSerialization.h"

#include "../external/SketchConnectivity/features/junction_features_impl.h"
#include "../stroke_strip_src/Serialization.h"
#include "../stroke_strip_src/SketchInfo.h"

#include "SketchConnectivity/classifier.h"
#include "glm/detail/type_vec.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

double end_end_angle = 20;
double end_end_threshold = 0.9;

void junction_end_merge(const Input &input, std::map<size_t, Cluster> &clusters,
                        std::vector<Junction> &junctions) {
  std::vector<FittedCurve> fits;
  prepare_fits(clusters, fits);
  size_t knn = 1;
  sketching::StrokeBVH bvh;
  construct_bvh(fits, bvh);

  std::map<size_t, sketching::Stroke const *> sid2s;
  for (auto const &s : bvh.strokes()) {
    sid2s[s.index] = &s;
  }

  std::map<int, FittedCurve> src_fits;
  std::map<int, FittedCurve> vis_fits;
  std::vector<Junction> cur_junctions;
  std::vector<Junction> cur_no_t_junctions;
  for (auto const &fit : fits) {
    cur_junctions.clear();
    cur_no_t_junctions.clear();
    bool is_head = true;
    predict_junctions(fit, is_head, knn, bvh, cur_junctions);
    for (const auto &junc : cur_junctions) {
      if (junc.type == sketching::JunctionType::X &&
          !(is_roughly_end(*sid2s[junc.from.first], junc.from.second) &&
            is_roughly_end(*sid2s[junc.to.first], junc.to.second)))
        continue;
      if (junc.type != sketching::JunctionType::T ||
          (is_roughly_end(*sid2s[junc.from.first], junc.from.second) &&
           is_roughly_end(*sid2s[junc.to.first], junc.to.second))) {
        cur_no_t_junctions.emplace_back(junc);
        if (cur_no_t_junctions.back().type == sketching::JunctionType::T ||
            cur_no_t_junctions.back().type == sketching::JunctionType::X) {
          // cur_no_t_junctions.back().type = sketching::JunctionType::R;
          cur_no_t_junctions.back().from.second =
            (cur_no_t_junctions.back().from.second > 0.5) ? 1 : 0;
          cur_no_t_junctions.back().to.second =
            (cur_no_t_junctions.back().to.second > 0.5) ? 1 : 0;
          auto v1 =
            sid2s[junc.from.first]->pos(cur_no_t_junctions.back().from.second *
                                        sid2s[junc.from.first]->length());
          cur_no_t_junctions.back().from_pos = glm::dvec2(v1.x(), v1.y());
          auto v2 =
            sid2s[junc.to.first]->pos(cur_no_t_junctions.back().to.second *
                                      sid2s[junc.to.first]->length());
          cur_no_t_junctions.back().to_pos = glm::dvec2(v2.x(), v2.y());
        }
      }
    }
    std::sort(cur_no_t_junctions.begin(), cur_no_t_junctions.end(),
              [](const Junction &a, const Junction &b) {
                if (a.type == sketching::JunctionType::X ||
                    b.type == sketching::JunctionType::X)
                  return glm::distance(a.from_pos, a.to_pos) <
                         glm::distance(b.from_pos, b.to_pos);
                return a.probability > b.probability;
              });

    for (auto const &junc : cur_no_t_junctions) {
      if (junc.probability >= end_end_threshold &&
          (junc.type == sketching::JunctionType::R ||
           glm::distance(junc.from_pos, junc.to_pos) <
             3 * std::min(sid2s[junc.from.first]->pen_width(),
                          sid2s[junc.to.first]->pen_width()))) {
        junctions.emplace_back(junc);
      }
    }

    cur_junctions.clear();
    cur_no_t_junctions.clear();
    is_head = false;
    predict_junctions(fit, is_head, knn, bvh, cur_junctions);
    for (const auto &junc : cur_junctions) {
      if (junc.type == sketching::JunctionType::X &&
          !(is_roughly_end(*sid2s[junc.from.first], junc.from.second) &&
            is_roughly_end(*sid2s[junc.to.first], junc.to.second)))
        continue;
      if (junc.type != sketching::JunctionType::T ||
          (is_roughly_end(*sid2s[junc.from.first], junc.from.second) &&
           is_roughly_end(*sid2s[junc.to.first], junc.to.second))) {
        cur_no_t_junctions.emplace_back(junc);
        if (cur_no_t_junctions.back().type == sketching::JunctionType::T ||
            cur_no_t_junctions.back().type == sketching::JunctionType::X) {
          cur_no_t_junctions.back().from.second =
            (cur_no_t_junctions.back().from.second > 0.5) ? 1 : 0;
          cur_no_t_junctions.back().to.second =
            (cur_no_t_junctions.back().to.second > 0.5) ? 1 : 0;
          auto v1 =
            sid2s[junc.from.first]->pos(cur_no_t_junctions.back().from.second *
                                        sid2s[junc.from.first]->length());
          cur_no_t_junctions.back().from_pos = glm::dvec2(v1.x(), v1.y());
          auto v2 =
            sid2s[junc.to.first]->pos(cur_no_t_junctions.back().to.second *
                                      sid2s[junc.to.first]->length());
          cur_no_t_junctions.back().to_pos = glm::dvec2(v2.x(), v2.y());
        }
      }
    }
    std::sort(cur_no_t_junctions.begin(), cur_no_t_junctions.end(),
              [](const Junction &a, const Junction &b) {
                if (a.type == sketching::JunctionType::X ||
                    b.type == sketching::JunctionType::X)
                  return glm::distance(a.from_pos, a.to_pos) <
                         glm::distance(b.from_pos, b.to_pos);
                return a.probability > b.probability;
              });
    for (auto const &junc : cur_no_t_junctions) {
      if (junc.probability >= end_end_threshold &&
          (junc.type == sketching::JunctionType::R ||
           glm::distance(junc.from_pos, junc.to_pos) <
             3 * std::min(sid2s[junc.from.first]->pen_width(),
                          sid2s[junc.to.first]->pen_width()))) {
        junctions.emplace_back(junc);
      }
    }
  }
}

Cluster::Stroke fit2stroke(const FittedCurve &fit) {
  Cluster::Stroke s;
  s.cluster_ind = fit.cluster_idx;
  s.stroke_ind = fit.cluster_idx;

  Sketch reparam_s;
  reparam_s.from_fits(fit, 0);
  {
    // double rate = 2;
    double rate = stroke_sampling_size;
    // reparam_s.reparameterize(std::min(rate, reparam_s.totalLen() / 5));
    reparam_s.reparameterize(
      std::min(rate, reparam_s.totalLen() / stroke_sampling_min_num));
  }

  s.points.reserve(reparam_s.points.size());
  for (const auto &p : reparam_s.points) {
    s.points.emplace_back(p.first.x, p.first.y);
  }

  return s;
}

bool is_continuation(const FittedCurve &fit1, const FittedCurve &fit2,
                     const Junction &junc, double &angle) {
  sketching::Stroke s1, s2;
  convert_stroke(fit1, s1);
  convert_stroke(fit2, s2);
  s1.ensure_arclengths();
  s2.ensure_arclengths();

  double w1 = s1.width_at((junc.from.first == fit1.cluster_idx)
                            ? junc.from.second * s1.length()
                            : junc.to.second * s1.length());
  double w2 = s2.width_at((junc.from.first == fit2.cluster_idx)
                            ? junc.from.second * s2.length()
                            : junc.to.second * s2.length());
  if (std::max(w1, w2) / std::min(w1, w2) > 4)
    return false;

  Cluster::Stroke ss1, ss2;
  ss1 = fit2stroke(fit1);
  ss2 = fit2stroke(fit2);

  double curv1 = std::abs(stroke_total_signed_curvature(ss1));
  double curv2 = std::abs(stroke_total_signed_curvature(ss2));
  if (curv1 >= (360. - 30) / 180 * M_PI || curv2 >= (360. - 30) / 180 * M_PI)
    return false;

  double endend_dist = std::max(std::max(s1.pen_width(), s2.pen_width()),
                                glm::distance(junc.from_pos, junc.to_pos));
  const auto tangent1 = sketching::features::stepaway_tangent(
    s1,
    (fit1.cluster_idx == junc.from.first) ? junc.from.second < 0.5
                                          : junc.to.second < 0.5,
    endend_dist);
  const auto tangent2 = sketching::features::stepaway_tangent(
    s2,
    (fit2.cluster_idx == junc.from.first) ? junc.from.second < 0.5
                                          : junc.to.second < 0.5,
    endend_dist);
  if (tangent1.dot(tangent2) > 0)
    return false;
  angle = std::acos(std::clamp(std::abs(tangent1.dot(tangent2)), 0.0, 1.0));
  angle *= 180 / M_PI;
  return angle < end_end_angle;
}

void resolve_continuation(const std::vector<double> &junc_angles,
                          std::vector<Junction> &in_junctions) {
  auto is_duplicate = [](const Junction &in_junc,
                         const std::vector<Junction> &results) {
    for (auto junc : results) {
      if ((junc.from == in_junc.from && junc.to == in_junc.to) ||
          (junc.from == in_junc.to && junc.to == in_junc.from))
        return true;
    }
    return false;
  };
  std::vector<Junction> junctions;
  for (auto const &junc : in_junctions)
    if (!is_duplicate(junc, junctions))
      junctions.emplace_back(junc);

  std::vector<Junction> results;
  results.reserve(junctions.size());

  auto find_branch = [&](const sketching::StrokeTime &st,
                         std::vector<size_t> &branch_junc) {
    std::vector<Junction> dedup;
    for (size_t j = 0; j < junctions.size(); ++j) {
      if (is_duplicate(junctions[j], dedup))
        continue;
      if (junctions[j].from == st || junctions[j].to == st) {
        branch_junc.emplace_back(j);
        dedup.emplace_back(junctions[j]);
      }
    }
  };
  std::set<size_t> removed_i;
  for (size_t i = 0; i < junctions.size(); ++i) {
    if (removed_i.count(i))
      continue;
    // Find branching
    std::vector<size_t> branch_junc;
    find_branch(junctions[i].from, branch_junc);
    assert(branch_junc.size() >= 1);
    double min_angle = std::numeric_limits<double>::infinity();
    size_t min_i = i;
    for (auto j : branch_junc) {
      if (!removed_i.count(j) && junc_angles[j] < min_angle) {
        min_angle = junc_angles[j];
        min_i = j;
      }
    }
    for (auto j : branch_junc) {
      if (j != min_i) {
        removed_i.emplace(j);
      }
    }

    if (!is_duplicate(junctions[min_i], results)) {
      results.emplace_back(junctions[min_i]);
    }

    branch_junc.clear();
    find_branch(junctions[i].to, branch_junc);
    assert(branch_junc.size() >= 1);
    min_angle = std::numeric_limits<double>::infinity();
    min_i = i;
    for (auto j : branch_junc) {
      if (!removed_i.count(j) && junc_angles[j] < min_angle) {
        min_angle = junc_angles[j];
        min_i = j;
      }
    }
    for (auto j : branch_junc) {
      if (j != min_i) {
        removed_i.emplace(j);
      }
    }

    if (!is_duplicate(junctions[min_i], results))
      results.emplace_back(junctions[min_i]);
  }

  in_junctions = results;
}

void merge_ends(const Input &input, std::map<size_t, Cluster> &clusters,
                std::map<size_t, Cluster> &merged_clusters,
                const std::vector<Junction> &junctions) {
  // 1. Into connected components
  std::map<size_t, size_t> components;
  for (auto const &cc : clusters) {
    components[cc.first] = components.size();
  }

  auto color_comp = [&components](size_t from, size_t to) {
    for (auto &cc : components) {
      if (cc.second == from)
        cc.second = to;
    }
  };
  bool changed = false;
  do {
    changed = false;
    for (auto const &junc : junctions) {
      if (components[junc.from.first] != components[junc.to.first]) {
        changed = true;
        size_t from =
          std::max(components[junc.from.first], components[junc.to.first]);
        size_t to =
          std::min(components[junc.from.first], components[junc.to.first]);
        color_comp(from, to);
      }
    }
  } while (changed);
  std::map<size_t, std::set<size_t>> c_components;
  for (auto &cc : components) {
    c_components[cc.second].emplace(cc.first);
  }

  // 2. Try merging
  auto try_merging = [&clusters, &input](const std::set<size_t> &c_indices,
                                         Cluster &cluster) -> bool {
    assert(!c_indices.empty());
    if (c_indices.size() == 1) {
      cluster = clusters[*c_indices.begin()];
      return true;
    }
    for (auto cid : c_indices) {
      cluster.strokes.insert(cluster.strokes.end(),
                             clusters[cid].strokes.begin(),
                             clusters[cid].strokes.end());
      cluster.original_input_strokes.insert(
        cluster.original_input_strokes.end(),
        clusters[cid].original_input_strokes.begin(),
        clusters[cid].original_input_strokes.end());
      cluster.periodic |= clusters[cid].periodic;
    }

    bool to_orient = true;
    bool test_time = true;
    fit_cluster(solve_context.norm_thickness, solve_context.width,
                solve_context.height, cluster, to_orient, test_time,
                &solve_context.matching_pair);
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
    for (auto cid : c_indices) {
      double w = *std::max_element(clusters[cid].fit.widths.begin(),
                                   clusters[cid].fit.widths.end());
      double d = hausdorff_distance_one_sided(clusters[cid].fit, cluster.fit);
      hausdorff = std::max(hausdorff, d);
      max_w = std::max(max_w, w);
      input_length += total_length(clusters[cid].fit.centerline);
    }

    double new_length = total_length(cluster.fit.centerline);
    double new_w =
      *std::max_element(cluster.fit.widths.begin(), cluster.fit.widths.end());
    double width_change = std::max(new_w, max_w) / std::min(new_w, max_w);

    return hausdorff < max_w && width_change < 2 &&
           (new_length - input_length) < 3;
  };
  std::map<size_t, Cluster> results;
  {
    results.clear();
    merged_clusters.clear();
    for (auto &cc : c_components) {
      // For debug output
      if (cc.second.size() == 1) {
        size_t new_c = results.size();
        results[new_c] = clusters[*cc.second.begin()];
        results[new_c].fit.cluster_idx = new_c;
        continue;
      }
      Cluster cluster;
      bool succeeded = try_merging(cc.second, cluster);
      if (!succeeded && cc.second.size() > 2) {
        for (auto c_remove : cc.second) {
          Cluster tmp_cluster;
          std::set<size_t> c_indices;

          for (auto c : cc.second) {
            if (c == c_remove)
              continue;
            c_indices.emplace(c);
          }

          succeeded = try_merging(c_indices, tmp_cluster);
          if (succeeded) {
            size_t new_c = results.size();
            results[new_c] = clusters[c_remove];
            results[new_c].fit.cluster_idx = new_c;

            size_t new_c2 = merged_clusters.size();
            merged_clusters[new_c2] = results[new_c];
            merged_clusters[new_c].fit.cluster_idx = new_c2;

            cluster = tmp_cluster;
            break;
          }
        }
      }
      if (succeeded) {
        size_t new_c = results.size();
        results[new_c] = cluster;
        results[new_c].fit.cluster_idx = new_c;

        size_t new_c2 = merged_clusters.size();
        merged_clusters[new_c2] = results[new_c];
        merged_clusters[new_c].fit.cluster_idx = new_c2;
      } else {
        for (auto cid : cc.second) {
          size_t new_c = results.size();
          results[new_c] = clusters[cid];
          results[new_c].fit.cluster_idx = new_c;

          size_t new_c2 = merged_clusters.size();
          merged_clusters[new_c2] = results[new_c];
          merged_clusters[new_c].fit.cluster_idx = new_c2;
        }
      }
    }
  }

  clusters = results;
}
