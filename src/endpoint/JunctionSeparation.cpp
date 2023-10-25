#include "JunctionSeparation.h"

#include "../Logger.h"
#include "../Util.h"
#include "../feature/FeatureIO.h"
#include "../solve/SolveEndMerge.h"
#include "../solve/SolveEndSplit.h"
#include "../solve/SolveLazy.h"
#include "../solve/SolveUtil.h"
#include "../stroke_strip_src/Serialization.h"
#include "Junction.h"
#include "JunctionMatching.h"
#include "JunctionSerialization.h"

bool separate_sub_cluster_probability(Cluster &c1, Cluster &c2,
                                      Cluster &merged_cluster,
                                      std::string prefix) {
  double norm_thickness = 1;
  std::vector<std::pair<double, PredictionInfo>> scores;
  bool later_cluster_only = false;
  bool in_place = false;
  std::map<size_t, Cluster> clusters;
  clusters[c1.fit.cluster_idx] = c1;
  clusters[c2.fit.cluster_idx] = c2;
  compare_clusters(c1, c1.fit.cluster_idx, clusters, later_cluster_only,
                   in_place, &merged_cluster, solve_context.width,
                   solve_context.height, norm_thickness, scores, "");

  std::sort(scores.begin(), scores.end(),
            [](const std::pair<double, PredictionInfo> &a,
               const std::pair<double, PredictionInfo> &b) {
              return std::abs(a.first) > std::abs(b.first);
            });
  if (scores.empty()) {
    logger().info(prefix + " Try merging {} and {} with p={}",
                  std::min(c1.fit.cluster_idx, c2.fit.cluster_idx),
                  std::max(c1.fit.cluster_idx, c2.fit.cluster_idx), 0);
    return false;
  } else {
    logger().info(prefix + " Try merging {} and {} with p={}",
                  std::min(c1.fit.cluster_idx, c2.fit.cluster_idx),
                  std::max(c1.fit.cluster_idx, c2.fit.cluster_idx),
                  std::abs(scores.front().first));

    return std::abs(scores.front().first) < solve_context.junc_cutoff;
  }
}

bool are_separable_clusters(const Input &input,
                            std::map<size_t, Cluster> clusters,
                            std::map<size_t, Cluster> &sep_clusters,
                            Cluster merged_cluster,
                            const std::string &vis_folder) {
  assert(sep_clusters.size() == 2);

  // 1. Prepare fit curves
  auto fit_width_cluster = [&input](std::map<size_t, Cluster> &clusters) {
    context.tighter_fit = true;
    for (auto &cc : clusters) {
      cc.second.fit.cluster_idx = cc.first;

      if (!cc.second.fit.widths.empty() &&
          cc.second.fit.widths.size() == cc.second.fit.centerline.size())
        continue;
      cc.second.fit.widths.clear();
      cc.second.fit.centerline.clear();
      cc.second.fit.fit_sample_xsecs.clear();
      bool to_orient = true;
      bool test_time = true;
      fit_cluster(solve_context.norm_thickness, solve_context.width,
                  solve_context.height, cc.second, to_orient, test_time,
                  &solve_context.matching_pair);

      // Fit with thickness
      bool in_width = context.widths;
      Input input_final;
      context.widths = true;
      input_final.clusters[0] = cc.second;
      cc.second.fit = fea_fitting.fit(&input_final).begin()->second;

      transform_fit_curve(cc.second.fit, input.orig_center,
                          solve_context.input_thickness);

      context.widths = in_width;
    }
    context.tighter_fit = false;
  };
  fit_width_cluster(clusters);
  fit_width_cluster(sep_clusters);
  std::map<size_t, Cluster> tmp_merged_cluster;
  tmp_merged_cluster[merged_cluster.fit.cluster_idx] = merged_cluster;
  merged_cluster.fit.centerline.clear();
  merged_cluster.fit.widths.clear();
  merged_cluster.fit.fit_sample_xsecs.clear();
  fit_width_cluster(tmp_merged_cluster);
  merged_cluster = tmp_merged_cluster[merged_cluster.fit.cluster_idx];

  std::vector<FittedCurve> fits;
  prepare_fits(clusters, fits);
  size_t knn = 1;
  sketching::StrokeBVH bvh;
  construct_bvh(fits, bvh);
  std::vector<FittedCurve> sep_fits;
  prepare_fits(sep_clusters, sep_fits);
  std::vector<Junction> junctions;

  std::vector<FittedCurve> all_fits = fits;
  std::map<int, FittedCurve> vis_fits;
  for (auto const &fit : sep_fits) {
    vis_fits[fit.cluster_idx] = fit;
    all_fits.emplace_back(fit);
    bool is_head = true;
    predict_junctions(fit, is_head, knn, bvh, junctions);
    is_head = false;
    predict_junctions(fit, is_head, knn, bvh, junctions);
  }

  // Filter out the junctions between the two
  {
    std::vector<Junction> junctions_tmp;
    for (auto const &junc : junctions) {
      if (!((junc.from.first == sep_fits[0].cluster_idx &&
             junc.to.first == sep_fits[1].cluster_idx) ||
            (junc.from.first == sep_fits[1].cluster_idx &&
             junc.to.first == sep_fits[0].cluster_idx)))
        junctions_tmp.emplace_back(junc);
    }
    junctions = std::move(junctions_tmp);
  }

  bool to_separate = false;
  {
    std::vector<Junction> decision_junctions;
    bool is_separate = is_separate_forward_check(
      sep_clusters.begin()->second, sep_clusters.rbegin()->second,
      merged_cluster, all_fits, junctions, decision_junctions);
    if (!is_separate)
      is_separate = is_separate_forward_check(
        sep_clusters.rbegin()->second, sep_clusters.begin()->second,
        merged_cluster, all_fits, junctions, decision_junctions);
    if (is_separate) {
      bool is_separate_prob = separate_sub_cluster_probability(
        sep_clusters.begin()->second, sep_clusters.rbegin()->second,
        merged_cluster, "Sep final:");
      if (is_separate_prob) {
        to_separate = true;
      }
    }
    if (!vis_folder.empty()) {
      context.widths = true;
      Input out_input = input;
      std::ofstream refit_svg(vis_folder + "/sep2_fit_width_final.svg");
      out_input.thickness = solve_context.input_thickness;
      fea_fitting.fit_colored_svg(
        refit_svg, out_input, vis_fits,
        solve_context.width * solve_context.input_thickness,
        solve_context.height * solve_context.input_thickness);
      context.widths = false;
    }
  }

  return to_separate;
}

void separable_clusters(const Input &input,
                        std::map<size_t, Cluster> &in_init_clusters,
                        std::map<size_t, Cluster> &in_sec_clusters,
                        std::vector<std::pair<size_t, size_t>> &sep_cid,
                        const std::string &vis_folder) {
  std::map<size_t, Cluster> init_clusters = in_init_clusters,
                            sec_clusters = in_sec_clusters;
  // 1. Prepare fit curves
  auto fit_width_cluster = [&input](std::map<size_t, Cluster> &clusters) {
    context.tighter_fit = true;
    for (auto &cc : clusters) {
      bool to_orient = true;
      bool test_time = true;
      fit_cluster(solve_context.norm_thickness, solve_context.width,
                  solve_context.height, cc.second, to_orient, test_time,
                  &solve_context.matching_pair);

      // Fit with thickness
      bool in_width = context.widths;
      Input input_final;
      context.widths = true;
      input_final.clusters[0] = cc.second;
      cc.second.fit = fea_fitting.fit(&input_final).begin()->second;
      cc.second.fit.cluster_idx = cc.first;

      transform_fit_curve(cc.second.fit, input.orig_center,
                          solve_context.input_thickness);

      context.widths = in_width;
    }
    context.tighter_fit = false;
  };
  fit_width_cluster(init_clusters);
  fit_width_cluster(sec_clusters);

  // 2. Prepare the input clusters (the sub-clusters split from the init
  // clusters)
  std::map<size_t, size_t> split2init;
  std::map<size_t, std::set<size_t>> init2split;
  {
    std::map<size_t, size_t> s2c_init;
    for (auto const &cc : init_clusters) {
      for (auto const &s : cc.second.strokes) {
        s2c_init[s.stroke_ind] = cc.first;
      }
    }
    for (auto const &cc : sec_clusters) {
      for (auto const &s : cc.second.strokes) {
        split2init[cc.first] = s2c_init[s.stroke_ind];
        init2split[s2c_init[s.stroke_ind]].emplace(cc.first);
      }
    }
  }

  std::vector<FittedCurve> fits;
  prepare_fits(sec_clusters, fits);
  size_t knn = 1;
  sketching::StrokeBVH bvh;
  construct_bvh(fits, bvh);

  std::vector<Junction> junctions;
  std::map<int, FittedCurve> src_fits;
  std::map<int, FittedCurve> vis_fits;
  for (auto const &fit : fits) {
    vis_fits[fit.cluster_idx] = fit;

    size_t init_cid = split2init[fit.cluster_idx];
    if (init2split[init_cid].size() <= 1)
      continue;
    std::vector<FittedCurve> other_fits;
    for (auto const &fit2 : fits) {
      size_t init_c = split2init[fit2.cluster_idx];
      if (init_c != init_cid)
        other_fits.emplace_back(fit2);
    }
    src_fits[fit.cluster_idx] = fit;

    bool is_head = true;
    predict_junctions(fit, is_head, knn, bvh, junctions);
    is_head = false;
    predict_junctions(fit, is_head, knn, bvh, junctions);
  }

  // Filter out the junctions within the same init cluster
  {
    std::vector<Junction> junctions_tmp;
    for (auto const &junc : junctions) {
      if (split2init[junc.to.first] != split2init[junc.from.first])
        junctions_tmp.emplace_back(junc);
    }
    junctions = std::move(junctions_tmp);
  }

  if (!vis_folder.empty()) {
    context.widths = true;
    Input out_input = input;
    std::ofstream refit_svg(vis_folder + "/junction_split_fit_width.svg");
    out_input.thickness = solve_context.input_thickness;
    fea_fitting.fit_svg(refit_svg, out_input, vis_fits,
                        solve_context.width * solve_context.input_thickness,
                        solve_context.height * solve_context.input_thickness);

    std::ofstream fit_svg(vis_folder + "/junction_source_fit_width.svg");
    out_input.thickness = solve_context.input_thickness;
    fea_fitting.fit_colored_svg(
      fit_svg, out_input, src_fits,
      solve_context.width * solve_context.input_thickness,
      solve_context.height * solve_context.input_thickness);
    context.widths = false;
  }

  {
    std::map<int, FittedCurve> sep_fits;
    std::vector<Junction> decision_junctions;
    for (auto const &fit : fits) {
      size_t init_cid = split2init[fit.cluster_idx];
      for (auto other_cid : init2split[init_cid]) {
        if (other_cid == fit.cluster_idx)
          continue;
        bool is_separate = is_separate_forward_check(
          sec_clusters[fit.cluster_idx], sec_clusters[other_cid],
          init_clusters[init_cid], fits, junctions, decision_junctions);
        if (is_separate) {
          bool is_separate_prob = separate_sub_cluster_probability(
            in_sec_clusters[fit.cluster_idx], in_sec_clusters[other_cid],
            in_init_clusters[init_cid], "Sep2:");
          if (is_separate_prob) {
            sep_cid.emplace_back(std::min(fit.cluster_idx, other_cid),
                                 std::max(fit.cluster_idx, other_cid));
            sep_fits[fit.cluster_idx] = vis_fits[fit.cluster_idx];
            sep_fits[other_cid] = vis_fits[other_cid];
          }
          break;
        }
      }
    }
    if (!vis_folder.empty()) {
      context.widths = true;
      Input out_input = input;
      std::ofstream refit_svg(vis_folder + "/sep2_fit_width.svg");
      out_input.thickness = solve_context.input_thickness;
      fea_fitting.fit_colored_svg(
        refit_svg, out_input, sep_fits,
        solve_context.width * solve_context.input_thickness,
        solve_context.height * solve_context.input_thickness);
      context.widths = false;
    }
  }
}
