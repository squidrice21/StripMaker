#include "SolveUtil.h"

#include <algorithm>
#include <unordered_map>

FeatureVector inf_feature(int width, int height) {
  bool has_global = !solve_context.cluster_features.empty();

  FeatureVector fea(width, height);
  const auto feature_calculators = get_stroke_cluster_features();
  const auto global_feature_calculators = get_secondary_cluster_features();
  fill_inf(feature_calculators, fea.descriptions_, fea.features_);
  if (has_global)
    fill_inf(global_feature_calculators, fea.descriptions_, fea.features_);
  fea.descriptions_.emplace_back("velocity");
  fea.features_.emplace_back(-1);
  fea.descriptions_.emplace_back("alignment");
  fea.features_.emplace_back(-1);

  return fea;
}

std::vector<double> pick_features(const FeatureVector &fea,
                                  const std::vector<std::string> &feature_set,
                                  std::vector<std::string> &feature_set_out) {
  std::vector<double> fea_vec;
  std::unordered_map<std::string, double> fea_dict;
  for (size_t i = 0; i < fea.features_.size(); ++i) {
    assert(!fea_dict.count(fea.descriptions_[i]));
    fea_dict[fea.descriptions_[i]] = fea.features_[i];
  }

  assert(!feature_set.empty());
  for (const auto &fea : feature_set) {
    assert(fea_dict.count(fea));
    fea_vec.emplace_back(fea_dict[fea]);
    feature_set_out.emplace_back(fea);
  }

  return fea_vec;
}

FeatureVector compute_feature(Cluster &cluster1, Cluster &cluster2,
                              Cluster &merged_cluster, int width, int height,
                              double norm_thickness, std::string vis_folder) {
  bool has_global = !solve_context.cluster_features.empty();

  FeatureVector fea(width, height);
  // bool reorient = false;
  bool reorient = true;
  if (!has_global) {
    fea.compute(cluster1, cluster2, norm_thickness, merged_cluster, true,
                &solve_context.matching_pair, reorient, vis_folder);
  } else {
    fea.compute_secondary(cluster1, cluster2, norm_thickness, merged_cluster,
                          solve_context.cluster_features,
                          &solve_context.matching_pair, reorient, vis_folder);
  }
  // Add the parameterization info as well.
  for (const auto &k_v : merged_cluster.obj_term_values) {
    fea.descriptions_.emplace_back(k_v.first);
    fea.features_.emplace_back(k_v.second);
  }

  return fea;
}

void print_csv_header(const std::vector<std::string> &feature_set, bool typical,
                      std::string &csv_str) {
  for (const auto &f_str : feature_set) {
    auto f = f_str;
    std::replace(f.begin(), f.end(), ' ', '_');
    csv_str += f + ",";
  }
  csv_str += "score,";
  if (typical)
    csv_str += "typical,";
  csv_str += "base_number,";
  csv_str += "target_number,";
  csv_str += "is_merging,";
  csv_str += "decision,";
  csv_str += "fit_output,";
  csv_str += "isoline";
  csv_str += "\n";
}

void print_csv_row(const std::vector<std::string> &feature_set,
                   const std::vector<double> &feature_vec, double score,
                   double typical, int base_number, int target_number,
                   bool is_merging, bool decision, std::string fit_output,
                   std::string isoline, std::string &csv_str) {
  if (csv_str.empty())
    print_csv_header(feature_set, (typical >= 0), csv_str);
  for (auto v : feature_vec) {
    csv_str += std::to_string(v) + ",";
  }
  csv_str += std::to_string(score) + ",";
  if (typical >= 0)
    csv_str += std::to_string(typical) + ",";
  csv_str += std::to_string(base_number) + ",";
  csv_str += std::to_string(target_number) + ",";
  csv_str += std::to_string(is_merging) + ",";
  csv_str += std::to_string(decision) + ",";
  csv_str += fit_output + ",";
  csv_str += isoline;
  csv_str += "\n";
}

bool to_merge(const Cluster &cluster1, const Cluster &cluster2, double score) {
  auto contains_single = [](const Cluster &cluster) {
    for (auto const &s : cluster.strokes)
      if (solve_context.final_single_sid.count(s.stroke_ind))
        return true;
    return false;
  };
  bool contains_single1 = contains_single(cluster1);
  bool contains_single2 = contains_single(cluster2);

  if (solve_context.stage == SolveContext::SolveStage::Split_Merging) {
    if (cluster1.strokes.size() == 1 && cluster2.strokes.size() == 1) {
      return score >= solve_context.split_merge_cutoffs[2];
    } else if ((cluster1.strokes.size() == 1 && cluster2.strokes.size() > 1) ||
               (cluster1.strokes.size() > 1 && cluster2.strokes.size() == 1)) {
      return score >= solve_context.split_merge_cutoffs[0];
    } else {
      return score >= solve_context.split_merge_cutoffs[1];
    }
  }

  if (solve_context.stage == SolveContext::SolveStage::Merging &&
      ((contains_single1 && !contains_single2) ||
       (!contains_single1 && contains_single2))) {
    return score >= solve_context.single_final_cutoff;
  } else if (solve_context.stage == SolveContext::SolveStage::Merging &&
             ((cluster1.strokes.size() == 1 && cluster2.strokes.size() == 1))) {
    return score >= solve_context.single_single_final_cutoff;
  } else {
    return score >= solve_context.cutoff;
  }
}
