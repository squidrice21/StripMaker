#include "SecondaryFeatureImplementations.h"

#include "../classifier/forest.h"
#include "../classifier/forest_sec.h"
#include "../solve/SolveUtil.h"

double get_feature_percentile(std::vector<FeatureVector> const &features,
                              double percentile, const std::string &fea_name) {
  assert(!features.empty());
  std::map<std::string, size_t> fea_map;
  for (size_t i = 0; i < features.front().descriptions_.size(); ++i) {
    fea_map[features.front().descriptions_[i]] = i;
  }

  std::vector<double> measure_vec;
  measure_vec.reserve(features.size());
  for (auto &fea : features) {
    assert(fea_map.count(fea_name));
    size_t i = fea_map[fea_name];
    const auto measures = fea.features_[i];
    measure_vec.emplace_back(measures);
  }
  if (measure_vec.empty())
    return std::numeric_limits<double>::infinity();

  assert(percentile >= 0 && percentile <= 1);
  int pos = (measure_vec.size() - 1) * percentile;
  pos = std::max(pos, 0);
  auto p = measure_vec.begin() + pos;
  std::nth_element(measure_vec.begin(), p, measure_vec.end());
  return measure_vec[pos];
}

double get_feature_mean(std::vector<FeatureVector> const &features,
                        const std::string &fea_name) {
  assert(!features.empty());
  std::map<std::string, size_t> fea_map;
  for (size_t i = 0; i < features.front().descriptions_.size(); ++i) {
    fea_map[features.front().descriptions_[i]] = i;
  }

  std::vector<double> measure_vec;
  measure_vec.reserve(features.size());
  for (auto &fea : features) {
    assert(fea_map.count(fea_name));
    size_t i = fea_map[fea_name];
    const auto measures = fea.features_[i];
    measure_vec.emplace_back(measures);
  }
  if (measure_vec.empty())
    return std::numeric_limits<double>::infinity();
  double sum = 0;
  for (const auto v : measure_vec)
    sum += v;
  return sum / measure_vec.size();
}

std::vector<double> InitialProbabilityStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  Cluster cluster1, cluster2;
  for (size_t i = 0; i < merged_cluster.strokes.size(); ++i) {
    if (merged_cluster.strokes[i].cluster_ind == 0) {
      cluster1.strokes.emplace_back(merged_cluster.strokes[i]);
      cluster1.original_input_strokes.emplace_back(
        merged_cluster.original_input_strokes[i]);
    } else {
      cluster2.strokes.emplace_back(merged_cluster.strokes[i]);
      cluster2.original_input_strokes.emplace_back(
        merged_cluster.original_input_strokes[i]);
    }
  }
  double norm_thickness = 1;
  Cluster merged_cluster_tmp = merged_cluster;
  FeatureVector fea = compute_feature(cluster1, cluster2, merged_cluster_tmp,
                                      -1, -1, norm_thickness, "");

  double alignment = -1;
  if (merged_cluster.obj_term_values.count("alignment"))
    alignment = merged_cluster.obj_term_values.at("alignment");

  // Must be successful since we provide the common parameterization
  assert(alignment >= 0);

  std::vector<std::string> feature_set;
  double prob = -1;
  if (solve_context.stage == SolveContext::SolveStage::Initial) {
    std::vector<double> fea_vec =
      pick_features(fea, solve_context.feature_set, feature_set);
    normalize(fea_vec);
    prob = clf(fea_vec);
  } else if (solve_context.stage == SolveContext::SolveStage::Splitting) {
    std::vector<double> fea_vec =
      pick_features(fea, solve_context.sec_feature_set, feature_set);
    normalize_sec(fea_vec);
    prob = clf_sec(fea_vec);
  }

  assert(prob >= 0);
  return std::vector<double>{prob};
}

std::vector<double> MeanGlobalClusterFeature::operator()(
  const FeatureVector &base_feature,
  const std::vector<FeatureVector> &features) const {
  std::vector<double> fea_vec;
  for (const auto &fea : feas_) {
    double angle = get_feature_mean(features, fea);
    fea_vec.emplace_back(angle);
  }
  return fea_vec;
}

std::vector<double> PercentileGlobalClusterFeature::operator()(
  const FeatureVector &base_feature,
  const std::vector<FeatureVector> &features) const {
  std::vector<double> fea_vec;
  for (const auto &fea : feas_) {
    double angle = get_feature_percentile(features, percentile_, fea);
    fea_vec.emplace_back(angle);
  }
  return fea_vec;
}

std::vector<double> MeanGlobalRatioClusterFeature::operator()(
  const FeatureVector &base_feature,
  const std::vector<FeatureVector> &features) const {
  std::map<std::string, double> fea_map;
  for (size_t i = 0; i < base_feature.descriptions_.size(); ++i) {
    fea_map[base_feature.descriptions_[i]] = base_feature.features_[i];
  }

  std::vector<double> fea_vec;
  for (const auto &b_fea : base_feas_) {
    assert(fea_map.count(b_fea));
    double numerator = fea_map[b_fea];
    for (const auto &fea : feas_) {
      double angle = get_feature_mean(features, fea);
      fea_vec.emplace_back(numerator / (angle + epsilon));
    }
  }
  return fea_vec;
}

std::vector<double> PercentileGlobalRatioClusterFeature::operator()(
  const FeatureVector &base_feature,
  const std::vector<FeatureVector> &features) const {
  std::map<std::string, double> fea_map;
  for (size_t i = 0; i < base_feature.descriptions_.size(); ++i) {
    fea_map[base_feature.descriptions_[i]] = base_feature.features_[i];
  }

  std::vector<double> fea_vec;
  for (const auto &b_fea : base_feas_) {
    assert(fea_map.count(b_fea));
    double numerator = fea_map[b_fea];
    for (const auto &fea : feas_) {
      double angle = get_feature_percentile(features, percentile_, fea);
      fea_vec.emplace_back(numerator / (angle + epsilon));
    }
  }
  return fea_vec;
}
