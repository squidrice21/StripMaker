#include "Features.h"

#include "../../stroke_strip_src/Cluster.h"
#include "../../stroke_strip_src/Fitting.h"
#include "../../stroke_strip_src/Serialization.h"
#include "../../stroke_strip_src/Utils.h"
#include "../Capture2Input.h"
#include "../Logger.h"
#include "../Util.h"
#include "../solve/SolveUtil.h"
#include "../stroke_strip_src/StrokeCutting.h"
#include "FeatureImplementations.h"
#include "SecondaryFeatureImplementations.h"
#include "SingleFeatureImplementations.h"

#include "glm/detail/type_vec.hpp"
#include <oneapi/tbb/parallel_for.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

Context context;
Capture2Input capture2input;

StrokeOrientation fea_orientation(context);
Parameterization fea_param(context);
FittingEigenSparse fea_fitting(context);

size_t significant_feature_bins = 3;
size_t feature_bins = 1;

std::vector<std::unique_ptr<StrokeClusterFeature>>
get_stroke_cluster_features() {
  auto features = std::vector<std::unique_ptr<StrokeClusterFeature>>();

  features.push_back(std::make_unique<TypeStrokeClusterFeature>());
  if (significant_feature_bins > 1) {
    features.push_back(std::make_unique<MeanAngleStrokeClusterFeature>(1));
    features.push_back(
      std::make_unique<PercentileAngleStrokeClusterFeature>(0.5, 1));
    features.push_back(std::make_unique<MeanDistanceStrokeClusterFeature>(1));
    features.push_back(
      std::make_unique<PercentileDistanceStrokeClusterFeature>(0.5, 1));
  }
  features.push_back(std::make_unique<OverlappingRatioStrokeClusterFeature>());

  features.push_back(
    std::make_unique<MeanWidthLengthMaxRatioStrokeClusterFeature>(
      feature_bins));
  features.push_back(
    std::make_unique<PercentileWidthLengthMaxRatioStrokeClusterFeature>(
      0.5, feature_bins));
  features.push_back(
    std::make_unique<MeanWidthLengthMinRatioStrokeClusterFeature>(
      feature_bins));
  features.push_back(
    std::make_unique<PercentileWidthLengthMinRatioStrokeClusterFeature>(
      0.5, feature_bins));
  features.push_back(
    std::make_unique<MeanWidthLengthRatioOverallStrokeClusterFeature>(
      feature_bins));
  features.push_back(
    std::make_unique<PercentileWidthLengthRatioOverallStrokeClusterFeature>(
      0.5, feature_bins));
  features.push_back(
    std::make_unique<MeanClusterwiseDistanceMaxWidthRatioStrokeClusterFeature>(
      feature_bins));
  features.push_back(
    std::make_unique<
      PercentileClusterwiseDistanceMaxWidthRatioStrokeClusterFeature>(
      0.5, feature_bins));
  features.push_back(
    std::make_unique<MeanClusterwiseDistanceMinWidthRatioStrokeClusterFeature>(
      feature_bins));
  features.push_back(
    std::make_unique<
      PercentileClusterwiseDistanceMinWidthRatioStrokeClusterFeature>(
      0.5, feature_bins));

  features.push_back(
    std::make_unique<MeanClusterwiseDistanceRatioStrokeClusterFeature>(
      feature_bins));
  features.push_back(
    std::make_unique<PercentileClusterwiseDistanceRatioStrokeClusterFeature>(
      0.5, feature_bins));

  features.push_back(
    std::make_unique<
      PercentileGapClusterwiseDistanceMaxRatioStrokeClusterFeature>(
      0.9, feature_bins));
  features.push_back(
    std::make_unique<
      PercentileGapClusterwiseDistanceMinRatioStrokeClusterFeature>(
      0.9, feature_bins));
  features.push_back(
    std::make_unique<
      PercentileGapClusterwiseDistanceMaxRatioStrokeClusterFeature>(
      0.1, feature_bins));
  features.push_back(
    std::make_unique<
      PercentileGapClusterwiseDistanceMinRatioStrokeClusterFeature>(
      0.1, feature_bins));

  features.push_back(
    std::make_unique<MaxWidthChangeRatioStrokeClusterFeature>());
  features.push_back(
    std::make_unique<MinWidthChangeRatioStrokeClusterFeature>());
  features.push_back(
    std::make_unique<MinAvgWidthChangeRatioStrokeClusterFeature>());
  features.push_back(
    std::make_unique<MaxAvgWidthChangeRatioStrokeClusterFeature>());

  features.push_back(std::make_unique<OverlappingLengthStrokeClusterFeature>());

  features.push_back(
    std::make_unique<MinPointwiseDistanceStrokeClusterFeature>());

  features.push_back(
    std::make_unique<OverlappingIsolineNumberStrokeClusterFeature>());

  features.push_back(
    std::make_unique<InvMaxAvgWidthChangeRatioStrokeClusterFeature>());
  features.push_back(
    std::make_unique<InvMinAvgWidthChangeRatioStrokeClusterFeature>());

  return features;
}

static std::vector<std::unique_ptr<StrokeClusterFeature>>
get_stroke_cluster_angle_features() {
  auto features = std::vector<std::unique_ptr<StrokeClusterFeature>>();

  features.push_back(std::make_unique<MeanAngleStrokeClusterFeature>(1));
  features.push_back(
    std::make_unique<PercentileAngleStrokeClusterFeature>(0.5, 1));
  features.push_back(
    std::make_unique<PercentileAngleStrokeClusterFeature>(1.0, 1));

  return features;
}

static std::vector<std::unique_ptr<StrokeClusterFeature>>
get_stroke_cluster_distance_features() {
  auto features = std::vector<std::unique_ptr<StrokeClusterFeature>>();

  features.push_back(std::make_unique<MeanDistanceStrokeClusterFeature>(1));
  features.push_back(
    std::make_unique<PercentileDistanceStrokeClusterFeature>(0.5, 1));

  return features;
}

std::vector<std::unique_ptr<StrokeClusterFeature>>
get_stroke_cluster_probability_features() {
  auto features = std::vector<std::unique_ptr<StrokeClusterFeature>>();

  features.push_back(
    std::make_unique<InitialProbabilityStrokeClusterFeature>());

  return features;
}

std::vector<std::unique_ptr<SecondaryClusterFeature>>
get_secondary_cluster_features() {
  auto features = std::vector<std::unique_ptr<SecondaryClusterFeature>>();

  features.push_back(std::make_unique<MeanGlobalRatioClusterFeature>(
    std::vector<std::string>{"average_distance", "median_distance"},
    std::vector<std::string>{"average_distance", "median_distance"}));
  features.push_back(std::make_unique<PercentileGlobalRatioClusterFeature>(
    std::vector<std::string>{"average_distance", "median_distance"},
    std::vector<std::string>{"average_distance", "median_distance"}, 0.5));
  features.push_back(std::make_unique<PercentileGlobalRatioClusterFeature>(
    std::vector<std::string>{"average_distance", "median_distance"},
    std::vector<std::string>{"average_distance", "median_distance"}, 0.9));

  return features;
}

std::vector<std::unique_ptr<TypicalClusterFeature>>
get_typical_cluster_features() {
  auto features = std::vector<std::unique_ptr<TypicalClusterFeature>>();

  features.push_back(
    std::make_unique<MeanAngleTypicalClusterFeature>(feature_bins));
  features.push_back(
    std::make_unique<PercentileAngleTypicalClusterFeature>(0.5, feature_bins));
  features.push_back(
    std::make_unique<MeanDistanceTypicalClusterFeature>(feature_bins));
  features.push_back(std::make_unique<PercentileDistanceTypicalClusterFeature>(
    0.5, feature_bins));
  features.push_back(std::make_unique<PercentileDistanceTypicalClusterFeature>(
    1.0, feature_bins));
  features.push_back(
    std::make_unique<StrokeCountTypicalClusterFeature>(feature_bins));
  features.push_back(
    std::make_unique<MeanXsecStrokeCountTypicalClusterFeature>(feature_bins));
  features.push_back(
    std::make_unique<PercentileXsecStrokeCountTypicalClusterFeature>(
      0.5, feature_bins));
  features.push_back(
    std::make_unique<MeanWidthTypicalClusterFeature>(feature_bins));
  features.push_back(
    std::make_unique<PercentileWidthTypicalClusterFeature>(0.5, feature_bins));
  features.push_back(
    std::make_unique<PercentileWidthTypicalClusterFeature>(1.0, feature_bins));
  features.push_back(
    std::make_unique<MeanGapRatioTypicalClusterFeature>(feature_bins));
  features.push_back(std::make_unique<PercentileGapRatioTypicalClusterFeature>(
    0.5, feature_bins));
  features.push_back(std::make_unique<PercentileGapRatioTypicalClusterFeature>(
    1.0, feature_bins));
  features.push_back(
    std::make_unique<MeanWidthLengthRatioTypicalClusterFeature>(feature_bins));
  features.push_back(
    std::make_unique<PercentileWidthLengthRatioTypicalClusterFeature>(
      0.5, feature_bins));

  return features;
}

std::string example_hash(const std::vector<size_t> &stroke_indices) {
  std::string stroke_str;
  for (auto const &sid : stroke_indices) {
    if (!stroke_str.empty())
      stroke_str += ' ';
    stroke_str += std::to_string(sid);
  }
  std::size_t str_hash = std::hash<std::string>{}(stroke_str);
  return std::to_string(str_hash);
}

std::string example_hash(const std::vector<int> &stroke_indices) {
  std::string stroke_str;
  for (auto const &sid : stroke_indices) {
    if (!stroke_str.empty())
      stroke_str += ' ';
    stroke_str += std::to_string(sid);
  }
  std::size_t str_hash = std::hash<std::string>{}(stroke_str);
  return std::to_string(str_hash);
}

std::string example_hash(const std::vector<size_t> &stroke_indices1,
                         const std::vector<size_t> &stroke_indices2) {
  std::string stroke_str;
  for (auto const &sid : stroke_indices1) {
    if (!stroke_str.empty())
      stroke_str += ' ';
    stroke_str += std::to_string(sid);
  }
  stroke_str += ',';
  for (auto const &sid : stroke_indices2) {
    if (!stroke_str.empty())
      stroke_str += ' ';
    stroke_str += std::to_string(sid);
  }
  std::size_t str_hash = std::hash<std::string>{}(stroke_str);
  return std::to_string(str_hash);
}

std::string get_hash_str(const Cluster &cluster) {
  std::vector<size_t> stroke_indices_comb;
  for (auto const &s : cluster.strokes) {
    stroke_indices_comb.emplace_back(s.stroke_ind);
  }
  std::string hash_str_comb = example_hash(stroke_indices_comb);

  return hash_str_comb;
}

std::string get_hash_str_comb(const Cluster &cluster1,
                              const Cluster &cluster2) {
  std::vector<size_t> stroke_indices_comb;
  for (auto const &s : cluster1.strokes) {
    stroke_indices_comb.emplace_back(s.stroke_ind);
  }
  for (auto const &s : cluster2.strokes) {
    stroke_indices_comb.emplace_back(s.stroke_ind);
  }
  std::sort(stroke_indices_comb.begin(), stroke_indices_comb.end());
  std::string hash_str_comb = example_hash(stroke_indices_comb);

  return hash_str_comb;
}

void fill_inf(
  const std::vector<std::unique_ptr<StrokeClusterFeature>> &feature_calculators,
  std::vector<std::string> &descriptions, std::vector<double> &features) {
  for (auto const &fea : feature_calculators) {
    if (fea->num_binned_ == 1)
      descriptions.emplace_back(fea->csv());
    else {
      for (size_t i = 0; i < fea->num_binned_; ++i) {
        descriptions.emplace_back(fea->csv() + "_b" + std::to_string(i));
      }
    }
    features.insert(features.end(), fea->num_binned_,
                    std::numeric_limits<double>::infinity());
  }
}

void fill_inf(const std::vector<std::unique_ptr<SecondaryClusterFeature>>
                &SecondaryClusterFeature,
              std::vector<std::string> &descriptions,
              std::vector<double> &features) {
  for (auto const &fea : SecondaryClusterFeature) {
    assert(fea->num_binned_ == 1);
    auto csvs = fea->csvs();
    descriptions.insert(descriptions.end(), csvs.begin(), csvs.end());
    features.insert(features.end(), csvs.size(),
                    std::numeric_limits<double>::infinity());
  }
}

static void compute_features(
  const std::vector<std::unique_ptr<StrokeClusterFeature>> &feature_calculators,
  const Cluster &cluster1, const Cluster &cluster2,
  const Cluster &merged_cluster, std::vector<std::string> &descriptions,
  std::vector<double> &features) {
  auto begin = std::chrono::high_resolution_clock::now();

  std::vector<size_t> start_indices;
  start_indices.emplace_back(0);
  size_t fea_size = 0;
  for (size_t i = 1; i < feature_calculators.size(); ++i) {
    auto const &fea = feature_calculators[i - 1];
    start_indices.emplace_back(start_indices[i - 1] + fea->num_binned_);
    fea_size += fea->num_binned_;
  }
  fea_size += feature_calculators.back()->num_binned_;

  descriptions.resize(fea_size);
  features.resize(fea_size);
  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<size_t>(0u, feature_calculators.size()),
    [&](const oneapi::tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        auto const &fea = feature_calculators[i];
        size_t start_i = start_indices[i];
        if (fea->num_binned_ == 1)
          descriptions[start_i] = fea->csv();
        else {
          for (size_t j = 0; j < fea->num_binned_; ++j) {
            descriptions[start_i + j] = fea->csv() + "_b" + std::to_string(j);
          }
        }
        auto fea_vec =
          fea->operator()(merged_cluster, cluster1.fit, cluster2.fit);
        for (size_t j = 0; j < fea_vec.size(); ++j) {
          features[start_i + j] = fea_vec[j];
        }
      }
    });

  auto end = std::chrono::high_resolution_clock::now();
  solve_context.feature_timer +=
    std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() /
    1000.0;
}

static bool merge_cached(
  Cluster &cluster1, Cluster &cluster2, double thickness,
  Cluster &merged_cluster, bool reorient, bool test_time, Input &input,
  std::map<int, FittedCurve> &fit_map, bool &pre_computed,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  assert(cluster1.strokes.size() == cluster1.original_input_strokes.size() &&
         cluster2.strokes.size() == cluster2.original_input_strokes.size());

  // Check if we want to avoid using fitting based orientation to avoid
  // complicated bookkeeping...
  bool no_fit_orient = false;
  if (matching_pair_ptr) {
    std::set<size_t> seen_cut_strokes;
    for (const auto &ss : *matching_pair_ptr) {
      seen_cut_strokes.emplace(ss.first.first);
      seen_cut_strokes.emplace(ss.second.first);
    }
    for (const auto &s : cluster1.strokes)
      if (seen_cut_strokes.count(s.stroke_ind)) {
        no_fit_orient = true;
        break;
      }
    for (const auto &s : cluster2.strokes)
      if (seen_cut_strokes.count(s.stroke_ind)) {
        no_fit_orient = true;
        break;
      }
  }

  reorient = no_fit_orient || test_time;

  Capture capture;
  std::map<size_t, size_t> index_mapping, s2c;
  for (const auto &s : (!no_fit_orient && !test_time)
                         ? cluster1.strokes
                         : cluster1.original_input_strokes) {
    capture.sketchedPolylines.emplace_back(stroke2polyline(s));
    capture.sketchedPolylines.back().stroke_ind =
      capture.sketchedPolylines.size() - 1;
    index_mapping[capture.sketchedPolylines.back().stroke_ind] = s.stroke_ind;
    s2c[s.stroke_ind] = s.cluster_ind;
    capture.sketchedPolylines.back().group_ind = 0;
  }
  for (const auto &s : (!no_fit_orient && !test_time)
                         ? cluster2.strokes
                         : cluster2.original_input_strokes) {
    capture.sketchedPolylines.emplace_back(stroke2polyline(s));
    capture.sketchedPolylines.back().stroke_ind =
      capture.sketchedPolylines.size() - 1;
    index_mapping[capture.sketchedPolylines.back().stroke_ind] = s.stroke_ind;
    s2c[s.stroke_ind] = s.cluster_ind;
    capture.sketchedPolylines.back().group_ind = 0;
  }
  capture.thickness = thickness;
  // Note that we assume the input is already preprocessed so the
  // preprocessing is skipped here.

  input = capture2input.from_capture(capture, false);
  for (auto &c : input.clusters) {
    for (auto &s : c.second.strokes) {
      s.stroke_ind = index_mapping[s.stroke_ind];
    }
    for (auto &s : c.second.original_input_strokes) {
      s.stroke_ind = index_mapping[s.stroke_ind];
      assert(s2c.count(s.stroke_ind));
      s.cluster_ind = s2c[s.stroke_ind];
    }
  }

  // Parameterize the merged cluster
  pre_computed = false;
  if (merged_cluster.strokes.empty()) {
    assert(input.clusters.begin()->second.strokes.size() ==
           capture.sketchedPolylines.size());
    for (size_t i = 0; i < cluster1.strokes.size(); ++i) {
      input.clusters.begin()->second.strokes[i].cluster_ind = 0;
    }
    for (size_t i = cluster1.strokes.size();
         i < cluster1.strokes.size() + cluster2.strokes.size(); ++i) {
      input.clusters.begin()->second.strokes[i].cluster_ind = 1;
    }

    std::map<size_t, size_t> stroke_reorder;
    if (test_time) {
      std::sort(input.clusters.begin()->second.strokes.begin(),
                input.clusters.begin()->second.strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });
      std::sort(input.clusters.begin()->second.original_input_strokes.begin(),
                input.clusters.begin()->second.original_input_strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });
    }
    for (size_t i = 0; i < input.clusters.begin()->second.strokes.size(); ++i) {
      const auto &s = input.clusters.begin()->second.strokes[i];
      stroke_reorder[s.stroke_ind] = i;
    }

    auto is_sparse_cluster =
      [](const std::vector<Cluster::Stroke> &strokes) -> bool {
      return strokes.size() > 1 && strokes.size() <= 3;
    };

    // Orient
    // When the fitting doesn't exist
    std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
      unmerge_mapping;
    if (cluster1.fit.centerline.empty() || cluster2.fit.centerline.empty() ||
        is_sparse_cluster(cluster1.strokes) ||
        is_sparse_cluster(cluster2.strokes)) {
      std::map<int, std::vector<int>> orientations =
        fea_orientation.orient_strokes(&input);
      // Merge cut point if possible (cut info provided; orientation is
      // consistent)
      if (matching_pair_ptr)
        fea_orientation.flip_cut_strokes(&input, orientations,
                                         *matching_pair_ptr, &unmerge_mapping);
      else
        fea_orientation.flip_strokes(&input, orientations);

      if (cluster1.fit.centerline.empty()) {
        for (size_t i = 0; i < cluster1.strokes.size(); ++i) {
          size_t input_i = stroke_reorder[cluster1.strokes[i].stroke_ind];
          cluster1.strokes[i].points =
            input.clusters.begin()->second.strokes[input_i].points;
          cluster1.strokes[i].u.resize(cluster1.strokes[i].points.size(), 0);
        }
        fit_cluster(thickness, input.width, input.height, cluster1, false,
                    test_time);
      }
      if (cluster2.fit.centerline.empty()) {
        for (size_t i = 0; i < cluster2.strokes.size(); ++i) {
          size_t input_i = stroke_reorder[cluster2.strokes[i].stroke_ind];
          cluster2.strokes[i].points =
            input.clusters.begin()->second.strokes[input_i].points;
          cluster2.strokes[i].u.resize(cluster2.strokes[i].points.size(), 0);
        }
        fit_cluster(thickness, input.width, input.height, cluster2, false,
                    test_time);
      }

      // Individual parameterization failed
      if (cluster1.fit.centerline.empty() || cluster2.fit.centerline.empty()) {
        return false;
      }
    } else if (reorient || no_fit_orient) {
      std::map<int, std::vector<int>> orientations =
        fea_orientation.orient_strokes(&input);
      // Merge cut point if possible (cut info provided; orientation is
      // consistent)
      if (matching_pair_ptr)
        fea_orientation.flip_cut_strokes(&input, orientations,
                                         *matching_pair_ptr, &unmerge_mapping);
      else
        fea_orientation.flip_strokes(&input, orientations);
    } else { // When the fitting exists
      Capture capture_fit;
      capture_fit.thickness = 1;
      {
        SketchUI::Polyline2D fit_poly;
        fit_poly.from_fits(cluster1.fit, 0);
        double rate = stroke_sampling_size;
        fit_poly.reparameterize(
          std::min(rate, fit_poly.totalLen() / stroke_sampling_min_num));

        capture_fit.sketchedPolylines.emplace_back(fit_poly);
        capture_fit.sketchedPolylines.back().stroke_ind = 0;
        capture_fit.sketchedPolylines.back().group_ind = 0;
      }
      {
        SketchUI::Polyline2D fit_poly;
        fit_poly.from_fits(cluster2.fit, 1);
        double rate = stroke_sampling_size;
        fit_poly.reparameterize(
          std::min(rate, fit_poly.totalLen() / stroke_sampling_min_num));

        capture_fit.sketchedPolylines.emplace_back(fit_poly);
        capture_fit.sketchedPolylines.back().stroke_ind = 1;
        capture_fit.sketchedPolylines.back().group_ind = 0;
      }
      Input input_fit = capture2input.from_capture(capture_fit);
      std::map<int, std::vector<int>> simple_orientations =
        fea_orientation.orient_strokes(&input_fit);
      // Flip the first cluster
      std::map<int, std::vector<int>> orientations;
      orientations[0] = std::vector<int>{};
      for (auto &s : input.clusters.begin()->second.strokes) {
        if (simple_orientations[0][0] != simple_orientations[0][1] &&
            s.cluster_ind == 0)
          orientations[0].emplace_back(-1);
        else
          orientations[0].emplace_back(1);
      }

      // It is checked that no merging actually happen here
      fea_orientation.flip_strokes(&input, orientations);
    }

    auto obj_term_values = fea_param.parameterize(&input);

    double alignment_threshold = (cluster1.strokes.front().cluster_ind ==
                                  cluster2.strokes.front().cluster_ind)
                                   ? pos_parameterization_alignment_threshold
                                   : neg_parameterization_alignment_threshold;

    // Parameterization failed, drop this sample
    if (input.clusters[0].strokes.empty() ||
        obj_term_values["alignment"] >= alignment_threshold) {
      return false;
    }

    fit_map = fea_fitting.fit(&input);
    const auto &fit_curve = fit_map.at(0);
    merged_cluster = input.clusters.begin()->second;
    merged_cluster.original_input_strokes =
      input.clusters.begin()->second.original_input_strokes;
    merged_cluster.fit = fit_curve;
    merged_cluster.obj_term_values = obj_term_values;

    // Sort sample points along the same cross section based on their V values
    for (auto &xsec : merged_cluster.xsecs) {
      sort_xsec(xsec);
      map_xsec(xsec, unmerge_mapping);
    }
    map_strokes(unmerge_mapping, merged_cluster.strokes);
    if (!unmerge_mapping.empty()) {
      for (auto &xsec : merged_cluster.xsecs) {
        update_xsec_tangent(xsec, merged_cluster.strokes);
      }
    }
  } else {
    input.clusters.clear();
    input.clusters[0] = merged_cluster;
    fit_map[0] = merged_cluster.fit;
    pre_computed = true;
  }

  // Assignment cluster indices of 1 and 2. So in the feature computation, we
  // can associate clusters and fitting curves
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
  assign_cluster_index(merged_cluster, std::numeric_limits<size_t>::max(),
                       merged_cluster);

  assign_cluster_index(cluster1, 0, merged_cluster);
  assign_cluster_index(cluster2, 1, merged_cluster);

  return true;
}

static void
save_visualization(const Input &input, const Cluster &cluster1,
                   const Cluster &cluster2, const std::string &vis_output,
                   const std::string &svg_name, bool pre_computed,
                   bool hash_only, const std::map<int, FittedCurve> &fit_map,
                   std::map<std::string, std::string> &vis_files_) {
  make_folder(vis_output, false);

  // File root
  vis_files_["root"] = vis_output;

  std::vector<size_t> stroke_indices1, stroke_indices2;
  for (auto const &s : cluster1.strokes) {
    stroke_indices1.emplace_back(s.stroke_ind);
  }
  for (auto const &s : cluster2.strokes) {
    stroke_indices2.emplace_back(s.stroke_ind);
  }
  std::string hash_str = (stroke_indices2.empty())
                           ? example_hash(stroke_indices1)
                           : example_hash(stroke_indices1, stroke_indices2);
  if (hash_only)
    hash_str = get_hash_str_comb(cluster1, cluster2);

  // Input
  std::string final_output_name =
    vis_output + "/input/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/input/" + hash_str;
  final_output_name += ".svg";
  vis_files_["input"] = final_output_name;
  {
    std::ofstream input_svg(final_output_name);
    fea_param.save_svg(input_svg, input);
  }

  // Input cluster
  final_output_name =
    vis_output + "/input_cluster/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/input_cluster/" + hash_str;
  final_output_name += ".svg";
  vis_files_["input_cluster"] = final_output_name;
  {
    std::ofstream input_svg(final_output_name);
    fea_param.save_non_isoline_svg(input_svg, input);
  }

  // Isoline
  final_output_name =
    vis_output + "/isoline/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/isoline/" + hash_str;
  final_output_name += "_isolines.svg";
  vis_files_["isoline"] = final_output_name;
  if (!pre_computed || (pre_computed && solve_context.stage >
                                          SolveContext::SolveStage::Initial)) {
    std::ofstream isolines_svg(final_output_name);
    Input filtered_input = input;
    filtered_input.clusters[0].xsecs =
      filter_xsecs(filtered_input.clusters[0].xsecs);
    fea_param.isolines_cluster_svg(isolines_svg, -1, filtered_input, 0);
  }

  // Fitting
  final_output_name =
    vis_output + "/fit/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/fit/" + hash_str;
  final_output_name += "_fit.svg";
  vis_files_["fit"] = final_output_name;
  if ((!pre_computed ||
       (pre_computed &&
        solve_context.stage > SolveContext::SolveStage::Initial)) &&
      !fit_map.empty()) {
    std::ofstream fit_svg(final_output_name);
    fea_fitting.fit_svg(fit_svg, input, fit_map);
  }

  // Parameterization
  final_output_name =
    vis_output + "/parameterization/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/parameterization/" + hash_str;
  final_output_name += "_parameterization.svg";
  vis_files_["parameterization"] = final_output_name;
  if (!pre_computed || (pre_computed && solve_context.stage >
                                          SolveContext::SolveStage::Initial)) {
    std::ofstream param_svg(final_output_name);
    input.param_svg(param_svg);
  }

  // Orientation
  final_output_name =
    vis_output + "/orientation/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/orientation/" + hash_str;
  final_output_name += "_orientation.svg";
  vis_files_["orientation"] = final_output_name;
  if (!pre_computed || (pre_computed && solve_context.stage >
                                          SolveContext::SolveStage::Initial)) {
    std::ofstream orientations_svg(final_output_name);
    fea_orientation.orientation_debug(orientations_svg, input);
  }
}

static void
save_visualization_paths(const Input &input, const Cluster &cluster1,
                         const Cluster &cluster2, const std::string &vis_output,
                         const std::string &svg_name, bool pre_computed,
                         bool hash_only,
                         const std::map<int, FittedCurve> &fit_map,
                         std::map<std::string, std::string> &vis_files_) {
  // File root
  vis_files_["root"] = vis_output;

  std::vector<size_t> stroke_indices1, stroke_indices2;
  for (auto const &s : cluster1.strokes) {
    stroke_indices1.emplace_back(s.stroke_ind);
  }
  for (auto const &s : cluster2.strokes) {
    stroke_indices2.emplace_back(s.stroke_ind);
  }
  std::string hash_str = (stroke_indices2.empty())
                           ? example_hash(stroke_indices1)
                           : example_hash(stroke_indices1, stroke_indices2);
  if (hash_only)
    hash_str = get_hash_str_comb(cluster1, cluster2);

  // Input
  std::string final_output_name =
    vis_output + "/input/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/input/" + hash_str;
  final_output_name += ".svg";
  vis_files_["input"] = final_output_name;

  // Input cluster
  final_output_name =
    vis_output + "/input_cluster/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/input_cluster/" + hash_str;
  final_output_name += ".svg";
  vis_files_["input_cluster"] = final_output_name;

  // Isoline
  final_output_name =
    vis_output + "/isoline/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/isoline/" + hash_str;
  final_output_name += "_isolines.svg";
  vis_files_["isoline"] = final_output_name;

  // Fitting
  final_output_name =
    vis_output + "/fit/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/fit/" + hash_str;
  final_output_name += "_fit.svg";
  vis_files_["fit"] = final_output_name;

  // Parameterization
  final_output_name =
    vis_output + "/parameterization/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/parameterization/" + hash_str;
  final_output_name += "_parameterization.svg";
  vis_files_["parameterization"] = final_output_name;

  // Orientation
  final_output_name =
    vis_output + "/orientation/" + svg_name +
    std::to_string(cluster1.strokes.front().cluster_ind) +
    ((cluster2.strokes.empty())
       ? ""
       : ("_" + std::to_string(cluster2.strokes.front().cluster_ind))) +
    "_" + hash_str;
  if (hash_only)
    final_output_name = vis_output + "/orientation/" + hash_str;
  final_output_name += "_orientation.svg";
  vis_files_["orientation"] = final_output_name;
}

void FeatureVector::compute(
  Cluster &cluster1, Cluster &cluster2, double thickness,
  Cluster &merged_cluster, bool test_time,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool reorient, std::string vis_output) {
  // Record input info
  cluster_info_vec_.resize(2);

  cluster_info_vec_[0].cluster_idx = cluster1.strokes.front().cluster_ind;
  for (const auto &s : cluster1.strokes)
    cluster_info_vec_[0].stroke_indices.emplace_back(s.stroke_ind);
  cluster_info_vec_[1].cluster_idx = cluster2.strokes.front().cluster_ind;
  for (const auto &s : cluster2.strokes)
    cluster_info_vec_[1].stroke_indices.emplace_back(s.stroke_ind);
  std::string svg_name = (cluster1.strokes.front().cluster_ind ==
                          cluster2.strokes.front().cluster_ind)
                           ? "pos_"
                           : "neg_";
  if (cluster1.strokes.front().cluster_ind < 0)
    svg_name = "pred_";

  const auto feature_calculators = get_stroke_cluster_features();
  features_.reserve(feature_calculators.size());
  descriptions_.reserve(feature_calculators.size());

  Input input;
  std::map<int, FittedCurve> fit_map;
  bool pre_computed;

  bool succeeded =
    merge_cached(cluster1, cluster2, thickness, merged_cluster, reorient,
                 test_time, input, fit_map, pre_computed, matching_pair_ptr);
  if (!succeeded) {
    fill_inf(feature_calculators, descriptions_, features_);
    if (!vis_output.empty()) {
      input.clusters[0] = merged_cluster;
      save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                               pre_computed, test_time, fit_map, vis_files_);
    }
    return;
  }
  obj_term_values_ = merged_cluster.obj_term_values;

  // Test if any side-by-side exists
  bool is_side_by_side = false;
  for (auto const &xsec : merged_cluster.xsecs)
    if (is_overlapping(xsec)) {
      is_side_by_side = true;
      break;
    }

  // Compute features
  if (!is_side_by_side) {
    // Save computatiion when there is no side-by-side section
    fill_inf(feature_calculators, descriptions_, features_);
  } else {
    compute_features(feature_calculators, cluster1, cluster2, merged_cluster,
                     descriptions_, features_);
  }

  if (!vis_output.empty()) {
    input.clusters[0] = merged_cluster;
    save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                             pre_computed, test_time, fit_map, vis_files_);
  }
}

void FeatureVector::compute_angle(
  Cluster &cluster1, Cluster &cluster2, double thickness,
  Cluster &merged_cluster,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool reorient, std::string vis_output) {
  // Record input info
  cluster_info_vec_.resize(2);

  cluster_info_vec_[0].cluster_idx = cluster1.strokes.front().cluster_ind;
  for (const auto &s : cluster1.strokes)
    cluster_info_vec_[0].stroke_indices.emplace_back(s.stroke_ind);
  cluster_info_vec_[1].cluster_idx = cluster2.strokes.front().cluster_ind;
  for (const auto &s : cluster2.strokes)
    cluster_info_vec_[1].stroke_indices.emplace_back(s.stroke_ind);
  std::string svg_name = (cluster1.strokes.front().cluster_ind ==
                          cluster2.strokes.front().cluster_ind)
                           ? "pos_"
                           : "neg_";
  if (cluster1.strokes.front().cluster_ind < 0)
    svg_name = "pred_";

  const auto feature_calculators = get_stroke_cluster_angle_features();
  features_.reserve(feature_calculators.size());
  descriptions_.reserve(feature_calculators.size());

  Input input;
  std::map<int, FittedCurve> fit_map;
  bool pre_computed;
  bool test_time = true;
  bool succeeded =
    merge_cached(cluster1, cluster2, thickness, merged_cluster, reorient,
                 test_time, input, fit_map, pre_computed, matching_pair_ptr);
  if (!succeeded) {
    fill_inf(feature_calculators, descriptions_, features_);
    return;
  }
  obj_term_values_ = merged_cluster.obj_term_values;

  // Test if any side-by-side exists
  bool is_side_by_side = false;
  for (auto const &xsec : merged_cluster.xsecs)
    if (is_overlapping(xsec)) {
      is_side_by_side = true;
      break;
    }

  // Compute features
  if (!is_side_by_side) {
    // Save computatiion when there is no side-by-side section
    fill_inf(feature_calculators, descriptions_, features_);
  } else {
    compute_features(feature_calculators, cluster1, cluster2, merged_cluster,
                     descriptions_, features_);
  }

  if (!vis_output.empty()) {
    input.clusters[0] = merged_cluster;
    save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                             pre_computed, true, fit_map, vis_files_);
  }
}

void FeatureVector::compute_distance(
  Cluster &cluster1, Cluster &cluster2, double thickness,
  Cluster &merged_cluster,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool reorient, std::string vis_output) {
  // Record input info
  cluster_info_vec_.resize(2);

  cluster_info_vec_[0].cluster_idx = cluster1.strokes.front().cluster_ind;
  for (const auto &s : cluster1.strokes)
    cluster_info_vec_[0].stroke_indices.emplace_back(s.stroke_ind);
  cluster_info_vec_[1].cluster_idx = cluster2.strokes.front().cluster_ind;
  for (const auto &s : cluster2.strokes)
    cluster_info_vec_[1].stroke_indices.emplace_back(s.stroke_ind);
  std::string svg_name = (cluster1.strokes.front().cluster_ind ==
                          cluster2.strokes.front().cluster_ind)
                           ? "pos_"
                           : "neg_";
  if (cluster1.strokes.front().cluster_ind < 0)
    svg_name = "pred_";

  const auto feature_calculators = get_stroke_cluster_distance_features();
  features_.reserve(feature_calculators.size());
  descriptions_.reserve(feature_calculators.size());

  Input input;
  std::map<int, FittedCurve> fit_map;
  bool pre_computed;
  bool test_time = true;
  bool succeeded =
    merge_cached(cluster1, cluster2, thickness, merged_cluster, reorient,
                 test_time, input, fit_map, pre_computed, matching_pair_ptr);
  if (!succeeded) {
    fill_inf(feature_calculators, descriptions_, features_);
    return;
  }
  obj_term_values_ = merged_cluster.obj_term_values;

  // Test if any side-by-side exists
  bool is_side_by_side = false;
  for (auto const &xsec : merged_cluster.xsecs)
    if (is_overlapping(xsec)) {
      is_side_by_side = true;
      break;
    }

  // Compute features
  if (!is_side_by_side) {
    // Save computatiion when there is no side-by-side section
    fill_inf(feature_calculators, descriptions_, features_);
  } else {
    compute_features(feature_calculators, cluster1, cluster2, merged_cluster,
                     descriptions_, features_);
  }

  if (!vis_output.empty()) {
    input.clusters[0] = merged_cluster;
    save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                             pre_computed, true, fit_map, vis_files_);
  }
}

std::vector<FeatureVector> compute_features(
  const std::vector<std::pair<Cluster, Cluster>> &in_clusters,
  const std::vector<Sample> &samples, const int width, const int height,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool test_time, const std::string &output_vis_folder,
  const std::string &cache_folder) {
  std::vector<std::pair<Cluster, Cluster>> clusters = in_clusters;

  double norm_thickness = 1;
  std::vector<FeatureVector> feas;
  const auto feature_calculators = get_stroke_cluster_features();

  for (auto &cc : clusters) {
    feas.emplace_back(FeatureVector(width, height));

    Cluster merged_cluster;
    std::string json_filename1 = "";
    std::string json_filename2 = "";
    std::string json_filename_all = "";
    if (!cache_folder.empty()) {
      std::vector<size_t> stroke_indices1, stroke_indices2, stroke_indices_comb;
      for (auto const &s : cc.first.strokes) {
        stroke_indices1.emplace_back(s.stroke_ind);
        stroke_indices_comb.emplace_back(s.stroke_ind);
      }
      for (auto const &s : cc.second.strokes) {
        stroke_indices2.emplace_back(s.stroke_ind);
        stroke_indices_comb.emplace_back(s.stroke_ind);
      }

      // Keep each cluster internally sorted
      std::sort(stroke_indices_comb.begin(), stroke_indices_comb.end());

      std::string hash_str1 = example_hash(stroke_indices1);
      std::string hash_str2 = example_hash(stroke_indices2);
      std::string hash_str_comb = example_hash(stroke_indices_comb);

      json_filename1 = cache_folder + "/" + hash_str1 + ".json";
      json_filename2 = cache_folder + "/" + hash_str2 + ".json";
      json_filename_all = cache_folder + "/" + hash_str_comb + ".json";
      bool hit1 = read_cache(json_filename1, cc.first);
      bool hit2 = read_cache(json_filename2, cc.second);

      // Cache exists but the individual parameterization failed
      if ((hit1 && cc.first.fit.centerline.empty()) ||
          (hit2 && cc.second.fit.centerline.empty())) {
        fill_inf(feature_calculators, feas.back().descriptions_,
                 feas.back().features_);
        continue;
      }
    }

    // Fit the two individual clusters
    // Write cache here instead of after the potential re-orientation since
    // the combined parameterization needs to be recomputed each time
    if (cc.first.fit.centerline.empty()) {
      bool to_orient = false;
      fit_cluster(norm_thickness, width, height, cc.first, to_orient, test_time,
                  matching_pair_ptr);
      if (!cache_folder.empty()) {
        write_json(json_filename1, cc.first);
      }
    }
    if (cc.second.fit.centerline.empty()) {
      bool to_orient = false;
      fit_cluster(norm_thickness, width, height, cc.second, to_orient,
                  test_time, matching_pair_ptr);
      if (!cache_folder.empty()) {
        write_json(json_filename2, cc.second);
      }
    }

    bool reorient = false;
    feas.back().compute(cc.first, cc.second, norm_thickness, merged_cluster,
                        false, matching_pair_ptr, reorient, output_vis_folder);

    // Add the parameterization info as well.
    for (const auto &k_v : merged_cluster.obj_term_values) {
      feas.back().descriptions_.emplace_back(k_v.first);
      feas.back().features_.emplace_back(k_v.second);
    }

    if (!cache_folder.empty()) {
      write_json(json_filename_all, merged_cluster);
    }
  }

  return feas;
}

std::vector<FeatureVector> compute_features(
  const Input &input, const std::vector<Sample> &samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool test_time, const std::string &output_vis_folder,
  const std::string &cache_folder) {
  std::vector<std::pair<Cluster, Cluster>> clusters;
  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : input.clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      sid2s[s.stroke_ind] = s;
    }
  }

  for (const auto &ss : samples) {
    clusters.emplace_back();
    for (auto s_idx : ss.first) {
      clusters.back().first.strokes.emplace_back(sid2s[s_idx]);
      clusters.back().first.strokes.back().cluster_ind = -1;
    }
    for (auto s_idx : ss.second) {
      clusters.back().second.strokes.emplace_back(sid2s[s_idx]);
      clusters.back().second.strokes.back().cluster_ind = -1;
    }
  }

  return compute_features(clusters, samples, input.width, input.height,
                          matching_pair_ptr, test_time, output_vis_folder,
                          cache_folder);
}

std::vector<FeatureVector> compute_features(
  const std::string &scap_filename, const std::vector<Sample> &samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool test_time, const std::string &output_vis_folder,
  const std::string &cache_folder) {
  int width, height;
  bool to_preprocess = false;
  std::vector<std::pair<Cluster, Cluster>> clusters;
  read_example(scap_filename, samples, width, height, clusters, to_preprocess);

  return compute_features(clusters, samples, width, height, matching_pair_ptr,
                          test_time, output_vis_folder, cache_folder);
}

////////////////////////////////////////

void FeatureVector::compute_typical(
  Cluster &cluster1, double thickness,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool reorient, std::string vis_output) {
  auto begin = std::chrono::high_resolution_clock::now();

  assert(!cluster1.strokes.empty());

  // Record input info
  cluster_info_vec_.emplace_back(ClusterInfo());
  cluster_info_vec_[0].cluster_idx = cluster1.strokes.front().cluster_ind;
  for (const auto &s : cluster1.strokes)
    cluster_info_vec_[0].stroke_indices.emplace_back(s.stroke_ind);

  int cluster_num = cluster1.strokes.front().cluster_ind;
  bool is_positive = true;
  for (const auto &s : cluster1.strokes) {
    if (s.cluster_ind != cluster_num) {
      is_positive = false;
      break;
    }
  }

  std::string svg_name = is_positive ? "pos_" : "neg_";
  if (cluster1.strokes.front().cluster_ind < 0)
    svg_name = "pred_";

  const auto feature_calculators = get_typical_cluster_features();
  features_.reserve(feature_calculators.size());
  descriptions_.reserve(feature_calculators.size());

  Input input;
  input.thickness = 1;
  input.clusters[0] = cluster1;
  obj_term_values_ = cluster1.obj_term_values;

  std::map<int, FittedCurve> fit_map;
  fit_map[0] = cluster1.fit;

  // Compute features
  for (auto const &fea : feature_calculators) {
    descriptions_.emplace_back(fea->csv());
    // Save computatiion when there is no side-by-side section
    if (!features_.empty() &&
        features_.front() == std::numeric_limits<double>::infinity())
      features_.emplace_back(std::numeric_limits<double>::infinity());
    else {
      auto fea_vec = fea->operator()(input.clusters[0]);
      features_.insert(features_.end(), fea_vec.begin(), fea_vec.end());
    }
  }

  // Add the parameterization info as well.
  for (const auto &k_v : obj_term_values_) {
    descriptions_.emplace_back(k_v.first);
    features_.emplace_back(k_v.second);
  }

  if (!vis_output.empty()) {
    bool pre_computed = false;
    bool test_time = false;
    Cluster cluster2;
    save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                             pre_computed, test_time, fit_map, vis_files_);
  }

  auto end = std::chrono::high_resolution_clock::now();
  solve_context.feature_timer +=
    std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() /
    1000.0;
}

////////////////////////////////////////

void FeatureVector::compute_probability(
  Cluster &cluster1, Cluster &cluster2, double thickness,
  Cluster &merged_cluster,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool reorient, std::string vis_output) {
  // Record input info
  cluster_info_vec_.resize(2);

  cluster_info_vec_[0].cluster_idx = cluster1.strokes.front().cluster_ind;
  for (const auto &s : cluster1.strokes)
    cluster_info_vec_[0].stroke_indices.emplace_back(s.stroke_ind);
  cluster_info_vec_[1].cluster_idx = cluster2.strokes.front().cluster_ind;
  for (const auto &s : cluster2.strokes)
    cluster_info_vec_[1].stroke_indices.emplace_back(s.stroke_ind);
  std::string svg_name = (cluster1.strokes.front().cluster_ind ==
                          cluster2.strokes.front().cluster_ind)
                           ? "pos_"
                           : "neg_";
  if (cluster1.strokes.front().cluster_ind < 0)
    svg_name = "pred_";

  const auto feature_calculators = get_stroke_cluster_probability_features();
  features_.reserve(feature_calculators.size());
  descriptions_.reserve(feature_calculators.size());

  Input input;
  std::map<int, FittedCurve> fit_map;
  bool pre_computed;
  bool test_time = false;
  bool succeeded =
    merge_cached(cluster1, cluster2, thickness, merged_cluster, reorient,
                 test_time, input, fit_map, pre_computed, matching_pair_ptr);
  if (!succeeded) {
    fill_inf(feature_calculators, descriptions_, features_);
    return;
  }
  obj_term_values_ = merged_cluster.obj_term_values;

  // Test if any side-by-side exists
  bool is_side_by_side = false;
  for (auto const &xsec : merged_cluster.xsecs)
    if (is_overlapping(xsec)) {
      is_side_by_side = true;
      break;
    }

  // Compute features
  if (!is_side_by_side) {
    // Save computatiion when there is no side-by-side section
    fill_inf(feature_calculators, descriptions_, features_);
  } else {
    compute_features(feature_calculators, cluster1, cluster2, merged_cluster,
                     descriptions_, features_);
  }

  if (!vis_output.empty()) {
    input.clusters[0] = merged_cluster;
    bool test_time = false;
    save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                             pre_computed, test_time, fit_map, vis_files_);
  }
}

void FeatureVector::compute_secondary(
  Cluster &cluster1, Cluster &cluster2, double thickness,
  Cluster &merged_cluster,
  const std::map<size_t, FeatureVector> &cluster_features,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool reorient, std::string vis_output) {
  // Record input info
  cluster_info_vec_.resize(2);

  cluster_info_vec_[0].cluster_idx = cluster1.strokes.front().cluster_ind;
  for (const auto &s : cluster1.strokes)
    cluster_info_vec_[0].stroke_indices.emplace_back(s.stroke_ind);
  cluster_info_vec_[1].cluster_idx = cluster2.strokes.front().cluster_ind;
  for (const auto &s : cluster2.strokes)
    cluster_info_vec_[1].stroke_indices.emplace_back(s.stroke_ind);
  std::string svg_name = (cluster1.strokes.front().cluster_ind ==
                          cluster2.strokes.front().cluster_ind)
                           ? "pos_"
                           : "neg_";
  if (cluster1.strokes.front().cluster_ind < 0)
    svg_name = "pred_";

  const auto feature_calculators = get_stroke_cluster_features();
  const auto global_feature_calculators = get_secondary_cluster_features();
  size_t global_size = 0;
  for (const auto &fea : global_feature_calculators)
    global_size += fea->csvs().size();
  features_.reserve(feature_calculators.size() + global_size);
  descriptions_.reserve(feature_calculators.size() + global_size);

  Input input;
  std::map<int, FittedCurve> fit_map;
  bool pre_computed;
  bool test_time = false;
  bool succeeded =
    merge_cached(cluster1, cluster2, thickness, merged_cluster, reorient,
                 test_time, input, fit_map, pre_computed, matching_pair_ptr);
  if (solve_context.stage == SolveContext::SolveStage::Merging && !succeeded) {
    bool seen_cut = false;
    if (matching_pair_ptr) {
      std::set<size_t> seen_cut_strokes;
      for (const auto &ss : *matching_pair_ptr) {
        seen_cut_strokes.emplace(ss.first.first);
        seen_cut_strokes.emplace(ss.second.first);
      }
      for (const auto &s : cluster1.strokes)
        if (seen_cut_strokes.count(s.stroke_ind)) {
          seen_cut = true;
          break;
        }
      for (const auto &s : cluster2.strokes)
        if (seen_cut_strokes.count(s.stroke_ind)) {
          seen_cut = true;
          break;
        }
    }
    if (seen_cut) {
      context.ignore_endpoint = true;
      succeeded = merge_cached(cluster1, cluster2, thickness, merged_cluster,
                               reorient, test_time, input, fit_map,
                               pre_computed, matching_pair_ptr);
      if (!succeeded) {
        succeeded =
          merge_cached(cluster1, cluster2, thickness, merged_cluster, reorient,
                       test_time, input, fit_map, pre_computed, nullptr);
      }
      context.ignore_endpoint = false;
    }
  }
  if (!succeeded) {
    fill_inf(feature_calculators, descriptions_, features_);
    fill_inf(global_feature_calculators, descriptions_, features_);
    return;
  }
  obj_term_values_ = merged_cluster.obj_term_values;

  // Test if any side-by-side exists
  bool is_side_by_side = false;
  for (auto const &xsec : merged_cluster.xsecs)
    if (is_overlapping(xsec)) {
      is_side_by_side = true;
      break;
    }

  // Compute features
  if (!is_side_by_side) {
    // Save computatiion when there is no side-by-side section
    fill_inf(feature_calculators, descriptions_, features_);
    fill_inf(global_feature_calculators, descriptions_, features_);
    return;
  } else {
    compute_features(feature_calculators, cluster1, cluster2, merged_cluster,
                     descriptions_, features_);
  }

  auto begin = std::chrono::high_resolution_clock::now();

  // Add the global features
  {
    FeatureVector base_fea;
    for (size_t i = 0; i < descriptions_.size(); ++i) {
      base_fea.descriptions_.emplace_back(descriptions_[i]);
      base_fea.features_.emplace_back(features_[i]);
    }
    for (auto const &ov : obj_term_values_) {
      base_fea.descriptions_.emplace_back(ov.first);
      base_fea.features_.emplace_back(ov.second);
    }

    std::set<size_t> seen_strokes;
    for (auto const &s : cluster1.original_input_strokes) {
      seen_strokes.emplace(s.stroke_ind);
    }
    for (auto const &s : cluster2.original_input_strokes) {
      seen_strokes.emplace(s.stroke_ind);
    }
    std::vector<FeatureVector> other_cluster_fea;
    for (const auto &c_fea : cluster_features) {
      bool seen_intersect = false;
      bool seen_non_intersect = false;
      for (auto sid : c_fea.second.cluster_info_vec_[0].stroke_indices) {
        if (seen_strokes.count(sid)) {
          seen_intersect = true;
          break;
        }
      }
      if (!seen_intersect) {
        other_cluster_fea.emplace_back(c_fea.second);
      }
    }

    if (other_cluster_fea.empty()) {
      for (const auto &c_fea : cluster_features) {
        other_cluster_fea.emplace_back(c_fea.second);
      }
    }

    for (auto const &fea : global_feature_calculators) {
      auto csvs = fea->csvs();
      descriptions_.insert(descriptions_.end(), csvs.begin(), csvs.end());
      // Save computatiion when there is no side-by-side section
      if (!features_.empty() &&
          features_.front() == std::numeric_limits<double>::infinity())
        features_.emplace_back(std::numeric_limits<double>::infinity());
      else {
        auto fea_vec = fea->operator()(base_fea, other_cluster_fea);
        features_.insert(features_.end(), fea_vec.begin(), fea_vec.end());
      }
    }
  }

  if (!vis_output.empty()) {
    input.clusters[0] = merged_cluster;
    bool test_time = false;
    save_visualization_paths(input, cluster1, cluster2, vis_output, svg_name,
                             pre_computed, test_time, fit_map, vis_files_);
  }

  auto end = std::chrono::high_resolution_clock::now();
  solve_context.feature_timer +=
    std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() /
    1000.0;
}
