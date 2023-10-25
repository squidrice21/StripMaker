/**
 * This file contains implementations of various global features used in
 * StripMaker. It includes functions to calculate ratios defined between
 * inner-strip measurements.
 */
#pragma once

#include "../stroke_strip_src/Cluster.h"
#include "Features.h"

#include <algorithm>
#include <string>
#include <vector>

struct InitialProbabilityStrokeClusterFeature final : StrokeClusterFeature {
  InitialProbabilityStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Probability from the initial classifier";
  }
  std::string csv() const { return "init_prob"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MeanGlobalClusterFeature final : SecondaryClusterFeature {
  MeanGlobalClusterFeature(const std::vector<std::string> &feas)
    : SecondaryClusterFeature()
    , feas_(feas) {}
  std::string description() const { return "Mean angular diff"; }
  std::vector<std::string> csvs() const {
    std::string prefix = "global_average_";
    std::vector<std::string> csv_str;
    for (const auto &fea : feas_) {
      csv_str.emplace_back(prefix + fea);
    }
    return csv_str;
  }
  std::vector<double>
  operator()(const FeatureVector &base_feature,
             const std::vector<FeatureVector> &features) const;

  std::vector<std::string> feas_;
};

struct PercentileGlobalClusterFeature final : SecondaryClusterFeature {
  PercentileGlobalClusterFeature(const std::vector<std::string> &feas,
                                 double percentile)
    : SecondaryClusterFeature()
    , percentile_(percentile)
    , feas_(feas) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) + " angular diff";
    else
      return "Min angular diff";
  }
  std::vector<std::string> csvs() const {
    std::string prefix;
    if (percentile_ == 0.5)
      prefix = "global_median_";
    else if (percentile_ != 0 && percentile_ != 1)
      prefix = "global_p" + std::to_string(int(percentile_ * 100)) + "th_";
    else if (percentile_ == 0)
      prefix = "global_min_";
    else
      prefix = "global_max_";
    std::vector<std::string> csv_str;
    for (const auto &fea : feas_) {
      csv_str.emplace_back(prefix + fea);
    }
    return csv_str;
  }
  std::vector<double>
  operator()(const FeatureVector &base_feature,
             const std::vector<FeatureVector> &features) const;

  double percentile_;
  std::vector<std::string> feas_;
};

struct MeanGlobalRatioClusterFeature final : SecondaryClusterFeature {
  MeanGlobalRatioClusterFeature(const std::vector<std::string> &base_feas,
                                const std::vector<std::string> &feas)
    : SecondaryClusterFeature()
    , base_feas_(base_feas)
    , feas_(feas) {}
  std::string description() const { return "Mean angular diff"; }
  std::vector<std::string> csvs() const {
    std::string prefix = "global_average_";
    std::vector<std::string> csv_str;
    for (const auto &b_fea : base_feas_) {
      for (const auto &fea : feas_) {
        csv_str.emplace_back(b_fea + "2" + prefix + fea);
      }
    }
    return csv_str;
  }
  std::vector<double>
  operator()(const FeatureVector &base_feature,
             const std::vector<FeatureVector> &features) const;

  std::vector<std::string> base_feas_;
  std::vector<std::string> feas_;

  double epsilon = 1e-3;
};

struct PercentileGlobalRatioClusterFeature final : SecondaryClusterFeature {
  PercentileGlobalRatioClusterFeature(const std::vector<std::string> &base_feas,
                                      const std::vector<std::string> &feas,
                                      double percentile)
    : SecondaryClusterFeature()
    , percentile_(percentile)
    , base_feas_(base_feas)
    , feas_(feas) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) + " angular diff";
    else
      return "Min angular diff";
  }
  std::vector<std::string> csvs() const {
    std::string prefix;
    if (percentile_ == 0.5)
      prefix = "global_median_";
    else if (percentile_ != 0 && percentile_ != 1)
      prefix = "global_p" + std::to_string(int(percentile_ * 100)) + "th_";
    else if (percentile_ == 0)
      prefix = "global_min_";
    else
      prefix = "global_max_";
    std::vector<std::string> csv_str;
    for (const auto &b_fea : base_feas_) {
      for (const auto &fea : feas_) {
        csv_str.emplace_back(b_fea + "2" + prefix + fea);
      }
    }
    return csv_str;
  }
  std::vector<double>
  operator()(const FeatureVector &base_feature,
             const std::vector<FeatureVector> &features) const;

  double percentile_;
  std::vector<std::string> base_feas_;
  std::vector<std::string> feas_;

  double epsilon = 1e-3;
};
