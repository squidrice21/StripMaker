/**
 * This file contains implementations of various per-strip features that may be
 * used for global feature computation.
 */
#pragma once

#include "../stroke_strip_src/Cluster.h"
#include "Features.h"

#include <algorithm>
#include <string>
#include <vector>

struct MeanAngleTypicalClusterFeature final : TypicalClusterFeature {
  MeanAngleTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean angular diff"; }
  std::string csv() const { return "average_angle"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileAngleTypicalClusterFeature final : TypicalClusterFeature {
  PercentileAngleTypicalClusterFeature(double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) + " angular diff";
    else
      return "Min angular diff";
  }
  std::string csv() const {
    if (percentile_ == 0.5)
      return "median_angle";
    else if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) + "th_angle";
    else if (percentile_ == 0)
      return "min_angle";
    else
      return "max_angle";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};

struct MeanDistanceTypicalClusterFeature final : TypicalClusterFeature {
  MeanDistanceTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean max local gap"; }
  std::string csv() const { return "average_distance"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileDistanceTypicalClusterFeature final : TypicalClusterFeature {
  PercentileDistanceTypicalClusterFeature(double percentile,
                                          size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " max local distance";
    else
      return "Min max local gap";
  }
  std::string csv() const {
    if (percentile_ == 0.5)
      return "median_distance";
    else if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) +
             "th_distance";
    else if (percentile_ == 0)
      return "min_distance";
    else
      return "max_distance";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};

struct StrokeCountTypicalClusterFeature final : TypicalClusterFeature {
  StrokeCountTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Stroke count"; }
  std::string csv() const { return "stroke_count"; }
  std::vector<double> operator()(const Cluster &cluster) const {
    return std::vector<double>{(double)cluster.strokes.size()};
  }
};

struct MeanXsecStrokeCountTypicalClusterFeature final : TypicalClusterFeature {
  MeanXsecStrokeCountTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean stroke count per xsec"; }
  std::string csv() const { return "average_xsec_stroke_count"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileXsecStrokeCountTypicalClusterFeature final
  : TypicalClusterFeature {
  PercentileXsecStrokeCountTypicalClusterFeature(double percentile,
                                                 size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke count per xsec";
    else
      return "Min stroke count per xsec";
  }
  std::string csv() const {
    if (percentile_ == 0.5)
      return "median_xsec_stroke_count";
    else if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) +
             "th_xsec_stroke_count";
    else if (percentile_ == 0)
      return "min_xsec_stroke_count";
    else
      return "max_xsec_stroke_count";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};

struct MeanWidthTypicalClusterFeature final : TypicalClusterFeature {
  MeanWidthTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean width"; }
  std::string csv() const { return "average_width"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileWidthTypicalClusterFeature final : TypicalClusterFeature {
  PercentileWidthTypicalClusterFeature(double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) + " width";
    else
      return "Min width";
  }
  std::string csv() const {
    if (percentile_ == 0.5)
      return "median_width";
    else if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) + "th_width";
    else if (percentile_ == 0)
      return "min_width";
    else
      return "max_width";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};

struct MeanGapRatioTypicalClusterFeature final : TypicalClusterFeature {
  MeanGapRatioTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean max local gap / width"; }
  std::string csv() const { return "average_distance2width_combined"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileGapRatioTypicalClusterFeature final : TypicalClusterFeature {
  PercentileGapRatioTypicalClusterFeature(double percentile,
                                          size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " max local gap / width";
    else
      return "Min max local gap / width";
  }
  std::string csv() const {
    if (percentile_ == 0.5)
      return "median_distance2width_combined";
    else if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) +
             "th_distance2width_combined";
    else if (percentile_ == 0)
      return "min_distance2width_combined";
    else
      return "max_distance2width_combined";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};

struct MeanWidthLengthRatioTypicalClusterFeature final : TypicalClusterFeature {
  MeanWidthLengthRatioTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean width / length"; }
  std::string csv() const { return "average_narrowness"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileWidthLengthRatioTypicalClusterFeature final
  : TypicalClusterFeature {
  PercentileWidthLengthRatioTypicalClusterFeature(double percentile,
                                                  size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " max width / length";
    else
      return "Min max width / length";
  }
  std::string csv() const {
    if (percentile_ == 0.5)
      return "median_narrowness";
    else if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) +
             "th_narrowness";
    else if (percentile_ == 0)
      return "min_narrowness";
    else
      return "max_narrowness";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};

struct MeanWidthChangeRatioTypicalClusterFeature final : TypicalClusterFeature {
  MeanWidthChangeRatioTypicalClusterFeature(size_t num_binned = 1)
    : TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean width change"; }
  std::string csv() const { return "avg_WidthChange"; }
  std::vector<double> operator()(const Cluster &cluster) const;
};

struct PercentileWidthChangeRatioTypicalClusterFeature final
  : TypicalClusterFeature {
  PercentileWidthChangeRatioTypicalClusterFeature(double percentile,
                                                  size_t num_binned = 1)
    : percentile_(percentile)
    , TypicalClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) + " max width change";
    else
      return "Min max width change";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "global_p" + std::to_string(int(percentile_ * 100)) +
             "th_WidthChange";
    else if (percentile_ == 0)
      return "min_WidthChange";
    else
      return "max_WidthChange";
  }
  std::vector<double> operator()(const Cluster &cluster) const;

  double percentile_;
};
