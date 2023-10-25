/**
 * This file contains implementations of various local features used in
 * StripMaker. It includes functions to calculate means and percentiles of
 * various point measures for overlapping cross sections. It also defines
 * several structs that represent stroke strip features, which are used to
 * extract features from strips of strokes.
 */
#pragma once

#include "../Util.h"
#include "../stroke_strip_src/Cluster.h"
#include "Features.h"

#include <string>
#include <vector>

/**
 * Sets all infs in the given vector to the specified value.
 *
 * @param vec The vector to modify.
 * @param value The value to set all elements to.
 */
void set_inf(std::vector<double> &vec, const double value);

/**
 * Calculates the mean of a point measure for overlapping cross sections, binned
 * into a specified number of bins.
 *
 * @param xsecs The vector of cross sections to calculate the mean for.
 * @param num_binned The number of bins to divide the cross sections into.
 * @param point_measure A function that takes a cross section and returns a
 * vector of point measures.
 * @return A vector of means for each bin.
 */
std::vector<double> get_overlapping_xsecs_binned_mean(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure);
/**
 * Calculates the binned percentile of overlapping cross sections based on a
 * given point measure.
 *
 * @param xsecs The vector of cross sections to calculate the binned percentile
 * for.
 * @param percentile The percentile to calculate.
 * @param num_binned The number of bins to use.
 * @param point_measure The function to use to calculate the point measure for
 * each cross section.
 * @return A vector of binned percentiles for the overlapping cross sections.
 */
std::vector<double> get_overlapping_xsecs_binned_percentile(
  const std::vector<Cluster::XSec> &xsecs, double percentile, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure);

/**
 * Calculates the mean of a point measure for each of the `num_binned` bins
 * of cross sections in `xsecs`.
 *
 * @param xsecs The cross sections to bin and calculate the mean of.
 * @param num_binned The number of bins to divide the cross sections into.
 * @param point_measure A function that takes a cross section and returns a
 *                      vector of point measures for that cross section.
 * @return A vector of means, one for each bin.
 */
std::vector<double> get_all_xsecs_binned_mean(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure);
/**
 * Calculates the binned percentile of a given point measure for all cross
 * sections in a vector.
 *
 * @param xsecs The vector of cross sections to calculate the binned percentile
 * for.
 * @param percentile The percentile to calculate.
 * @param num_binned The number of bins to use for the calculation.
 * @param point_measure A function that takes a cross section and returns a
 * vector of point measures.
 * @return A vector of binned percentiles for each cross section.
 */
std::vector<double> get_all_xsecs_binned_percentile(
  const std::vector<Cluster::XSec> &xsecs, double percentile, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure);

/**
 * Returns a vector of function returned values of bins in the overlapping
 * section.
 *
 * @param xsecs The vector of cross sections to bin.
 * @param num_binned The number of bins to use.
 * @param func_binned The function to use to calculate the value for each bin.
 * @return A vector of returned values of a given function for bins in the
 * overlapping section.
 */
std::vector<double> get_all_xsecs_binned(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<double(span<Cluster::XSec const> const)> func_binned);

/**
 * Calculates the Euclidean distance between two strips in the overlapping
 * section.
 *
 * @param xsec The cross sections (slices) to calculate distances for.
 * @return A vector containing the Euclidean distances between each pair of
 * cross sections.
 */
std::vector<double> get_xsec_euclidean_clusterwise(const Cluster::XSec &xsec);

////////////////////////////////////////////////////////

struct MeanAngleStrokeClusterFeature final : StrokeClusterFeature {
  MeanAngleStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean stroke-cluster angle"; }
  std::string csv() const { return "average_angle"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileAngleStrokeClusterFeature final : StrokeClusterFeature {
  PercentileAngleStrokeClusterFeature(double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Percentile " + std::to_string(percentile_) +
           " stroke-cluster angle";
  }
  std::string csv() const {
    if (percentile_ == 1)
      return "max_angle";
    return "median_angle";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanWidthLengthRatio1StrokeClusterFeature final : StrokeClusterFeature {
  MeanWidthLengthRatio1StrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster width / length of cluster 1";
  }
  std::string csv() const { return "avg_narrowness_cluster1"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileWidthLengthRatio1StrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileWidthLengthRatio1StrokeClusterFeature(double percentile,
                                                  size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster width / length of cluster 1";
    else
      return "Min stroke-cluster width / length of cluster 1";
  }
  std::string csv() const { return "median_narrowness_cluster1"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanWidthLengthRatio2StrokeClusterFeature final : StrokeClusterFeature {
  MeanWidthLengthRatio2StrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster width / length of cluster 2";
  }
  std::string csv() const { return "avg_narrowness_cluster2"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileWidthLengthRatio2StrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileWidthLengthRatio2StrokeClusterFeature(double percentile,
                                                  size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster width / length of cluster 2";
    else
      return "Min stroke-cluster width / length of cluster 2";
  }
  std::string csv() const { return "median_narrowness_cluster2"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanWidthLengthMaxRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanWidthLengthMaxRatioStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Max Mean stroke-cluster width / length of cluster";
  }
  std::string csv() const { return "maxC12_avg_narrowness"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileWidthLengthMaxRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileWidthLengthMaxRatioStrokeClusterFeature(double percentile,
                                                    size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Max Percentile " + std::to_string(percentile_) +
             " stroke-cluster width / length of cluster";
    else
      return "Max Min stroke-cluster width / length of cluster";
  }
  std::string csv() const { return "maxC12_median_narrowness"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanWidthLengthMinRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanWidthLengthMinRatioStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Min Mean stroke-cluster width / length of cluster";
  }
  std::string csv() const { return "minC12_avg_narrowness"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileWidthLengthMinRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileWidthLengthMinRatioStrokeClusterFeature(double percentile,
                                                    size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Min Percentile " + std::to_string(percentile_) +
             " stroke-cluster width / length of cluster";
    else
      return "Min Min stroke-cluster width / length of cluster";
  }
  std::string csv() const { return "minC12_median_narrowness"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanWidthLengthRatioOverallStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanWidthLengthRatioOverallStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster width / length of cluster 1+2";
  }
  std::string csv() const { return "avg_narrowness_combined"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileWidthLengthRatioOverallStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileWidthLengthRatioOverallStrokeClusterFeature(double percentile,
                                                        size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster width / length of cluster 1+2";
    else
      return "Min stroke-cluster width / length of cluster 1+2";
  }
  std::string csv() const { return "median_narrowness_combined"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanWidthStrokeClusterFeature final : StrokeClusterFeature {
  MeanWidthStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean overlapping combined width"; }
  std::string csv() const { return "average_width_combined"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileWidthStrokeClusterFeature final : StrokeClusterFeature {
  PercentileWidthStrokeClusterFeature(double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " overlapping combined width";
    else
      return "Min overlapping combined width";
  }
  std::string csv() const { return "median_width_combined"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanClusterwiseDistanceRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistanceRatioStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster clusterwise distance / merged cluster width";
  }
  std::string csv() const { return "avg_distance2width_combined"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileClusterwiseDistanceRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileClusterwiseDistanceRatioStrokeClusterFeature(double percentile,
                                                         size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster clusterwise distance / merged cluster width";
    else
      return "Min stroke-cluster clusterwise distance / merged cluster width";
  }
  std::string csv() const { return "median_distance2width_combined"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanClusterwiseDistanceMaxWidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistanceMaxWidthRatioStrokeClusterFeature(
    size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster clusterwise distance / max(width1, width2)";
  }
  std::string csv() const { return "avg_distance2LocalMaxWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileClusterwiseDistanceMaxWidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileClusterwiseDistanceMaxWidthRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster clusterwise distance / max(width1, width2)";
    else
      return "Min stroke-cluster clusterwise distance / max(width1, width2)";
  }
  std::string csv() const { return "median_distance2LocalMaxWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanClusterwiseDistanceMinWidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistanceMinWidthRatioStrokeClusterFeature(
    size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster clusterwise distance / min(width1, width2)";
  }
  std::string csv() const { return "avg_distance2LocalMinWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileClusterwiseDistanceMinWidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileClusterwiseDistanceMinWidthRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster clusterwise distance / min(width1, width2)";
    else
      return "Min stroke-cluster clusterwise distance / min(width1, width2)";
  }
  std::string csv() const { return "median_distance2LocalMinWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanClusterwiseDistanceMaxC12WidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistanceMaxC12WidthRatioStrokeClusterFeature(
    size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster clusterwise distance / max(width1, width2)";
  }
  std::string csv() const { return "avg_distance2MaxC12Width"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileClusterwiseDistanceMaxC12WidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileClusterwiseDistanceMaxC12WidthRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster clusterwise distance / max(width1, width2)";
    else
      return "Min stroke-cluster clusterwise distance / max(width1, width2)";
  }
  std::string csv() const { return "median_distance2MaxC12Width"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanClusterwiseDistanceMinC12WidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistanceMinC12WidthRatioStrokeClusterFeature(
    size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster clusterwise distance / min(width1, width2)";
  }
  std::string csv() const { return "avg_distance2MinC12Width"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileClusterwiseDistanceMinC12WidthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileClusterwiseDistanceMinC12WidthRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster clusterwise distance / min(width1, width2)";
    else
      return "Min stroke-cluster clusterwise distance / min(width1, width2)";
  }
  std::string csv() const { return "median_distance2MinC12Width"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

// Warn: May be buggy
struct PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster per-cluster gap / clusterwise distance";
    else if (percentile_ == 0)
      return "Min stroke-cluster per-cluster gap / clusterwise distance";
    else
      return "Max stroke-cluster per-cluster gap / clusterwise distance";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "p" + std::to_string(int(percentile_ * 100)) +
             "th_LocalMaxGap2width_cluster1";
    else if (percentile_ == 0)
      return "min_LocalMaxGap2width_cluster1";
    else
      return "max_LocalMaxGap2width_cluster1";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

// Warn: May be buggy
struct PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster per-cluster gap / clusterwise distance";
    else if (percentile_ == 0)
      return "Min stroke-cluster per-cluster gap / clusterwise distance";
    else
      return "Max stroke-cluster per-cluster gap / clusterwise distance";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "p" + std::to_string(int(percentile_ * 100)) +
             "th_LocalMaxGap2width_cluster2";
    else if (percentile_ == 0)
      return "min_LocalMaxGap2width_cluster2";
    else
      return "max_LocalMaxGap2width_cluster2";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct PercentileGapClusterwiseDistanceMaxRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileGapClusterwiseDistanceMaxRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "Max Percentile " + std::to_string(percentile_) +
             " stroke-cluster per-cluster gap / clusterwise distance";
    else if (percentile_ == 0)
      return "Max Min stroke-cluster per-cluster gap / clusterwise distance";
    else
      return "Max Max stroke-cluster per-cluster gap / clusterwise distance";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "maxC12_" + std::to_string(int(percentile_ * 100)) +
             "th_LocalMedianGap2width";
    else if (percentile_ == 0)
      return "maxC12_min_LocalMedianGap2width";
    else
      return "maxC12_max_LocalMedianGap2width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct PercentileGapClusterwiseDistanceMinRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileGapClusterwiseDistanceMinRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "Min Percentile " + std::to_string(percentile_) +
             " stroke-cluster per-cluster gap / clusterwise distance";
    else if (percentile_ == 0)
      return "Min Min stroke-cluster per-cluster gap / clusterwise distance";
    else
      return "Min Max stroke-cluster per-cluster gap / clusterwise distance";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "minC12_" + std::to_string(int(percentile_ * 100)) +
             "th_LocalMedianGap2width";
    else if (percentile_ == 0)
      return "minC12_min_LocalMedianGap2width";
    else
      return "minC12_max_LocalMedianGap2width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct FirstWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  FirstWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "First stroke-cluster width change ratio";
  }
  std::string csv() const { return "first_local_nonoverlapping2overlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};
struct LastWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  LastWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Last stroke-cluster width change ratio";
  }
  std::string csv() const { return "last_local_nonoverlapping2overlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MeanIndividualWidth1StrokeClusterFeature final : StrokeClusterFeature {
  MeanIndividualWidth1StrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster cluster 1 width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileIndividualWidth1StrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileIndividualWidth1StrokeClusterFeature(double percentile,
                                                 size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster cluster 1 width";
    else
      return "Min stroke-cluster cluster 1 width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanIndividualWidth2StrokeClusterFeature final : StrokeClusterFeature {
  MeanIndividualWidth2StrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster cluster 2 width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileIndividualWidth2StrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileIndividualWidth2StrokeClusterFeature(double percentile,
                                                 size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster cluster 2 width";
    else
      return "Min stroke-cluster cluster 2 width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanMergedWidthStrokeClusterFeature final : StrokeClusterFeature {
  MeanMergedWidthStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean stroke-cluster merged width"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileMergedWidthStrokeClusterFeature final : StrokeClusterFeature {
  PercentileMergedWidthStrokeClusterFeature(double percentile,
                                            size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster merged width";
    else
      return "Min stroke-cluster merged width";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MinPointwiseDistanceStrokeClusterFeature final : StrokeClusterFeature {
  MinPointwiseDistanceStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Min stroke-cluster pointwise distance";
  }
  std::string csv() const { return "min_pointwise_distance"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MeanDistanceStrokeClusterFeature final : StrokeClusterFeature {
  MeanDistanceStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const { return "Mean stroke-cluster distance"; }
  std::string csv() const { return "average_distance"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileDistanceStrokeClusterFeature final : StrokeClusterFeature {
  PercentileDistanceStrokeClusterFeature(double percentile,
                                         size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster distance";
    else
      return "Min stroke-cluster distance";
  }
  std::string csv() const { return "median_distance"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct OverlappingRatioStrokeClusterFeature final : StrokeClusterFeature {
  OverlappingRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Side-by-side length / length of cluster 1+2";
  }
  std::string csv() const { return "overlapping_length2combined_length"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  MinWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Min stroke-cluster width change ratio";
  }
  std::string csv() const { return "minC12_local_nonoverlapping2overlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};
struct MaxWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  MaxWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Max stroke-cluster width change ratio";
  }
  std::string csv() const { return "maxC12_local_nonoverlapping2overlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinAvgWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  MinAvgWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Min average stroke-cluster width change ratio";
  }
  std::string csv() const { return "maxC12_AvgOverlapping2AvgNonoverlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};
struct MaxAvgWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  MaxAvgWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Max average stroke-cluster width change ratio";
  }
  std::string csv() const { return "minC12_AvgOverlapping2AvgNonoverlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct InvMinWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  InvMinWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Min stroke-cluster width change ratio inv";
  }
  std::string csv() const { return "minC12_local_overlapping2nonoverlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};
struct InvMaxWidthChangeRatioStrokeClusterFeature final : StrokeClusterFeature {
  InvMaxWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Max stroke-cluster width change ratio inv";
  }
  std::string csv() const { return "maxC12_local_overlapping2nonoverlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct InvMaxAvgWidthChangeRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  InvMaxAvgWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Max average stroke-cluster width change ratio inv";
  }
  std::string csv() const { return "maxC12_AvgNonoverlapping2AvgOverlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};
struct InvMinAvgWidthChangeRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  InvMinAvgWidthChangeRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Min average stroke-cluster width change ratio inv";
  }
  std::string csv() const { return "minC12_AvgNonoverlapping2AvgOverlapping"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct TypeStrokeClusterFeature final : StrokeClusterFeature {
  TypeStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Example type"; }
  std::string csv() const { return "type"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MeanClusterwiseDistaceMaxWidthStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistaceMaxWidthStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean clusterwise distance / max(c1, c2 average cluster width)";
  }
  std::string csv() const { return "AvgDistance2MaxC12AvgWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MedianClusterwiseDistaceMaxWidthStrokeClusterFeature final
  : StrokeClusterFeature {
  MedianClusterwiseDistaceMaxWidthStrokeClusterFeature(double percentile,
                                                       size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Median clusterwise distance / max(c1, c2 average cluster width)";
  }
  std::string csv() const { return "MedianDistance2MaxC12AvgWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
  double percentile_;
};

struct MeanClusterwiseDistaceMinWidthStrokeClusterFeature final
  : StrokeClusterFeature {
  MeanClusterwiseDistaceMinWidthStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean clusterwise distance / min(c1, c2 average cluster width)";
  }
  std::string csv() const { return "AvgDistance2MinC12AvgWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MedianClusterwiseDistaceMinWidthStrokeClusterFeature final
  : StrokeClusterFeature {
  MedianClusterwiseDistaceMinWidthStrokeClusterFeature(double percentile,
                                                       size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Median clusterwise distance / min(c1, c2 average cluster width)";
  }
  std::string csv() const { return "MedianDistance2MinC12AvgWidth"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
  double percentile_;
};

struct PersistenceStrokeClusterFeature final : StrokeClusterFeature {
  PersistenceStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Gap persistence ratio"; }
  std::string csv() const { return "persistence_ratio"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

//

struct MaxOverlappingRatioStrokeClusterFeature final : StrokeClusterFeature {
  MaxOverlappingRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Side-by-side length / min(length of cluster 1, 2)";
  }
  std::string csv() const { return "maxC12_overlapping_length2ind_length"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinOverlappingRatioStrokeClusterFeature final : StrokeClusterFeature {
  MinOverlappingRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "Side-by-side length / max(length of cluster 1, 2)";
  }
  std::string csv() const { return "minC12_overlapping_length2ind_length"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct OverlappingLengthStrokeClusterFeature final : StrokeClusterFeature {
  OverlappingLengthStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Side-by-side length normalized"; }
  std::string csv() const { return "overlapping_length"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileGapOverlappingLengthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileGapOverlappingLengthRatioStrokeClusterFeature(double percentile,
                                                          size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster per-cluster gap / overlapping length";
    else if (percentile_ == 0)
      return "Min stroke-cluster per-cluster gap / overlapping length";
    else
      return "Max stroke-cluster per-cluster gap / overlapping length";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "p" + std::to_string(int(percentile_ * 100)) +
             "th_LocalMaxGap2overlapping_length";
    else if (percentile_ == 0)
      return "min_LocalMaxGap2overlapping_length";
    else
      return "max_LocalMaxGap2overlapping_length";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct PercentileWidthOverlappingLengthRatioStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileWidthOverlappingLengthRatioStrokeClusterFeature(
    double percentile, size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster per-cluster width / overlapping length";
    else if (percentile_ == 0)
      return "Min stroke-cluster per-cluster width / overlapping length";
    else
      return "Max stroke-cluster per-cluster width / overlapping length";
  }
  std::string csv() const {
    if (percentile_ != 0 && percentile_ != 1)
      return "p" + std::to_string(int(percentile_ * 100)) +
             "th_LocalWidth2overlapping_length";
    else if (percentile_ == 0)
      return "min_LocalWidth2overlapping_length";
    else
      return "max_LocalWidth2overlapping_length";
  }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct OverlappingIsolineNumberStrokeClusterFeature final
  : StrokeClusterFeature {
  OverlappingIsolineNumberStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const {
    return "the number of isolines in overlapping";
  }
  std::string csv() const { return "overlapping_isolines_number"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

////////////////////////////////////////////////////////////////////////////////

struct MaxStepawayRatioStrokeClusterFeature final : StrokeClusterFeature {
  MaxStepawayRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Max step-away ratio"; }
  std::string csv() const { return "max_stepawayratio"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinStepawayRatioStrokeClusterFeature final : StrokeClusterFeature {
  MinStepawayRatioStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Min step-away ratio"; }
  std::string csv() const { return "min_stepawayratio"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MaxStepawayRatio2StrokeClusterFeature final : StrokeClusterFeature {
  MaxStepawayRatio2StrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Max step-away ratio 2"; }
  std::string csv() const { return "max_stepawayratio2"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinStepawayRatio2StrokeClusterFeature final : StrokeClusterFeature {
  MinStepawayRatio2StrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Min step-away ratio 2"; }
  std::string csv() const { return "min_stepawayratio2"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MaxStepawayRatio3StrokeClusterFeature final : StrokeClusterFeature {
  MaxStepawayRatio3StrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Max step-away ratio 3"; }
  std::string csv() const { return "max_stepawayratio3"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinStepawayRatio3StrokeClusterFeature final : StrokeClusterFeature {
  MinStepawayRatio3StrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Min step-away ratio 3"; }
  std::string csv() const { return "min_stepawayratio3"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MaxStepawayAngleStrokeClusterFeature final : StrokeClusterFeature {
  MaxStepawayAngleStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Max step-away angle"; }
  std::string csv() const { return "max_stepawayangle"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MinStepawayAngleStrokeClusterFeature final : StrokeClusterFeature {
  MinStepawayAngleStrokeClusterFeature()
    : StrokeClusterFeature(1) {}
  std::string description() const { return "Min step-away angle"; }
  std::string csv() const { return "min_stepawayangle"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct MeanOutsideAngleStrokeClusterFeature final : StrokeClusterFeature {
  MeanOutsideAngleStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster angle where clusters are separated";
  }
  std::string csv() const { return "average_outside_angle"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileOutsideAngleStrokeClusterFeature final : StrokeClusterFeature {
  PercentileOutsideAngleStrokeClusterFeature(double percentile,
                                             size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Percentile " + std::to_string(percentile_) +
           " stroke-cluster angle where clusters are separated";
  }
  std::string csv() const { return "median_outside_angle"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};

struct MeanOutsideDistanceStrokeClusterFeature final : StrokeClusterFeature {
  MeanOutsideDistanceStrokeClusterFeature(size_t num_binned = 1)
    : StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    return "Mean stroke-cluster distance where clusters are separated";
  }
  std::string csv() const { return "average_outside_distance"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;
};

struct PercentileOutsideDistanceStrokeClusterFeature final
  : StrokeClusterFeature {
  PercentileOutsideDistanceStrokeClusterFeature(double percentile,
                                                size_t num_binned = 1)
    : percentile_(percentile)
    , StrokeClusterFeature(std::max((size_t)1, num_binned)) {}
  std::string description() const {
    if (percentile_ != 0)
      return "Percentile " + std::to_string(percentile_) +
             " stroke-cluster distance where clusters are separated";
    else
      return "Min stroke-cluster distance where clusters are separated";
  }
  std::string csv() const { return "median_outside_distance"; }
  std::vector<double> operator()(const Cluster &merged_cluster,
                                 const FittedCurve &fit1,
                                 const FittedCurve &fit2) const;

  double percentile_;
};
