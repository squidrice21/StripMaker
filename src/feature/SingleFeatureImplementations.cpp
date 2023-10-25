#include "SingleFeatureImplementations.h"
#include "FeatureImplementations.h"

#include <algorithm>
#include <vector>

static std::vector<double>
get_angle_diff_typical_cluster(const Cluster::XSec &xsec) {
  std::vector<double> angle_vec;

  if (!xsec.points.empty()) {
    double max_gap = -1;
    size_t div_i = 0;
    std::vector<size_t> sides;
    sides.resize(xsec.points.size(), 0);

    for (size_t i = 1; i < xsec.points.size(); ++i) {
      double g = std::max(
        glm::distance(xsec.points[i].point, xsec.points[i - 1].point) - 1.0,
        0.0);
      if (g > max_gap) {
        div_i = i;
      }
    }

    for (size_t i = 1; i < xsec.points.size(); ++i) {
      if (i >= div_i) {
        sides[i] = 1;
      }
    }

    for (auto &connection : xsec.connections) {
      if (sides[connection.a_idx] != sides[connection.b_idx]) {
        double dot = glm::dot(xsec.points[connection.a_idx].tangent,
                              xsec.points[connection.b_idx].tangent);
        dot = std::min(std::max(dot, -1.0), 1.0);
        double angle = glm::acos(dot) * 180.0 / M_PI;
        angle_vec.emplace_back(angle);
      }
    }
  }

  return angle_vec;
}

static std::vector<double> get_gap_typical_cluster(const Cluster::XSec &xsec) {
  std::vector<double> merged_width;
  if (!xsec.points.empty()) {
    double max_gap = -1;
    for (size_t i = 1; i < xsec.points.size(); ++i) {
      max_gap = std::max(
        glm::distance(xsec.points[i].point, xsec.points[i - 1].point) - 1.0,
        0.0);
    }
    if (max_gap >= 0)
      merged_width.emplace_back(max_gap);
  }
  return merged_width;
}

static std::vector<double>
get_gap_ratio_typical_cluster(const Cluster::XSec &xsec) {
  std::vector<double> merged_width;
  if (!xsec.points.empty()) {
    double max_gap = -1;
    double w =
      glm::distance(xsec.points.front().point, xsec.points.back().point) + 1;
    for (size_t i = 1; i < xsec.points.size(); ++i) {
      max_gap = std::max(
        glm::distance(xsec.points[i].point, xsec.points[i - 1].point) - 1.0,
        0.0);
    }
    if (max_gap >= 0)
      merged_width.emplace_back(max_gap / w);
  }
  return merged_width;
}

static std::vector<double>
get_width_typical_cluster(const Cluster::XSec &xsec) {
  std::vector<double> merged_width;
  if (!xsec.points.empty()) {
    double w =
      glm::distance(xsec.points.front().point, xsec.points.back().point) + 1;
    merged_width.emplace_back(w);
  }
  return merged_width;
}

static std::vector<double>
get_stroke_count_typical_cluster(const Cluster::XSec &xsec) {
  return std::vector<double>{(double)xsec.points.size()};
}

std::vector<double>
MeanAngleTypicalClusterFeature::operator()(const Cluster &cluster) const {
  std::vector<double> p_local_gap = get_all_xsecs_binned_mean(
    cluster.xsecs, Feature::num_binned_, get_angle_diff_typical_cluster);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double>
PercentileAngleTypicalClusterFeature::operator()(const Cluster &cluster) const {
  std::vector<double> p_local_gap = get_all_xsecs_binned_percentile(
    cluster.xsecs, percentile_, Feature::num_binned_,
    get_angle_diff_typical_cluster);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double>
MeanDistanceTypicalClusterFeature::operator()(const Cluster &cluster) const {
  std::vector<double> p_local_gap = get_all_xsecs_binned_mean(
    cluster.xsecs, Feature::num_binned_, get_gap_typical_cluster);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double> PercentileDistanceTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  std::vector<double> p_local_gap = get_all_xsecs_binned_percentile(
    cluster.xsecs, percentile_, Feature::num_binned_, get_gap_typical_cluster);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double> MeanXsecStrokeCountTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  std::vector<double> p_local_count = get_all_xsecs_binned_mean(
    cluster.xsecs, Feature::num_binned_, get_stroke_count_typical_cluster);
  set_inf(p_local_count, 0);

  return p_local_count;
}

std::vector<double> PercentileXsecStrokeCountTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  std::vector<double> p_local_count = get_all_xsecs_binned_percentile(
    cluster.xsecs, percentile_, Feature::num_binned_,
    get_stroke_count_typical_cluster);
  set_inf(p_local_count, 0);

  return p_local_count;
}

std::vector<double>
MeanWidthTypicalClusterFeature::operator()(const Cluster &cluster) const {
  std::vector<double> p_local_width = get_all_xsecs_binned_mean(
    cluster.xsecs, Feature::num_binned_, get_width_typical_cluster);
  set_inf(p_local_width, 0);

  return p_local_width;
}

std::vector<double>
PercentileWidthTypicalClusterFeature::operator()(const Cluster &cluster) const {
  std::vector<double> p_local_width = get_all_xsecs_binned_percentile(
    cluster.xsecs, percentile_, Feature::num_binned_,
    get_width_typical_cluster);
  set_inf(p_local_width, 0);

  return p_local_width;
}

std::vector<double>
MeanGapRatioTypicalClusterFeature::operator()(const Cluster &cluster) const {
  std::vector<double> p_local_gap = get_all_xsecs_binned_mean(
    cluster.xsecs, Feature::num_binned_, get_gap_ratio_typical_cluster);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double> PercentileGapRatioTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  std::vector<double> p_local_gap = get_all_xsecs_binned_percentile(
    cluster.xsecs, percentile_, Feature::num_binned_,
    get_gap_ratio_typical_cluster);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double> MeanWidthLengthRatioTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  std::vector<double> p_local_width = get_all_xsecs_binned_mean(
    cluster.xsecs, Feature::num_binned_, get_width_typical_cluster);
  double l = total_length(cluster.fit.centerline);
  for (auto &w : p_local_width)
    w /= l;
  set_inf(p_local_width, 0);

  return p_local_width;
}

std::vector<double> PercentileWidthLengthRatioTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  std::vector<double> p_local_width = get_all_xsecs_binned_percentile(
    cluster.xsecs, percentile_, Feature::num_binned_,
    get_width_typical_cluster);
  double l = total_length(cluster.fit.centerline);
  for (auto &w : p_local_width)
    w /= l;
  set_inf(p_local_width, 0);

  return p_local_width;
}

std::vector<double> MeanWidthChangeRatioTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  auto aggregated_func = [](span<Cluster::XSec const> const xsecs) -> double {
    std::vector<double> measure_vec;
    measure_vec.reserve(xsecs.size());
    for (size_t i = 0; i + 1 < xsecs.size(); ++i) {
      auto &xsec = xsecs[i];
      const auto measures = get_width_typical_cluster(xsec);
      auto &next = xsecs[i + 1];
      const auto measures_next = get_width_typical_cluster(next);

      if (measures.empty() || measures_next.empty())
        continue;

      double r1 = measures.front() / measures_next.front();
      double r2 = measures_next.front() / measures.front();

      if (!measures.empty())
        measure_vec.emplace_back(std::max(r1, r2));
    }
    if (measure_vec.empty())
      return 1;

    double sum = 0;
    for (const auto v : measure_vec)
      sum += v;
    return sum / measure_vec.size();
  };
  std::vector<double> p_local_gap =
    get_all_xsecs_binned(cluster.xsecs, Feature::num_binned_, aggregated_func);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}

std::vector<double> PercentileWidthChangeRatioTypicalClusterFeature::operator()(
  const Cluster &cluster) const {
  double percentile = percentile_;
  auto aggregated_func =
    [&percentile](span<Cluster::XSec const> const xsecs) -> double {
    std::vector<double> measure_vec;
    measure_vec.reserve(xsecs.size());
    for (size_t i = 0; i + 1 < xsecs.size(); ++i) {
      auto &xsec = xsecs[i];
      const auto measures = get_width_typical_cluster(xsec);
      auto &next = xsecs[i + 1];
      const auto measures_next = get_width_typical_cluster(next);

      if (measures.empty() || measures_next.empty())
        continue;

      double r1 = measures.front() / measures_next.front();
      double r2 = measures_next.front() / measures.front();

      if (!measures.empty())
        measure_vec.emplace_back(std::max(r1, r2));
    }
    if (measure_vec.empty())
      return 1;

    assert(percentile >= 0 && percentile <= 1);
    int pos = (measure_vec.size() - 1) * percentile;
    pos = std::max(pos, 0);
    auto p = measure_vec.begin() + pos;
    std::nth_element(measure_vec.begin(), p, measure_vec.end());
    return measure_vec[pos];
  };
  std::vector<double> p_local_gap =
    get_all_xsecs_binned(cluster.xsecs, Feature::num_binned_, aggregated_func);
  set_inf(p_local_gap, 0);

  return p_local_gap;
}
