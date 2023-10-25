#include "FeatureImplementations.h"

#include "../Logger.h"
#include "../Util.h"
#include "glm/detail/type_vec.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <set>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

double sandwiched_weight = 10;
size_t end_noisy_range = 0;

void set_inf(std::vector<double> &vec, const double value) {
  for (auto &e : vec)
    if (e == std::numeric_limits<double>::infinity())
      e = value;
}

static bool is_point_noisy(const Cluster::XSecPoint &p) {
  return (p.stroke_within_idx + 1) <= end_noisy_range ||
         (p.stroke_xsec_count - p.stroke_within_idx) <= end_noisy_range;
}

static double total_length(const Cluster &merged_cluster,
                           const FittedCurve &fit, int c) {
  if (!fit.centerline.empty())
    return total_length(fit.centerline);
  auto min_u = [](const Cluster &cluster, int c) {
    double result = std::numeric_limits<double>::infinity();
    for (auto &stroke : cluster.strokes) {
      if (c < 0 && stroke.cluster_ind >= std::numeric_limits<size_t>::max())
        continue;
      if (c >= 0 && stroke.cluster_ind != c)
        continue;
      for (double u : stroke.u) {
        result = std::min(result, u);
      }
    }
    return result;
  };
  auto max_u = [](const Cluster &cluster, int c) {
    double result = -std::numeric_limits<double>::infinity();
    for (auto &stroke : cluster.strokes) {
      if (c < 0 && stroke.cluster_ind >= std::numeric_limits<size_t>::max())
        continue;
      if (c >= 0 && stroke.cluster_ind != c)
        continue;
      for (double u : stroke.u) {
        result = std::max(result, u);
      }
    }
    return result;
  };

  return std::abs(max_u(merged_cluster, c) - min_u(merged_cluster, c));
}

static std::vector<double> get_xsec_angle(const Cluster::XSec &xsec) {
  std::vector<double> angle_vec;
  for (auto &connection : xsec.connections) {
    assert(connection.a_idx < xsec.points.size() &&
           connection.b_idx < xsec.points.size());
    if (is_point_noisy(xsec.points[connection.a_idx]) ||
        is_point_noisy(xsec.points[connection.b_idx]))
      continue;
    if (xsec.points[connection.a_idx].cluster_ind !=
        xsec.points[connection.b_idx].cluster_ind) {
      double dot = glm::dot(xsec.points[connection.a_idx].tangent,
                            xsec.points[connection.b_idx].tangent);
      dot = std::min(std::max(dot, -1.0), 1.0);
      double angle = glm::acos(dot) * 180.0 / M_PI;
      angle_vec.emplace_back(angle);
    }
  }
  return angle_vec;
}
static std::vector<double> get_xsec_euclidean(const Cluster::XSec &xsec) {
  std::vector<double> euclidean_clusterwise;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]) || is_point_noisy(xsec.points[i + 1]))
      continue;
    glm::dvec2 dir = xsec.points[i].point - xsec.points[i + 1].point;
    dir = glm::normalize(dir);
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      double dis =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;

      dis = std::max(dis, 0.0);
      euclidean_clusterwise.emplace_back(dis);
      // The two clusters are not clearly divided
      if (euclidean_clusterwise.size() > 1) {
        euclidean_clusterwise.clear();
        euclidean_clusterwise.emplace_back(0);
        break;
      }
    }
  }

  return euclidean_clusterwise;
}
static std::vector<double> get_xsec_merged_width(const Cluster::XSec &xsec) {
  std::vector<double> merged_width;
  std::vector<size_t> non_noisy_indices;
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      continue;
    non_noisy_indices.emplace_back(i);
  }
  if (!non_noisy_indices.empty()) {
    double xsec_width =
      glm::distance(xsec.points[non_noisy_indices.front()].point,
                    xsec.points[non_noisy_indices.back()].point) +
      1.0;
    merged_width.emplace_back(xsec_width);
  }
  return merged_width;
}
static bool is_xsec_noisy(const Cluster::XSec &xsec) {
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      return true;
  }
  return false;
}
static std::vector<double>
get_xsec_euclidean_clusterwise_ratio(const Cluster::XSec &xsec) {
  std::vector<double> euclidean_clusterwise;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]) || is_point_noisy(xsec.points[i + 1]))
      continue;
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      double dis =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      dis = std::max(dis, 0.0);
      euclidean_clusterwise.emplace_back(dis);
      // The two clusters are not clearly divided
      if (euclidean_clusterwise.size() > 1) {
        euclidean_clusterwise.clear();
        euclidean_clusterwise.emplace_back(0);
        break;
      }
    }
  }
  if (!xsec.points.empty() && !euclidean_clusterwise.empty()) {
    double xsec_width = get_xsec_merged_width(xsec).front();
    euclidean_clusterwise[0] /= xsec_width;
  }
  return euclidean_clusterwise;
}
std::vector<double> get_xsec_euclidean_clusterwise(const Cluster::XSec &xsec) {
  std::vector<double> euclidean_clusterwise;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]) || is_point_noisy(xsec.points[i + 1]))
      continue;
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      double dis =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      dis = std::max(dis, 0.0);
      euclidean_clusterwise.emplace_back(dis);
      // The two clusters are not clearly divided
      if (euclidean_clusterwise.size() > 1) {
        euclidean_clusterwise.clear();
        euclidean_clusterwise.emplace_back(0);
        break;
      }
    }
  }
  return euclidean_clusterwise;
}
static std::vector<double> get_xsec_width(const Cluster::XSec &xsec) {
  std::vector<double> euclidean_clusterwise;
  auto xsec_widths = get_xsec_merged_width(xsec);
  if (!xsec_widths.empty()) {
    double xsec_width = xsec_widths.front();
    euclidean_clusterwise.emplace_back(1 / xsec_width);
  }
  return euclidean_clusterwise;
}

static std::vector<double> get_xsec_individual_widths(const Cluster::XSec &xsec,
                                                      int cid) {
  std::map<size_t, std::vector<glm::dvec2>> cluster_samples;
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      continue;
    cluster_samples[xsec.points[i].cluster_ind].emplace_back(
      xsec.points[i].point);
  }

  std::vector<double> cluster_width;
  for (const auto &c_samples : cluster_samples) {
    if (cid < 0 || c_samples.first == cid) {
      double dist =
        glm::distance(c_samples.second.front(), c_samples.second.back()) + 1.0;
      cluster_width.emplace_back(dist);
    }
  }
  if (cid < 0)
    std::sort(cluster_width.begin(), cluster_width.end());
  return cluster_width;
}

static std::vector<double> get_xsec_width_ratio(const Cluster::XSec &xsec) {
  std::vector<double> cluster_width = get_xsec_individual_widths(xsec, -1);
  std::vector<double> euclidean_clusterwise;
  auto xsec_widths = get_xsec_merged_width(xsec);
  if (!xsec_widths.empty()) {
    double xsec_width = xsec_widths.front();
    euclidean_clusterwise.emplace_back(cluster_width.back() / xsec_width);
  }
  return euclidean_clusterwise;
}
static std::vector<double>
get_xsec_width_clusterwise_euclidean_ratio(const Cluster::XSec &xsec) {
  const auto width = get_xsec_width_ratio(xsec);
  const auto cwise_dist = get_xsec_euclidean_clusterwise_ratio(xsec);
  assert(width.size() == cwise_dist.size());
  std::vector<double> ratio;
  for (size_t i = 0; i < width.size(); ++i) {
    if (width[i] != 0)
      ratio.emplace_back(cwise_dist[i] / width[i]);
    else
      ratio.emplace_back(0);
  }
  return ratio;
}
static std::vector<double>
get_xsec_clusterwise_euclidean_max_width_ratio(const Cluster::XSec &xsec) {
  std::vector<double> cluster_width = get_xsec_individual_widths(xsec, -1);
  if (cluster_width.empty())
    return std::vector<double>();
  double max_width =
    *std::max_element(cluster_width.begin(), cluster_width.end());
  max_width = std::max(max_width, 1.0);

  const auto cwise_dist = get_xsec_euclidean_clusterwise(xsec);
  if (cwise_dist.empty())
    return std::vector<double>();

  return std::vector<double>{cwise_dist.front() / max_width};
}
static std::vector<double>
get_xsec_clusterwise_euclidean_min_width_ratio(const Cluster::XSec &xsec) {
  std::vector<double> cluster_width = get_xsec_individual_widths(xsec, -1);
  if (cluster_width.empty())
    return std::vector<double>();
  double min_width =
    *std::min_element(cluster_width.begin(), cluster_width.end());
  min_width = std::max(min_width, 1.0);

  const auto cwise_dist = get_xsec_euclidean_clusterwise(xsec);
  if (cwise_dist.empty())
    return std::vector<double>();

  return std::vector<double>{cwise_dist.front() / min_width};
}

static std::vector<double>
get_xsec_clusterwise_euclidean_cluster_width_ratio(const Cluster::XSec &xsec,
                                                   int cid) {
  std::vector<double> cluster_width = get_xsec_individual_widths(xsec, cid);
  if (cluster_width.empty())
    return std::vector<double>();
  double c_width = std::max(cluster_width[0], 1.0);

  const auto cwise_dist = get_xsec_euclidean_clusterwise(xsec);
  if (cwise_dist.empty())
    return std::vector<double>();

  return std::vector<double>{cwise_dist.front() / c_width};
}

static std::vector<double> get_xsec_max_gap_ratio(const Cluster::XSec &xsec,
                                                  int cluster_idx) {
  std::map<size_t, std::vector<glm::dvec2>> cluster_samples;
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      continue;
    cluster_samples[xsec.points[i].cluster_ind].emplace_back(
      xsec.points[i].point);
  }

  std::vector<double> cluster_max_gap;
  for (const auto &c_samples : cluster_samples) {
    if (cluster_idx >= 0 && c_samples.first != cluster_idx)
      continue;
    double max_gap = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i + 1 < c_samples.second.size(); ++i) {
      double gap =
        glm::distance(c_samples.second[i], c_samples.second[i + 1]) - 1.0;
      gap = std::max(gap, 0.0);
      max_gap = std::max(gap, max_gap);
    }
    if (max_gap > -std::numeric_limits<double>::infinity()) {
      max_gap = std::max(max_gap, 0.0);
      cluster_max_gap.emplace_back(max_gap);
    }
  }
  std::sort(cluster_max_gap.begin(), cluster_max_gap.end());

  std::vector<double> euclidean_clusterwise;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]) || is_point_noisy(xsec.points[i + 1]))
      continue;
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      double dis =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      dis = std::max(dis, 0.0);
      euclidean_clusterwise.emplace_back(dis);
      // The two clusters are not clearly divided
      if (euclidean_clusterwise.size() > 1) {
        euclidean_clusterwise.clear();
        euclidean_clusterwise.emplace_back(0);
        break;
      }
    }
  }
  if (!euclidean_clusterwise.empty()) {
    if (euclidean_clusterwise[0] != 0) {
      if (!cluster_max_gap.empty())
        euclidean_clusterwise[0] =
          cluster_max_gap.back() / euclidean_clusterwise[0];
    } else {
      euclidean_clusterwise[0] = sandwiched_weight;
    }

    // The corresponding cluster (or both) only have a single sample locally
    if (cluster_max_gap.empty())
      euclidean_clusterwise[0] = -1;
  }
  return euclidean_clusterwise;
}
static std::vector<double> get_xsec_min_gap_ratio(const Cluster::XSec &xsec,
                                                  int cluster_idx) {
  std::map<size_t, std::vector<glm::dvec2>> cluster_samples;
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      continue;
    cluster_samples[xsec.points[i].cluster_ind].emplace_back(
      xsec.points[i].point);
  }

  std::vector<double> cluster_max_gap;
  for (const auto &c_samples : cluster_samples) {
    if (cluster_idx >= 0 && c_samples.first != cluster_idx)
      continue;
    double max_gap = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i + 1 < c_samples.second.size(); ++i) {
      double gap =
        glm::distance(c_samples.second[i], c_samples.second[i + 1]) - 1.0;
      gap = std::max(gap, 0.0);
      max_gap = std::max(gap, max_gap);
    }
    if (max_gap > -std::numeric_limits<double>::infinity()) {
      max_gap = std::max(max_gap, 0.0);
      cluster_max_gap.emplace_back(max_gap);
    }
  }
  std::sort(cluster_max_gap.begin(), cluster_max_gap.end());

  std::vector<double> euclidean_clusterwise;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]) || is_point_noisy(xsec.points[i + 1]))
      continue;
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      double dis =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      dis = std::max(dis, 0.0);
      euclidean_clusterwise.emplace_back(dis);
      // The two clusters are not clearly divided
      if (euclidean_clusterwise.size() > 1) {
        euclidean_clusterwise.clear();
        euclidean_clusterwise.emplace_back(0);
        break;
      }
    }
  }
  if (!euclidean_clusterwise.empty()) {
    if (euclidean_clusterwise[0] != 0) {
      if (!cluster_max_gap.empty())
        euclidean_clusterwise[0] =
          cluster_max_gap.front() / euclidean_clusterwise[0];
    } else {
      euclidean_clusterwise[0] = sandwiched_weight;
    }

    // The corresponding cluster (or both) only have a single sample locally
    if (cluster_max_gap.empty())
      euclidean_clusterwise[0] = std::numeric_limits<double>::infinity();
  }
  return euclidean_clusterwise;
}

static std::vector<double> get_xsec_max_gap(const Cluster::XSec &xsec,
                                            int cluster_idx) {
  std::map<size_t, std::vector<glm::dvec2>> cluster_samples;
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      continue;
    if (cluster_idx >= 0)
      cluster_samples[xsec.points[i].cluster_ind].emplace_back(
        xsec.points[i].point);
    else
      cluster_samples[0].emplace_back(xsec.points[i].point);
  }

  std::vector<double> cluster_max_gap;
  for (const auto &c_samples : cluster_samples) {
    if (cluster_idx >= 0 && c_samples.first != cluster_idx)
      continue;
    double max_gap = -std::numeric_limits<double>::infinity();
    for (size_t i = 0; i + 1 < c_samples.second.size(); ++i) {
      double gap =
        glm::distance(c_samples.second[i], c_samples.second[i + 1]) - 1.0;
      gap = std::max(gap, 0.0);
      max_gap = std::max(gap, max_gap);
    }
    if (max_gap > -std::numeric_limits<double>::infinity()) {
      max_gap = std::max(max_gap, 0.0);
      cluster_max_gap.emplace_back(max_gap);
    }
  }
  std::sort(cluster_max_gap.begin(), cluster_max_gap.end());

  if (cluster_max_gap.empty())
    return std::vector<double>{std::numeric_limits<double>::infinity()};

  return std::vector<double>{cluster_max_gap.back()};
}

static std::vector<double> get_xsec_median_gap_ratio(const Cluster::XSec &xsec,
                                                     int cluster_idx) {
  std::map<size_t, std::vector<glm::dvec2>> cluster_samples;
  for (size_t i = 0; i < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]))
      continue;
    cluster_samples[xsec.points[i].cluster_ind].emplace_back(
      xsec.points[i].point);
  }

  std::vector<double> cluster_max_gap;
  for (const auto &c_samples : cluster_samples) {
    if (cluster_idx >= 0 && c_samples.first != cluster_idx)
      continue;
    for (size_t i = 0; i + 1 < c_samples.second.size(); ++i) {
      double gap =
        glm::distance(c_samples.second[i], c_samples.second[i + 1]) - 1.0;
      gap = std::max(gap, 0.0);
      cluster_max_gap.emplace_back(gap);
    }
  }
  if (cluster_max_gap.empty())
    return std::vector<double>();
  auto m = cluster_max_gap.begin() + cluster_max_gap.size() / 2;
  std::nth_element(cluster_max_gap.begin(), m, cluster_max_gap.end());
  double gap_size = cluster_max_gap[cluster_max_gap.size() / 2];

  std::vector<double> euclidean_clusterwise;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (is_point_noisy(xsec.points[i]) || is_point_noisy(xsec.points[i + 1]))
      continue;
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      double dis =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      dis = std::max(dis, 0.0);
      euclidean_clusterwise.emplace_back(dis);
      // The two clusters are not clearly divided
      if (euclidean_clusterwise.size() > 1) {
        euclidean_clusterwise.clear();
        euclidean_clusterwise.emplace_back(0);
        break;
      }
    }
  }
  if (!euclidean_clusterwise.empty()) {
    if (euclidean_clusterwise[0] != 0) {
      if (!cluster_max_gap.empty())
        euclidean_clusterwise[0] = gap_size / euclidean_clusterwise[0];
    } else {
      euclidean_clusterwise[0] = sandwiched_weight;
    }

    // The corresponding cluster (or both) only have a single sample locally
    if (cluster_max_gap.empty())
      euclidean_clusterwise.clear();
  }
  return euclidean_clusterwise;
}

static std::vector<double>
get_xsec_max_individual_width_local(const Cluster::XSec &xsec) {
  std::vector<double> cluster_width = get_xsec_individual_widths(xsec, -1);
  if (cluster_width.empty())
    return std::vector<double>{1.0};
  return std::vector<double>{
    *std::max_element(cluster_width.begin(), cluster_width.end())};
}
static std::vector<double>
get_xsec_min_individual_width_local(const Cluster::XSec &xsec) {
  std::vector<double> cluster_width = get_xsec_individual_widths(xsec, -1);
  if (cluster_width.empty())
    return std::vector<double>{1.0};
  return std::vector<double>{
    *std::min_element(cluster_width.begin(), cluster_width.end())};
}

static double get_xsecs_mean(
  span<Cluster::XSec const> const xsecs,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  std::vector<double> measure_vec;
  measure_vec.reserve(xsecs.size());
  for (auto &xsec : xsecs) {
    const auto measures = point_measure(xsec);
    if (!measures.empty())
      measure_vec.insert(measure_vec.end(), measures.begin(), measures.end());
  }
  if (measure_vec.empty())
    return std::numeric_limits<double>::infinity();
  double sum = 0;
  for (const auto v : measure_vec)
    sum += v;
  return sum / measure_vec.size();
}
static double get_xsecs_mean(
  const std::vector<Cluster::XSec> &xsecs,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  return get_xsecs_mean(span<Cluster::XSec const>{xsecs.data(), xsecs.size()},
                        point_measure);
}
static double get_xsecs_percentile(
  span<Cluster::XSec const> const xsecs, double percentile,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  std::vector<double> measure_vec;
  measure_vec.reserve(xsecs.size());
  for (auto &xsec : xsecs) {
    const auto measures = point_measure(xsec);
    if (!measures.empty())
      measure_vec.insert(measure_vec.end(), measures.begin(), measures.end());
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
static double get_xsecs_percentile(
  const std::vector<Cluster::XSec> &xsecs, double percentile,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  return get_xsecs_percentile(
    span<Cluster::XSec const>{xsecs.data(), xsecs.size()}, percentile,
    point_measure);
}

std::vector<double> get_all_xsecs_binned(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<double(span<Cluster::XSec const> const)> func_binned) {
  std::vector<double> bin_vec;
  std::vector<double> bin_w_vec;
  bin_vec.resize(num_binned, std::numeric_limits<double>::infinity());
  bin_w_vec.resize(num_binned, std::numeric_limits<double>::infinity());

  if (xsecs.empty())
    return bin_vec;

  std::vector<std::pair<double, double>> bin_ranges;
  bin_ranges.reserve(num_binned);
  double step_size = 1.0 / num_binned;
  for (size_t i = 0; i < num_binned; ++i) {
    bin_ranges.emplace_back(i * step_size, (i + 1) * step_size);
  }

  auto aggregated_w_func = [](span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_mean(xsecs, get_xsec_merged_width);
  };

  for (size_t i = 0; i < num_binned; ++i) {
    size_t start = std::floor(xsecs.size() * bin_ranges[i].first);
    size_t end = std::floor(xsecs.size() * bin_ranges[i].second);
    bin_vec[i] =
      func_binned(span<Cluster::XSec const>{xsecs.data() + start, end - start});
    bin_w_vec[i] = aggregated_w_func(
      span<Cluster::XSec const>{xsecs.data() + start, end - start});
  }

  // Order the bins in the same way, so it's independent to the U values.
  double first_v = -1, last_v = -1;
  bool seen_first = false;
  for (size_t i = 0; i < bin_w_vec.size(); ++i) {
    if (bin_w_vec[i] != std::numeric_limits<double>::infinity()) {
      if (!seen_first) {
        first_v = bin_w_vec[i];
        seen_first = true;
      }
      last_v = bin_w_vec[i];
    }
  }

  if (first_v > last_v)
    std::reverse(bin_vec.begin(), bin_vec.end());

  return bin_vec;
}
std::vector<double> get_all_xsecs_binned_mean(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  auto aggregated_func =
    [&point_measure](span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_mean(xsecs, point_measure);
  };
  return get_all_xsecs_binned(xsecs, num_binned, aggregated_func);
}
std::vector<double> get_all_xsecs_binned_percentile(
  const std::vector<Cluster::XSec> &xsecs, double percentile, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  auto aggregated_func = [&point_measure, &percentile](
                           span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_percentile(xsecs, percentile, point_measure);
  };
  return get_all_xsecs_binned(xsecs, num_binned, aggregated_func);
}

static std::vector<double> get_overlapping_xsecs_binned(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<double(span<Cluster::XSec const> const)> func_binned) {
  std::vector<Cluster::XSec> overlapping_xsecs;
  overlapping_xsecs.reserve(xsecs.size());
  for (const auto &xsec : xsecs)
    if (is_overlapping(xsec))
      overlapping_xsecs.emplace_back(xsec);

  std::vector<double> bin_vec;
  std::vector<double> bin_w_vec;
  bin_vec.resize(num_binned, std::numeric_limits<double>::infinity());
  bin_w_vec.resize(num_binned, std::numeric_limits<double>::infinity());

  if (overlapping_xsecs.empty())
    return bin_vec;

  std::vector<std::pair<double, double>> bin_ranges;
  bin_ranges.reserve(num_binned);
  double step_size = 1.0 / num_binned;
  for (size_t i = 0; i < num_binned; ++i) {
    bin_ranges.emplace_back(i * step_size, (i + 1) * step_size);
  }

  auto aggregated_w_func = [](span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_mean(xsecs, get_xsec_merged_width);
  };

  for (size_t i = 0; i < num_binned; ++i) {
    size_t start = std::floor(overlapping_xsecs.size() * bin_ranges[i].first);
    size_t end = std::floor(overlapping_xsecs.size() * bin_ranges[i].second);
    if (end == start)
      end++;
    bin_vec[i] = func_binned(
      span<Cluster::XSec const>{overlapping_xsecs.data() + start, end - start});
    bin_w_vec[i] = aggregated_w_func(
      span<Cluster::XSec const>{overlapping_xsecs.data() + start, end - start});
  }

  // Order the bins in the same way, so it's independent to the U values.
  double first_v = -1, last_v = -1;
  bool seen_first = false;
  for (size_t i = 0; i < bin_w_vec.size(); ++i) {
    if (bin_w_vec[i] != std::numeric_limits<double>::infinity()) {
      if (!seen_first) {
        first_v = bin_w_vec[i];
        seen_first = true;
      }
      last_v = bin_w_vec[i];
    }
  }

  if (first_v > last_v)
    std::reverse(bin_vec.begin(), bin_vec.end());

  return bin_vec;
}
std::vector<double> get_overlapping_xsecs_binned_mean(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  auto aggregated_func =
    [&point_measure](span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_mean(xsecs, point_measure);
  };
  return get_overlapping_xsecs_binned(xsecs, num_binned, aggregated_func);
}

std::vector<double> get_overlapping_xsecs_binned_percentile(
  const std::vector<Cluster::XSec> &xsecs, double percentile, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  auto aggregated_func = [&point_measure, &percentile](
                           span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_percentile(xsecs, percentile, point_measure);
  };
  return get_overlapping_xsecs_binned(xsecs, num_binned, aggregated_func);
}

static std::vector<double> get_sep_overlapping_xsecs_binned(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<double(span<Cluster::XSec const> const)> func_binned) {
  std::vector<Cluster::XSec> overlapping_xsecs;
  overlapping_xsecs.reserve(xsecs.size());
  for (const auto &xsec : xsecs)
    if (is_overlapping(xsec) && is_separatible(xsec))
      overlapping_xsecs.emplace_back(xsec);

  std::vector<double> bin_vec;
  std::vector<double> bin_w_vec;
  bin_vec.resize(num_binned, std::numeric_limits<double>::infinity());
  bin_w_vec.resize(num_binned, std::numeric_limits<double>::infinity());

  if (overlapping_xsecs.empty())
    return bin_vec;

  std::vector<std::pair<double, double>> bin_ranges;
  bin_ranges.reserve(num_binned);
  double step_size = 1.0 / num_binned;
  for (size_t i = 0; i < num_binned; ++i) {
    bin_ranges.emplace_back(i * step_size, (i + 1) * step_size);
  }

  auto aggregated_w_func = [](span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_mean(xsecs, get_xsec_merged_width);
  };

  for (size_t i = 0; i < num_binned; ++i) {
    size_t start = std::floor(overlapping_xsecs.size() * bin_ranges[i].first);
    size_t end = std::floor(overlapping_xsecs.size() * bin_ranges[i].second);
    if (end == start)
      end++;
    bin_vec[i] = func_binned(
      span<Cluster::XSec const>{overlapping_xsecs.data() + start, end - start});
    bin_w_vec[i] = aggregated_w_func(
      span<Cluster::XSec const>{overlapping_xsecs.data() + start, end - start});
  }

  // Order the bins in the same way, so it's independent to the U values.
  double first_v = -1, last_v = -1;
  bool seen_first = false;
  for (size_t i = 0; i < bin_w_vec.size(); ++i) {
    if (bin_w_vec[i] != std::numeric_limits<double>::infinity()) {
      if (!seen_first) {
        first_v = bin_w_vec[i];
        seen_first = true;
      }
      last_v = bin_w_vec[i];
    }
  }

  if (first_v > last_v)
    std::reverse(bin_vec.begin(), bin_vec.end());

  return bin_vec;
}
std::vector<double> get_sep_overlapping_xsecs_binned_mean(
  const std::vector<Cluster::XSec> &xsecs, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  auto aggregated_func =
    [&point_measure](span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_mean(xsecs, point_measure);
  };
  return get_sep_overlapping_xsecs_binned(xsecs, num_binned, aggregated_func);
}

std::vector<double> get_sep_overlapping_xsecs_binned_percentile(
  const std::vector<Cluster::XSec> &xsecs, double percentile, size_t num_binned,
  std::function<std::vector<double>(const Cluster::XSec &)> point_measure) {
  auto aggregated_func = [&point_measure, &percentile](
                           span<Cluster::XSec const> const xsecs) -> double {
    return get_xsecs_percentile(xsecs, percentile, point_measure);
  };
  return get_sep_overlapping_xsecs_binned(xsecs, num_binned, aggregated_func);
}

static std::pair<double, double>
get_xsecs_width_changes(const std::vector<Cluster::XSec> &xsecs) {
  std::vector<double> change_ratios;
  std::pair<double, double> ret_change_ratios;

  if (end_noisy_range > 0) {
    for (size_t i = 0; i + 2 < xsecs.size(); ++i) {
      size_t j = i + 2;
      // For now, skip the boundary xsecs of a side-by-side section
      if (is_overlapping(xsecs[i + 1]) &&
          is_overlapping(xsecs[i]) != is_overlapping(xsecs[j])) {
        auto xsec_widths1 = get_xsec_merged_width(xsecs[i]);
        auto xsec_widths2 = get_xsec_merged_width(xsecs[j]);
        assert(!xsec_widths1.empty() && !xsec_widths2.empty());
        double xsec_width1 = xsec_widths1.front();
        double xsec_width2 = xsec_widths2.front();
        change_ratios.emplace_back((is_overlapping(xsecs[j])) ? xsec_width2
                                                              : xsec_width1);
      }
    }
  } else {
    for (size_t i = 0; i + 1 < xsecs.size(); ++i) {
      if (is_overlapping(xsecs[i]) != is_overlapping(xsecs[i + 1])) {
        double xsec_width1 = glm::distance(xsecs[i].points.front().point,
                                           xsecs[i].points.back().point) +
                             1.0;
        double xsec_width2 = glm::distance(xsecs[i + 1].points.front().point,
                                           xsecs[i + 1].points.back().point) +
                             1.0;
        change_ratios.emplace_back(
          (is_overlapping(xsecs[i + 1])) ? xsec_width2 : xsec_width1);
      }
    }
  }
  if (!change_ratios.empty()) {
    ret_change_ratios.first = change_ratios.front();
    if (change_ratios.size() > 1)
      ret_change_ratios.second = change_ratios.back();
    else
      ret_change_ratios.second = 0;
  } else {
    ret_change_ratios.first = ret_change_ratios.second = 0;
  }

  return ret_change_ratios;
}

static std::pair<double, double>
get_xsecs_width_change_ratios(const std::vector<Cluster::XSec> &xsecs) {
  std::vector<double> change_ratios;
  std::pair<double, double> ret_change_ratios;
  if (end_noisy_range > 0) {
    for (size_t i = 0; i + 2 < xsecs.size(); ++i) {
      size_t j = i + 2;
      // For now, skip the boundary xsecs of a side-by-side section
      if (is_overlapping(xsecs[i + 1]) &&
          is_overlapping(xsecs[i]) != is_overlapping(xsecs[j])) {
        auto xsec_widths1 = get_xsec_merged_width(xsecs[i]);
        auto xsec_widths2 = get_xsec_merged_width(xsecs[j]);
        assert(!xsec_widths1.empty() && !xsec_widths2.empty());
        double xsec_width1 = xsec_widths1.front();
        double xsec_width2 = xsec_widths2.front();
        change_ratios.emplace_back((is_overlapping(xsecs[j]))
                                     ? xsec_width1 / xsec_width2
                                     : xsec_width2 / xsec_width1);
      }
    }
  } else {
    for (size_t i = 0; i + 1 < xsecs.size(); ++i) {
      if (is_overlapping(xsecs[i]) != is_overlapping(xsecs[i + 1])) {
        double xsec_width1 = glm::distance(xsecs[i].points.front().point,
                                           xsecs[i].points.back().point) +
                             1.0;
        double xsec_width2 = glm::distance(xsecs[i + 1].points.front().point,
                                           xsecs[i + 1].points.back().point) +
                             1.0;
        change_ratios.emplace_back((is_overlapping(xsecs[i + 1]))
                                     ? xsec_width1 / xsec_width2
                                     : xsec_width2 / xsec_width1);
      }
    }
  }
  if (!change_ratios.empty()) {
    ret_change_ratios.first = change_ratios.front();
    if (change_ratios.size() > 1)
      ret_change_ratios.second = change_ratios.back();
    else
      ret_change_ratios.second = 1;
  } else {
    ret_change_ratios.first = ret_change_ratios.second = 1;
  }

  return ret_change_ratios;
}

static std::vector<double>
get_overlapping_xsecs_width(const Cluster::XSec &xsec) {
  std::vector<double> overlapping_xsecs_widths;

  if (is_overlapping(xsec)) {
    std::vector<size_t> non_noisy_indices;
    for (size_t i = 0; i < xsec.points.size(); ++i) {
      if (is_point_noisy(xsec.points[i]))
        continue;
      non_noisy_indices.emplace_back(i);
    }
    if (!non_noisy_indices.empty()) {
      double xsec_width =
        glm::distance(xsec.points[non_noisy_indices.front()].point,
                      xsec.points[non_noisy_indices.back()].point) +
        1.0;
      overlapping_xsecs_widths.emplace_back(xsec_width);
    }
  }

  return overlapping_xsecs_widths;
}

static std::pair<std::pair<double, double>, std::pair<double, double>>
get_xsecs_stepaway_ratios(const std::vector<Cluster::XSec> &xsecs) {
  std::pair<std::pair<double, double>, std::pair<double, double>>
    stepaway_ratios;
  std::vector<std::pair<double, double>> overlapping_xsecs_widths_distances;

  for (auto const &xsec : xsecs) {
    if (is_overlapping(xsec)) {
      std::vector<size_t> non_noisy_indices;
      for (size_t i = 0; i < xsec.points.size(); ++i) {
        if (is_point_noisy(xsec.points[i]))
          continue;
        non_noisy_indices.emplace_back(i);
      }
      if (!non_noisy_indices.empty()) {
        double xsec_width =
          glm::distance(xsec.points[non_noisy_indices.front()].point,
                        xsec.points[non_noisy_indices.back()].point) +
          1.0;
        std::vector<double> dist = get_xsec_euclidean(xsec);
        if (dist.empty())
          continue;
        overlapping_xsecs_widths_distances.emplace_back(xsec_width, dist[0]);
      }
    }
  }

  if (overlapping_xsecs_widths_distances.size() < 2)
    return std::pair<std::pair<double, double>, std::pair<double, double>>(
      std::make_pair(1, 1), std::make_pair(1, 1));

  // Assume the overlapping section is continuous
  auto get_head_stepaway_ratio =
    [&overlapping_xsecs_widths_distances]() -> double {
    double head_dist = overlapping_xsecs_widths_distances.front().second;
    size_t stepaway_steps = std::ceil(
      overlapping_xsecs_widths_distances.front().first / stroke_sampling_size);
    double stepaway_dist =
      overlapping_xsecs_widths_distances
        [std::min(stepaway_steps,
                  overlapping_xsecs_widths_distances.size() - 1)]
          .second;

    return (head_dist > 0) ? std::min(1.0, stepaway_dist / head_dist) : 1.0;
  };
  stepaway_ratios.first =
    std::make_pair(overlapping_xsecs_widths_distances.front().first,
                   get_head_stepaway_ratio());
  std::reverse(overlapping_xsecs_widths_distances.begin(),
               overlapping_xsecs_widths_distances.end());
  stepaway_ratios.second =
    std::make_pair(overlapping_xsecs_widths_distances.front().first,
                   get_head_stepaway_ratio());

  return stepaway_ratios;
}

static std::pair<std::pair<double, double>, std::pair<double, double>>
get_xsecs_stepaway_ratios2(const std::vector<Cluster::XSec> &xsecs) {
  std::pair<std::pair<double, double>, std::pair<double, double>>
    stepaway_ratios;
  std::vector<std::pair<double, double>> overlapping_xsecs_widths_distances;

  for (auto const &xsec : xsecs) {
    if (is_overlapping(xsec)) {
      std::vector<size_t> non_noisy_indices;
      for (size_t i = 0; i < xsec.points.size(); ++i) {
        if (is_point_noisy(xsec.points[i]))
          continue;
        non_noisy_indices.emplace_back(i);
      }
      if (!non_noisy_indices.empty()) {
        double xsec_width =
          glm::distance(xsec.points[non_noisy_indices.front()].point,
                        xsec.points[non_noisy_indices.back()].point) +
          1.0;
        std::vector<double> dist = get_xsec_euclidean(xsec);
        if (dist.empty())
          continue;
        overlapping_xsecs_widths_distances.emplace_back(xsec_width, dist[0]);
      }
    }
  }

  if (overlapping_xsecs_widths_distances.size() < 2)
    return std::pair<std::pair<double, double>, std::pair<double, double>>(
      std::make_pair(1, 1), std::make_pair(1, 1));

  // Assume the overlapping section is continuous
  auto get_head_stepaway_ratio =
    [&overlapping_xsecs_widths_distances]() -> double {
    double head_dist = overlapping_xsecs_widths_distances.front().first;
    size_t stepaway_steps = std::ceil(
      overlapping_xsecs_widths_distances.front().first / stroke_sampling_size);
    double stepaway_dist =
      overlapping_xsecs_widths_distances
        [std::min(stepaway_steps,
                  overlapping_xsecs_widths_distances.size() - 1)]
          .first;

    return (head_dist > 0) ? std::min(1.0, stepaway_dist / head_dist) : 1.0;
  };
  stepaway_ratios.first =
    std::make_pair(overlapping_xsecs_widths_distances.front().first,
                   get_head_stepaway_ratio());
  std::reverse(overlapping_xsecs_widths_distances.begin(),
               overlapping_xsecs_widths_distances.end());
  stepaway_ratios.second =
    std::make_pair(overlapping_xsecs_widths_distances.front().first,
                   get_head_stepaway_ratio());

  return stepaway_ratios;
}

static std::pair<std::pair<double, double>, std::pair<double, double>>
get_xsecs_stepaway_angles(const std::vector<Cluster::XSec> &xsecs) {
  std::pair<std::pair<double, double>, std::pair<double, double>>
    stepaway_angles;
  std::vector<std::pair<double, Cluster::XSec>> overlapping_xsecs_widths;

  for (auto const &xsec : xsecs) {
    if (is_overlapping(xsec)) {
      std::vector<size_t> non_noisy_indices;
      for (size_t i = 0; i < xsec.points.size(); ++i) {
        if (is_point_noisy(xsec.points[i]))
          continue;
        non_noisy_indices.emplace_back(i);
      }
      if (!non_noisy_indices.empty()) {
        double xsec_width =
          glm::distance(xsec.points[non_noisy_indices.front()].point,
                        xsec.points[non_noisy_indices.back()].point) +
          1.0;
        overlapping_xsecs_widths.emplace_back(xsec_width, xsec);
      }
    }
  }

  if (overlapping_xsecs_widths.size() < 2)
    return std::pair<std::pair<double, double>, std::pair<double, double>>(
      std::make_pair(1, 0), std::make_pair(1, 0));

  auto get_xsec_avg_positions =
    [](const Cluster::XSec &xsec) -> std::map<size_t, Cluster::XSecPoint> {
    std::map<size_t, Cluster::XSecPoint> c2avg_pos;

    c2avg_pos[0].stroke_idx = 0;
    c2avg_pos[0].point = glm::dvec2(0, 0);
    c2avg_pos[1].stroke_idx = 0;
    c2avg_pos[1].point = glm::dvec2(0, 0);
    for (const auto &p : xsec.points) {
      c2avg_pos[p.cluster_ind].stroke_idx++;
      c2avg_pos[p.cluster_ind].point += p.point;
    }
    c2avg_pos[0].point /= c2avg_pos[0].stroke_idx;
    c2avg_pos[1].point /= c2avg_pos[1].stroke_idx;
    return c2avg_pos;
  };

  // Assume the overlapping section is continuous
  auto get_head_stepaway_angle = [&overlapping_xsecs_widths,
                                  &get_xsec_avg_positions]() -> double {
    auto head_xsec = overlapping_xsecs_widths.front().second;
    size_t stepaway_steps =
      std::ceil(overlapping_xsecs_widths.front().first / stroke_sampling_size);
    auto stepaway_xsec =
      overlapping_xsecs_widths[std::min(stepaway_steps,
                                        overlapping_xsecs_widths.size() - 1)]
        .second;
    auto head_pos = get_xsec_avg_positions(head_xsec);
    auto stepaway_pos = get_xsec_avg_positions(stepaway_xsec);
    double sign = (glm::distance(head_pos[0].point, head_pos[1].point) >
                   glm::distance(stepaway_pos[0].point, stepaway_pos[1].point))
                    ? 1
                    : -1;
    glm::dvec2 tangent0 =
      glm::normalize(head_pos[0].point - stepaway_pos[0].point);
    glm::dvec2 tangent1 =
      glm::normalize(head_pos[1].point - stepaway_pos[1].point);
    double dot = glm::dot(tangent0, tangent1);
    dot = std::min(std::max(dot, -1.0), 1.0);
    double angle = glm::acos(dot) * 180.0 / M_PI;
    return sign * angle;
  };
  stepaway_angles.first = std::make_pair(overlapping_xsecs_widths.front().first,
                                         get_head_stepaway_angle());
  std::reverse(overlapping_xsecs_widths.begin(),
               overlapping_xsecs_widths.end());
  stepaway_angles.second = std::make_pair(
    overlapping_xsecs_widths.front().first, get_head_stepaway_angle());

  return stepaway_angles;
}

static std::pair<double, double>
get_xsecs_avg_width_non_overlapping(const std::vector<Cluster::XSec> &xsecs,
                                    bool periodic) {
  std::map<size_t, int> sort_indices;
  // 0: non-side-by-side is in the descending direction;
  // 1: non-side-by-side is in the ascending direction.
  sort_indices[0] = -1;
  sort_indices[1] = -1;

  std::vector<size_t> overlapping_interval;
  for (size_t i = 0;; ++i) {
    size_t next = (i + 1) % xsecs.size();
    if (is_overlapping(xsecs[i]) != is_overlapping(xsecs[next])) {
      if (!is_overlapping(xsecs[i]))
        sort_indices[0] = i;
      else {
        sort_indices[1] = next;
      }
    }
    if (sort_indices[0] >= 0 && sort_indices[1] >= 0)
      break;
    if (!periodic && i + 2 >= xsecs.size())
      break;
    else if (periodic && i + 1 >= xsecs.size())
      break;
  }
  overlapping_interval.emplace_back(sort_indices[0]);
  overlapping_interval.emplace_back(sort_indices[1]);

  if (overlapping_interval[0] < 0 && overlapping_interval[1] < 0)
    return std::make_pair(-1, -1);
  assert(!periodic || overlapping_interval.size() == 2);

  std::vector<double> widths;
  std::map<size_t, double> c2w;
  int cid1 = -1;
  for (int i = overlapping_interval[0]; i >= 0;) {
    if (!xsecs[i].points.empty() && !is_overlapping(xsecs[i])) {
      cid1 = xsecs[i].points.front().cluster_ind;
      std::vector<double> local_width =
        get_xsec_individual_widths(xsecs[i], cid1);
      if (!local_width.empty()) {
        widths.emplace_back(local_width.front());
      }
    }
    if (i == 0 && !periodic) {
      break;
    }
    if (periodic && i == overlapping_interval[1])
      break;

    i = (i - 1 + xsecs.size()) % xsecs.size();
  }
  // Non-overlapping exists on one side
  if (cid1 >= 0) {
    if (!widths.empty())
      c2w[cid1] =
        std::accumulate(widths.begin(), widths.end(), 0.0) / widths.size();
    else
      c2w[cid1] = -1;
  }

  if (overlapping_interval.size() == 2) {
    widths.clear();
    int cid2 = -1;
    for (int i = overlapping_interval[1]; i >= 0;) {
      if (!xsecs[i].points.empty() && !is_overlapping(xsecs[i])) {
        cid2 = xsecs[i].points.front().cluster_ind;
        std::vector<double> local_width =
          get_xsec_individual_widths(xsecs[i], cid2);
        if (!local_width.empty()) {
          widths.emplace_back(local_width.front());
        }
      }
      if (i + 1 == xsecs.size() && !periodic) {
        break;
      }
      if (periodic && i == overlapping_interval[0])
        break;

      i = (i + 1) % xsecs.size();
    }
    // For the case that one subcluster is fully inside the other, we
    // distinguish the left and right non-side-by-side by assigning it a
    // different cluster index
    if (cid2 >= 0) {
      if (c2w.count(cid2))
        cid2 = 1 - cid2;
      if (!widths.empty())
        c2w[cid2] =
          std::accumulate(widths.begin(), widths.end(), 0.0) / widths.size();
      else
        c2w[cid2] = -1;
    }
  }

  // Fully overlapping
  if (!c2w.count(0) && !c2w.count(1))
    return std::pair<double, double>(-1, -1);
  return std::make_pair((c2w.count(0)) ? c2w[0] : c2w[1],
                        (c2w.count(1)) ? c2w[1] : c2w[0]);
}

static std::vector<double>
get_xsecs_all_width_changes(const std::vector<Cluster::XSec> &xsecs) {
  std::vector<double> change_values;
  for (size_t i = 0; i + 1 < xsecs.size(); ++i) {
    auto xsec_widths1 = get_xsec_merged_width(xsecs[i]);
    auto xsec_widths2 = get_xsec_merged_width(xsecs[i + 1]);
    if (!xsec_widths1.empty() && !xsec_widths2.empty()) {
      double xsec_width1 = xsec_widths1.front();
      double xsec_width2 = xsec_widths2.front();
      change_values.emplace_back(std::abs(xsec_width1 - xsec_width2));
    }
  }

  return change_values;
}

std::vector<double>
MeanAngleStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                          const FittedCurve &fit1,
                                          const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(filtered_xsecs, Feature::num_binned_,
                                           get_xsec_angle);
}

std::vector<double>
PercentileAngleStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                                const FittedCurve &fit1,
                                                const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_angle);
}

std::vector<double>
MeanDistanceStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                             const FittedCurve &fit1,
                                             const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(filtered_xsecs, Feature::num_binned_,
                                           get_xsec_euclidean);
}

std::vector<double> PercentileDistanceStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_euclidean);
}

std::vector<double>
MeanWidthStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                          const FittedCurve &fit1,
                                          const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(filtered_xsecs, Feature::num_binned_,
                                           get_xsec_merged_width);
}

std::vector<double>
PercentileWidthStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                                const FittedCurve &fit1,
                                                const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_merged_width);
}

std::vector<double>
MeanClusterwiseDistanceRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_, get_xsec_euclidean_clusterwise_ratio);
}

std::vector<double>
PercentileClusterwiseDistanceRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_euclidean_clusterwise_ratio);
}

std::vector<double>
MeanClusterwiseDistanceMaxWidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_max_width_ratio);
}

std::vector<double>
PercentileClusterwiseDistanceMaxWidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_max_width_ratio);
}

std::vector<double>
MeanClusterwiseDistanceMinWidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_min_width_ratio);
}

std::vector<double>
PercentileClusterwiseDistanceMinWidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_min_width_ratio);
}

std::vector<double>
MeanClusterwiseDistanceMaxC12WidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileIndividualWidth1StrokeClusterFeature w1_fea(0.5, 1);
  double w1 = w1_fea.operator()(merged_cluster, fit1, fit2)[0];
  PercentileIndividualWidth2StrokeClusterFeature w2_fea(0.5, 1);
  double w2 = w2_fea.operator()(merged_cluster, fit1, fit2)[0];

  bool max_c1 = w1 > w2;

  int c = (max_c1) ? 0 : 1;
  auto get_xsec_clusterwise_euclidean_width_ratio =
    [&c](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width =
      get_xsec_clusterwise_euclidean_cluster_width_ratio(xsec, c);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_width_ratio);
}

std::vector<double>
PercentileClusterwiseDistanceMaxC12WidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileIndividualWidth1StrokeClusterFeature w1_fea(0.5, 1);
  double w1 = w1_fea.operator()(merged_cluster, fit1, fit2)[0];
  PercentileIndividualWidth2StrokeClusterFeature w2_fea(0.5, 1);
  double w2 = w2_fea.operator()(merged_cluster, fit1, fit2)[0];

  bool max_c1 = w1 > w2;

  int c = (max_c1) ? 0 : 1;
  auto get_xsec_clusterwise_euclidean_width_ratio =
    [&c](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width =
      get_xsec_clusterwise_euclidean_cluster_width_ratio(xsec, c);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_width_ratio);
}

std::vector<double>
MeanClusterwiseDistanceMinC12WidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileIndividualWidth1StrokeClusterFeature w1_fea(0.5, 1);
  double w1 = w1_fea.operator()(merged_cluster, fit1, fit2)[0];
  PercentileIndividualWidth2StrokeClusterFeature w2_fea(0.5, 1);
  double w2 = w2_fea.operator()(merged_cluster, fit1, fit2)[0];

  bool min_c1 = w1 < w2;

  int c = (min_c1) ? 0 : 1;
  auto get_xsec_clusterwise_euclidean_width_ratio =
    [&c](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width =
      get_xsec_clusterwise_euclidean_cluster_width_ratio(xsec, c);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_width_ratio);
}

std::vector<double>
PercentileClusterwiseDistanceMinC12WidthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileIndividualWidth1StrokeClusterFeature w1_fea(0.5, 1);
  double w1 = w1_fea.operator()(merged_cluster, fit1, fit2)[0];
  PercentileIndividualWidth2StrokeClusterFeature w2_fea(0.5, 1);
  double w2 = w2_fea.operator()(merged_cluster, fit1, fit2)[0];

  bool min_c1 = w1 < w2;

  int c = (min_c1) ? 0 : 1;
  auto get_xsec_clusterwise_euclidean_width_ratio =
    [&c](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width =
      get_xsec_clusterwise_euclidean_cluster_width_ratio(xsec, c);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_clusterwise_euclidean_width_ratio);
}

std::vector<double>
PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  int c1 = 0;
  auto get_xsec_max_gap_ratio_c1 = [&c1](const Cluster::XSec &xsec) {
    return get_xsec_max_gap_ratio(xsec, c1);
  };
  auto get_xsec_min_gap_ratio_c1 = [&c1](const Cluster::XSec &xsec) {
    return get_xsec_min_gap_ratio(xsec, c1);
  };
  auto get_xsec_median_gap_ratio_c1 = [&c1](const Cluster::XSec &xsec) {
    return get_xsec_median_gap_ratio(xsec, c1);
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  if (percentile_ == 1.0)
    return get_overlapping_xsecs_binned_percentile(filtered_xsecs, percentile_,
                                                   Feature::num_binned_,
                                                   get_xsec_max_gap_ratio_c1);
  else if (percentile_ == 0.0) {
    auto ratio = get_overlapping_xsecs_binned_percentile(
      filtered_xsecs, percentile_, Feature::num_binned_,
      get_xsec_min_gap_ratio_c1);
    // infinity means all cross-sections only contain a single sample
    set_inf(ratio, -1);
    return ratio;
  } else {
    auto ratio = get_overlapping_xsecs_binned_percentile(
      filtered_xsecs, percentile_, Feature::num_binned_,
      get_xsec_median_gap_ratio_c1);
    // infinity means all cross-sections only contain a single sample
    set_inf(ratio, -1);
    return ratio;
  }
}

std::vector<double>
PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  int c2 = 1;
  auto get_xsec_max_gap_ratio_c2 = [&c2](const Cluster::XSec &xsec) {
    return get_xsec_max_gap_ratio(xsec, c2);
  };
  auto get_xsec_min_gap_ratio_c2 = [&c2](const Cluster::XSec &xsec) {
    return get_xsec_min_gap_ratio(xsec, c2);
  };
  auto get_xsec_median_gap_ratio_c2 = [&c2](const Cluster::XSec &xsec) {
    return get_xsec_median_gap_ratio(xsec, c2);
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  if (percentile_ == 1.0)
    return get_overlapping_xsecs_binned_percentile(filtered_xsecs, percentile_,
                                                   Feature::num_binned_,
                                                   get_xsec_max_gap_ratio_c2);
  else if (percentile_ == 0.0) {
    auto ratio = get_overlapping_xsecs_binned_percentile(
      filtered_xsecs, percentile_, Feature::num_binned_,
      get_xsec_min_gap_ratio_c2);
    // infinity means all cross-sections only contain a single sample
    set_inf(ratio, -1);
    return ratio;
  } else {
    auto ratio = get_overlapping_xsecs_binned_percentile(
      filtered_xsecs, percentile_, Feature::num_binned_,
      get_xsec_median_gap_ratio_c2);
    // infinity means all cross-sections only contain a single sample
    set_inf(ratio, -1);
    return ratio;
  }
}

std::vector<double>
PercentileGapClusterwiseDistanceMaxRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature rr1(percentile_,
                                                                 1);
  PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature rr2(percentile_,
                                                                 1);

  PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature r1(
    percentile_, Feature::num_binned_);
  PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature r2(
    percentile_, Feature::num_binned_);
  if (rr1.operator()(merged_cluster, fit1, fit2)[0] >
      rr2.operator()(merged_cluster, fit1, fit2)[0])
    return r1.operator()(merged_cluster, fit1, fit2);
  else
    return r2.operator()(merged_cluster, fit1, fit2);
}
std::vector<double>
PercentileGapClusterwiseDistanceMinRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature rr1(percentile_,
                                                                 1);
  PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature rr2(percentile_,
                                                                 1);

  PercentileGapClusterwiseDistanceRatio1StrokeClusterFeature r1(
    percentile_, Feature::num_binned_);
  PercentileGapClusterwiseDistanceRatio2StrokeClusterFeature r2(
    percentile_, Feature::num_binned_);
  if (rr1.operator()(merged_cluster, fit1, fit2)[0] <
      rr2.operator()(merged_cluster, fit1, fit2)[0])
    return r1.operator()(merged_cluster, fit1, fit2);
  else
    return r2.operator()(merged_cluster, fit1, fit2);
}

std::vector<double> FirstWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  const auto ratios = get_xsecs_width_change_ratios(filtered_xsecs);
  return std::vector<double>{ratios.first};
}
std::vector<double> LastWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  const auto ratios = get_xsecs_width_change_ratios(filtered_xsecs);
  return std::vector<double>{ratios.second};
}

std::vector<double> MeanIndividualWidth1StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  int c1 = 0;
  auto get_width1 = [&c1](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c1);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_all_xsecs_binned_mean(filtered_xsecs, Feature::num_binned_,
                                   get_width1);
}

std::vector<double> PercentileIndividualWidth1StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  int c1 = 0;
  auto get_width1 = [&c1](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c1);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_all_xsecs_binned_percentile(filtered_xsecs, percentile_,
                                         Feature::num_binned_, get_width1);
}

std::vector<double> MeanIndividualWidth2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  int c2 = 1;
  auto get_width2 = [&c2](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c2);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_all_xsecs_binned_mean(filtered_xsecs, Feature::num_binned_,
                                   get_width2);
}

std::vector<double> PercentileIndividualWidth2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  int c2 = 1;
  auto get_width2 = [&c2](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c2);
    return cluster_width;
  };

  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_all_xsecs_binned_percentile(filtered_xsecs, percentile_,
                                         Feature::num_binned_, get_width2);
}

std::vector<double>
MeanMergedWidthStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                                const FittedCurve &fit1,
                                                const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_all_xsecs_binned_mean(filtered_xsecs, Feature::num_binned_,
                                   get_xsec_merged_width);
}

std::vector<double> PercentileMergedWidthStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  return get_all_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_merged_width);
}

std::vector<double> MeanWidthLengthRatio1StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  MeanIndividualWidth1StrokeClusterFeature w1_fea(Feature::num_binned_);
  auto w1 = w1_fea.operator()(merged_cluster, fit1, fit2);
  double l1 = total_length(merged_cluster, fit1, 0);
  for (auto &w : w1)
    w /= l1;
  return w1;
}

std::vector<double> PercentileWidthLengthRatio1StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileIndividualWidth1StrokeClusterFeature w1_fea(percentile_,
                                                        Feature::num_binned_);
  auto w1 = w1_fea.operator()(merged_cluster, fit1, fit2);
  double l1 = total_length(merged_cluster, fit1, 0);
  for (auto &w : w1)
    w /= l1;
  return w1;
}

std::vector<double> MeanWidthLengthRatio2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  MeanIndividualWidth2StrokeClusterFeature w2_fea(Feature::num_binned_);
  auto w2 = w2_fea.operator()(merged_cluster, fit1, fit2);
  double l2 = total_length(merged_cluster, fit2, 1);
  for (auto &w : w2)
    w /= l2;
  return w2;
}

std::vector<double> PercentileWidthLengthRatio2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileIndividualWidth2StrokeClusterFeature w2_fea(percentile_,
                                                        Feature::num_binned_);
  auto w2 = w2_fea.operator()(merged_cluster, fit1, fit2);
  double l2 = total_length(merged_cluster, fit2, 1);
  for (auto &w : w2)
    w /= l2;
  return w2;
}

std::vector<double> MeanWidthLengthMaxRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  MeanWidthLengthRatio1StrokeClusterFeature rr1(1);
  MeanWidthLengthRatio2StrokeClusterFeature rr2(1);
  double mean_r1 = rr1.operator()(merged_cluster, fit1, fit2)[0];
  double mean_r2 = rr2.operator()(merged_cluster, fit1, fit2)[0];

  MeanWidthLengthRatio1StrokeClusterFeature r1(Feature::num_binned_);
  MeanWidthLengthRatio2StrokeClusterFeature r2(Feature::num_binned_);
  if (mean_r1 > mean_r2)
    return r1.operator()(merged_cluster, fit1, fit2);
  else
    return r2.operator()(merged_cluster, fit1, fit2);
}

std::vector<double>
PercentileWidthLengthMaxRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileWidthLengthRatio1StrokeClusterFeature rr1(percentile_, 1);
  PercentileWidthLengthRatio2StrokeClusterFeature rr2(percentile_, 1);

  PercentileWidthLengthRatio1StrokeClusterFeature r1(percentile_,
                                                     Feature::num_binned_);
  PercentileWidthLengthRatio2StrokeClusterFeature r2(percentile_,
                                                     Feature::num_binned_);
  if (rr1.operator()(merged_cluster, fit1, fit2)[0] >
      rr2.operator()(merged_cluster, fit1, fit2)[0])
    return r1.operator()(merged_cluster, fit1, fit2);
  else
    return r2.operator()(merged_cluster, fit1, fit2);
}
std::vector<double> MeanWidthLengthMinRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  MeanWidthLengthRatio1StrokeClusterFeature rr1(1);
  MeanWidthLengthRatio2StrokeClusterFeature rr2(1);
  double mean_r1 = rr1.operator()(merged_cluster, fit1, fit2)[0];
  double mean_r2 = rr2.operator()(merged_cluster, fit1, fit2)[0];

  MeanWidthLengthRatio1StrokeClusterFeature r1(Feature::num_binned_);
  MeanWidthLengthRatio2StrokeClusterFeature r2(Feature::num_binned_);

  auto bin1 = r1.operator()(merged_cluster, fit1, fit2);
  auto bin2 = r2.operator()(merged_cluster, fit1, fit2);

  if (mean_r1 < mean_r2)
    return r1.operator()(merged_cluster, fit1, fit2);
  else
    return r2.operator()(merged_cluster, fit1, fit2);
}

std::vector<double>
PercentileWidthLengthMinRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileWidthLengthRatio1StrokeClusterFeature rr1(percentile_, 1);
  PercentileWidthLengthRatio2StrokeClusterFeature rr2(percentile_, 1);

  PercentileWidthLengthRatio1StrokeClusterFeature r1(percentile_,
                                                     Feature::num_binned_);
  PercentileWidthLengthRatio2StrokeClusterFeature r2(percentile_,
                                                     Feature::num_binned_);
  if (rr1.operator()(merged_cluster, fit1, fit2)[0] <
      rr2.operator()(merged_cluster, fit1, fit2)[0])
    return r1.operator()(merged_cluster, fit1, fit2);
  else
    return r2.operator()(merged_cluster, fit1, fit2);
}

std::vector<double> MeanWidthLengthRatioOverallStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  MeanMergedWidthStrokeClusterFeature w_fea(Feature::num_binned_);
  auto w = w_fea.operator()(merged_cluster, fit1, fit2);
  double l = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);
  for (auto &ww : w)
    ww /= l;
  return w;
}

std::vector<double>
PercentileWidthLengthRatioOverallStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  PercentileMergedWidthStrokeClusterFeature w_fea(percentile_,
                                                  Feature::num_binned_);
  auto w = w_fea.operator()(merged_cluster, fit1, fit2);
  double l = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);
  for (auto &ww : w)
    ww /= l;
  return w;
}

std::vector<double> MinPointwiseDistanceStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::unordered_map<size_t, std::vector<Cluster::Stroke const *>> c2strokes;
  for (const auto &s : merged_cluster.strokes) {
    c2strokes[s.cluster_ind].emplace_back(&s);
  }

  double min_dist = std::numeric_limits<double>::infinity();
  for (const auto &c_s1 : c2strokes) {
    for (const auto &c_s2 : c2strokes) {
      if (c_s1.first == c_s2.first)
        continue;

      for (const auto s1 : c_s1.second) {
        for (const auto s2 : c_s2.second) {
          double dist = stroke_sample_distance(*s1, *s2);
          min_dist = std::min(min_dist, dist);
        }
      }
    }
  }

  return std::vector<double>{min_dist};
}

std::vector<double> OverlappingRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::set<size_t> c_indices;
  std::vector<std::pair<double, double>> u_secs;
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  for (size_t i = 0; i < filtered_xsecs.size(); ++i) {
    c_indices.clear();
    for (const auto &p : filtered_xsecs[i].points) {
      c_indices.emplace(p.cluster_ind);
      if (c_indices.size() == 2)
        break;
    }

    if (c_indices.size() == 2) {
      // The individual xsec length is computed as half of the length before and
      // half of the one after
      int prev = (int)i - 1;
      int next = i + 1;
      if (prev >= 0) {
        u_secs.emplace_back(filtered_xsecs[prev].u, filtered_xsecs[i].u);
      }
      if (next < filtered_xsecs.size()) {
        u_secs.emplace_back(filtered_xsecs[i].u, filtered_xsecs[next].u);
      }
    }
  }

  // Sum length using the u values on the fitting
  double overlapping_length = 0;
  if (merged_cluster.max_u() == 0)
    return std::vector<double>{std::numeric_limits<double>::infinity()};

  for (const auto &sec : u_secs) {
    double ratio =
      std::abs(sec.second - sec.first) / std::abs(merged_cluster.max_u());
    overlapping_length += 0.5 * ratio;
  }

  return std::vector<double>{overlapping_length};
}

std::vector<double> MinWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  const auto ratios = get_xsecs_width_change_ratios(filtered_xsecs);
  return std::vector<double>{std::min(ratios.first, ratios.second)};
}
std::vector<double> MaxWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  const auto ratios = get_xsecs_width_change_ratios(filtered_xsecs);
  return std::vector<double>{std::max(ratios.first, ratios.second)};
}

std::vector<double> MinAvgWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  std::pair<double, double> non_overlapping_avg_width =
    get_xsecs_avg_width_non_overlapping(filtered_xsecs,
                                        merged_cluster.periodic);
  double mean_overlapping_width = get_overlapping_xsecs_binned_mean(
    filtered_xsecs, 1, get_overlapping_xsecs_width)[0];
  double min_width =
    std::min(non_overlapping_avg_width.first, non_overlapping_avg_width.second);
  return std::vector<double>{
    (min_width >= 0) ? mean_overlapping_width / min_width : 1};
}
std::vector<double> MaxAvgWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  std::pair<double, double> non_overlapping_avg_width =
    get_xsecs_avg_width_non_overlapping(filtered_xsecs,
                                        merged_cluster.periodic);
  double mean_overlapping_width = get_overlapping_xsecs_binned_mean(
    filtered_xsecs, 1, get_overlapping_xsecs_width)[0];
  double max_width =
    std::max(non_overlapping_avg_width.first, non_overlapping_avg_width.second);
  return std::vector<double>{
    (max_width >= 0) ? mean_overlapping_width / max_width : 1};
}

std::vector<double> InvMinWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto ratios = get_xsecs_width_change_ratios(filtered_xsecs);
  ratios.first = 1 / ratios.first;
  ratios.second = 1 / ratios.second;
  return std::vector<double>{std::min(ratios.first, ratios.second)};
}
std::vector<double> InvMaxWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto ratios = get_xsecs_width_change_ratios(filtered_xsecs);
  ratios.first = 1 / ratios.first;
  ratios.second = 1 / ratios.second;
  return std::vector<double>{std::max(ratios.first, ratios.second)};
}

std::vector<double> InvMinAvgWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  std::pair<double, double> non_overlapping_avg_width =
    get_xsecs_avg_width_non_overlapping(filtered_xsecs,
                                        merged_cluster.periodic);
  double mean_overlapping_width = get_overlapping_xsecs_binned_mean(
    filtered_xsecs, 1, get_overlapping_xsecs_width)[0];
  double min_width =
    std::min(non_overlapping_avg_width.first, non_overlapping_avg_width.second);
  return std::vector<double>{
    (min_width >= 0) ? min_width / mean_overlapping_width : 1};
}
std::vector<double> InvMaxAvgWidthChangeRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  std::pair<double, double> non_overlapping_avg_width =
    get_xsecs_avg_width_non_overlapping(filtered_xsecs,
                                        merged_cluster.periodic);
  double mean_overlapping_width = get_overlapping_xsecs_binned_mean(
    filtered_xsecs, 1, get_overlapping_xsecs_width)[0];
  double max_width =
    std::max(non_overlapping_avg_width.first, non_overlapping_avg_width.second);
  return std::vector<double>{
    (max_width >= 0) ? max_width / mean_overlapping_width : 1};
}

std::vector<double>
TypeStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                     const FittedCurve &fit1,
                                     const FittedCurve &fit2) const {
  std::unordered_map<size_t, std::vector<size_t>> c2strokes;
  for (const auto &s : merged_cluster.strokes) {
    c2strokes[s.cluster_ind].emplace_back(s.stroke_ind);
  }
  std::vector<size_t> sindices1, sindices2;
  sindices1 = c2strokes[0];
  sindices2 = c2strokes[1];
  if (sindices1.size() != 1)
    std::swap(sindices1, sindices2);

  // Type 1: Single-single
  if (sindices1.size() == 1 && sindices2.size() == 1)
    return std::vector<double>{1};
  else if (sindices1.size() == 1 && sindices2.size() > 1) {
    // Type 2: Single-multiple, single later than multiple
    if (sindices1[0] > *(std::max_element(sindices2.begin(), sindices2.end())))
      return std::vector<double>{2};
    // Type 3: Single-multiple, single not later than multiple
    return std::vector<double>{3};
  } else if (sindices1.size() > 1 && sindices2.size() > 1) {
    // Type 4: Multiple-multiple, two are clearly divided timewise
    size_t min1 = *(std::min_element(sindices1.begin(), sindices1.end()));
    size_t min2 = *(std::min_element(sindices2.begin(), sindices2.end()));
    size_t max1 = *(std::max_element(sindices1.begin(), sindices1.end()));
    size_t max2 = *(std::max_element(sindices2.begin(), sindices2.end()));
    if (min1 > max2 || min2 > max1)
      return std::vector<double>{4};
    // Type 5: Multiple-multiple, two are not clearly divided timewise
    return std::vector<double>{5};
  } else
    std::abort();
}

std::vector<double>
MeanClusterwiseDistaceMaxWidthStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  int c1 = 0;
  auto get_width1 = [&c1](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c1);
    return cluster_width;
  };
  double mean_width_1 = get_xsecs_mean(filtered_xsecs, get_width1);

  int c2 = 1;
  auto get_width2 = [&c2](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c2);
    return cluster_width;
  };
  double mean_width_2 = get_xsecs_mean(filtered_xsecs, get_width2);

  double max_mean_width = std::max(mean_width_1, mean_width_2);

  auto mean_clusterwise = get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_, get_xsec_euclidean_clusterwise);

  for (auto &w : mean_clusterwise)
    w /= max_mean_width;

  return mean_clusterwise;
}

std::vector<double>
MedianClusterwiseDistaceMaxWidthStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  int c1 = 0;
  auto get_width1 = [&c1](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c1);
    return cluster_width;
  };
  double mean_width_1 = get_xsecs_mean(filtered_xsecs, get_width1);

  int c2 = 1;
  auto get_width2 = [&c2](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c2);
    return cluster_width;
  };
  double mean_width_2 = get_xsecs_mean(filtered_xsecs, get_width2);

  double max_mean_width = std::max(mean_width_1, mean_width_2);

  auto median_clusterwise = get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_euclidean_clusterwise);

  for (auto &w : median_clusterwise)
    w /= max_mean_width;

  return median_clusterwise;
}

std::vector<double>
MeanClusterwiseDistaceMinWidthStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  int c1 = 0;
  auto get_width1 = [&c1](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c1);
    return cluster_width;
  };
  double mean_width_1 = get_xsecs_mean(filtered_xsecs, get_width1);

  int c2 = 1;
  auto get_width2 = [&c2](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c2);
    return cluster_width;
  };
  double mean_width_2 = get_xsecs_mean(filtered_xsecs, get_width2);

  double min_mean_width = std::min(mean_width_1, mean_width_2);

  auto mean_clusterwise = get_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_, get_xsec_euclidean_clusterwise);

  for (auto &w : mean_clusterwise)
    w /= min_mean_width;

  return mean_clusterwise;
}

std::vector<double>
MedianClusterwiseDistaceMinWidthStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  int c1 = 0;
  auto get_width1 = [&c1](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c1);
    return cluster_width;
  };
  double mean_width_1 = get_xsecs_mean(filtered_xsecs, get_width1);

  int c2 = 1;
  auto get_width2 = [&c2](const Cluster::XSec &xsec) -> std::vector<double> {
    std::vector<double> cluster_width = get_xsec_individual_widths(xsec, c2);
    return cluster_width;
  };
  double mean_width_2 = get_xsecs_mean(filtered_xsecs, get_width2);

  double min_mean_width = std::min(mean_width_1, mean_width_2);

  auto median_clusterwise = get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_,
    get_xsec_euclidean_clusterwise);

  for (auto &w : median_clusterwise)
    w /= min_mean_width;

  return median_clusterwise;
}

std::vector<double>
PersistenceStrokeClusterFeature::operator()(const Cluster &merged_cluster,
                                            const FittedCurve &fit1,
                                            const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  // Find the largest distance between the two parts
  std::vector<double> local_cwise_distances;
  local_cwise_distances.reserve(filtered_xsecs.size());
  for (const auto &xsec : filtered_xsecs) {
    auto dist = get_xsec_euclidean_clusterwise(xsec);
    if (!dist.empty())
      local_cwise_distances.emplace_back(dist[0]);
    else
      local_cwise_distances.emplace_back(-1);
  }
  double max_dist = -1;
  size_t max_i = 0;
  for (size_t i = 0; i < local_cwise_distances.size(); ++i) {
    if (local_cwise_distances[i] > max_dist) {
      max_dist = local_cwise_distances[i];
      max_i = i;
    }
  }

  if (max_dist <= 0)
    return std::vector<double>{0};

  std::pair<size_t, size_t> left_right{
    filtered_xsecs[max_i].points.front().cluster_ind,
    filtered_xsecs[max_i].points.back().cluster_ind};
  int start_i;
  for (start_i = max_i; start_i >= 0; --start_i) {
    int j = start_i - 1;
    if (j <= 0 || local_cwise_distances[j] <= 0)
      break;

    // Check if they are flipped
    std::pair<size_t, size_t> left_right_j{
      filtered_xsecs[j].points.front().cluster_ind,
      filtered_xsecs[j].points.back().cluster_ind};
    if (left_right_j != left_right)
      break;
  }
  assert(start_i >= 0);

  double cwise_count = 0;
  for (size_t i = start_i; i < local_cwise_distances.size(); ++i) {
    cwise_count++;
    size_t j = i + 1;
    if (j >= local_cwise_distances.size() || local_cwise_distances[j] <= 0)
      break;

    // Check if they are flipped
    std::pair<size_t, size_t> left_right_j{
      filtered_xsecs[j].points.front().cluster_ind,
      filtered_xsecs[j].points.back().cluster_ind};
    if (left_right_j != left_right)
      break;
  }

  return std::vector<double>{cwise_count / filtered_xsecs.size()};
}

std::vector<double> MaxOverlappingRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  double l1 = total_length(merged_cluster, fit1, 0);
  double l2 = total_length(merged_cluster, fit2, 1);
  double l12 = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);

  double min_length = std::min(l1, l2);

  OverlappingRatioStrokeClusterFeature fea;
  auto overlapping_ratio = fea.operator()(merged_cluster, fit1, fit2).front();
  return std::vector<double>{overlapping_ratio * l12 / min_length};
}

std::vector<double> MinOverlappingRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  double l1 = total_length(merged_cluster, fit1, 0);
  double l2 = total_length(merged_cluster, fit2, 1);
  double l12 = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);

  double max_length = std::max(l1, l2);

  OverlappingRatioStrokeClusterFeature fea;
  auto overlapping_ratio = fea.operator()(merged_cluster, fit1, fit2).front();
  return std::vector<double>{overlapping_ratio * l12 / max_length};
}

std::vector<double> OverlappingLengthStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  double l12 = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);
  OverlappingRatioStrokeClusterFeature fea;
  auto overlapping_ratio = fea.operator()(merged_cluster, fit1, fit2).front();
  return std::vector<double>{overlapping_ratio * l12};
}

std::vector<double>
PercentileGapOverlappingLengthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);

  auto get_xsec_max_gap_merged = [](const Cluster::XSec &xsec) {
    return get_xsec_max_gap(xsec, -1);
  };
  std::vector<double> max_gap = get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_max_gap_merged);

  double l12 = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);
  OverlappingRatioStrokeClusterFeature fea;
  auto overlapping_ratio = fea.operator()(merged_cluster, fit1, fit2).front();
  return std::vector<double>{max_gap.front() / (overlapping_ratio * l12)};
}

std::vector<double>
PercentileWidthOverlappingLengthRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  std::vector<double> merged_width = get_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_merged_width);

  double l12 = total_length(
    merged_cluster, (!fit1.centerline.empty()) ? merged_cluster.fit : fit1, -1);
  OverlappingRatioStrokeClusterFeature fea;
  auto overlapping_ratio = fea.operator()(merged_cluster, fit1, fit2).front();
  return std::vector<double>{merged_width.front() / (overlapping_ratio * l12)};
}

std::vector<double> OverlappingIsolineNumberStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  double overlapping_xsecs_num = 0;
  for (const auto &xsec : filtered_xsecs)
    if (is_overlapping(xsec))
      overlapping_xsecs_num++;

  return std::vector<double>{overlapping_xsecs_num};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> MaxStepawayRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto stepaway_ratios = get_xsecs_stepaway_ratios(filtered_xsecs);

  return (stepaway_ratios.first.first > stepaway_ratios.second.first)
           ? std::vector<double>{stepaway_ratios.first.second}
           : std::vector<double>{stepaway_ratios.second.second};
}

std::vector<double> MinStepawayRatioStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto stepaway_ratios = get_xsecs_stepaway_ratios(filtered_xsecs);

  return (stepaway_ratios.first.first <= stepaway_ratios.second.first)
           ? std::vector<double>{stepaway_ratios.first.second}
           : std::vector<double>{stepaway_ratios.second.second};
}

std::vector<double> MaxStepawayRatio2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto stepaway_ratios = get_xsecs_stepaway_ratios2(filtered_xsecs);

  return (stepaway_ratios.first.first > stepaway_ratios.second.first)
           ? std::vector<double>{stepaway_ratios.first.second}
           : std::vector<double>{stepaway_ratios.second.second};
}

std::vector<double> MinStepawayRatio2StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto stepaway_ratios = get_xsecs_stepaway_ratios2(filtered_xsecs);

  return (stepaway_ratios.first.first <= stepaway_ratios.second.first)
           ? std::vector<double>{stepaway_ratios.first.second}
           : std::vector<double>{stepaway_ratios.second.second};
}

std::vector<double> MaxStepawayRatio3StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  Cluster cluster1, cluster2, fit_cluster1, fit_cluster2, fit_merged_cluster;
  cluster1.strokes.emplace_back();
  cluster1.strokes.front().cluster_ind = 0;
  cluster1.strokes.front().stroke_ind = 0;
  cluster1.fit = fit1;
  cluster2.strokes.emplace_back();
  cluster2.strokes.front().cluster_ind = 1;
  cluster2.strokes.front().stroke_ind = 1;
  cluster2.fit = fit2;

  reuse_cluster(cluster1, fit_cluster1, true);
  reuse_cluster(cluster2, fit_cluster2, true);

  FeatureVector fea(0, 0);
  bool reorient = true;
  fea.compute_angle(fit_cluster1, fit_cluster2, 1.0, fit_merged_cluster,
                    nullptr, reorient, "");

  auto stepaway_ratios = get_xsecs_stepaway_ratios(fit_merged_cluster.xsecs);

  return (stepaway_ratios.first.first > stepaway_ratios.second.first)
           ? std::vector<double>{stepaway_ratios.first.second}
           : std::vector<double>{stepaway_ratios.second.second};
}

std::vector<double> MinStepawayRatio3StrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  Cluster cluster1, cluster2, fit_cluster1, fit_cluster2, fit_merged_cluster;
  cluster1.strokes.emplace_back();
  cluster1.strokes.front().cluster_ind = 0;
  cluster1.strokes.front().stroke_ind = 0;
  cluster1.fit = fit1;
  cluster2.strokes.emplace_back();
  cluster2.strokes.front().cluster_ind = 1;
  cluster2.strokes.front().stroke_ind = 1;
  cluster2.fit = fit2;

  reuse_cluster(cluster1, fit_cluster1, true);
  reuse_cluster(cluster2, fit_cluster2, true);

  FeatureVector fea(0, 0);
  bool reorient = true;
  fea.compute_angle(fit_cluster1, fit_cluster2, 1.0, fit_merged_cluster,
                    nullptr, reorient, "");

  auto stepaway_ratios = get_xsecs_stepaway_ratios(fit_merged_cluster.xsecs);

  return (stepaway_ratios.first.first <= stepaway_ratios.second.first)
           ? std::vector<double>{stepaway_ratios.first.second}
           : std::vector<double>{stepaway_ratios.second.second};
}

std::vector<double> MaxStepawayAngleStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto stepaway_angles = get_xsecs_stepaway_angles(filtered_xsecs);

  return (stepaway_angles.first.first > stepaway_angles.second.first)
           ? std::vector<double>{stepaway_angles.first.second}
           : std::vector<double>{stepaway_angles.second.second};
}

std::vector<double> MinStepawayAngleStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto stepaway_angles = get_xsecs_stepaway_angles(filtered_xsecs);

  return (stepaway_angles.first.first <= stepaway_angles.second.first)
           ? std::vector<double>{stepaway_angles.first.second}
           : std::vector<double>{stepaway_angles.second.second};
}

std::vector<double> MeanOutsideAngleStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto angle_vec = get_sep_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_, get_xsec_angle);
  for (auto &v : angle_vec)
    if (v >= std::numeric_limits<double>::infinity())
      v = 0;
  return angle_vec;
}

std::vector<double> PercentileOutsideAngleStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto angle_vec = get_sep_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_angle);
  for (auto &v : angle_vec)
    if (v >= std::numeric_limits<double>::infinity())
      v = 0;
  return angle_vec;
}

std::vector<double> MeanOutsideDistanceStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto dist_vec = get_sep_overlapping_xsecs_binned_mean(
    filtered_xsecs, Feature::num_binned_, get_xsec_euclidean);
  for (auto &v : dist_vec)
    if (v >= std::numeric_limits<double>::infinity())
      v = 0;
  return dist_vec;
}

std::vector<double> PercentileOutsideDistanceStrokeClusterFeature::operator()(
  const Cluster &merged_cluster, const FittedCurve &fit1,
  const FittedCurve &fit2) const {
  std::vector<Cluster::XSec> filtered_xsecs =
    filter_xsecs(merged_cluster.xsecs);
  auto dist_vec = get_sep_overlapping_xsecs_binned_percentile(
    filtered_xsecs, percentile_, Feature::num_binned_, get_xsec_euclidean);
  for (auto &v : dist_vec)
    if (v >= std::numeric_limits<double>::infinity())
      v = 0;
  return dist_vec;
}
