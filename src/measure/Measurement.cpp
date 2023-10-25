#include "Measurement.h"
#include "../../stroke_strip_src/Fitting.h"
#include "../Logger.h"
#include "../Util.h"
#include "../uncutting/Uncutting.h"
#include "glm/detail/type_vec.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

double spiral_jump_threshold = M_PI / 6;

// Note this function contains parameters
void detect_spiral_stroke(Cluster::Stroke &s) {
  double curv = stroke_total_signed_curvature(s);
  curv = std::abs(curv);
  if (curv < 2 * M_PI)
    return;

  Cluster cluster;
  Cluster::Stroke original_s = s;
  set_spiral_cut_angle(s, &cluster);

  double max_gap = -std::numeric_limits<double>::infinity();
  for (const auto &xsec : cluster.xsecs) {
    for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
      double gap =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      gap = std::max(gap, 0.0);
      max_gap = std::max(gap, max_gap);
    }
  }

  if (max_gap >= 50 || max_gap < 0) {
    s = original_s;
    return;
  }

  double spiral_length = total_length(cluster.fit.centerline);
  double radius = spiral_length / (2 * M_PI);
  std::vector<glm::dvec2> centers = spiral_mass_centers(s);
  assert(centers.size() == 2);
  double center_diff = glm::distance(centers[0], centers[1]);
  double diff_ratio = center_diff / radius;

  if (diff_ratio < 0.25 ||
      ((curv > 2 * M_PI * 3 || max_gap <= 3) && diff_ratio < 0.45)) {
    s.spiral = true;
    return;
  }

  s = cluster.original_input_strokes.front();
}

double stroke_total_signed_curvature(const Cluster::Stroke &s) {
  double curv = 0;
  for (size_t i = 1; i + 1 < s.points.size(); ++i) {
    glm::dvec2 prev = glm::normalize(s.points[i] - s.points[i - 1]);
    glm::dvec2 cur = glm::normalize(s.points[i + 1] - s.points[i]);
    glm::dvec2 prev_norm(-prev.y, prev.x);

    double sign = (glm::dot(prev_norm, cur) > 0) ? 1 : -1;
    double dot = glm::dot(prev, cur);
    double angle = sign * std::acos(std::min(std::max(-1., dot), 1.));
    curv += angle;
  }

  return curv;
}

double max_spiral_gap(const Cluster::Stroke &s) {
  // Fit assuming stroke is a spiral
  Cluster::Stroke ss = s;
  ss.spiral = true;
  Input input;
  set_spiral_cut_angle(ss, &(input.clusters[0]));

  double max_gap = -std::numeric_limits<double>::infinity();
  for (const auto &xsec : input.clusters[0].xsecs) {
    for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
      double gap =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      gap = std::max(gap, 0.0);
      max_gap = std::max(gap, max_gap);
    }
  }

  return max_gap;
}

double min_spiral_gap(const Cluster::Stroke &s) {
  // Fit assuming stroke is a spiral
  Cluster::Stroke ss = s;
  ss.spiral = true;
  Input input;
  set_spiral_cut_angle(ss, &(input.clusters[0]));

  double min_gap = std::numeric_limits<double>::infinity();
  for (const auto &xsec : input.clusters[0].xsecs) {
    for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
      double gap =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      gap = std::max(gap, 0.0);
      min_gap = std::min(gap, min_gap);
    }
  }

  return min_gap;
}

double spiral_closed_length(const Cluster::Stroke &s) {
  // Fit assuming stroke is a spiral
  Cluster::Stroke ss = s;
  ss.spiral = true;
  Input input;
  set_spiral_cut_angle(ss, &(input.clusters[0]));

  double max_gap = -std::numeric_limits<double>::infinity();
  for (const auto &xsec : input.clusters[0].xsecs) {
    for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
      double gap =
        glm::distance(xsec.points[i].point, xsec.points[i + 1].point) - 1.0;
      gap = std::max(gap, 0.0);
      max_gap = std::max(gap, max_gap);
    }
  }

  if (max_gap < 0)
    return -1;

  return total_length(input.clusters[0].fit.centerline);
}

std::vector<glm::dvec2> spiral_mass_centers(const Cluster::Stroke &s) {
  std::vector<glm::dvec2> mass_centers;

  auto get_mass_center =
    [](const std::vector<glm::dvec2> &points) -> glm::dvec2 {
    double curv = 0;
    glm::dvec2 center(0, 0);
    if (points.empty())
      return center;
    center = points.front();

    size_t i = 1;
    for (; i + 1 < points.size(); ++i) {
      glm::dvec2 prev = glm::normalize(points[i] - points[i - 1]);
      glm::dvec2 cur = glm::normalize(points[i + 1] - points[i]);
      glm::dvec2 prev_norm(-prev.y, prev.x);

      double sign = (glm::dot(prev_norm, cur) > 0) ? 1 : -1;
      double dot = glm::dot(prev, cur);
      double angle = sign * std::acos(std::min(std::max(-1., dot), 1.));
      curv += angle;

      center += points[i];

      if (std::abs(curv) > 2 * M_PI) {
        center.x /= (i + 1);
        center.y /= (i + 1);
        return center;
      }
    }

    center.x /= (i + 1);
    center.y /= (i + 1);
    return center;
  };

  auto center1 = get_mass_center(s.points);
  std::vector<glm::dvec2> points = s.points;
  std::reverse(points.begin(), points.end());
  auto center2 = get_mass_center(points);

  mass_centers.emplace_back(center1);
  mass_centers.emplace_back(center2);

  return mass_centers;
}

void set_spiral_cut_angle(Cluster::Stroke &s, Cluster *spiral_cluster_ptr) {
  if (std::abs(stroke_total_signed_curvature(s)) < 2 * M_PI)
    return;

  std::vector<double> cut_angles{2 * M_PI, M_PI, M_PI / 2};
  bool pass = true;

  for (auto cut_angle : cut_angles) {
    if (cut_angle > s.spiral_cut_angle)
      continue;

    pass = true;
    Cluster::Stroke ss = s;
    ss.spiral = true;
    ss.spiral_cut_angle = cut_angle;

    Input input;
    input.clusters[0].original_input_strokes.emplace_back(s);
    input.clusters[0].strokes.emplace_back(ss);
    fea_param.parameterize(&input);

    const auto &fit_map = fea_fitting.fit(&input);
    const auto &fit_curve = fit_map.at(0);
    input.clusters[0].fit = fit_curve;

    for (size_t i = 1; i + 1 < input.clusters[0].fit.centerline.size(); ++i) {
      glm::dvec2 prev = glm::normalize(input.clusters[0].fit.centerline[i] -
                                       input.clusters[0].fit.centerline[i - 1]);
      glm::dvec2 cur = glm::normalize(input.clusters[0].fit.centerline[i + 1] -
                                      input.clusters[0].fit.centerline[i]);
      glm::dvec2 prev_norm(-prev.y, prev.x);

      double dot = glm::dot(prev, cur);
      double angle = std::acos(std::min(std::max(-1., dot), 1.));

      if (angle > spiral_jump_threshold) {
        pass = false;
        break;
      }
    }

    if (pass) {
      if (spiral_cluster_ptr)
        *spiral_cluster_ptr = std::move(input.clusters[0]);
      s.spiral_cut_angle = ss.spiral_cut_angle;
      break;
    }
  }
}
