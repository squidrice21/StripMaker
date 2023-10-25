#pragma once

#include <glm/glm.hpp>

#include <functional>
#include <map>
#include <ostream>
#include <sstream>
#include <vector>

struct Cluster {
  struct Stroke {
    size_t stroke_ind;
    size_t cluster_ind;
    std::vector<glm::dvec2> points;
    std::vector<double> u;
    bool spiral = false;
    std::vector<double> spiral_angles;
    double spiral_cut_angle = 2 * M_PI;
  };

  struct XSecPoint {
    size_t stroke_ind;
    size_t cluster_ind;
    size_t stroke_idx;
    double i;
    glm::dvec2 point;
    glm::dvec2 to_next;
    glm::dvec2 tangent;

    // Indexing within a stroke
    size_t stroke_within_idx;
    size_t stroke_xsec_count;
  };

  struct XSecConnection {
    size_t a_idx;
    size_t b_idx;
    double weight;
  };

  struct XSec {
    std::vector<XSecPoint> points;
    std::vector<XSecConnection> connections;
    size_t center_idx;
    double u;
    bool connector;

    double distance_weight(size_t i) const;
    glm::dvec2 avg_tangent() const;
    glm::dvec2 avg_point() const;
    double xsec_width() const;
  };
  struct FittedCurve {
    std::vector<glm::dvec2> centerline;
    std::vector<double> widths;
    size_t cluster_idx;
    std::vector<XSec> fit_sample_xsecs;
  };

  std::vector<Stroke> strokes;
  std::vector<Stroke> original_input_strokes;
  std::vector<XSec> xsecs;
  std::vector<std::vector<XSec>> orthogonal_xsecs_list;
  std::vector<std::pair<int, bool>> potential_strokes;
  bool periodic = false;
  FittedCurve fit;
  std::map<std::string, double> obj_term_values;
  bool param_failed = false;

  double min_u() const;
  double max_u() const;
};

glm::dvec2 point(const std::vector<glm::dvec2> &points, double i);
glm::dvec2 tangent(const std::vector<glm::dvec2> &points, double i);
double average_u(const std::vector<double> &points, double i);
glm::dvec2 discrete_tangent(const std::vector<glm::dvec2> &points, size_t i);
glm::dvec2 normal(const glm::dvec2 &v);
glm::dvec2 midpoint(const std::vector<glm::dvec2> &points);
double spiral_angle(const std::vector<double> &spiral_angles, double i);
double total_length(const std::vector<glm::dvec2> &points);
double total_u_length(const std::vector<double> &u);

struct Input {
  double thickness;
  double width;
  double height;
  std::map<int, Cluster> clusters;

  double orig_thickness;
  glm::dvec2 orig_center;

  std::string repr() const {
    size_t s_count = 0;
    for (auto const &c : clusters)
      s_count += c.second.strokes.size();

    auto ss = std::stringstream();
    ss.precision(6);
    ss << "Input(" << width << "x" << height;
    ss << "; len: " << s_count;
    ss << ")";
    return ss.str();
  }

  void param_svg_gradient(std::ostream &os, bool rainbow = false) const;
  void param_svg(std::ostream &os, bool rainbow = false) const;
  void orientation_svg(
    std::ostream &os,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;
  void cluster_svg(
    std::ostream &os, int new_stroke_ind,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;
  void input_svg(
    std::ostream &os,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;
  void non_isoline_svg(
    std::ostream &os,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;

  bool vector_finder(std::vector<int> vec, int number) const;
  void cluster_svg_stroke_isoline(
    std::ostream &os, int new_stroke_ind, int cluster_num,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;
  void only_cluster_svg(
    std::ostream &os, int new_stroke_ind, int cluster_num,
    std::map<int, int> stroke_to_cluster,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;
  void only_divided_cluster_svg(
    std::ostream &os, int new_stroke_ind, int cluster_num,
    std::vector<int> base_stroke_ind, std::vector<int> target_stroke_ind,
    std::function<void(std::ostream &)> cb = [](std::ostream &) {}) const;
};

using FittedCurve = Cluster::FittedCurve;
