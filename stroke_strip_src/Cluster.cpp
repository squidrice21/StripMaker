#include <cmath>

#include "Cluster.h"
#include "SvgUtils.h"
#include "Utils.h"

#include "glm/detail/type_vec.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::string> colors = {
  "#00a391", "#a52a2a", "#ffd700", "#45769e", "#f77ff2", "#ff0000",
  "#db7093", "#1e90ff", "#00ee00", "#ff8c00", "#9046cc", "#6a7b53",
  "#ff1493", "#7fffd4", "#9acd32", "#ba55d3", "#ff00ff", "#d2691e",
  "#458c45", "#954595", "#f0e68c", "#00bfff", "#ffa07a",
};

glm::dvec2 simple_discrete_tangent(const std::vector<glm::dvec2> &points,
                                   size_t i) {
  glm::dvec2 result(0.0, 0.0);
  if (i > 0) {
    result += glm::normalize(points[i] - points[i - 1]);
  } else {
    result += glm::normalize(points[i + 1] - points[i]);
  }
  return glm::normalize(result);
}

glm::dvec2 discrete_tangent(const std::vector<glm::dvec2> &points, size_t i) {
  glm::dvec2 result(0.0, 0.0);
  if (i > 0) {
    result += glm::normalize(points[i] - points[i - 1]);
  }
  if (i < points.size() - 1) {
    result += glm::normalize(points[i + 1] - points[i]);
  }
  return glm::normalize(result);
}

glm::dvec2 tangent(const std::vector<glm::dvec2> &points, double i) {
  glm::dvec2 tangent(0.0, 0.0);
  size_t a = std::floor(i);
  if (double(a) == i) {
    tangent = discrete_tangent(points, (size_t)i);
  } else {
    size_t b = std::ceil(i);
    double mix = i - double(a);
    tangent = glm::normalize((1 - mix) * discrete_tangent(points, a) +
                             mix * discrete_tangent(points, b));
  }

  if (!glm::any(glm::isnan(tangent))) {
    return tangent;
  }
  return simple_discrete_tangent(points, (size_t)std::floor(i));
}

glm::dvec2 point(const std::vector<glm::dvec2> &points, double i) {
  size_t a = std::floor(i);
  size_t b = std::ceil(i);
  double mix = i - double(a);
  return (1 - mix) * points[a] + mix * points[b];
}

double average_u(const std::vector<double> &points, double i) {
  size_t a = std::floor(i);
  size_t b = std::ceil(i);
  double mix = i - double(a);
  return (1 - mix) * points[a] + mix * points[b];
}

double spiral_angle(const std::vector<double> &spiral_angles, double i) {
  size_t a = std::floor(i);
  size_t b = std::ceil(i);
  double mix = i - double(a);
  return (1 - mix) * spiral_angles[a] + mix * spiral_angles[b];
}

glm::dvec2 normal(const glm::dvec2 &v) { return glm::dvec2(-v.y, v.x); }

double total_length(const std::vector<glm::dvec2> &points) {
  double len = 0;
  for (size_t i = 1; i < points.size(); ++i) {
    len += glm::distance(points[i - 1], points[i]);
  }
  return len;
}

double total_u_length(const std::vector<double> &u) {
  double len = 0;
  for (size_t i = 1; i < u.size(); ++i) {
    len += std::abs(u[i] - u[i - 1]);
  }
  return len;
}

glm::dvec2 midpoint(const std::vector<glm::dvec2> &points) {
  return points[points.size() / 2];
}

void Input::param_svg_gradient(std::ostream &os, bool rainbow) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  glm::dvec2 gradient_start(-263.453, -1378.15);
  glm::dvec2 gradient_end(758.577, 281.46);

  SVG::begin(os, x, y, w, h);
  for (auto &kv : clusters) {
    double max_u = kv.second.max_u();

    double gradient_length = glm::distance(gradient_start, gradient_end);

    // Make a group
    os << "<g>" << std::endl;

    glm::dvec2 p_start;
    glm::dvec2 p_end;
    double u_min = std::numeric_limits<double>::infinity(), u_max = -1;
    for (auto &stroke : kv.second.strokes) {
      for (size_t i = 0; i < stroke.points.size(); ++i) {
        if (stroke.u[i] > u_max) {
          u_max = stroke.u[i];
          p_end = (stroke.points[i] - center) * thickness;
        }
        if (stroke.u[i] < u_min) {
          u_min = stroke.u[i];
          p_start = (stroke.points[i] - center) * thickness;
        }
      }
    }

    glm::dvec2 far_p = (glm::distance(gradient_start, p_start) <
                        glm::distance(gradient_start, p_end))
                         ? p_end
                         : p_start;
    bool is_left = false;
    if (far_p[0] < gradient_start[0]) {
      glm::dvec2 gradient_end2(-651.237, -616.191);
      gradient_length = glm::distance(gradient_start, gradient_end2);
      is_left = true;
    }
    double start_alpha =
      (glm::distance(gradient_start, p_start)) / gradient_length;
    double end_alpha =
      (glm::distance(gradient_start, p_start) + max_u) / gradient_length;

    bool is_start = false;
    if ((!is_left && glm::distance(gradient_start, p_start) <
                       glm::distance(gradient_start, p_end)) ||
        (is_left && p_start[1] < p_end[1])) {
      start_alpha = (glm::distance(gradient_start, p_start)) / gradient_length;
      end_alpha =
        (glm::distance(gradient_start, p_start) + max_u) / gradient_length;
      is_start = true;
    } else {
      end_alpha = (glm::distance(gradient_start, p_end)) / gradient_length;
      start_alpha =
        (glm::distance(gradient_start, p_end + max_u)) / gradient_length;
      is_start = false;
    }

    for (auto &stroke : kv.second.strokes) {
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }

      for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
        std::stringstream ss;
        if (rainbow) {
          ss << "hsl(";
          ss << int(stroke.u[i] / max_u * 360) << ", ";
          ss << "90%, 50%)";
        } else {
          double alpha = stroke.u[i] / max_u;
          alpha = (1 - alpha) * start_alpha + alpha * end_alpha;
          double alpha_alpha = stroke.u[i] / max_u;
          if (!is_start)
            alpha_alpha = 1 - alpha_alpha;

          glm::dvec3 start_color(0, 7, 55);
          glm::dvec3 end_color(13, 255, 239);
          // glm::dvec3 color = alpha * (end_color - start_color) + start_color;
          glm::dvec3 color =
            alpha_alpha * (end_color - start_color) + start_color;

          glm::dvec3 white(255, 255, 255);
          color = 0.5 * alpha_alpha * white + (1 - 0.5 * alpha_alpha) * color;
          color[0] = glm::clamp(color[0], 0.0, 255.0);
          color[1] = glm::clamp(color[1], 0.0, 255.0);
          color[2] = glm::clamp(color[2], 0.0, 255.0);
          alpha = std::min(alpha, 1.0);

          assert(alpha <= 1);
          ss << "#";
          ss << std::setfill('0') << std::setw(2);
          ss << std::hex << int(color[0]); // red
          // ss << "00"; // green
          ss << std::setfill('0') << std::setw(2);
          ss << std::hex << int(color[1]); // red
          // ss << std::hex << int(255 - alpha * 255); // blue
          ss << std::setfill('0') << std::setw(2);
          ss << std::hex << int(color[2]); // red
          // ss << std::setfill('0') << std::setw(2);
          // ss << std::hex << int(255 - (alpha_alpha * 155 + 100)); // alpha
        }
        SVG::line(os, scaled[i].x, scaled[i].y, scaled[i + 1].x,
                  scaled[i + 1].y, thickness, ss.str());
      }
    }

    os << "</g>" << std::endl;
  }
  SVG::end(os);
}

void Input::param_svg(std::ostream &os, bool rainbow) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  for (auto &kv : clusters) {
    double max_u = kv.second.max_u();

    for (auto &stroke : kv.second.strokes) {
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
        std::stringstream ss;
        if (rainbow) {
          ss << "hsl(";
          ss << int(stroke.u[i] / max_u * 360) << ", ";
          ss << "90%, 50%)";
        } else {
          ss << "#";
          ss << std::setfill('0') << std::setw(2);
          ss << std::hex << int(stroke.u[i] / max_u * 255); // red
          ss << "00"; // green
          ss << std::setfill('0') << std::setw(2);
          ss << std::hex << int(255 - stroke.u[i] / max_u * 255); // blue
        }
        SVG::line(os, scaled[i].x, scaled[i].y, scaled[i + 1].x,
                  scaled[i + 1].y, thickness, ss.str());
      }
    }
  }
  SVG::end(os);
}

void Input::orientation_svg(std::ostream &os,
                            std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  size_t stroke_num = 0;
  for (auto &kv : clusters) {
    for (auto &stroke : kv.second.strokes) {
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      std::string color(colors[stroke_num % colors.size()]);
      SVG::polyline(os, scaled, thickness, color);

      for (size_t i = 2; i < stroke.points.size(); i += 5) {
        glm::dvec2 origin = scaled[i];
        glm::dvec2 next = origin + normal(tangent(scaled, i)) * 4. * thickness;
        SVG::line(os, origin.x, origin.y, next.x, next.y, thickness, color);
      }
      ++stroke_num;
    }
  }
  cb(os);
  SVG::end(os);
}

void Input::cluster_svg(std::ostream &os, int new_stroke_ind,
                        std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  size_t cluster_num = 0;
  for (auto &kv : clusters) {
    // std::string color(colors[cluster_num % colors.size()]);
    for (auto &stroke : kv.second.strokes) {
      if (stroke.cluster_ind != std::numeric_limits<size_t>::max())
        continue;
      std::string color = "#e0e0e0";
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      // SVG::text(os, scaled[0], stroke.stroke_ind, 10);
      if (new_stroke_ind >= 0 && stroke.stroke_ind == new_stroke_ind)
        color = "#006f00";
      SVG::polyline(os, scaled, thickness, color);
    }
    for (auto &stroke : kv.second.strokes) {
      if (stroke.cluster_ind == std::numeric_limits<size_t>::max())
        continue;
      std::string color(colors[stroke.cluster_ind % colors.size()]);
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      // SVG::text(os, scaled[0], stroke.stroke_ind, 10);
      if (new_stroke_ind >= 0 && stroke.stroke_ind == new_stroke_ind)
        color = "#006f00";
      SVG::polyline(os, scaled, thickness, color);
    }
    ++cluster_num;
  }
  cb(os);
  SVG::end(os);
}

void Input::input_svg(std::ostream &os,
                      std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  size_t cluster_num = 2;
  for (auto &kv : clusters) {
    std::string color(colors[cluster_num % colors.size()]);
    for (auto &stroke : kv.second.strokes) {
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      // SVG::text(os, scaled[0], stroke.stroke_ind, 10);
      SVG::polyline(os, scaled, thickness, color);
    }
  }
  cb(os);
  SVG::end(os);
}

void Input::non_isoline_svg(std::ostream &os,
                            std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  size_t cluster_num = 0;
  for (auto &kv : clusters) {
    // std::string color(colors[cluster_num % colors.size()]);
    for (auto &stroke : kv.second.strokes) {
      if (stroke.cluster_ind != std::numeric_limits<size_t>::max())
        continue;
      std::string color = "#e0e0e0";
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      // SVG::text(os, scaled[0], stroke.stroke_ind, 10);
      SVG::polyline(os, scaled, thickness, color);
    }
    for (auto &stroke : kv.second.strokes) {
      if (stroke.cluster_ind == std::numeric_limits<size_t>::max())
        continue;
      std::string color(colors[stroke.cluster_ind % colors.size()]);
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      // SVG::text(os, scaled[0], stroke.stroke_ind, 10);
      SVG::polyline(os, scaled, thickness, color);
    }
    ++cluster_num;
  }
  cb(os);
  SVG::end(os);
}
void Input::cluster_svg_stroke_isoline(
  std::ostream &os, int new_stroke_ind, int cluster_num,
  std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  for (auto &kv : clusters) {
    if (kv.first == cluster_num) {
      std::string color(colors[kv.first % colors.size()]);
      for (auto &stroke : kv.second.strokes) {
        std::vector<glm::dvec2> scaled;
        scaled.reserve(stroke.points.size());
        for (auto &pt : stroke.points) {
          scaled.push_back((pt - center) * thickness);
        }
        // SVG::text(os, scaled[0], kv.first, 10);
        if (stroke.stroke_ind == new_stroke_ind)
          color = "#006F00";
        SVG::polyline(os, scaled, thickness, color);
      }
    }
  }
  cb(os);
  SVG::end(os);
}

void Input::only_cluster_svg(std::ostream &os, int new_stroke_ind,
                             int cluster_num,
                             std::map<int, int> stroke_to_cluster,
                             std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  auto &kv = clusters.at(cluster_num);
  std::string color;
  int first_cluster_num = -2;
  for (auto &stroke : kv.strokes) {
    std::vector<glm::dvec2> scaled;
    scaled.reserve(stroke.points.size());
    for (auto &pt : stroke.points) {
      scaled.push_back((pt - center) * thickness);
    }
    // SVG::text(os, scaled[0], stroke_to_cluster[stroke.stroke_ind], 10);
    if (first_cluster_num == -2) {
      first_cluster_num = stroke_to_cluster[stroke.stroke_ind];
    }
    if (first_cluster_num == stroke_to_cluster[stroke.stroke_ind]) {
      color = "#FF0000";
    } else {
      color = "#0000FF";
    }
    SVG::polyline(os, scaled, thickness, color);
  }

  cb(os);
  SVG::end(os);
}

void Input::only_divided_cluster_svg(
  std::ostream &os, int new_stroke_ind, int cluster_num,
  std::vector<int> base_stroke_ind, std::vector<int> target_stroke_ind,
  std::function<void(std::ostream &)> cb) const {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(*this, w, h, x, y, center);

  SVG::begin(os, x, y, w, h);
  std::string color;
  for (auto &kv : clusters) {
    for (auto &stroke : kv.second.strokes) {
      std::vector<glm::dvec2> scaled;
      scaled.reserve(stroke.points.size());
      for (auto &pt : stroke.points) {
        scaled.push_back((pt - center) * thickness);
      }
      // SVG::text(os, scaled[0], stroke_to_cluster[stroke.stroke_ind], 10);
      if (vector_finder(base_stroke_ind, stroke.stroke_ind)) {
        color = "#FF0000";
      } else {
        color = "#0000FF";
      }
      SVG::polyline(os, scaled, thickness, color);
    }
  }

  cb(os);
  SVG::end(os);
}

bool Input::vector_finder(std::vector<int> vec, int number) const {
  auto itr = std::find(vec.begin(), vec.end(), number);
  size_t index = std::distance(vec.begin(), itr);
  if (index != vec.size()) {
    return true;
  } else {
    return false;
  }
}

double Cluster::XSec::distance_weight(size_t i) const {
  if (points.size() == 1)
    return 1.0;

  double weight = 0.0;
  if (i == 0) {
    weight += glm::distance(points[0].point, points[1].point);
  } else {
    weight += glm::distance(points[i - 1].point, points[i].point);
  }

  if (i == points.size() - 1) {
    weight += glm::distance(points[points.size() - 2].point,
                            points[points.size() - 1].point);
  } else {
    weight += glm::distance(points[i].point, points[i + 1].point);
  }

  return weight;
}

glm::dvec2 Cluster::XSec::avg_point() const {
  glm::dvec2 sum(0.0, 0.0);
  for (auto &p : points) {
    sum += p.point;
  }
  if (points.size() > 0)
    sum /= double(points.size());
  return sum;
}

glm::dvec2 Cluster::XSec::avg_tangent() const {
  if (points.size() == 1)
    return points.front().tangent;

  glm::dvec2 sum(0.0, 0.0);

  std::vector<double> dists(points.size(), 0.0);
  for (size_t i = 0; i < points.size() - 1; ++i) {
    dists[i] = glm::distance(points[i].point, points[i + 1].point);
  }
  for (size_t i = 0; i < points.size(); ++i) {
    double weight = 0.0;
    if (i == 0) {
      weight += dists[0];
    } else {
      weight += dists[i - 1];
    }
    if (i == points.size() - 1) {
      weight += dists.back();
    } else {
      weight += dists[i + 1];
    }
    weight += 0.1; // Regularize to avoid zeros

    sum += weight * points[i].tangent;
  }

  sum = glm::normalize(sum);
  if (glm::any(glm::isnan(sum))) {
    sum = glm::dvec2();
  }
  return sum;
}

double Cluster::XSec::xsec_width() const {
  double max_w = 0;
  for (auto &connection : connections) {
    auto &pt_a = points[connection.a_idx];
    auto &pt_b = points[connection.b_idx];
    double sq_d = glm::dot(pt_a.point - pt_b.point, pt_a.point - pt_b.point);
    max_w = std::max(max_w, sq_d);
  }
  max_w = std::sqrt(max_w);
  return max_w;
}

double Cluster::max_u() const {
  double result = -std::numeric_limits<double>::infinity();
  for (auto &stroke : strokes) {
    for (double u : stroke.u) {
      result = std::max(result, u);
    }
  }
  return result;
}

double Cluster::min_u() const {
  double result = std::numeric_limits<double>::infinity();
  for (auto &stroke : strokes) {
    for (double u : stroke.u) {
      result = std::min(result, u);
    }
  }
  return result;
}
