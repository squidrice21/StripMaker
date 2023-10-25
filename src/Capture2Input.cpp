#include <cmath>

#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

#include <algorithm>
#include <deque>
#include <iostream>

#include <glm/glm.hpp>

#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Fitting.h"
#include "../stroke_strip_src/Parameterization.h"
#include "../stroke_strip_src/StrokeCutting.h"
#include "../stroke_strip_src/StrokeOrientation.h"
#include "Capture2Input.h"
#include "Util.h"
#include "measure/Measurement.h"

double stroke_sampling_size = 1.2; // double stroke_sampling_size = 1;
double stroke_sampling_size_ref = 1.2;
size_t stroke_sampling_min_num = 5;
size_t stroke_sampling_total_ref = 21000;

// double stroke_sampling_size = 2.75;
// size_t stroke_sampling_min_num = 3;

Input Capture2Input::from_capture(Capture capture, bool recenter) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  if (!recenter)
    center = glm::dvec2(0, 0);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    clusters[polyline.group_ind].strokes.emplace_back();
    auto &stroke = clusters[polyline.group_ind].strokes.back();
    stroke.points.reserve(polyline.points.size());
    stroke.u.reserve(polyline.points.size());
    stroke.stroke_ind = polyline.stroke_ind;
    stroke.cluster_ind = polyline.group_ind;
    for (size_t i = 0; i < polyline.points.size(); ++i) {
      auto &point = polyline.points[i];

      // Recenter and normalize to stroke width
      stroke.points.push_back(
        (glm::dvec2(point.first.x, point.first.y) - center) /
        capture.thickness);
      stroke.u.push_back(0);
    }
  }

  for (auto &c : clusters) {
    // Label spiral strokes
    if (context.to_spiral) {
      for (auto &s : c.second.strokes) {
        detect_spiral_stroke(s);
      }
    }
    c.second.original_input_strokes = c.second.strokes;
  }

  return input4capture;
}

Input Capture2Input::from_capture4set(Capture capture,
                                      std::map<int, bool> &rev_list,
                                      glm::dvec2 &ct, int stroke_num,
                                      int num_cluster) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  ct = center;
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    for (auto &potential_group : polyline.potential_groups) {
      clusters[potential_group.first].potential_strokes.emplace_back(
        std::make_pair(polyline.stroke_ind, potential_group.second));
    }
    if (polyline.points.empty())
      continue;
    if (polyline.stroke_ind <= stroke_num) {
      clusters[polyline.group_ind].strokes.emplace_back();
      auto &stroke = clusters[polyline.group_ind].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      stroke.cluster_ind = polyline.group_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];

        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  return input4capture;
}

Input Capture2Input::from_capture4cluster(Capture capture, int base_num,
                                          std::map<int, bool> &rev_list,
                                          int target_num, bool flip) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    if (polyline.group_ind == base_num || polyline.group_ind == target_num) {
      int used_group_ind = base_num;
      clusters[used_group_ind].strokes.emplace_back();
      auto &stroke = clusters[used_group_ind].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];
        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
      if (flip && polyline.group_ind == target_num) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  return input4capture;
}

void Capture2Input::add_stroke(Capture capture, Input &input, int stroke_num,
                               int group_num) {

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input.thickness = capture.thickness;

  double rate = stroke_sampling_size;
  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.stroke_ind == stroke_num) {
      input.clusters[group_num].strokes.emplace_back();
      auto &stroke = input.clusters[group_num].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      stroke.cluster_ind = polyline.group_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];
        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
    }
  }
}

Input Capture2Input::from_capture4orientation(Capture capture,
                                              std::map<int, FittedCurve> fits,
                                              glm::dvec2 ct, int stroke_num) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  std::map<int, SketchUI::Polyline2D> fits_polyline;

  for (auto &curve : fits) {
    SketchUI::Polyline2D line;
    line.from_fits(curve.second, curve.first);
    for (auto &point : line.points) {
      glm::dvec2 pt =
        (glm::dvec2(point.first.x, point.first.y) * capture.thickness + ct);
      point.first.x = pt.x;
      point.first.y = pt.y;
    }
    fits_polyline[curve.first] = line;
  }

  capture.sketchedPolylines[stroke_num].reparameterize(
    std::min(rate * capture.thickness,
             capture.sketchedPolylines[stroke_num].totalLen() / 3));

  for (auto &point : capture.sketchedPolylines[stroke_num].points) {
    min_x = std::min(min_x, point.first.x);
    max_x = std::max(max_x, point.first.x);
    min_y = std::min(min_y, point.first.y);
    max_y = std::max(max_y, point.first.y);
  }

  for (auto &polyline : fits_polyline) {
    polyline.second.reparameterize(
      std::min(rate * capture.thickness, polyline.second.totalLen() / 3));

    for (auto &point : polyline.second.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }

  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : fits_polyline) {
    if (polyline.second.points.empty())
      continue;
    {
      clusters[polyline.first].strokes.emplace_back();
      auto &stroke = clusters[polyline.first].strokes.back();
      stroke.points.reserve(polyline.second.points.size());
      stroke.u.reserve(polyline.second.points.size());
      stroke.stroke_ind = polyline.second.stroke_ind;
      stroke.cluster_ind = polyline.second.group_ind;
      for (size_t i = 0; i < polyline.second.points.size(); ++i) {
        auto &point = polyline.second.points[i];
        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
    }

    if (capture.sketchedPolylines[stroke_num].points.empty())
      continue;
    {
      clusters[polyline.first].strokes.emplace_back();
      auto &stroke = clusters[polyline.first].strokes.back();
      stroke.points.reserve(
        capture.sketchedPolylines[stroke_num].points.size());
      stroke.u.reserve(capture.sketchedPolylines[stroke_num].points.size());
      stroke.stroke_ind = capture.sketchedPolylines[stroke_num].stroke_ind;
      stroke.cluster_ind = capture.sketchedPolylines[stroke_num].group_ind;
      for (size_t i = 0;
           i < capture.sketchedPolylines[stroke_num].points.size(); ++i) {
        auto &point = capture.sketchedPolylines[stroke_num].points[i];
        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
    }
  }

  return input4capture;
}

Input Capture2Input::from_capture4input(Capture capture,
                                        std::map<int, bool> &rev_list,
                                        std::map<int, bool> crossing_cluster,
                                        std::map<int, bool> rev_clusters,
                                        int stroke_num, int num_cluster) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));
    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  std::map<int, bool> cls;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    // std::cout << "polyline cluster is " << polyline.group_ind <<std::endl;
    if (polyline.stroke_ind < stroke_num) {
      for (auto &potential_group : polyline.potential_groups) {
        clusters[potential_group.first].potential_strokes.emplace_back(
          std::make_pair(polyline.stroke_ind, potential_group.second));
      }
      cls[polyline.group_ind] = true;
      clusters[polyline.group_ind].strokes.emplace_back();
      auto &stroke = clusters[polyline.group_ind].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];

        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.stroke_ind == stroke_num) {
      for (auto &c : cls) {
        // std::cout << "num: " << num << std::endl;
        clusters[c.first].strokes.emplace_back();
        auto &stroke = clusters[c.first].strokes.back();
        stroke.points.reserve(polyline.points.size());
        stroke.u.reserve(polyline.points.size());
        stroke.stroke_ind = polyline.stroke_ind;
        for (size_t i = 0; i < polyline.points.size(); ++i) {
          auto &point = polyline.points[i];

          // Recenter and normalize to stroke width
          stroke.points.push_back(
            (glm::dvec2(point.first.x, point.first.y) - center) /
            capture.thickness);
          stroke.u.push_back(0);
        }
        if (rev_clusters[c.first]) {
          std::reverse(stroke.points.begin(), stroke.points.end());
        }
      }
    }
  }

  return input4capture;
}

Input Capture2Input::from_capture4two_cluster(Capture capture,
                                              std::map<int, bool> &rev_list,
                                              int base_num, int base_stroke,
                                              int upper_ind) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    if (polyline.group_ind == base_num && polyline.stroke_ind <= upper_ind) {
      int used_group_ind = polyline.group_ind;
      if (polyline.stroke_ind > base_stroke)
        used_group_ind = -1;
      clusters[used_group_ind].strokes.emplace_back();
      auto &stroke = clusters[used_group_ind].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];
        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  return input4capture;
}

Input Capture2Input::from_capture4set_merge(Capture capture,
                                            std::map<int, bool> &rev_list,
                                            int base_num, int base_stroke,
                                            int upper_ind) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    if (polyline.group_ind == base_num && upper_ind) {
      clusters[polyline.group_ind].strokes.emplace_back();
      auto &stroke = clusters[polyline.group_ind].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];

        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  return input4capture;
}

Input Capture2Input::from_capture4set_merge_negative(
  Capture capture, std::map<int, bool> &rev_list, int base_num,
  int target_num) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    if (polyline.group_ind == base_num || polyline.group_ind == target_num) {
      clusters[base_num].strokes.emplace_back();
      auto &stroke = clusters[base_num].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];

        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  return input4capture;
}

Input Capture2Input::from_capture4two_cluster_negative(
  Capture capture, std::map<int, bool> &rev_list, int base_num,
  int target_num) {
  Input input4capture;
  auto &clusters = input4capture.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input4capture.thickness = capture.thickness;

  double rate = stroke_sampling_size;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(
      rate * capture.thickness, polyline.totalLen() / stroke_sampling_min_num));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input4capture.orig_center = center;
  input4capture.width = (max_x - min_x) / capture.thickness;
  input4capture.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    if (polyline.group_ind == base_num || polyline.group_ind == target_num) {
      // int group_num = base_num;
      // if (polyline.group_ind == target_num)
      // 	group_num = -1;
      clusters[polyline.group_ind].strokes.emplace_back();
      auto &stroke = clusters[polyline.group_ind].strokes.back();
      stroke.points.reserve(polyline.points.size());
      stroke.u.reserve(polyline.points.size());
      stroke.stroke_ind = polyline.stroke_ind;
      for (size_t i = 0; i < polyline.points.size(); ++i) {
        auto &point = polyline.points[i];

        // Recenter and normalize to stroke width
        stroke.points.push_back(
          (glm::dvec2(point.first.x, point.first.y) - center) /
          capture.thickness);
        stroke.u.push_back(0);
      }
      if (rev_list[polyline.stroke_ind]) {
        std::reverse(stroke.points.begin(), stroke.points.end());
      }
    }
  }

  return input4capture;
}
