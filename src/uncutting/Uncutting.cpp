#include "Uncutting.h"

#include "../Logger.h"
#include "../Util.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

double match_dist_threshold = 5;
double endpoint_match_dist_threshold = 1e-4;
// double endpoint_match_dist_threshold = 0.1;
// double angle_threshold = 120;
double angle_threshold = 25;

glm::dvec2 stepaway_tangent(const Cluster::Stroke &stroke, const bool is_end,
                            const double stepaway_amount) {
  double furthest_sq_dist = stepaway_amount;
  if (!is_end) {
    for (size_t i = 1; i < stroke.points.size(); ++i) {
      const double sq_dist =
        glm::length(stroke.points[i] - stroke.points[i - 1]);
      if (sq_dist > furthest_sq_dist) {
        glm::dvec2 target_point =
          stroke.points[i - 1] +
          glm::normalize(stroke.points[i] - stroke.points[i - 1]) *
            furthest_sq_dist;
        return glm::normalize(stroke.points[0] - target_point);
      } else {
        furthest_sq_dist -= sq_dist;
      }
    }
    return glm::normalize(stroke.points[0] -
                          stroke.points[stroke.points.size() - 1]);
  } else {
    for (size_t i = stroke.points.size() - 2; i >= 0; --i) {
      const double sq_dist =
        glm::length(stroke.points[i] - stroke.points[i + 1]);
      if (sq_dist > furthest_sq_dist) {
        glm::dvec2 target_point =
          stroke.points[i + 1] +
          glm::normalize(stroke.points[i] - stroke.points[i + 1]) *
            furthest_sq_dist;
        return glm::normalize(stroke.points[stroke.points.size() - 1] -
                              target_point);
      } else {
        furthest_sq_dist -= sq_dist;
      }
    }
    return glm::normalize(stroke.points[stroke.points.size() - 1] -
                          stroke.points[0]);
  }
}

double hausdorff_distance_one_sided(const Sketch &s1, const Sketch &s2) {
  double hausdorff_dist = -1;
  for (auto const &p : s1.points) {
    double c1_dist = stroke_sample_distance(to_glm_dvec2(p.first), s2);
    hausdorff_dist = std::max(hausdorff_dist, c1_dist);
  }

  return hausdorff_dist;
}

double hausdorff_distance(const Sketch &s1, const Sketch &s2) {
  double hausdorff_dist = -1;
  for (auto const &p : s1.points) {
    double c1_dist = stroke_sample_distance(to_glm_dvec2(p.first), s2);
    hausdorff_dist = std::max(hausdorff_dist, c1_dist);
  }
  for (auto const &p : s2.points) {
    double c2_dist = stroke_sample_distance(to_glm_dvec2(p.first), s1);
    hausdorff_dist = std::max(hausdorff_dist, c2_dist);
  }

  return hausdorff_dist;
}

double hausdorff_distance(const FittedCurve &s1, const FittedCurve &s2) {
  double hausdorff_dist = -1;
  for (auto const &p : s1.centerline) {
    double c1_dist = stroke_sample_distance(p, s2);
    hausdorff_dist = std::max(hausdorff_dist, c1_dist);
  }
  for (auto const &p : s2.centerline) {
    double c2_dist = stroke_sample_distance(p, s1);
    hausdorff_dist = std::max(hausdorff_dist, c2_dist);
  }

  return hausdorff_dist;
}

double hausdorff_distance_one_sided(const FittedCurve &s1,
                                    const FittedCurve &s2) {
  double hausdorff_dist = -1;
  for (auto const &p : s1.centerline) {
    double c1_dist = stroke_sample_distance(p, s2);
    hausdorff_dist = std::max(hausdorff_dist, c1_dist);
  }

  return hausdorff_dist;
}

void map_to_ref(const Capture &capture, const Capture &capture_ref,
                std::map<size_t, size_t> &in2ref) {
  for (const auto &s1 : capture.sketchedPolylines) {
    double min_dist = std::numeric_limits<double>::infinity();
    int s2_idx = -1;
    for (const auto &s2 : capture_ref.sketchedPolylines) {
      double dist12 = hausdorff_distance_one_sided(s1, s2);
      if (dist12 < min_dist) {
        min_dist = dist12;
        assert(s2.stroke_ind >= 0);
        s2_idx = s2.stroke_ind;
      }
    }

    if (s2_idx >= 0) {
      in2ref[s1.stroke_ind] = s2_idx;
      if (min_dist > match_dist_threshold) {
        logger().warn("Matched {} to {} with distance {}", s1.stroke_ind,
                      s2_idx, min_dist);
      }
    }
  }
}

glm::dvec2 capture_center(const Capture &cluster, double &width,
                          double &height) {
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &s : cluster.sketchedPolylines) {
    for (const auto &point : s.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }

  width = max_x - min_x;
  height = max_y - min_y;

  return glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
}

void matching_point(
  Input &input,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair) {
  for (auto &cluster : input.clusters) {
    for (auto &stroke : cluster.second.strokes) {
      for (auto &cluster_sub : input.clusters) {
        for (auto &stroke_sub : cluster_sub.second.strokes) {
          if (stroke.stroke_ind >= stroke_sub.stroke_ind)
            continue;
          if (glm::distance(stroke.points[0], stroke_sub.points[0]) <
              endpoint_match_dist_threshold)
            matching_pair.emplace_back(
              std::make_pair(std::make_pair(stroke, false),
                             std::make_pair(stroke_sub, false)));
          else if (glm::distance(
                     stroke.points[0],
                     stroke_sub.points[stroke_sub.points.size() - 1]) <
                   endpoint_match_dist_threshold)
            matching_pair.emplace_back(std::make_pair(
              std::make_pair(stroke, false), std::make_pair(stroke_sub, true)));
          else if (glm::distance(stroke.points[stroke.points.size() - 1],
                                 stroke_sub.points[0]) <
                   endpoint_match_dist_threshold)
            matching_pair.emplace_back(std::make_pair(
              std::make_pair(stroke, true), std::make_pair(stroke_sub, false)));
          else if (glm::distance(
                     stroke.points[stroke.points.size() - 1],
                     stroke_sub.points[stroke_sub.points.size() - 1]) <
                   endpoint_match_dist_threshold)
            matching_pair.emplace_back(std::make_pair(
              std::make_pair(stroke, true), std::make_pair(stroke_sub, true)));
        }
      }
    }
  }
}

void check_GT_cluster(
  Input &input,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair) {
  std::vector<int> delete_list;
  for (size_t i = 0; i < matching_pair.size(); i++) {
    Cluster::Stroke &stroke = matching_pair[i].first.first;
    Cluster::Stroke &stroke_sub = matching_pair[i].second.first;
    if (stroke.cluster_ind != stroke_sub.cluster_ind)
      delete_list.emplace_back(i);
  }
  std::sort(delete_list.begin(), delete_list.end(),
            [](int a, int b) { return a > b; });
  for (size_t i = 0; i < delete_list.size(); i++)
    matching_pair.erase(matching_pair.begin() + delete_list[i]);
}

void check_tangent(
  Input &input,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair) {
  std::vector<int> delete_list;
  for (size_t i = 0; i < matching_pair.size(); i++) {
    Cluster::Stroke &stroke = matching_pair[i].first.first;
    Cluster::Stroke &stroke_sub = matching_pair[i].second.first;
    glm::dvec2 vec =
      stepaway_tangent(stroke, matching_pair[i].first.second, 2.0);
    glm::dvec2 vec_sub =
      stepaway_tangent(stroke_sub, matching_pair[i].second.second, 2.0);
    // if (!matching_pair[i].first.second)
    //   vec = glm::normalize(stroke.points[0] - stroke.points[1]);
    // else
    //   vec = glm::normalize(stroke.points[stroke.points.size() - 1] -
    //                        stroke.points[stroke.points.size() - 2]);
    // if (!matching_pair[i].second.second)
    //   vec_sub = glm::normalize(stroke_sub.points[0] - stroke_sub.points[1]);
    // else
    //   vec_sub = glm::normalize(stroke_sub.points[stroke_sub.points.size() -
    //   1] -
    //                            stroke_sub.points[stroke_sub.points.size() -
    //                            2]);
    std::cout << vec.x << ", " << vec.y << std::endl;
    std::cout << vec_sub.x << ", " << vec_sub.y << std::endl;
    double dot = glm::dot(vec, -vec_sub);
    std::cout << "dot: " << dot << std::endl;
    double angle = std::acos(glm::clamp(dot, -1.0, 1.0));
    std::cout << "angle: " << angle << std::endl;
    if (angle > M_PI / 180 * angle_threshold)
      delete_list.emplace_back(i);
  }
  std::sort(delete_list.begin(), delete_list.end(),
            [](int a, int b) { return a > b; });

  for (size_t i = 0; i < delete_list.size(); i++)
    matching_pair.erase(matching_pair.begin() + delete_list[i]);
}

void check_original_input(
  const Capture &capture, const Capture &capture_orig,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair) {
  std::map<size_t, size_t> in2ref;
  map_to_ref(capture, capture_orig, in2ref);

  std::vector<int> delete_list;
  for (size_t i = 0; i < matching_pair.size(); i++) {
    Cluster::Stroke &stroke = matching_pair[i].first.first;
    Cluster::Stroke &stroke_sub = matching_pair[i].second.first;
    assert(in2ref.count(stroke.stroke_ind) &&
           in2ref.count(stroke_sub.stroke_ind));
    // Mapping to the same stroke in the original input means that the stroke
    // pair here is caused by cutting in preprocessing
    if (in2ref[stroke.stroke_ind] != in2ref[stroke_sub.stroke_ind])
      delete_list.emplace_back(i);
  }
  std::sort(delete_list.begin(), delete_list.end(),
            [](int a, int b) { return a > b; });
  for (size_t i = 0; i < delete_list.size(); i++)
    matching_pair.erase(matching_pair.begin() + delete_list[i]);
}

void connect_stroke(
  Capture &capture,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair,
  const std::string &csv_filename) {
  std::vector<int> delete_stroke_list;

  if (!csv_filename.empty()) {
    std::ofstream file(csv_filename);
    file << "stroke_ind_a,is_end_a,stroke_ind_b,is_end_b\n";
    for (auto &pair : matching_pair) {
      int base_stroke_ind = pair.first.first.stroke_ind;
      int target_stroke_ind = pair.second.first.stroke_ind;
      file << base_stroke_ind << "," << pair.first.second << ",";
      file << target_stroke_ind << "," << pair.second.second << "\n";
    }
  }

  for (size_t i = 0; i < matching_pair.size(); i++) {
    auto &pair = matching_pair[i];
    int base_stroke_ind = pair.first.first.stroke_ind;
    int target_stroke_ind = pair.second.first.stroke_ind;
    std::cout << base_stroke_ind << pair.first.second << ","
              << target_stroke_ind << pair.second.second << std::endl;
  }

  std::map<int, std::pair<int, int>> belonged_strokes;
  for (size_t i = 0; i < matching_pair.size(); i++) {
    auto &pair = matching_pair[i];
    int base_stroke_ind = pair.first.first.stroke_ind;
    int target_stroke_ind = pair.second.first.stroke_ind;
    Cluster::Stroke base = pair.first.first;
    Cluster::Stroke target = pair.first.first;
    if (pair.first.second) {
      if (pair.second.second) {
        std::reverse(capture.sketchedPolylines[base_stroke_ind].points.begin(),
                     capture.sketchedPolylines[base_stroke_ind].points.end());
        capture.sketchedPolylines[target_stroke_ind].points.insert(
          capture.sketchedPolylines[target_stroke_ind].points.end(),
          capture.sketchedPolylines[base_stroke_ind].points.begin(),
          capture.sketchedPolylines[base_stroke_ind].points.end());
        // have to check
        capture.sketchedPolylines[target_stroke_ind].reparameterize(0.1);
        // capture
        //   .sketchedPolylines[target_stroke_ind]
        //   .group_ind = -1;
        delete_stroke_list.emplace_back(base_stroke_ind);
        for (size_t j = i + 1; j < matching_pair.size(); j++) {
          auto &pair_sub = matching_pair[j];
          if (pair_sub.first.first.stroke_ind == base_stroke_ind) {
            pair_sub.first.first = target;
            if (!pair_sub.first.second)
              pair_sub.first.second = true;
            else
              std::cout << "two options_0" << std::endl;
          }
          if (pair_sub.second.first.stroke_ind == base_stroke_ind) {
            pair_sub.second.first = target;
            if (!pair_sub.second.second)
              pair_sub.second.second = true;
            else
              std::cout << "two options_1" << std::endl;
          }
        }
      } else {
        capture.sketchedPolylines[base_stroke_ind].points.insert(
          capture.sketchedPolylines[base_stroke_ind].points.end(),
          capture.sketchedPolylines[target_stroke_ind].points.begin(),
          capture.sketchedPolylines[target_stroke_ind].points.end());
        capture.sketchedPolylines[base_stroke_ind].reparameterize(0.1);
        // capture.sketchedPolylines[base_stroke_ind]
        //   .group_ind = -1;
        delete_stroke_list.emplace_back(target_stroke_ind);
        for (size_t j = i + 1; j < matching_pair.size(); j++) {
          auto &pair_sub = matching_pair[j];
          if (pair_sub.first.first.stroke_ind == target_stroke_ind) {
            pair_sub.first.first = base;
          }
          if (pair_sub.second.first.stroke_ind == target_stroke_ind) {
            pair_sub.second.first = base;
          }
        }
      }
    } else {
      if (pair.second.second) {
        std::cout << "two options_5" << std::endl;
        capture.sketchedPolylines[target_stroke_ind].points.insert(
          capture.sketchedPolylines[target_stroke_ind].points.end(),
          capture.sketchedPolylines[base_stroke_ind].points.begin(),
          capture.sketchedPolylines[base_stroke_ind].points.end());
        std::cout << "two options_5" << std::endl;
        capture.sketchedPolylines[target_stroke_ind].reparameterize(0.1);
        // capture
        //   .sketchedPolylines[target_stroke_ind]
        //   .group_ind = -1;
        delete_stroke_list.emplace_back(base_stroke_ind);
        std::cout << "two options_5" << std::endl;
        for (size_t j = i + 1; j < matching_pair.size(); j++) {
          auto &pair_sub = matching_pair[j];
          if (pair_sub.first.first.stroke_ind == base_stroke_ind) {
            pair_sub.first.first = target;
          }
          if (pair_sub.second.first.stroke_ind == base_stroke_ind) {
            pair_sub.second.first = target;
          }
        }
      } else {
        std::reverse(capture.sketchedPolylines[base_stroke_ind].points.begin(),
                     capture.sketchedPolylines[base_stroke_ind].points.end());
        capture.sketchedPolylines[base_stroke_ind].points.insert(
          capture.sketchedPolylines[base_stroke_ind].points.end(),
          capture.sketchedPolylines[target_stroke_ind].points.begin(),
          capture.sketchedPolylines[target_stroke_ind].points.end());
        capture.sketchedPolylines[base_stroke_ind].reparameterize(0.1);
        // capture.sketchedPolylines[base_stroke_ind]
        //   .group_ind = -1;
        delete_stroke_list.emplace_back(target_stroke_ind);
        for (size_t j = i + 1; j < matching_pair.size(); j++) {
          auto &pair_sub = matching_pair[j];
          if (pair_sub.first.first.stroke_ind == target_stroke_ind) {
            pair_sub.first.first = base;
          }
          if (pair_sub.second.first.stroke_ind == target_stroke_ind) {
            pair_sub.second.first = base;
          }
          if (pair_sub.first.first.stroke_ind == base_stroke_ind) {
            if (pair_sub.first.second)
              pair_sub.first.second = false;
            else
              std::cout << "two options_2" << std::endl;
          }
          if (pair_sub.second.first.stroke_ind == base_stroke_ind) {
            if (pair_sub.second.second)
              pair_sub.second.second = false;
            else
              std::cout << "two options_3" << std::endl;
          }
        }
      }
    }
  }
  std::vector<int> delete_stroke_index_list;
  for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
    if (std::find(delete_stroke_list.begin(), delete_stroke_list.end(),
                  capture.sketchedPolylines[i].stroke_ind) !=
        delete_stroke_list.end())
      delete_stroke_index_list.emplace_back(i);
  }

  std::sort(delete_stroke_index_list.begin(), delete_stroke_index_list.end(),
            [](int a, int b) { return a > b; });

  for (auto &delete_index : delete_stroke_index_list) {
    std::cout << delete_index << std::endl;
    capture.sketchedPolylines.erase(capture.sketchedPolylines.begin() +
                                    delete_index);
  }
}

void connect_stroke(
  Capture &capture,
  const std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                              std::pair<Cluster::Stroke, bool>>> &matching_pair,
  const std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                              std::pair<Cluster::Stroke, bool>>> &cut_pair,
  const std::string &vis_csv_filename, const std::string &cut_csv_filename) {
  if (!vis_csv_filename.empty()) {
    std::ofstream file(vis_csv_filename);
    file << "stroke_ind_a,is_end_a,stroke_ind_b,is_end_b\n";
    for (const auto &pair : matching_pair) {
      int base_stroke_ind = pair.first.first.stroke_ind;
      int target_stroke_ind = pair.second.first.stroke_ind;
      file << base_stroke_ind << "," << pair.first.second << ",";
      file << target_stroke_ind << "," << pair.second.second << "\n";
    }
  }

  // Merge strokes
  std::map<size_t, SketchUI::Polyline2D *> s_map;
  for (auto &s : capture.sketchedPolylines) {
    s_map[s.stroke_ind] = &s;
  }

  std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    matching_pair_indices, cut_pair_indices;
  for (const auto &pair : matching_pair) {
    auto s_head =
      std::make_pair(pair.first.first.stroke_ind, pair.first.second);
    auto s_tail =
      std::make_pair(pair.second.first.stroke_ind, pair.second.second);
    auto head_find = std::find_if(
      matching_pair_indices.begin(), matching_pair_indices.end(),
      [&s_head](
        const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> &p) {
        return p.first == s_head || p.second == s_head;
      });
    auto tail_find = std::find_if(
      matching_pair_indices.begin(), matching_pair_indices.end(),
      [&s_tail](
        const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> &p) {
        return p.first == s_tail || p.second == s_tail;
      });
    // Check if we have cases where one endpoint is matched to multiple strokes
    if (head_find == matching_pair_indices.end() &&
        tail_find == matching_pair_indices.end())
      matching_pair_indices.emplace(std::make_pair(
        std::make_pair(pair.first.first.stroke_ind, pair.first.second),
        std::make_pair(pair.second.first.stroke_ind, pair.second.second)));
  }
  // Deduplicate
  for (auto &pair : cut_pair) {
    auto query = std::make_pair(
      std::make_pair(pair.first.first.stroke_ind, pair.first.second),
      std::make_pair(pair.second.first.stroke_ind, pair.second.second));
    auto s_head = query.first;
    auto s_tail = query.second;
    auto head_find = std::find_if(
      matching_pair_indices.begin(), matching_pair_indices.end(),
      [&s_head](
        const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> &p) {
        return p.first == s_head || p.second == s_head;
      });
    auto tail_find = std::find_if(
      matching_pair_indices.begin(), matching_pair_indices.end(),
      [&s_tail](
        const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> &p) {
        return p.first == s_tail || p.second == s_tail;
      });

    // Check if we have cases where one endpoint is matched to multiple strokes
    if (head_find == matching_pair_indices.end() &&
        tail_find == matching_pair_indices.end() &&
        !matching_pair_indices.count(query))
      cut_pair_indices.emplace(query);
  }

  std::vector<std::set<size_t>> merge_group;
  std::map<size_t, size_t> merge_group_index;
  for (const auto &pair : matching_pair_indices) {
    if (!merge_group_index.count(pair.first.first) &&
        !merge_group_index.count(pair.second.first)) { // Both not seen
      merge_group_index[pair.first.first] = merge_group.size();
      merge_group_index[pair.second.first] = merge_group.size();
      merge_group.emplace_back(
        std::set<size_t>{pair.first.first, pair.second.first});
    } else if (merge_group_index.count(pair.first.first) !=
               merge_group_index.count(pair.second.first)) { // One is seen
      auto &seen =
        (merge_group_index.count(pair.first.first)) ? pair.first : pair.second;
      auto &unseen =
        (!merge_group_index.count(pair.first.first)) ? pair.first : pair.second;
      size_t seen_idx = merge_group_index[seen.first];
      merge_group[seen_idx].emplace(unseen.first);
      merge_group_index[unseen.first] = seen_idx;
    } else { // Both seen
      size_t to_idx = merge_group_index[pair.first.first];
      size_t from_idx = merge_group_index[pair.second.first];
      merge_group[to_idx].insert(merge_group[from_idx].begin(),
                                 merge_group[from_idx].end());
      for (auto i : merge_group[from_idx]) {
        merge_group_index[i] = to_idx;
      }
      merge_group[from_idx].clear();
    }
  }

  // Merge strokes
  Capture merged_capture = capture;
  merged_capture.sketchedPolylines.clear();

  auto merge_strokes = [](SketchUI::Polyline2D &to_s,
                          const SketchUI::Polyline2D &from_s, bool is_to_tail,
                          bool is_from_tail, std::pair<size_t, bool> &s_head,
                          std::pair<size_t, bool> &s_tail) {
    // Move one point away from the endpoint to avoid segment with 0-length.
    if (from_s.points.size() <= 1)
      return;
    if (is_to_tail && !is_from_tail) {
      to_s.points.insert(to_s.points.end(), from_s.points.begin() + 1,
                         from_s.points.end());
      s_tail.first = from_s.stroke_ind;
      s_tail.second = true;
    } else if (is_to_tail && is_from_tail) {
      to_s.points.insert(to_s.points.end(), from_s.points.rbegin() + 1,
                         from_s.points.rend());
      s_tail.first = from_s.stroke_ind;
      s_tail.second = false;
    } else if (!is_to_tail && is_from_tail) {
      to_s.points.erase(to_s.points.begin());
      to_s.points.insert(to_s.points.begin(), from_s.points.begin(),
                         from_s.points.end());
      s_head.first = from_s.stroke_ind;
      s_head.second = false;
    } else if (!is_to_tail && !is_from_tail) {
      to_s.points.erase(to_s.points.begin());
      to_s.points.insert(to_s.points.begin(), from_s.points.rbegin(),
                         from_s.points.rend());
      s_head.first = from_s.stroke_ind;
      s_head.second = true;
    }
  };

  std::map<size_t, size_t> cut_map;
  std::set<std::pair<size_t, bool>> valid_cut_points;
  for (const auto &g : merge_group) {
    if (g.empty())
      continue;
    SketchUI::Polyline2D s = *s_map[*g.begin()];
    std::pair<size_t, bool> s_head(s.stroke_ind, false),
      s_tail(s.stroke_ind, true);

    bool merged = false;
    do {
      // Find merging pairs
      auto head_find = std::find_if(
        matching_pair_indices.begin(), matching_pair_indices.end(),
        [&s_head](
          const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>
            &p) { return p.first == s_head || p.second == s_head; });
      auto tail_find = std::find_if(
        matching_pair_indices.begin(), matching_pair_indices.end(),
        [&s_tail](
          const std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>
            &p) { return p.first == s_tail || p.second == s_tail; });
      auto merge_p =
        (head_find != matching_pair_indices.end()) ? head_find : tail_find;
      if (merge_p != matching_pair_indices.end()) {
        auto from_s_p = (merge_p->first == s_head || merge_p->first == s_tail)
                          ? merge_p->second
                          : merge_p->first;
        auto to_s_p = (merge_p->first == s_head || merge_p->first == s_tail)
                        ? merge_p->first
                        : merge_p->second;
        bool is_to_tail = to_s_p == s_tail;
        // Merge
        merge_strokes(s, *s_map[from_s_p.first], is_to_tail, from_s_p.second,
                      s_head, s_tail);
        matching_pair_indices.erase(merge_p);
        merged = true;
      } else {
        merged = false;
      }
    } while (merged);
    // Set the time stamp to be the last one in the group
    s.stroke_ind = *g.rbegin();
    cut_map[s_head.first * 10 + s_head.second] = s.stroke_ind * 10 + 0;
    cut_map[s_tail.first * 10 + s_tail.second] = s.stroke_ind * 10 + 1;
    valid_cut_points.emplace(s_head);
    valid_cut_points.emplace(s_tail);
    merged_capture.sketchedPolylines.emplace_back(s);
  }

  // Add not merged strokes
  for (const auto &s : capture.sketchedPolylines) {
    if (!merge_group_index.count(s.stroke_ind)) {
      merged_capture.sketchedPolylines.emplace_back(s);
      valid_cut_points.emplace(s.stroke_ind, false);
      valid_cut_points.emplace(s.stroke_ind, true);
    }
  }
  std::sort(merged_capture.sketchedPolylines.begin(),
            merged_capture.sketchedPolylines.end(),
            [](const SketchUI::Polyline2D &a, const SketchUI::Polyline2D &b) {
              return a.stroke_ind < b.stroke_ind;
            });
  capture = merged_capture;

  // Output the final cut points
  {
    std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
      final_cuts;
    for (const auto &cut : cut_pair_indices) {
      std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> c = cut;
      if (valid_cut_points.count(c.first) && valid_cut_points.count(c.second)) {
        size_t k1 = c.first.first * 10 + c.first.second;
        if (cut_map.count(k1)) {
          c.first.first = cut_map[k1] / 10;
          c.first.second = cut_map[k1] % 10;
        }
        size_t k2 = c.second.first * 10 + c.second.second;
        if (cut_map.count(k2)) {
          c.second.first = cut_map[k2] / 10;
          c.second.second = cut_map[k2] % 10;
        }
        final_cuts.emplace(c);
      }
    }

    // Check if there are exact duplications
    s_map.clear();
    for (auto &s : capture.sketchedPolylines) {
      s_map[s.stroke_ind] = &s;
    }
    std::map<size_t, size_t> updated_map;

    for (size_t i = 0; i < capture.sketchedPolylines.size(); ++i) {
      updated_map[capture.sketchedPolylines[i].stroke_ind] = i;
      capture.sketchedPolylines[i].stroke_ind = i;
    }

    auto map_s_index = [&updated_map](size_t s_idx) -> size_t {
      assert(updated_map.count(s_idx));
      return updated_map[s_idx];
    };
    if (!cut_csv_filename.empty()) {
      std::ofstream file(cut_csv_filename);
      // file << "stroke_ind_a,is_end_a,stroke_ind_b,is_end_b\n";
      for (const auto &f_pair : final_cuts) {
        std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>> pair =
          f_pair;
        pair.first.first = map_s_index(pair.first.first);
        pair.second.first = map_s_index(pair.second.first);
        assert(pair.first.first != pair.second.first);
        int base_stroke_ind = pair.first.first;
        int target_stroke_ind = pair.second.first;
        file << base_stroke_ind << "," << pair.first.second << ",";
        file << target_stroke_ind << "," << pair.second.second << "\n";
      }
    }
  }
}
