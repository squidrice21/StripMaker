#pragma once

#include <map>
#include <mutex>
#include <ostream>
#include <vector>

#include "Cluster.h"
#include "Context.h"

class StrokeOrientation {
public:
  StrokeOrientation(Context &context);
  std::map<int, std::vector<int>> orient_strokes(Input *input);
  void flip_strokes(Input *input,
                    const std::map<int, std::vector<int>> &orientations);
  void flip_cut_strokes(
    Input *input, const std::map<int, std::vector<int>> &orientations,
    const std::vector<std::pair<std::pair<size_t, bool>,
                                std::pair<size_t, bool>>> &matching_pair,
    std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
      *unmerge_mapping_ptr = nullptr);
  void orientation_debug(std::ostream &os, const Input &input);

  int MAX_VIOLATIONS = 30;
  // int MAX_VIOLATIONS = 10;

private:
  Context &context;

  std::mutex debug_lock;

  struct DebugLine {
    glm::dvec2 from;
    glm::dvec2 to;
    std::string color;
  };
  std::vector<DebugLine> debug_lines;
  void add_debug_line(DebugLine line);

  std::vector<int> orient_cluster_strokes(const Cluster &cluster);

  struct PairOrientation {
    int orientation;
    double weight;
  };
  PairOrientation orient_stroke_pair(const Cluster::Stroke &a,
                                     const Cluster::Stroke &b);

  struct PolicyResult {
    std::vector<std::pair<glm::dvec2, glm::dvec2>> violations;
    std::vector<double> connection_angles;
    std::vector<double> connection_dists;
    double shortest_cut;
  };
  PolicyResult evaluate_policy(std::vector<int> overlaps_a,
                               std::vector<int> overlaps_b,
                               const Cluster::Stroke &a,
                               const Cluster::Stroke &b, int policy);

  double weight_for_angle(double angle);
  int closest_ortho_idx_on_curve(const Cluster::Stroke &from, size_t from_i,
                                 const Cluster::Stroke &to);
  size_t closest_idx_on_curve(const Cluster::Stroke &from, size_t from_i,
                              const Cluster::Stroke &to);
  void fill_overlaps(const Cluster::Stroke &from, const Cluster::Stroke &to,
                     std::vector<int> &overlaps);
};
