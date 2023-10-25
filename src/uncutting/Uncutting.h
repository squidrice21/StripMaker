#pragma once

#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/SketchInfo.h"

#include <functional>
#include <map>
#include <ostream>
#include <sstream>
#include <vector>

double hausdorff_distance_one_sided(const Sketch &s1, const Sketch &s2);
double hausdorff_distance(const Sketch &s1, const Sketch &s2);
double hausdorff_distance(const FittedCurve &s1, const FittedCurve &s2);
double hausdorff_distance_one_sided(const FittedCurve &s1,
                                    const FittedCurve &s2);

void map_to_ref(const Capture &capture, const Capture &capture_ref,
                std::map<size_t, size_t> &in2ref);

glm::dvec2 capture_center(const Capture &cluster, double &width,
                          double &height);

void matching_point(
  Input &input,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair);
void check_GT_cluster(
  Input &input,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair);

void check_tangent(
  Input &input,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair);

void check_original_input(
  const Capture &capture, const Capture &capture_orig,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair);

void connect_stroke(
  Capture &capture,
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>> &matching_pair,
  const std::string &csv_filename);
void connect_stroke(
  Capture &capture,
  const std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                              std::pair<Cluster::Stroke, bool>>> &matching_pair,
  const std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                              std::pair<Cluster::Stroke, bool>>> &cut_pair,
  const std::string &vis_csv_filename, const std::string &cut_csv_filename);

glm::dvec2 stepaway_tangent(const Cluster::Stroke &stroke, const bool is_end,
                            const double stepaway_amount);
