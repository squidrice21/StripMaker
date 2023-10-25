/**
 * This file contains helper functions (I/O, file system operations, simple
 * functions shared by various parts).
 */
#pragma once

#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Context.h"
#include "../stroke_strip_src/FittingEigenSparse.h"
#include "../stroke_strip_src/Parameterization.h"
#include "../stroke_strip_src/SketchInfo.h"
#include "../stroke_strip_src/StrokeOrientation.h"

#include "glm/detail/type_vec.hpp"
#include "nonstd/span.hpp"

#include <fstream>
#include <set>
#include <vector>

template <typename T>
using span = ::nonstd::span<T>;

extern Context context;
extern StrokeOrientation fea_orientation;
extern Parameterization fea_param;
extern FittingEigenSparse fea_fitting;

extern double match_dist_threshold;
extern double stroke_sampling_size;
extern double stroke_sampling_size_ref;
extern size_t stroke_sampling_min_num;
extern size_t stroke_sampling_total_ref;

extern double spiral_total_signed_curvature_threshold;
extern double spiral_jump_threshold;

// Train time filter to auto detect incorrect examples due to branching or
// continuity ambiguity, and partial information
extern double fillin_gap_change_threshold;
extern double continuity_overlapping_threshold;
extern double continuity_width_ratio;
extern double large_angle_threshold;

// Test time filter to save computation
extern double fit_angle_threshold;
extern double filter_min_length;
extern double filter_max_overlapping_ratio;
extern double filter_min_overlapping_ratio;

extern size_t min_overlapping_xsec_count;
extern double alignment_threshold;
extern double final_merge_ratio_threshold;

extern double pos_parameterization_alignment_threshold;
extern double neg_parameterization_alignment_threshold;
extern double pointwise_distance_threshold;
extern double sandwiched_weight;
extern size_t feature_bins;
extern size_t end_noisy_range;

extern double end_end_angle;

// Typical cluster classifier
extern bool generate_typical_partial_pos_samples;

struct Input;

inline SketchUI::Polyline2D stroke2polyline(const Cluster::Stroke &stroke) {
  SketchUI::Polyline2D poly;
  for (const auto &p : stroke.points) {
    SketchUI::Point2D pp(p.x, p.y);
    poly.points.emplace_back(std::pair<SketchUI::Point2D, std::int64_t>(pp, 0));
  }

  return poly;
}

void read_input(const std::string &filename, Capture &capture, int &width,
                int &height, bool to_preprocess = true);
void read_input(const std::string &filename, Input &input, int &width,
                int &height, double &input_thickness,
                bool to_preprocess = true);
void read_example(
  const std::string &scap_filename,
  const std::vector<std::pair<std::vector<int>, std::vector<int>>> &samples,
  int &width, int &height, std::vector<std::pair<Cluster, Cluster>> &clusters,
  bool to_preprocess = true);

void make_folder(std::string output_name, bool to_erase = true);
void make_cache_folder(std::string output_name, bool to_erase = true);

inline bool file_exists(const std::string &file_name) {
  std::ifstream f(file_name.c_str());
  return f.good();
}

void reorient_param_input(
  Input &input, bool to_orient, bool test_time,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    &matching_pair);
void fit_cluster(
  double thickness, double width, double height, Cluster &cluster,
  bool to_orient, bool test_time,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr = nullptr);

void transform_fit_curve(FittedCurve &fit, glm::dvec2 orig_center,
                         double orig_thickness);

inline double stroke_sample_distance(const Cluster::Stroke &s1,
                                     const Cluster::Stroke &s2) {
  double dist = std::numeric_limits<double>::infinity();

  auto segments_distance = [](glm::dvec2 v1, glm::dvec2 v2, glm::dvec2 v3,
                              glm::dvec2 v4) -> double {
    auto project_to_seg = [](glm::dvec2 v1, glm::dvec2 v2,
                             glm::dvec2 v3) -> double {
      // Project v3, v4 to v1v2
      glm::dvec2 v21 = v2 - v1;
      glm::dvec2 n_dir(-v21.y, v21.x);

      double a = v3.x - v1.x, b = v3.y - v1.y;
      double t =
        (n_dir.x * b - n_dir.y * a) / (n_dir.x * v21.y - n_dir.y * v21.x);

      // Check
      double d = glm::dot(v1 + t * v21 - v3, v21);

      return t;
    };

    auto dist_to_seg = [&project_to_seg](glm::dvec2 v1, glm::dvec2 v2,
                                         glm::dvec2 v3) -> double {
      double t3 = project_to_seg(v1, v2, v3);

      double d3 = std::min(glm::distance(v1, v3), glm::distance(v2, v3));
      glm::dvec2 v21 = v2 - v1;
      if (t3 > 0 && t3 < 1)
        d3 = std::min(d3, glm::distance(v1 + t3 * v21, v3));

      return d3;
    };

    return std::min(
      dist_to_seg(v1, v2, v3),
      std::min(dist_to_seg(v1, v2, v4),
               std::min(dist_to_seg(v3, v4, v1), dist_to_seg(v3, v4, v2))));
  };

  for (size_t i = 0; i + 1 < s1.points.size(); i++) {
    for (size_t j = 0; j + 1 < s2.points.size(); j++) {
      double seg_d = segments_distance(s1.points[i], s1.points[i + 1],
                                       s2.points[j], s2.points[j + 1]);
      dist = std::min(seg_d, dist);
    }
  }

  return dist;
}

inline glm::dvec2 to_glm_dvec2(const SketchPoint &p) {
  return glm::dvec2(p.x, p.y);
}

inline double stroke_sample_distance(const Sketch &s1, const Sketch &s2) {
  double dist = std::numeric_limits<double>::infinity();

  auto segments_distance = [](glm::dvec2 v1, glm::dvec2 v2, glm::dvec2 v3,
                              glm::dvec2 v4) -> double {
    auto project_to_seg = [](glm::dvec2 v1, glm::dvec2 v2,
                             glm::dvec2 v3) -> double {
      // Project v3, v4 to v1v2
      glm::dvec2 v21 = v2 - v1;
      glm::dvec2 n_dir(-v21.y, v21.x);

      double a = v3.x - v1.x, b = v3.y - v1.y;
      double t =
        (n_dir.x * b - n_dir.y * a) / (n_dir.x * v21.y - n_dir.y * v21.x);

      // Check
      double d = glm::dot(v1 + t * v21 - v3, v21);

      return t;
    };

    auto dist_to_seg = [&project_to_seg](glm::dvec2 v1, glm::dvec2 v2,
                                         glm::dvec2 v3) -> double {
      double t3 = project_to_seg(v1, v2, v3);

      double d3 = std::min(glm::distance(v1, v3), glm::distance(v2, v3));
      glm::dvec2 v21 = v2 - v1;
      if (t3 > 0 && t3 < 1)
        d3 = std::min(d3, glm::distance(v1 + t3 * v21, v3));

      return d3;
    };

    return std::min(
      dist_to_seg(v1, v2, v3),
      std::min(dist_to_seg(v1, v2, v4),
               std::min(dist_to_seg(v3, v4, v1), dist_to_seg(v3, v4, v2))));
  };

  for (size_t i = 0; i + 1 < s1.points.size(); i++) {
    for (size_t j = 0; j + 1 < s2.points.size(); j++) {
      double seg_d = segments_distance(
        to_glm_dvec2(s1.points[i].first), to_glm_dvec2(s1.points[i + 1].first),
        to_glm_dvec2(s2.points[j].first), to_glm_dvec2(s2.points[j + 1].first));
      dist = std::min(seg_d, dist);
    }
  }

  return dist;
}

inline double stroke_sample_distance(glm::dvec2 v, const Cluster::Stroke &s) {
  double dist = std::numeric_limits<double>::infinity();

  auto segments_distance = [](glm::dvec2 v1, glm::dvec2 v2,
                              glm::dvec2 v3) -> double {
    auto project_to_seg = [](glm::dvec2 v1, glm::dvec2 v2,
                             glm::dvec2 v3) -> double {
      // Project v3, v4 to v1v2
      glm::dvec2 v21 = v2 - v1;
      glm::dvec2 n_dir(-v21.y, v21.x);

      double a = v3.x - v1.x, b = v3.y - v1.y;
      double t =
        (n_dir.x * b - n_dir.y * a) / (n_dir.x * v21.y - n_dir.y * v21.x);

      // Check
      double d = glm::dot(v1 + t * v21 - v3, v21);

      return t;
    };

    auto dist_to_seg = [&project_to_seg](glm::dvec2 v1, glm::dvec2 v2,
                                         glm::dvec2 v3) -> double {
      double t3 = project_to_seg(v1, v2, v3);

      double d3 = std::min(glm::distance(v1, v3), glm::distance(v2, v3));
      glm::dvec2 v21 = v2 - v1;
      if (t3 > 0 && t3 < 1)
        d3 = std::min(d3, glm::distance(v1 + t3 * v21, v3));

      return d3;
    };

    return dist_to_seg(v1, v2, v3);
  };

  for (size_t i = 0; i + 1 < s.points.size(); i++) {
    double seg_d = segments_distance(s.points[i], s.points[i + 1], v);
    dist = std::min(seg_d, dist);
  }

  return dist;
}

inline double stroke_sample_distance(glm::dvec2 v, const Sketch &s) {
  Cluster::Stroke ss;
  for (const auto &p : s.points)
    ss.points.emplace_back(to_glm_dvec2(p.first));

  return stroke_sample_distance(v, ss);
}

inline double stroke_sample_distance(glm::dvec2 v, const FittedCurve &s) {
  Cluster::Stroke ss;
  for (const auto &p : s.centerline)
    ss.points.emplace_back(p);

  return stroke_sample_distance(v, ss);
}

inline double stroke_sample_projection_time(glm::dvec2 v, const Sketch &s) {
  double dist = std::numeric_limits<double>::infinity();

  auto segments_distance = [](glm::dvec2 v1, glm::dvec2 v2, glm::dvec2 v3,
                              double &seg_t) -> double {
    auto project_to_seg = [](glm::dvec2 v1, glm::dvec2 v2,
                             glm::dvec2 v3) -> double {
      // Project v3, v4 to v1v2
      glm::dvec2 v21 = v2 - v1;
      glm::dvec2 n_dir(-v21.y, v21.x);

      double a = v3.x - v1.x, b = v3.y - v1.y;
      double t =
        (n_dir.x * b - n_dir.y * a) / (n_dir.x * v21.y - n_dir.y * v21.x);

      // Check
      double d = glm::dot(v1 + t * v21 - v3, v21);

      return t;
    };

    auto dist_to_seg = [&project_to_seg](glm::dvec2 v1, glm::dvec2 v2,
                                         glm::dvec2 v3,
                                         double &seg3_t) -> double {
      double t3 = project_to_seg(v1, v2, v3);

      double d13 = glm::distance(v1, v3);
      double d23 = glm::distance(v2, v3);
      double d3 = 0;
      if (d13 < d23) {
        d3 = d13;
        seg3_t = 0;
      } else {
        d3 = d23;
        seg3_t = 1;
      }

      glm::dvec2 v21 = v2 - v1;
      if (t3 > 0 && t3 < 1) {
        double d = glm::distance(v1 + t3 * v21, v3);
        if (d < d3) {
          d3 = d;
          seg3_t = t3;
        }
      }

      return d3;
    };

    return dist_to_seg(v1, v2, v3, seg_t);
  };

  size_t min_i = 0;
  double min_t = 0;
  for (size_t i = 0; i + 1 < s.points.size(); i++) {
    double t;
    double seg_d = segments_distance(to_glm_dvec2(s.points[i].first),
                                     to_glm_dvec2(s.points[i + 1].first), v, t);
    if (seg_d < dist) {
      dist = seg_d;
      min_i = i;
      min_t = t;
    }
  }

  double s_t = s.get_t(min_i, min_t);

  return s_t;
}
inline size_t stroke_sample_projection_index(glm::dvec2 v,
                                             const Cluster::Stroke &s,
                                             double &dist) {
  dist = std::numeric_limits<double>::infinity();

  auto segments_distance = [](glm::dvec2 v1, glm::dvec2 v2, glm::dvec2 v3,
                              double &seg_t) -> double {
    auto project_to_seg = [](glm::dvec2 v1, glm::dvec2 v2,
                             glm::dvec2 v3) -> double {
      // Project v3, v4 to v1v2
      glm::dvec2 v21 = v2 - v1;
      glm::dvec2 n_dir(-v21.y, v21.x);

      double a = v3.x - v1.x, b = v3.y - v1.y;
      double t =
        (n_dir.x * b - n_dir.y * a) / (n_dir.x * v21.y - n_dir.y * v21.x);

      // Check
      double d = glm::dot(v1 + t * v21 - v3, v21);

      return t;
    };

    auto dist_to_seg = [&project_to_seg](glm::dvec2 v1, glm::dvec2 v2,
                                         glm::dvec2 v3,
                                         double &seg3_t) -> double {
      double t3 = project_to_seg(v1, v2, v3);

      double d13 = glm::distance(v1, v3);
      double d23 = glm::distance(v2, v3);
      double d3 = 0;
      if (d13 < d23) {
        d3 = d13;
        seg3_t = 0;
      } else {
        d3 = d23;
        seg3_t = 1;
      }

      glm::dvec2 v21 = v2 - v1;
      if (t3 > 0 && t3 < 1) {
        double d = glm::distance(v1 + t3 * v21, v3);
        if (d < d3) {
          d3 = d;
          seg3_t = t3;
        }
      }

      return d3;
    };

    return dist_to_seg(v1, v2, v3, seg_t);
  };

  size_t min_i = 0;
  double min_t = 0;
  for (size_t i = 0; i + 1 < s.points.size(); i++) {
    double t;
    double seg_d = segments_distance(s.points[i], s.points[i + 1], v, t);
    if (seg_d < dist) {
      dist = seg_d;
      min_i = i;
      min_t = t;
    }
  }

  return (min_t < 0.5) ? min_i : min_i + 1;
}

void map_xsec(Cluster::XSec &xsec,
              const std::map<std::pair<size_t, size_t>,
                             std::pair<size_t, size_t>> &unmerge_mapping);
void map_strokes(const std::map<std::pair<size_t, size_t>,
                                std::pair<size_t, size_t>> &unmerge_mapping,
                 std::vector<Cluster::Stroke> &strokes);
void update_xsec_tangent(Cluster::XSec &xsec,
                         const std::vector<Cluster::Stroke> &strokes);
std::vector<Cluster::XSec>
filter_xsecs(const std::vector<Cluster::XSec> &xsecs);

bool is_overlapping(
  const Cluster::XSec &xsecs,
  const std::map<size_t, size_t> &s2c = std::map<size_t, size_t>());
bool is_overlapping_stroke(const Cluster &cluster, size_t sid1, size_t sid2);
bool is_overlapping_cluster(const Cluster &cluster, size_t sid1,
                            const Cluster &cluster2, int key_sid_ind = -1);

bool is_separatible(const Cluster::XSec &xsec);

double u_distance(const Cluster &cluster1, const Cluster &cluster2,
                  const Cluster &merged_cluster);

void reuse_cluster(const Cluster &fit_cluster, Cluster &reuse_cluster,
                   bool high_reso = false);

bool read_cache(const std::string &file_name, Cluster &cluster);
void dump_clusters(Capture &in_capture, std::map<size_t, Cluster> &clusters,
                   Input &input, int width, int height, double input_thickness,
                   std::string output_folder, std::string out_name,
                   std::string fit_name);
