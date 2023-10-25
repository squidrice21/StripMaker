#include "Util.h"

#include "../stroke_strip_src/Serialization.h"
#include "../stroke_strip_src/StrokeCutting.h"
#include "../stroke_strip_src/Utils.h"
#include "Capture2Input.h"
#include "glm/detail/type_vec.hpp"
#include "measure/Measurement.h"

#ifndef __APPLE__
#include <direct.h>
#else
#include <sys/stat.h>
#endif

#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_map>
#include <utility>

double pos_parameterization_alignment_threshold = 3;
double neg_parameterization_alignment_threshold = 3;

double pointwise_distance_threshold = 10;

size_t shading_flag = 3;

double spiral_total_signed_curvature_threshold = (360. + 30) / 180 * M_PI;

void read_input(const std::string &filename, Capture &capture, int &width,
                int &height, bool to_preprocess) {
  std::ifstream scap_ifs(filename);
  std::stringstream buffer;
  buffer << scap_ifs.rdbuf();
  scap_ifs.close();

  // Capture capture;
  capture.from_string(buffer.str());

  std::string canvas_size_line =
    buffer.str().substr(0, buffer.str().find_first_of("\n"));
  sscanf(canvas_size_line.c_str(), "#%d\t%d", &width, &height);

  // Assign missing ids
  int max_ind = -1;
  int max_cind = -1;
  for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
    max_ind = std::max(max_ind, capture.getSketchedPolyline(i).stroke_ind);
    max_cind = std::max(max_cind, capture.getSketchedPolyline(i).group_ind);
  }

  if (to_preprocess) {
    bool to_reindex = false;
    preprocess_cluster(1, width, height, &capture, to_reindex);
  }
  std::sort(capture.sketchedPolylines.begin(), capture.sketchedPolylines.end(),
            [](const SketchUI::Polyline2D &a, const SketchUI::Polyline2D &b) {
              return a.stroke_ind < b.stroke_ind;
            });
}

void read_input(const std::string &filename, Input &input, int &width,
                int &height, double &input_thickness, bool to_preprocess) {
  Capture capture;
  read_input(filename, capture, width, height, to_preprocess);

  Capture2Input capture2input;
  input = capture2input.from_capture(capture);
  input_thickness = capture.thickness;
  input.thickness = 1;
  input.orig_thickness = input_thickness;
}

void make_folder(std::string output_name, bool to_erase) {
  if (to_erase)
    output_name.erase(output_name.length() - 5, 5);
  char *output_folder = strdup(output_name.c_str());
#ifndef __APPLE__
  int ret = _mkdir(output_folder);
#else
  int ret = mkdir(output_folder, S_IRWXU);
#endif

  std::string output_isoline_folder = output_name + "/input";
  output_folder = strdup(output_isoline_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  output_isoline_folder = output_name + "/input_cluster";
  output_folder = strdup(output_isoline_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  output_isoline_folder = output_name + "/isoline";
  output_folder = strdup(output_isoline_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  output_isoline_folder = output_name + "/parameterization";
  output_folder = strdup(output_isoline_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  std::string output_fit_folder = output_name + "/fit";
  output_folder = strdup(output_fit_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  std::string output_ori_folder = output_name + "/orientation";
  output_folder = strdup(output_ori_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  std::string output_cache_folder = output_name + "/cache";
  output_folder = strdup(output_cache_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif

  std::string output_steps_folder = output_name + "/steps";
  output_folder = strdup(output_steps_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif
}

void make_cache_folder(std::string output_name, bool to_erase) {
  if (to_erase)
    output_name.erase(output_name.length() - 5, 5);
  char *output_folder = strdup(output_name.c_str());
#ifndef __APPLE__
  int ret = _mkdir(output_folder);
#else
  int ret = mkdir(output_folder, S_IRWXU);
#endif

  std::string output_cache_folder = output_name + "/cache";
  output_folder = strdup(output_cache_folder.c_str());
#ifndef __APPLE__
  ret = _mkdir(output_folder);
#else
  ret = mkdir(output_folder, S_IRWXU);
#endif
}

void read_example(
  const std::string &scap_filename,
  const std::vector<std::pair<std::vector<int>, std::vector<int>>> &samples,
  int &width, int &height, std::vector<std::pair<Cluster, Cluster>> &clusters,
  bool to_preprocess) {
  Input input;
  double input_thickness;
  read_input(scap_filename, input, width, height, input_thickness,
             to_preprocess);

  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : input.clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      sid2s[s.stroke_ind] = s;
    }
  }

  for (const auto &ss : samples) {
    clusters.emplace_back();
    for (auto s_idx : ss.first) {
      clusters.back().first.strokes.emplace_back(sid2s[s_idx]);
      clusters.back().first.strokes.back().cluster_ind = -1;
    }
    for (auto s_idx : ss.second) {
      clusters.back().second.strokes.emplace_back(sid2s[s_idx]);
      clusters.back().second.strokes.back().cluster_ind = -1;
    }
  }
}

void reorient_param_input(
  Input &in_input, bool to_orient, bool test_time,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    &matching_pair) {
  for (auto &idx_cluster : in_input.clusters) {
    auto &cluster = idx_cluster.second;

    // Reset the input stroke orientations
    assert(cluster.strokes.size() == cluster.original_input_strokes.size());
    cluster.strokes = cluster.original_input_strokes;

    Input input;
    input.thickness = in_input.thickness;
    input.width = in_input.width;
    input.height = in_input.height;
    input.clusters[0] = cluster;
    input.clusters[0].xsecs.clear();

    std::vector<size_t> ori_max_violation =
      (to_orient) ? std::vector<size_t>{30, 10} : std::vector<size_t>{30};
    Input input_init = input;
    int saved_MAX_VIOLATIONS = fea_orientation.MAX_VIOLATIONS;
    std::map<std::string, double> obj_term_values;
    std::map<std::pair<size_t, size_t>, std::pair<size_t, size_t>>
      unmerge_mapping_itr;
    for (auto max_vio : ori_max_violation) {
      input = input_init;
      if (to_orient) {
        fea_orientation.MAX_VIOLATIONS = max_vio;
        std::map<int, std::vector<int>> orientations =
          fea_orientation.orient_strokes(&input);
        unmerge_mapping_itr.clear();
        fea_orientation.flip_cut_strokes(&input, orientations, matching_pair,
                                         &unmerge_mapping_itr);
        if (input.clusters[0].strokes.empty()) {
          input = input_init;
          unmerge_mapping_itr.clear();
          fea_orientation.flip_cut_strokes(
            &input, orientations,
            std::vector<
              std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>(),
            &unmerge_mapping_itr);
        }
        fea_orientation.MAX_VIOLATIONS = saved_MAX_VIOLATIONS;
      }
      obj_term_values = fea_param.parameterize(&input);
      cluster.param_failed = input.clusters[0].param_failed;
      cluster.periodic = input.clusters[0].periodic;
      if (obj_term_values["alignment"] <
            pos_parameterization_alignment_threshold &&
          !input.clusters[0].strokes.empty())
        break;
    }
    // Sort sample points along the same cross section based on their V values
    for (auto &xsec : input.clusters[0].xsecs) {
      sort_xsec(xsec);
      map_xsec(xsec, unmerge_mapping_itr);
    }
    if (!input.clusters[0].strokes.empty() &&
        obj_term_values["alignment"] < pos_parameterization_alignment_threshold)
      cluster.obj_term_values = obj_term_values;
    cluster.orthogonal_xsecs_list = input.clusters[0].orthogonal_xsecs_list;
    cluster.strokes = input.clusters[0].strokes;
    map_strokes(unmerge_mapping_itr, cluster.strokes);
    cluster.xsecs = input.clusters[0].xsecs;

    if (!unmerge_mapping_itr.empty()) {
      for (auto &xsec : cluster.xsecs) {
        update_xsec_tangent(xsec, cluster.strokes);
      }
    }
  }
}

void fit_cluster(
  double thickness, double width, double height, Cluster &cluster,
  bool to_orient, bool test_time,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr) {
  assert(!cluster.strokes.empty());
  if (cluster.strokes.size() == 1) {
    cluster.fit.centerline = cluster.strokes.front().points;
    return;
  }
  Input input;
  input.thickness = thickness;
  input.width = width;
  input.height = height;
  input.clusters[0] = cluster;
  input.clusters[0].xsecs.clear();

  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    matching_pair;
  reorient_param_input(input, to_orient, test_time,
                       (matching_pair_ptr) ? *matching_pair_ptr
                                           : matching_pair);

  // If the default fails to orient the inputs, assign empty fit curve
  if (!input.clusters[0].strokes.empty() &&
      input.clusters[0].obj_term_values["alignment"] <
        pos_parameterization_alignment_threshold) {
    const auto &fit_map = fea_fitting.fit(&input);
    const auto &fit_curve = fit_map.at(0);
    cluster.fit = fit_curve;
    cluster.obj_term_values = input.clusters[0].obj_term_values;
    cluster.periodic = input.clusters[0].periodic;
  }

  cluster.strokes = input.clusters[0].strokes;
  cluster.xsecs = input.clusters[0].xsecs;
}

void transform_fit_curve(FittedCurve &fit, glm::dvec2 orig_center,
                         double orig_thickness) {
  // Scale
  for (auto &p : fit.centerline) {
    p *= orig_thickness;
  }
  for (auto &w : fit.widths) {
    w *= orig_thickness;
  }

  // Translate
  for (auto &p : fit.centerline) {
    p += orig_center;
  }
}

void map_xsec(Cluster::XSec &xsec,
              const std::map<std::pair<size_t, size_t>,
                             std::pair<size_t, size_t>> &unmerge_mapping) {
  if (unmerge_mapping.empty())
    return;
  for (auto &p : xsec.points) {
    auto q = std::make_pair(p.stroke_ind, size_t(std::ceil(p.i)));
    if (unmerge_mapping.count(q)) {
      auto record = unmerge_mapping.at(q);
      p.stroke_ind = record.first;
      if (record.second - (size_t(std::ceil(p.i)) - p.i) < 0)
        assert(p.i > 0);
      p.i = record.second - (size_t(std::ceil(p.i)) - p.i);
    }
  }
}

void map_strokes(const std::map<std::pair<size_t, size_t>,
                                std::pair<size_t, size_t>> &unmerge_mapping,
                 std::vector<Cluster::Stroke> &strokes) {
  std::vector<Cluster::Stroke> in_strokes = strokes;
  strokes.clear();
  for (const auto &s : in_strokes) {
    if (!unmerge_mapping.count(std::make_pair(s.stroke_ind, size_t(1)))) {
      strokes.emplace_back(s);
    } else {
      std::unordered_map<size_t, Cluster::Stroke> unmerged_strokes;
      for (size_t i = 0; i < s.points.size(); ++i) {
        auto q = std::make_pair(s.stroke_ind, i);
        assert(unmerge_mapping.count(q));
        auto record = unmerge_mapping.at(q);
        unmerged_strokes[record.first].stroke_ind = record.first;
        unmerged_strokes[record.first].cluster_ind = s.cluster_ind;
        if (unmerged_strokes[record.first].points.size() < record.second &&
            record.second == 1) {
          assert(i > 0);
          unmerged_strokes[record.first].points.emplace_back(s.points[i - 1]);
          unmerged_strokes[record.first].u.emplace_back(s.u[i - 1]);
        }
        unmerged_strokes[record.first].points.emplace_back(s.points[i]);
        unmerged_strokes[record.first].u.emplace_back(s.u[i]);
        if (i > 0 &&
            unmerge_mapping.at(std::make_pair(s.stroke_ind, i - 1)).first ==
              record.first &&
            unmerge_mapping.at(std::make_pair(s.stroke_ind, i - 1)).second >
              record.second &&
            record.second == 1 &&
            unmerge_mapping.at(std::make_pair(s.stroke_ind, i + 1)).second !=
              0) {
          unmerged_strokes[record.first].points.emplace_back(s.points[i + 1]);
          unmerged_strokes[record.first].u.emplace_back(s.u[i + 1]);
        }
      }
      for (auto const &ss : unmerged_strokes) {
        strokes.emplace_back(ss.second);
      }
    }
  }
  std::sort(strokes.begin(), strokes.end(),
            [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
              return a.stroke_ind < b.stroke_ind;
            });
}

void update_xsec_tangent(Cluster::XSec &xsec,
                         const std::vector<Cluster::Stroke> &strokes) {
  std::map<size_t, const Cluster::Stroke *> sid2s;
  for (const auto &s : strokes) {
    sid2s[s.stroke_ind] = &s;
  }
  for (auto &p : xsec.points) {
    if (sid2s[p.stroke_ind]->points.size()) {
      p.tangent = tangent(
        sid2s[p.stroke_ind]->points,
        double(std::min(std::max(0.0, p.i),
                        double(sid2s[p.stroke_ind]->points.size() - 1))));
    }
  }
}

std::vector<Cluster::XSec>
filter_xsecs(const std::vector<Cluster::XSec> &xsecs) {
  std::vector<Cluster::XSec> filtered_xsecs;
  filtered_xsecs.reserve(xsecs.size());
  for (auto &xsec : xsecs) {
    Cluster::XSec filtered_xsec;
    filtered_xsec.u = xsec.u;
    filtered_xsec.connector = xsec.connector;
    for (size_t i = 0; i < xsec.points.size(); ++i) {
      if (xsec.points[i].cluster_ind < std::numeric_limits<size_t>::max())
        filtered_xsec.points.emplace_back(xsec.points[i]);
    }
    if (filtered_xsec.points.empty())
      continue;
    for (size_t i = 0; i + 1 < filtered_xsec.points.size(); ++i) {
      for (size_t j = i + 1; j < filtered_xsec.points.size(); ++j) {
        filtered_xsec.connections.push_back({i, j, 1.0});
      }
    }

    filtered_xsec.center_idx = filtered_xsec.points.size() - 1;
    filtered_xsecs.emplace_back(filtered_xsec);
  }

  return filtered_xsecs;
}

bool is_overlapping(const Cluster::XSec &xsecs,
                    const std::map<size_t, size_t> &s2c) {
  std::unordered_set<size_t> c_indices;
  for (const auto &p : xsecs.points) {
    if (s2c.empty())
      c_indices.emplace(p.cluster_ind);
    else if (s2c.count(p.stroke_ind)) {
      c_indices.emplace(s2c.at(p.stroke_ind));
    }
  }
  return c_indices.size() > 1;
}

bool is_overlapping_stroke(const Cluster &cluster, size_t sid1, size_t sid2) {
  assert(sid1 != sid2);
  for (const auto &xsec : cluster.xsecs) {
    std::set<size_t> s_indices;
    for (const auto &p : xsec.points) {
      if (p.stroke_ind == sid1 || p.stroke_ind == sid2)
        s_indices.emplace(p.stroke_ind);
      if (s_indices.size() == 2)
        return true;
    }
  }

  return false;
}

bool is_overlapping_cluster(const Cluster &cluster, size_t sid1,
                            const Cluster &cluster2, int key_sid_ind) {
  std::set<size_t> c_sids;
  for (const auto &s : cluster2.strokes) {
    c_sids.emplace(s.stroke_ind);
  }

  assert(!c_sids.count(sid1));
  for (size_t i = 0; i < cluster.xsecs.size(); ++i) {
    if (key_sid_ind >= 0 && (int)i != key_sid_ind)
      continue;
    const auto &xsec = cluster.xsecs[i];
    bool seen_sid1 = false;
    bool seen_sid2 = false;
    for (const auto &p : xsec.points) {
      if (p.stroke_ind == sid1)
        seen_sid1 = true;
      if (c_sids.count(p.stroke_ind))
        seen_sid2 = true;
      if (seen_sid1 && seen_sid2)
        return true;
    }
  }

  return false;
}

bool is_separatible(const Cluster::XSec &xsec) {
  bool seen_cross = false;
  for (size_t i = 0; i + 1 < xsec.points.size(); ++i) {
    if (xsec.points[i].cluster_ind != xsec.points[i + 1].cluster_ind) {
      // The two clusters are not clearly divided
      if (seen_cross) {
        return false;
      }
      seen_cross = true;
    }
  }

  return true;
}

double u_distance(const Cluster &cluster1, const Cluster &cluster2,
                  const Cluster &merged_cluster) {
  auto get_u_range =
    [&merged_cluster](const Cluster &cluster) -> std::pair<double, double> {
    std::pair<double, double> range(std::numeric_limits<double>::infinity(),
                                    -std::numeric_limits<double>::infinity());
    std::set<size_t> sids;
    for (auto const &s : cluster.strokes) {
      sids.emplace(s.stroke_ind);
    }
    for (auto const &xsec : merged_cluster.xsecs) {
      for (auto const &p : xsec.points) {
        if (sids.count(p.stroke_ind)) {
          range.first = std::min(range.first, xsec.u);
          range.second = std::max(range.second, xsec.u);
          break;
        }
      }
    }

    return range;
  };
  auto range1 = get_u_range(cluster1);
  auto range2 = get_u_range(cluster2);
  if ((range1.second >= range2.first && range1.second <= range2.second) ||
      (range2.second >= range1.first && range2.second <= range1.second))
    return 0;
  if (range1.second < range2.first)
    return range2.first - range1.second;
  else
    return range1.first - range2.second;
}

void reuse_cluster(const Cluster &fit_cluster, Cluster &reuse_cluster,
                   bool high_reso) {
  reuse_cluster.fit = fit_cluster.fit;

  // Instead of using the actual stroke set, use the fitting curves to
  // check overlapping.
  Cluster::Stroke s;
  s.cluster_ind = fit_cluster.strokes.front().cluster_ind;
  s.stroke_ind = fit_cluster.strokes.front().stroke_ind;

  Sketch reparam_s;
  reparam_s.from_fits(fit_cluster.fit, 0);
  if (!high_reso) {
    double rate = stroke_sampling_size;

    reparam_s.reparameterize(
      std::min(rate, reparam_s.totalLen() / stroke_sampling_min_num));
  } else {
    double rate = 0.5;
    reparam_s.reparameterize(std::min(rate, reparam_s.totalLen() / 10));
  }

  s.points.reserve(reparam_s.points.size());
  for (const auto &p : reparam_s.points) {
    s.points.emplace_back(p.first.x, p.first.y);
  }

  // Copy the spiral flag
  for (const auto &ss : fit_cluster.strokes) {
    s.spiral |= ss.spiral;
  }

  reuse_cluster.strokes.emplace_back(s);
  reuse_cluster.original_input_strokes.emplace_back(s);
}

bool read_cache(const std::string &file_name, Cluster &cluster) {
  auto file_exists = [](const std::string &file_name) -> bool {
    std::ifstream f(file_name.c_str());
    return f.good();
  };
  if (file_exists(file_name)) {
    read_json(file_name, cluster);
    return true;
  }
  return false;
}

void dump_clusters(Capture &in_capture, std::map<size_t, Cluster> &clusters,
                   Input &input, int width, int height, double input_thickness,
                   std::string output_folder, std::string out_name,
                   std::string fit_name) {
  std::map<int, FittedCurve> fit_map_final;
  std::map<size_t, size_t> s2c;
  for (const auto &cid2c : clusters) {
    fit_map_final[cid2c.first] = cid2c.second.fit;
    for (const auto &s : cid2c.second.strokes) {
      s2c[s.stroke_ind] = cid2c.first;
    }
  }

  for (auto &s : in_capture.sketchedPolylines) {
    assert(s2c.count(s.stroke_ind));
    s.group_ind = s2c[s.stroke_ind];
  }
  std::ofstream scap_ofs(output_folder + "/" + out_name + "_out.scap");
  std::string out_buffer = in_capture.to_string();
  scap_ofs << "#" << width << "\t" << height << std::endl;
  scap_ofs.write(out_buffer.c_str(), out_buffer.size());
  scap_ofs.close();

  for (auto &c_fit : fit_map_final)
    transform_fit_curve(c_fit.second, input.orig_center, input_thickness);
  std::ofstream fit_svg(output_folder + "/" + fit_name + ".svg");
  input.thickness = input_thickness;
  fea_fitting.fit_svg(fit_svg, input, fit_map_final, width, height);
  input.thickness = 1;

  {
    Capture capture;
    for (const auto &c : clusters) {
      capture.sketchedPolylines.emplace_back();
      capture.sketchedPolylines.back().stroke_ind = c.first;
      capture.sketchedPolylines.back().group_ind = c.first;
      for (const auto &p : fit_map_final.at(c.first).centerline)
        capture.sketchedPolylines.back().points.emplace_back(
          SketchUI::Point2D(p.x, p.y), 0);
    }
    capture.thickness = input_thickness;

    {
      std::string scap_out =
        std::regex_replace(output_folder + "/" + fit_name + "_width.svg",
                           std::regex("svg"), "scap");
      std::ofstream scap_ofs(scap_out);
      std::string out_buffer = capture.to_string();
      scap_ofs << "#" << width << "\t" << height << std::endl;
      scap_ofs.write(out_buffer.c_str(), out_buffer.size());
      scap_ofs.close();
    }
    {
      std::string scap_out = std::regex_replace(
        output_folder + "/" + fit_name + ".svg", std::regex("svg"), "scap");
      std::ofstream scap_ofs(scap_out);
      std::string out_buffer = capture.to_string();
      scap_ofs << "#" << width << "\t" << height << std::endl;
      scap_ofs.write(out_buffer.c_str(), out_buffer.size());
      scap_ofs.close();
    }
  }

  fit_map_final.clear();
  bool in_width = context.widths;
  context.widths = true;
  Input input_final;
  double norm_thickness = 1;
  for (auto &cid2c : clusters) {
    input_final.clusters.clear();
    input_final.clusters[cid2c.first] = cid2c.second;
    if (input_final.clusters[cid2c.first].xsecs.empty())
      fea_param.parameterize(&input_final);
    fit_map_final[cid2c.first] = cid2c.second.fit;
    fit_map_final[cid2c.first] = fea_fitting.fit(&input_final).begin()->second;
    cid2c.second.fit = fit_map_final[cid2c.first];
    cid2c.second.fit.cluster_idx = cid2c.first;
  }
  std::ofstream fit_svg_width(output_folder + "/" + fit_name + "_width.svg");
  for (auto &c_fit : fit_map_final)
    transform_fit_curve(c_fit.second, input.orig_center, input_thickness);
  input.thickness = input_thickness;
  fea_fitting.fit_svg(fit_svg_width, input, fit_map_final, width, height);
  context.widths = in_width;
}
