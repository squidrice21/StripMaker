#include "Parameterization.h"
#include "../src/Util.h"
#include "SvgUtils.h"
#include "Utils.h"

#include <glm/detail/type_vec.hpp>
#include <glm/gtx/norm.hpp>
#include <limits>
#include <oneapi/tbb/parallel_for.h>

#include <chrono>
#include <cmath>
#include <complex>
#include <deque>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// For solve retrial
double lambda_init = 1e-11;
double lambda_init_l1 = 1e-7;
size_t max_retry = 5;
bool regularize_by_default = false;

// For periodic
double periodic_alignment_threshold = 2;

// For spiral
double periodic_offset_threshold = M_PI / 2;

static void index_xsecs(Cluster &cluster) {
  std::map<size_t, size_t> s_num_xsecs;
  std::map<size_t, std::vector<Cluster::XSec *>> xsecs_sort;

  for (const auto &s : cluster.strokes) {
    s_num_xsecs[s.stroke_ind] = 0;
    xsecs_sort[s.stroke_ind] = std::vector<Cluster::XSec *>();
  }

  std::unordered_set<size_t> seen_s;
  for (auto &xsec : cluster.xsecs) {
    seen_s.clear();

    for (const auto &p : xsec.points) {
      if (!seen_s.count(p.stroke_ind)) {
        // Count number of xsecs per stroke
        s_num_xsecs[p.stroke_ind]++;

        // Order xsecs into per-stroke bins
        xsecs_sort[p.stroke_ind].emplace_back(&xsec);

        seen_s.emplace(p.stroke_ind);
      }
    }
  }

  // Sort xsecs
  for (auto &s_xsecs : xsecs_sort) {
    std::sort(s_xsecs.second.begin(), s_xsecs.second.end(),
              [](const Cluster::XSec *a, const Cluster::XSec *b) {
                return a->u < b->u;
              });
    size_t xsec_i = 0;
    for (auto xsec : s_xsecs.second) {
      for (auto &p : xsec->points) {
        if (p.stroke_ind == s_xsecs.first) {
          p.stroke_within_idx = xsec_i;
          p.stroke_xsec_count = s_num_xsecs[s_xsecs.first];
        }
      }
      xsec_i++;
    }
  }
}

Parameterization::Parameterization(Context &context)
  : context(context) {}

std::map<std::string, double> Parameterization::parameterize(Input *input) {
  std::map<std::string, double> obj_term_values;
  map_clusters<int>(*input, [&](Cluster &c) {
    auto begin = std::chrono::high_resolution_clock::now();

    // Cut the spiral strokes
    std::vector<Cluster::Stroke> spiral_ss;
    for (auto &s : c.strokes) {
      std::vector<Cluster::Stroke> cut_ss;
      if (s.spiral) {
        // Compute spiral angles
        compute_spiral_angles(s);
        cut_spiral_stroke(s, cut_ss, s.spiral_cut_angle);
        spiral_ss.insert(spiral_ss.end(), cut_ss.begin(), cut_ss.end());
      } else {
        s.spiral_angles.resize(s.points.size(), -1);
        spiral_ss.emplace_back(s);
      }
    }
    // Don't reorder so we have the cut strokes together in the vector
    c.strokes = spiral_ss;

    // Recenter avoid numeric overflowing
    glm::dvec2 center;
    center_cluster(c, center);

    obj_term_values = parameterize_cluster(&c);

    if (!c.strokes.empty())
      index_xsecs(c);

    decenter_cluster(c, center);

    // Put the cut spiral strokes back
    spiral_ss.clear();
    assemble_spiral_stroke(c.strokes, spiral_ss, c.xsecs);
    c.strokes = spiral_ss;

    auto end = std::chrono::high_resolution_clock::now();
    context.parameterization_timer +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
        .count() /
      1000.0;

    return 0;
  });

  return obj_term_values;
}

void Parameterization::isolines_svg(std::ostream &os, int new_stroke_ind,
                                    const Input &input) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);
  input.cluster_svg(os, new_stroke_ind, [&](std::ostream &os) {
    for (auto &kv : input.clusters) {
      for (auto &xsec : kv.second.xsecs) {
        for (auto &connection : xsec.connections) {
          if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
            std::stringstream ss;
            ss << "#";
            if (connection.weight > 0.5) {
              ss << std::setfill('0') << std::setw(2) << std::hex
                 << int(map(connection.weight, 0.5, 1.0, 255.0, 0.0)); // red
              ss << std::setfill('0') << std::setw(1) << std::hex
                 << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0,
                            255.0)); // green
              ss << "00"; // blue
            } else {
              ss << "FF", // red
                ss << std::setfill('0') << std::setw(1) << std::hex
                   << int(map(connection.weight, 0.0, 0.5, 0.0,
                              255.0 / 2.0)); // green
              ss << "00"; // blue
            }
            glm::dvec2 p1 =
              (xsec.points[connection.a_idx].point - center) * input.thickness;
            glm::dvec2 p2 =
              (xsec.points[connection.b_idx].point - center) * input.thickness;
            SVG::line(os, p1.x, p1.y, p2.x, p2.y, 0.1, ss.str());
          }
        }
      }
    }
  });
}

void Parameterization::orthogonal_isolines_svg(std::ostream &os,
                                               int new_stroke_ind,
                                               const Input &input, size_t itr) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);
  input.cluster_svg(os, new_stroke_ind, [&](std::ostream &os) {
    for (auto &kv : input.clusters) {
      if (itr >= kv.second.orthogonal_xsecs_list.size())
        continue;
      for (auto &xsec : kv.second.orthogonal_xsecs_list[itr]) {
        for (auto &connection : xsec.connections) {
          if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
            std::stringstream ss;
            ss << "#";
            if (connection.weight > 0.5) {
              ss << std::setfill('0') << std::setw(2) << std::hex
                 << int(map(connection.weight, 0.5, 1.0, 255.0, 0.0)); // red
              ss << std::setfill('0') << std::setw(1) << std::hex
                 << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0,
                            255.0)); // green
              ss << "00"; // blue
            } else {
              ss << "FF", // red
                ss << std::setfill('0') << std::setw(1) << std::hex
                   << int(map(connection.weight, 0.0, 0.5, 0.0,
                              255.0 / 2.0)); // green
              ss << "00"; // blue
            }
            glm::dvec2 p1 =
              (xsec.points[connection.a_idx].point - center) * input.thickness;
            glm::dvec2 p2 =
              (xsec.points[connection.b_idx].point - center) * input.thickness;
            SVG::line(os, p1.x, p1.y, p2.x, p2.y, 0.1, ss.str());
          }
        }
      }
    }
  });
}

void Parameterization::isolines_cluster_svg(std::ostream &os,
                                            int new_stroke_ind,
                                            const Input &input,
                                            int cluster_num) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);

  input.cluster_svg(os, new_stroke_ind, [&](std::ostream &os) {
    auto &kv = input.clusters.at(cluster_num);
    {
      for (auto &xsec : kv.xsecs) {
        for (auto &connection : xsec.connections) {
          if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
            std::stringstream ss;
            ss << "#";
            if (connection.weight > 0.5) {
              ss << std::setfill('0') << std::setw(2) << std::hex
                 << int(map(connection.weight, 0.5, 1.0, 255.0, 0.0)); // red
              ss << std::setfill('0') << std::setw(1) << std::hex
                 << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0,
                            255.0)); // green
              ss << "00"; // blue
            } else {
              ss << "FF", // red
                ss << std::setfill('0') << std::setw(1) << std::hex
                   << int(map(connection.weight, 0.0, 0.5, 0.0,
                              255.0 / 2.0)); // green
              ss << "00"; // blue
            }
            glm::dvec2 p1 =
              (xsec.points[connection.a_idx].point - center) * input.thickness;
            glm::dvec2 p2 =
              (xsec.points[connection.b_idx].point - center) * input.thickness;
            SVG::line(os, p1.x, p1.y, p2.x, p2.y, 0.1, ss.str());
          }
        }
      }
    }
  });
}

void Parameterization::save_svg(std::ostream &os, const Input &input) {
  input.input_svg(os, [&](std::ostream &os) {});
}

void Parameterization::save_non_isoline_svg(std::ostream &os,
                                            const Input &input) {
  input.non_isoline_svg(os, [&](std::ostream &os) {});
}

void Parameterization::debug_svg(std::ostream &os, int new_stroke_ind,
                                 const Input &input) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);
  input.cluster_svg(os, new_stroke_ind, [&](std::ostream &os) {
    for (auto &line : debug_lines) {
      line.from -= center;
      line.to -= center;
      line.from *= input.thickness;
      line.to *= input.thickness;
      SVG::line(os, line.from.x, line.from.y, line.to.x, line.to.y, 1.0,
                line.color);
    }
  });
}

void Parameterization::isolines_only_cluster_svg(
  std::ostream &os, int new_stroke_ind, const Input &input, int cluster_num,
  std::map<int, int> stroke_to_cluster) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);

  input.only_cluster_svg(
    os, new_stroke_ind, cluster_num, stroke_to_cluster, [&](std::ostream &os) {
      auto &kv = input.clusters.at(cluster_num);
      {
        for (auto &xsec : kv.xsecs) {
          for (auto &connection : xsec.connections) {
            if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
              std::stringstream ss;
              ss << "#";
              if (connection.weight > 0.5) {
                ss << std::setfill('0') << std::setw(2) << std::hex
                   << int(map(connection.weight, 0.5, 1.0, 255.0,
                              0.0)); // red
                ss << std::setfill('0') << std::setw(1) << std::hex
                   << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0,
                              255.0)); // green
                ss << "00"; // blue
              } else {
                ss << "FF", // red
                  ss << std::setfill('0') << std::setw(1) << std::hex
                     << int(map(connection.weight, 0.0, 0.5, 0.0,
                                255.0 / 2.0)); // green
                ss << "00"; // blue
              }
              glm::dvec2 p1 = (xsec.points[connection.a_idx].point - center) *
                              input.thickness;
              glm::dvec2 p2 = (xsec.points[connection.b_idx].point - center) *
                              input.thickness;
              SVG::line(os, p1.x, p1.y, p2.x, p2.y, 0.1, ss.str());
            }
          }
        }
      }
    });
}

void Parameterization::isolines_only_divided_cluster_svg(
  std::ostream &os, int new_stroke_ind, const Input &input, int cluster_num,
  std::vector<int> base_stroke_ind, std::vector<int> target_stroke_ind) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);

  input.only_divided_cluster_svg(
    os, new_stroke_ind, cluster_num, base_stroke_ind, target_stroke_ind,
    [&](std::ostream &os) {
      auto &kv = input.clusters.at(cluster_num);
      {
        for (auto &xsec : kv.xsecs) {
          for (auto &connection : xsec.connections) {
            if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
              std::stringstream ss;
              ss << "#";
              if (connection.weight > 0.5) {
                ss << std::setfill('0') << std::setw(2) << std::hex
                   << int(map(connection.weight, 0.5, 1.0, 255.0,
                              0.0)); // red
                ss << std::setfill('0') << std::setw(1) << std::hex
                   << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0,
                              255.0)); // green
                ss << "00"; // blue
              } else {
                ss << "FF", // red
                  ss << std::setfill('0') << std::setw(1) << std::hex
                     << int(map(connection.weight, 0.0, 0.5, 0.0,
                                255.0 / 2.0)); // green
                ss << "00"; // blue
              }
              glm::dvec2 p1 = (xsec.points[connection.a_idx].point - center) *
                              input.thickness;
              glm::dvec2 p2 = (xsec.points[connection.b_idx].point - center) *
                              input.thickness;
              SVG::line(os, p1.x, p1.y, p2.x, p2.y, 0.1, ss.str());
            }
          }
        }
      }
    });
}

void Parameterization::debug_svg(std::ostream &os, const Input &input) {
  double w, h, x, y;
  glm::dvec2 center;
  center_input(input, w, h, x, y, center);
  input.cluster_svg(os, -1, [&](std::ostream &os) {
    for (auto &line : debug_lines) {
      line.from -= center;
      line.to -= center;
      line.from *= input.thickness;
      line.to *= input.thickness;
      SVG::line(os, line.from.x, line.from.y, line.to.x, line.to.y, 1.0,
                line.color);
    }
  });
}

void Parameterization::add_debug_line(Parameterization::DebugLine line) {
  std::lock_guard<std::mutex> lock(viz_lock);
  debug_lines.push_back(line);
}

std::map<std::string, double>
Parameterization::parameterize_cluster(Cluster *cluster) {
  std::map<std::string, double> obj_term_values;
  obj_term_values["velocity"] = 0;
  obj_term_values["alignment"] = 0;
  debug_lines.clear();

  // Clear the saved u values
  for (auto &stroke : cluster->strokes) {
    stroke.u.assign(stroke.u.size(), 0);
  }

  // use the existing parameterized value
  // if (!(cluster->xsecs.empty())) {
  //   std::vector<Cluster::XSec> xsecs = orthogonal_xsecs(*cluster);
  //   cluster->xsecs.insert(cluster->xsecs.end(), xsecs.begin(), xsecs.end());
  // } else {
  //   cluster->xsecs = orthogonal_xsecs(*cluster);
  // }

  // Recompute cross sections intead
  cluster->xsecs = orthogonal_xsecs(*cluster);
  cluster->orthogonal_xsecs_list.emplace_back(cluster->xsecs);

  std::vector<std::vector<double>> prev_u;

  auto record_current_u = [&]() -> void {
    prev_u.clear();
    for (auto &stroke : cluster->strokes) {
      prev_u.push_back(stroke.u);
    }
  };

  auto converged = [&]() -> bool {
    bool ok = true;
    for (size_t stroke = 0; ok && stroke < prev_u.size(); ++stroke) {
      for (size_t i = 0; ok && i < prev_u[stroke].size(); ++i) {
        double diff =
          std::abs(prev_u[stroke][i] - cluster->strokes[stroke].u[i]);
        if (i == prev_u[stroke].size() - 1) {
          // 変える場所　でかくする
          ok = diff * diff <
               0.5 * 0.5 *
                 glm::distance2(cluster->strokes[stroke].points[i - 1],
                                cluster->strokes[stroke].points[i]);
        } else {
          ok = diff * diff <
               0.5 * 0.5 *
                 glm::distance2(cluster->strokes[stroke].points[i],
                                cluster->strokes[stroke].points[i + 1]);
        }
      }
    }
    return ok;
  };

  for (int it = 0; it < 5; ++it) {
    record_current_u();
    ensure_connected(cluster);

    auto obj_values =
      (regularize_by_default)
        ? params_from_xsecs_regularized(cluster, true, it == 0, nullptr)
        : params_from_xsecs(cluster, true, it == 0, nullptr);

    if (cluster->param_failed) {
      double final_lambda;
      obj_values =
        params_from_xsecs_retry(cluster, final_lambda, true, it == 0, nullptr);
      std::cout << "Relaxed param retry: " << !cluster->param_failed
                << ", itr: " << it << ", lambda: " << final_lambda << std::endl;
    }

    for (auto k_v : obj_values) {
      obj_term_values[k_v.first] = k_v.second;
    }

    // Solve failed
    if (cluster->param_failed) {
      std::cout << "Relaxed param failed" << std::endl;
      cluster->strokes.clear();
      return obj_term_values;
    }

    if (converged()) {
      break;
    } else {
      // Adjust connection weights
      double total_len = cluster->max_u();
      for (auto &xsec : cluster->xsecs) {
        if (xsec.connector)
          continue;
        auto get_u = [&](size_t idx) -> double {
          auto &pt = xsec.points[idx];
          double mix = pt.i - std::floor(pt.i);
          return (1. - mix) *
                   cluster->strokes[pt.stroke_idx].u[std::floor(pt.i)] +
                 mix * cluster->strokes[pt.stroke_idx].u[std::ceil(pt.i)];
        };
        for (auto &connection : xsec.connections) {
          double diff = get_u(connection.a_idx) - get_u(connection.b_idx);
          connection.weight = gaussian(diff, total_len / 30., 0.);
        }
      }
    }
  }

  for (int it = 0; it < 5; ++it) {
    record_current_u();
    bool ensure_sample_reference = true;
    cluster->xsecs =
      xsecs_from_params(*cluster, ensure_sample_reference,
                        stroke_sampling_size / stroke_sampling_size_ref);
    ensure_connected(cluster);
    double lambda = 0;
    auto obj_values =
      (regularize_by_default)
        ? params_from_xsecs_regularized(cluster, false, false, nullptr)
        : params_from_xsecs(cluster, false, false, nullptr);
    if (cluster->param_failed) {
      obj_values =
        params_from_xsecs_retry(cluster, lambda, false, false, nullptr);
      std::cout << "Final param retry: " << !cluster->param_failed
                << ", itr: " << it << ", lambda: " << lambda << std::endl;
    }
    cluster->orthogonal_xsecs_list.emplace_back(cluster->xsecs);

    for (auto k_v : obj_values) {
      obj_term_values[k_v.first] = k_v.second;
    }

    // Solve failed
    if (cluster->param_failed) {
      cluster->strokes.clear();
      return obj_term_values;
    }

    if (converged() && lambda == 0) {
      break;
    }
  }

  record_current_u();
  auto cluster_saved = *cluster;
  check_periodic(cluster);
  if (cluster->periodic) {
    auto old_u = [&](size_t stroke, double i) {
      double mix = i - std::floor(i);
      return (1. - mix) * prev_u[stroke][std::floor(i)] +
             mix * prev_u[stroke][std::ceil(i)];
    };
    auto is_cut = [&](const Cluster::XSec &xsec) {
      for (auto &conn : xsec.connections) {
        if (std::abs(old_u(xsec.points[conn.a_idx].stroke_idx,
                           xsec.points[conn.a_idx].i) -
                     old_u(xsec.points[conn.b_idx].stroke_idx,
                           xsec.points[conn.b_idx].i)) > 5.) {
          return true;
        }
      }
      return false;
    };
    auto cut_it =
      std::find_if(cluster->xsecs.begin(), cluster->xsecs.end(), is_cut);
    if (cut_it != cluster->xsecs.end()) {
      auto cut = *cut_it;
      ensure_connected(cluster, &cut);
      auto obj_values =
        (regularize_by_default)
          ? params_from_xsecs_regularized(cluster, false, false, &cut)
          : params_from_xsecs(cluster, false, false, &cut);
      if (cluster->param_failed) {
        double final_lambda;
        obj_values =
          params_from_xsecs_retry(cluster, final_lambda, false, false, &cut);
        std::cout << "Periodic param retry: " << !cluster->param_failed
                  << ", lambda: " << final_lambda << std::endl;
      }

      // Solve failed
      if (cluster->param_failed) {
        cluster->strokes.clear();
        return obj_term_values;
      }

      // If the periodic parameterization is bad (and when there are more than
      // one stroke), use the non-periodic one
      if (cluster->strokes.size() > 1 &&
          obj_values["alignment"] >= periodic_alignment_threshold) {
        cluster->periodic = false;
        *cluster = cluster_saved;
        bool ensure_sample_reference = false;
        cluster->xsecs =
          xsecs_from_params(*cluster, ensure_sample_reference,
                            stroke_sampling_size / stroke_sampling_size_ref);
        return obj_term_values;
      }
    }
  }
  bool ensure_sample_reference = false;
  cluster->xsecs =
    xsecs_from_params(*cluster, ensure_sample_reference,
                      stroke_sampling_size / stroke_sampling_size_ref);

  // In case the periodic parameterization fails and it's not indicated by the
  // two terms
  if (cluster->periodic) {
    auto max_xsec_width = [](const Cluster &cluster) -> double {
      double max_w = 0;
      for (const auto &xsec : cluster.xsecs)
        max_w = std::max(max_w, xsec.xsec_width());
      return max_w;
    };
    double before_width = max_xsec_width(cluster_saved);
    double after_width = max_xsec_width(*cluster);
    if (std::max(1., after_width) / std::max(1., before_width) > 20 ||
        after_width > 20) {
      bool ensure_sample_reference = false;
      cluster->periodic = false;
      *cluster = cluster_saved;
      cluster->xsecs =
        xsecs_from_params(*cluster, ensure_sample_reference,
                          stroke_sampling_size / stroke_sampling_size_ref);
      return obj_term_values;
    }
  }

  return obj_term_values;
}

std::map<std::string, double>
Parameterization::params_from_xsecs_retry(Cluster *cluster,
                                          double &final_lambda, bool relaxed,
                                          bool initial, Cluster::XSec *cut) {
  double lambda = (initial) ? lambda_init_l1 : lambda_init;
  cluster->param_failed = false;
  std::map<std::string, double> obj_values;
  size_t retry_i = 0;
  do {
    obj_values = params_from_xsecs(cluster, relaxed, initial, cut, lambda);
    if (!cluster->param_failed)
      break;
    lambda *= 10;
    retry_i++;
  } while (retry_i < max_retry);
  final_lambda = lambda;
  return obj_values;
}

std::map<std::string, double>
Parameterization::params_from_xsecs(Cluster *cluster, bool relaxed,
                                    bool initial, Cluster::XSec *cut,
                                    double lambda) {
  std::map<std::string, double> obj_term_values;
  obj_term_values["velocity"] = 0;
  obj_term_values["alignment"] = 0;

  GRBModel model(context.grb);

  std::map<size_t, std::pair<size_t, size_t>> arr2point;
  std::map<std::pair<size_t, size_t>, size_t> point2arr;
  auto to_point_index =
    [&arr2point](const size_t idx_1d) -> std::pair<size_t, size_t> {
    // assert(arr2point.count(idx_1d));
    return arr2point[idx_1d];
  };
  auto to_array_index = [&point2arr](const size_t s_idx,
                                     const size_t p_idx) -> size_t {
    const auto p = std::make_pair(s_idx, p_idx);
    // assert(point2arr.count(p));
    return point2arr[p];
  };

  std::vector<GRBVar> param_vars;
  for (size_t j = 0; j < cluster->strokes.size(); ++j) {
    auto &stroke = cluster->strokes[j];
    for (size_t i = 0; i < stroke.points.size(); ++i) {
      std::string name =
        std::string("u_{") + std::to_string(j) + "," + std::to_string(i) + "}";
      const auto p = std::make_pair(j, i);

      arr2point[param_vars.size()] = p;
      point2arr[p] = param_vars.size();
      param_vars.push_back(
        model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, name));
    }
  }

  auto points_cross_cut = [&](size_t stroke, int a, int b) {
    for (auto &pt : cut->points) {
      // if (pt.stroke_idx == stroke && std::ceil(pt.i) >= a && std::floor(pt.i)
      // <= b) {
      if (pt.stroke_idx == stroke && pt.i >= a && pt.i <= b) {
        return true;
      }
    }
    return false;
  };

  std::vector<LinExpr> velocity_terms;
  std::vector<LinExpr> alignment_terms;
  std::vector<double> velocity_scales, alignment_scales;

  velocity_terms.reserve(cluster->xsecs.size());
  size_t alignment_count = 0;
  for (auto &xsec : cluster->xsecs) {
    alignment_count += xsec.connections.size();
  }
  alignment_terms.reserve(alignment_count);
  alignment_scales.reserve(alignment_count);

  for (auto &xsec : cluster->xsecs) {
    glm::dvec2 tangent = xsec.avg_tangent();

    auto connection_crosses_cut = [&](const Cluster::XSecConnection conn) {
      for (size_t idx : {conn.a_idx, conn.b_idx}) {
        size_t stroke = xsec.points[idx].stroke_idx;
        double i = xsec.points[idx].i;

        for (auto &pt : cut->points) {
          if (pt.stroke_idx == stroke && std::ceil(pt.i) >= std::ceil(i) &&
              std::floor(pt.i) <= std::floor(i)) {
            return true;
          }
        }
      }
      return false;
    };

    // 1. Add velocity term
    {
      velocity_terms.emplace_back();
      velocity_scales.emplace_back(0);
      auto &velocity_term = velocity_terms.back();

      double total_weight = 0.0;
      for (size_t i = 0; i < xsec.points.size(); ++i) {
        auto &point = xsec.points[i];
        double weight =
          xsec.distance_weight(i) +
          1; // Regularize to avoid zero weights on very close strokes
        double coefficient = weight * glm::dot(point.tangent, tangent) /
                             glm::length(point.to_next);
        if (std::isnan(coefficient)) {
          throw "coefficient NaN";
        }
        if (std::isinf(coefficient)) {
          // throw "coefficient inf";
          continue;
        }

        int a = std::floor(point.i);
        int b = std::ceil(point.i);
        double mix = point.i - a;
        if (a == b) {
          if (b == cluster->strokes[point.stroke_idx].points.size() - 1) {
            --a;
          } else {
            ++b;
          }
        }

        if (cut && points_cross_cut(point.stroke_idx, a, b))
          continue;

        total_weight += weight;
        /*velocity_term +=
          coefficient * (1 - mix) *
          (param_vars[point.stroke_idx][b] - param_vars[point.stroke_idx][a]);*/

        velocity_term.var_idx_coeff[to_array_index(point.stroke_idx, b)] +=
          coefficient * (1 - mix);
        velocity_term.var_idx_coeff[to_array_index(point.stroke_idx, a)] +=
          -coefficient * (1 - mix);
        velocity_scales.back() += coefficient * (1 - mix);
      }

      if (total_weight > 0) {
        for (auto &var_coef : velocity_term.var_idx_coeff)
          var_coef.second /= total_weight;
        velocity_scales.back() /= total_weight;
        velocity_term.const_coeff -= 1.0;
      } else {
        velocity_terms.pop_back();
        velocity_scales.pop_back();
      }
    }

    // 2. Add alignment terms
    if (!cut || !std::any_of(xsec.connections.begin(), xsec.connections.end(),
                             connection_crosses_cut)) {
      for (auto &connection : xsec.connections) {
        auto &pt_a = xsec.points[connection.a_idx];
        auto &pt_b = xsec.points[connection.b_idx];

        alignment_terms.emplace_back();
        auto &alignment_term = alignment_terms.back();

        double mix_a = pt_a.i - std::floor(pt_a.i);
        double mix_b = pt_b.i - std::floor(pt_b.i);
        double proj_dist = glm::dot(pt_a.point - pt_b.point, tangent);

        alignment_term
          .var_idx_coeff[to_array_index(pt_a.stroke_idx, std::floor(pt_a.i))] +=
          connection.weight / (xsec.connections.size()) * (1 - mix_a);
        alignment_term
          .var_idx_coeff[to_array_index(pt_a.stroke_idx, std::ceil(pt_a.i))] +=
          connection.weight / (xsec.connections.size()) * mix_a;
        alignment_term
          .var_idx_coeff[to_array_index(pt_b.stroke_idx, std::floor(pt_b.i))] +=
          -connection.weight / (xsec.connections.size()) * (1 - mix_b);
        alignment_term
          .var_idx_coeff[to_array_index(pt_b.stroke_idx, std::ceil(pt_b.i))] +=
          -connection.weight / (xsec.connections.size()) * mix_b;
        if (!relaxed)
          alignment_term.const_coeff +=
            -connection.weight / (xsec.connections.size()) * proj_dist;

        alignment_scales.emplace_back(connection.weight /
                                      (xsec.connections.size()));
      }
    }
  }

  // 3. Enforce monotonicity
  for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
       ++stroke_idx) {
    auto &stroke = cluster->strokes[stroke_idx];
    for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
      if (cut && points_cross_cut(stroke_idx, i, i + 1))
        continue;

      model.addConstr(param_vars[to_array_index(stroke_idx, i + 1)] -
                        param_vars[to_array_index(stroke_idx, i)] >=
                      0.5 *
                        glm::distance(stroke.points[i + 1], stroke.points[i]));
    }
  }

  // model.write("/Users/chenxil/projects/sketch_clustering/sketch_clustering/"
  //             "models/model_mono_const_" +
  //             std::to_string(model_number) + ".lp");

  // Find the stroke with the first point sample u == 0
  size_t s_min_u = 0;
  if (lambda > 0) {
    double min_u = std::numeric_limits<double>::infinity();
    for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
         ++stroke_idx) {
      auto &stroke = cluster->strokes[stroke_idx];
      assert(!stroke.u.empty());
      assert(stroke.u.front() <= stroke.u.back());
      if (stroke.u.front() < min_u) {
        min_u = stroke.u.front();
        s_min_u = stroke_idx;
      }
    }
  }

  // 4. Boundary
  if (cut) {
    for (auto &pt : cut->points) {
      double dist_from_cut =
        glm::distance(pt.point, point(cluster->strokes[pt.stroke_idx].points,
                                      std::ceil(pt.i)));
      model.addConstr(
        param_vars[to_array_index(pt.stroke_idx, std::ceil(pt.i))] ==
        dist_from_cut);
    }
  } else {
    if (lambda == 0)
      model.addConstr(param_vars[to_array_index(0, 0)] == 0.0);
    else // If regularization exists, need to pin to a specific u value
      model.addConstr(param_vars[to_array_index(s_min_u, 0)] ==
                      cluster->strokes[s_min_u].u.front());
  }

  // Expand quadratic terms
  // for (auto &t : velocity_terms) {
  //   t.expand();
  // }
  // for (auto &t : alignment_terms) {
  //   t.expand();
  // }
  // TBB:
  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<size_t>(0u, velocity_terms.size()),
    [&](const oneapi::tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        velocity_terms[i].expand();
      }
    });
  oneapi::tbb::parallel_for(
    oneapi::tbb::blocked_range<size_t>(0u, alignment_terms.size()),
    [&](const oneapi::tbb::blocked_range<size_t> &range) {
      for (size_t i = range.begin(); i != range.end(); ++i) {
        alignment_terms[i].expand();
      }
    });

  GRBQuadExpr objective;
  if (initial) {
    // double scale = 1e-5;
    double scale = 1e-3;
    if (alignment_terms.size() > 6 * velocity_terms.size()) {
      // scale = 1e-6;
      scale = 1e-4;
    }
    objective = l2_norm_sq(&model, velocity_terms, param_vars) +
                scale * l1_norm(&model, alignment_terms, param_vars);
  } else {
    objective = l2_norm_sq(&model, velocity_terms, param_vars) +
                l2_norm_sq(&model, alignment_terms, param_vars);
  }

  // Regularization
  if (lambda > 0) {
    for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
         ++stroke_idx) {
      auto &stroke = cluster->strokes[stroke_idx];
      for (size_t i = 0; i < stroke.points.size(); ++i) {
        // lambda * (u_i - u_i0) * (u_i - u_i0)
        objective += lambda * (param_vars[to_array_index(stroke_idx, i)] *
                               param_vars[to_array_index(stroke_idx, i)]) -
                     (lambda * 2 * stroke.u[i]) *
                       param_vars[to_array_index(stroke_idx, i)] +
                     lambda * stroke.u[i] * stroke.u[i];
      }
    }
  }

  try {
    model.setObjective(objective, GRB_MINIMIZE);
  } catch (GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    model.write("model.lp");

    throw e;
  }

  // {
  //   std::string model_name =
  //     "/Users/chenxil/projects/sketch_clustering/sketch_clustering/"
  //     "models/model_" +
  //     std::to_string(model_number) + ".lp";
  //   model.write(model_name);
  //   model_number++;
  //   std::cout << "Saved: " << model_name << std::endl;
  // }

  context.optimize_model(&model);

  // model.write("/Users/chenxil/projects/sketch_clustering/sketch_clustering/"
  //             "models/model_" +
  //             std::to_string(model_number) + ".sol");
  // model_number++;

  double min_u = std::numeric_limits<double>::infinity();
  for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
       ++stroke_idx) {
    auto &stroke = cluster->strokes[stroke_idx];
    for (size_t i = 0; i < stroke.points.size(); ++i) {
      try {
        stroke.u[i] =
          param_vars[to_array_index(stroke_idx, i)].get(GRB_DoubleAttr_X);
        min_u = std::min(min_u, stroke.u[i]);
      } catch (GRBException e) {
        // std::string model_name =
        //   "/Users/chenxil/projects/sketch_clustering/sketch_clustering/"
        //   "models/model_" +
        //   std::to_string(model_number) + ".lp";
        // model.write(model_name);
        // model_number++;

        // std::cout << "params_from_xsecs: Error code = " << e.getErrorCode()
        //           << "; " << std::to_string(i) << " / "
        //           << std::to_string(stroke.points.size())
        //           << "; Saved: " << model_name << std::endl;
        // std::cout << e.getMessage() << std::endl;

        // throw e;
        // cluster->strokes.clear();
        cluster->param_failed = true;
        return obj_term_values;
      }
    }
  }
  for (auto &stroke : cluster->strokes) {
    for (size_t i = 0; i < stroke.u.size(); ++i) {
      stroke.u[i] -= min_u;
    }
  }

  // Examine the solution
  double obj = model.get(GRB_DoubleAttr_ObjVal);
  std::vector<double> velocity_values, alignment_values;
  for (size_t i = 0; i < velocity_terms.size(); ++i) {
    const auto &v = velocity_terms[i];
    velocity_values.emplace_back(std::abs(v.getValue(param_vars)));
  }
  for (size_t i = 0; i < alignment_terms.size(); ++i) {
    const auto &a = alignment_terms[i];
    alignment_values.emplace_back(
      std::abs(a.getValue(param_vars) / alignment_scales[i]));
  }
  std::sort(velocity_values.begin(), velocity_values.end(), std::greater<>());
  std::sort(alignment_values.begin(), alignment_values.end(), std::greater<>());
  // std::cout << "Obj: " << obj << std::endl;

  if (!velocity_values.empty())
    obj_term_values["velocity"] = velocity_values.front();
  if (!alignment_values.empty()) {
    obj_term_values["alignment"] = alignment_values.front();
    // std::cout << obj_term_values["alignment"] << std::endl;
  }

  return obj_term_values;
}

std::map<std::string, double>
Parameterization::params_from_xsecs_regularized(Cluster *cluster, bool relaxed,
                                                bool initial,
                                                Cluster::XSec *cut) {
  double final_lambda;
  std::map<std::string, double> obj_values =
    params_from_xsecs_retry(cluster, final_lambda, relaxed, initial, cut);
  if (relaxed)
    std::cout << "Relaxed param regularized: " << !cluster->param_failed
              << ", lambda: " << final_lambda << std::endl;
  else
    std::cout << "Param regularized: " << !cluster->param_failed
              << ", lambda: " << final_lambda << std::endl;
  return obj_values;
}

std::vector<Cluster::XSec>
Parameterization::orthogonal_xsecs(const Cluster &cluster,
                                   bool periodic_check) {
  const double TARGET_ANGLE = 20. / 180. * M_PI;
  const double SPIRAL_TARGET_ANGLE = 40. / 180. * M_PI;
  // const double SPIRAL_TARGET_ANGLE = 20. / 180. * M_PI;
  double angle_tolerance;
  double MIN_GAP_CUTOFF;
  double MAX_GAP_CUTOFF;
  if (periodic_check) {
    angle_tolerance = 0.15; // 0.1; // 0.05;
    MIN_GAP_CUTOFF = 15.;
    MAX_GAP_CUTOFF = 30.;
  } else {
    angle_tolerance = 0.0;
    // MIN_GAP_CUTOFF = 5.;
    // MAX_GAP_CUTOFF = 20.;
    MIN_GAP_CUTOFF = 5.;
    MAX_GAP_CUTOFF = 15.;
  }

  std::vector<Cluster::XSec> xsecs;

  if (!(cluster.xsecs.empty())) {
    // for (size_t i = 0;
    //      i < cluster.strokes[cluster.strokes.size() - 1].points.size(); ++i)
    //      {
    //   xsecs.push_back(
    //     orthogonal_xsec_at(cluster, cluster.strokes.size() - 1, i, true));
    // }

    // This is for periodic check: Only find xsecs from the strokes containing
    // the min and max U values.
    auto pair = std::minmax_element(
      cluster.xsecs.begin(), cluster.xsecs.end(),
      [](const Cluster::XSec &a, const Cluster::XSec &b) { return a.u < b.u; });
    const Cluster::XSec &max_xsec = *pair.second;
    const Cluster::XSec &min_xsec = *pair.first;
    assert(!min_xsec.points.empty() && !max_xsec.points.empty());
    std::unordered_set<size_t> minmax_stroke_indices{
      max_xsec.points.front().stroke_ind, min_xsec.points.front().stroke_ind};
    for (size_t stroke = 0; stroke < cluster.strokes.size(); ++stroke) {
      if (minmax_stroke_indices.count(cluster.strokes[stroke].stroke_ind)) {
        for (size_t i = 0; i < cluster.strokes[stroke].points.size(); ++i) {
          xsecs.push_back(orthogonal_xsec_at(cluster, stroke, i, false));
        }
      }
    }
  } else {
    for (size_t stroke = 0; stroke < cluster.strokes.size(); ++stroke) {
      for (size_t i = 0; i < cluster.strokes[stroke].points.size(); ++i) {
        xsecs.push_back(orthogonal_xsec_at(cluster, stroke, i, false));
      }
    }
  }

  // Local filtering
  std::map<size_t, const Cluster::Stroke *> sid_s;
  for (size_t stroke = 0; stroke < cluster.strokes.size(); ++stroke) {
    sid_s[cluster.strokes[stroke].stroke_ind] = &cluster.strokes[stroke];
  }

  auto is_spiral_periodic = [&sid_s](const Cluster::XSecPoint &p1,
                                     const Cluster::XSecPoint &p2) -> bool {
    int orig_sid1 =
      (sid_s[p1.stroke_ind]->spiral) ? p1.stroke_ind % spiral_offset : -1;
    int orig_sid2 =
      (sid_s[p2.stroke_ind]->spiral) ? p2.stroke_ind % spiral_offset : -1;
    if (orig_sid1 == orig_sid2 && orig_sid2 >= 0) {
      double angle1 = spiral_angle(sid_s[p1.stroke_ind]->spiral_angles, p1.i);
      double angle2 = spiral_angle(sid_s[p2.stroke_ind]->spiral_angles, p2.i);
      double angle_diff = std::abs(angle1 - angle2);
      if (angle_diff <= 2 * M_PI - periodic_offset_threshold)
        return false;
      double periodic_offset = std::fmod(angle_diff, (2 * M_PI));
      periodic_offset = std::min(periodic_offset, 2 * M_PI - periodic_offset);
      if (periodic_offset <= periodic_offset_threshold)
        return true;
    }
    return false;
  };

  // 1. Angle filtering
  {
    struct AngleSample {
      glm::dvec2 mid;
      double angle;
    };
    std::vector<AngleSample> samples;
    for (auto &xsec : xsecs) {
      for (size_t i = 1; i < xsec.points.size(); ++i) {
        samples.push_back({
          (xsec.points[i - 1].point + xsec.points[i].point) * 0.5,
          std::acos(
            glm::dot(xsec.points[i - 1].tangent, xsec.points[i].tangent)),
        });
      }
    }
    for (auto &xsec : xsecs) {
      std::deque<AngleSample> neighbourhood;
      double r_sq =
        std::pow(std::max(30.0, glm::distance(xsec.points.front().point,
                                              xsec.points.back().point)),
                 2);

      auto sample_within_radius = [&](const AngleSample &sample) {
        for (auto &point : xsec.points) {
          if (glm::distance2(point.point, sample.mid) <= r_sq) {
            return true;
          }
        }
        return false;
      };

      // TODO k-d tree lookup?
      for (auto &sample : samples) {
        if (sample_within_radius(sample)) {
          neighbourhood.push_back(sample);
        }
      }

      double max_angle = TARGET_ANGLE;
      if (!neighbourhood.empty()) {

        // Sort by angle descending
        std::sort(neighbourhood.begin(), neighbourhood.end(),
                  [](const AngleSample &a, const AngleSample &b) {
                    return a.angle > b.angle;
                  });

        // Remove samples until the `angle_tolerence` percentile is below the
        // target value
        if (neighbourhood.back().angle <= TARGET_ANGLE) {
          max_angle = neighbourhood.front().angle;
          while (
            !neighbourhood.empty() &&
            neighbourhood[int(angle_tolerance * (neighbourhood.size() - 1))]
                .angle > TARGET_ANGLE) {
            neighbourhood.pop_front();
            max_angle = neighbourhood.front().angle;
          }
        }
      }

      // Remove connections above the threshold to the right of the center
      for (int i = xsec.center_idx + 1; i < xsec.points.size(); ++i) {
        double angle_diff = std::acos(
          glm::dot(xsec.points[i - 1].tangent, xsec.points[i].tangent));
        if ((angle_diff > max_angle &&
             !is_spiral_periodic(xsec.points[i - 1], xsec.points[i])) ||
            (angle_diff > SPIRAL_TARGET_ANGLE &&
             is_spiral_periodic(xsec.points[i - 1], xsec.points[i]))) {
          xsec.points = std::vector<Cluster::XSecPoint>(
            xsec.points.begin(), xsec.points.begin() + i);
          break;
        }
      }
      // Remove connections above the threshold to the left of the center
      for (int i = xsec.center_idx - 1; i >= 0; --i) {
        double angle_diff = std::acos(
          glm::dot(xsec.points[i + 1].tangent, xsec.points[i].tangent));
        if ((angle_diff > max_angle &&
             !is_spiral_periodic(xsec.points[i + 1], xsec.points[i])) ||
            (angle_diff > SPIRAL_TARGET_ANGLE &&
             is_spiral_periodic(xsec.points[i + 1], xsec.points[i]))) {
          xsec.points = std::vector<Cluster::XSecPoint>(
            xsec.points.begin() + i + 1, xsec.points.end());
          xsec.center_idx -= i + 1;
          break;
        }
      }
    }
  }

  // 2. Gap filtering
  {
    struct GapSample {
      glm::dvec2 left;
      glm::dvec2 right;
      double gap_sq;
    };
    std::vector<GapSample> samples;
    for (auto &xsec : xsecs) {
      for (size_t i = 1; i < xsec.points.size(); ++i) {
        samples.push_back({
          xsec.points[i - 1].point,
          xsec.points[i].point,
          glm::distance2(xsec.points.front().point, xsec.points.back().point),
        });
      }
    }
    for (auto &xsec : xsecs) {
      if (xsec.points.size() == 1)
        continue;

      std::unordered_set<double> neighbourhood;
      double r_sq =
        std::max(200.0 * 200.0, glm::distance2(xsec.points.front().point,
                                               xsec.points.back().point));

      auto sample_within_radius = [&](const GapSample &sample) {
        for (auto &point : xsec.points) {
          for (auto &other : {sample.left, sample.right}) {
            if (glm::distance2(point.point, other) <= r_sq) {
              return true;
            }
          }
        }
        return false;
      };

      // TODO k-d tree lookup?
      for (auto &sample : samples) {
        if (sample_within_radius(sample)) {
          neighbourhood.insert(sample.gap_sq);
        }
      }

      if (neighbourhood.empty())
        continue;

      // Get median gap
      std::vector<double> gaps(neighbourhood.begin(), neighbourhood.end());
      size_t median_offset = gaps.size() * 0.75;
      std::nth_element(gaps.begin(), gaps.begin() + median_offset, gaps.end());
      double median_gap_sq = gaps[median_offset];
      double gap_cutoff = std::min(
        MAX_GAP_CUTOFF * MAX_GAP_CUTOFF,
        std::max(MIN_GAP_CUTOFF * MIN_GAP_CUTOFF, 1.2 * 1.2 * median_gap_sq));

      // Remove connections above the threshold to the right of the center
      for (int i = xsec.center_idx + 1; i < xsec.points.size(); ++i) {
        if (glm::distance2(xsec.points[i - 1].point, xsec.points[i].point) >
              gap_cutoff &&
            !is_spiral_periodic(xsec.points[i - 1], xsec.points[i])) {
          xsec.points = std::vector<Cluster::XSecPoint>(
            xsec.points.begin(), xsec.points.begin() + i);
          break;
        }
      }
      // Remove connections above the threshold to the left of the center
      for (int i = xsec.center_idx - 1; i >= 0; --i) {
        if (glm::distance2(xsec.points[i + 1].point, xsec.points[i].point) >
              gap_cutoff &&
            !is_spiral_periodic(xsec.points[i + 1], xsec.points[i])) {
          xsec.points = std::vector<Cluster::XSecPoint>(
            xsec.points.begin() + i + 1, xsec.points.end());
          xsec.center_idx -= i + 1;
          break;
        }
      }
    }
  }

  // 3. Spiral filtering
  bool has_spiral = false;
  for (auto const &s : cluster.strokes) {
    if (s.spiral) {
      has_spiral = true;
      break;
    }
  }
  if (has_spiral) {
    struct SpiralSample {
      glm::dvec2 left;
      glm::dvec2 right;
      double angle_diff;
    };

    for (auto &xsec : xsecs) {
      if (xsec.points.size() == 1)
        continue;

      std::vector<SpiralSample> samples;

      // Decide which spiral stroke as the reference
      int spiral_ref = -1;
      std::map<size_t, size_t> spiral_occurences;
      for (size_t i = 0; i < xsec.points.size(); ++i) {
        if (sid_s[xsec.points[i].stroke_ind]->spiral) {
          size_t orig_sid = xsec.points[i].stroke_ind % spiral_offset;
          spiral_occurences[orig_sid]++;
        }
      }
      if (!spiral_occurences.empty())
        spiral_ref = spiral_occurences.begin()->first;
      for (const auto &so : spiral_occurences) {
        if (so.second > 1) {
          spiral_ref = so.first;
          break;
        }
      }
      for (size_t i = 1; i < xsec.points.size(); ++i) {
        double angle_diff = 0;
        double orig_left_angle =
          spiral_angle(sid_s[xsec.points[i - 1].stroke_ind]->spiral_angles,
                       xsec.points[i - 1].i);
        double left_angle = -1;
        int left_i = -1;
        for (int j = i - 1; j >= 0; --j) {
          if (sid_s[xsec.points[j].stroke_ind]->spiral) {
            size_t orig_sid = xsec.points[j].stroke_ind % spiral_offset;
            if (spiral_ref != orig_sid)
              continue;
          }
          double angle = spiral_angle(
            sid_s[xsec.points[j].stroke_ind]->spiral_angles, xsec.points[j].i);
          if (angle >= 0) {
            left_angle = angle;
            left_i = j;
            break;
          }
        }

        double orig_right_angle = spiral_angle(
          sid_s[xsec.points[i].stroke_ind]->spiral_angles, xsec.points[i].i);
        double right_angle = -1;
        int right_i = -1;
        for (size_t j = i; j < xsec.points.size(); ++j) {
          if (sid_s[xsec.points[j].stroke_ind]->spiral) {
            size_t orig_sid = xsec.points[j].stroke_ind % spiral_offset;
            if (spiral_ref != orig_sid)
              continue;
          }
          double angle = spiral_angle(
            sid_s[xsec.points[j].stroke_ind]->spiral_angles, xsec.points[j].i);
          if (angle >= 0) {
            right_angle = angle;
            right_i = j;
            break;
          }
        }

        if (orig_left_angle >= 0 && orig_right_angle >= 0)
          angle_diff = std::abs(left_angle - right_angle);
        else if (left_angle >= 0 && right_angle >= 0) {
          double total_dist = glm::distance(xsec.points[left_i].point,
                                            xsec.points[right_i].point);
          double left_dist1 =
            glm::distance(xsec.points[left_i].point, xsec.points[i - 1].point);
          double left_dist2 =
            glm::distance(xsec.points[left_i].point, xsec.points[i].point);
          if (left_dist1 < 0.5 * total_dist && left_dist2 > 0.5 * total_dist) {
            angle_diff = std::abs(left_angle - right_angle);
          }
        }
        samples.push_back({
          xsec.points[i - 1].point,
          xsec.points[i].point,
          angle_diff,
        });
      }

      auto sample_misalign = [&](const SpiralSample &sample) {
        double periodic_offset = std::fmod(sample.angle_diff, (2 * M_PI));
        if (sample.angle_diff <= 2 * M_PI - periodic_offset_threshold)
          return false;
        periodic_offset = std::min(periodic_offset, 2 * M_PI - periodic_offset);
        if (periodic_offset > periodic_offset_threshold)
          return true;
        return false;
      };

      int cut_i = -1;
      for (size_t i = 0; i < samples.size(); ++i) {
        if (sample_misalign(samples[i])) {
          cut_i = i + 1;
          break;
        }
      }

      if (cut_i < 0)
        continue;

      xsec.points = std::vector<Cluster::XSecPoint>(
        xsec.points.begin(), xsec.points.begin() + cut_i);
      xsec.center_idx = xsec.points.size() - 1;
    }
  }

  for (auto &xsec : xsecs) {
    for (size_t i = 0; i < xsec.points.size(); ++i) {
      for (size_t j = i + 1; j < xsec.points.size(); ++j) {
        xsec.connections.push_back({i, j, 1.0});
      }
    }
  }

  return xsecs;
}

Cluster::XSec Parameterization::orthogonal_xsec_at(const Cluster &cluster,
                                                   size_t stroke, double i,
                                                   bool is_fast) {
  const double SCALE = 100.;

  auto &stroke_points = cluster.strokes[stroke].points;

  size_t a = std::floor(i);
  double mix = i - double(a);
  glm::dvec2 origin = point(stroke_points, i);
  glm::dvec2 tan = tangent(stroke_points, i);
  glm::dvec2 ortho = normal(tan);
  glm::dvec2 to_next(0.0, 0.0);
  if (a < stroke_points.size() - 1) {
    to_next = (1 - mix) * (stroke_points[a + 1] - stroke_points[a]);
  } else {
    to_next = stroke_points[a] - stroke_points[a - 1];
  }

  glm::dvec2 pt1 = origin - SCALE * ortho;
  glm::dvec2 pt2 = origin + SCALE * ortho;

  struct PotentialInt {
    Cluster::XSecPoint point;
    double signed_dist;
  };
  std::vector<PotentialInt> potential_ints = {
    {{cluster.strokes[stroke].stroke_ind, cluster.strokes[stroke].cluster_ind,
      stroke, i, origin, to_next, tan},
     0.0}};

  if (is_fast) {
    for (size_t j = 0; j < cluster.strokes.size(); ++j) {
      auto ints = intersections(cluster.strokes[j].points, pt1, pt2);
      double min_distance = std::numeric_limits<double>::infinity();
      for (auto &intersection : ints) {

        // Ignore origin point
        if (stroke == j && std::abs(intersection.i - i) < 1e-1)
          continue;

        glm::dvec2 to_next(0.0, 0.0);
        double other_i = intersection.i;
        size_t floor_other_i = std::floor(intersection.i);
        if (floor_other_i < cluster.strokes[j].points.size() - 1) {
          to_next = (1 - (other_i - double(floor_other_i))) *
                    (cluster.strokes[j].points[floor_other_i + 1] -
                     cluster.strokes[j].points[floor_other_i]);
        } else if (floor_other_i == other_i) {
          to_next = cluster.strokes[j].points[floor_other_i] -
                    cluster.strokes[j].points[floor_other_i - 1];
        }

        double distance = glm::distance(origin, intersection.pt);

        if (potential_ints.size() == 1) {
          potential_ints.push_back(
            {{cluster.strokes[j].stroke_ind, cluster.strokes[j].cluster_ind, j,
              intersection.i, intersection.pt, to_next,
              tangent(cluster.strokes[j].points, intersection.i)},
             glm::dot(intersection.pt - origin, ortho)});
          if (glm::any(glm::isnan(potential_ints.back().point.tangent))) {
            potential_ints.pop_back();
          }
        } else if (min_distance > distance &&
                   !(glm::any(
                     glm::isnan(potential_ints.back().point.tangent)))) {
          potential_ints.pop_back();
          potential_ints.push_back(
            {{cluster.strokes[j].stroke_ind, cluster.strokes[j].cluster_ind, j,
              intersection.i, intersection.pt, to_next,
              tangent(cluster.strokes[j].points, intersection.i)},
             glm::dot(intersection.pt - origin, ortho)});
          min_distance = distance;
        }
      }
    }
  } else {
    for (size_t j = 0; j < cluster.strokes.size(); ++j) {
      auto ints = intersections(cluster.strokes[j].points, pt1, pt2);
      for (auto &intersection : ints) {

        // Ignore origin point
        if (stroke == j && std::abs(intersection.i - i) < 1e-1)
          continue;

        glm::dvec2 to_next(0.0, 0.0);
        double other_i = intersection.i;
        size_t floor_other_i = std::floor(intersection.i);
        if (floor_other_i < cluster.strokes[j].points.size() - 1) {
          to_next = (1 - (other_i - double(floor_other_i))) *
                    (cluster.strokes[j].points[floor_other_i + 1] -
                     cluster.strokes[j].points[floor_other_i]);
        } else if (floor_other_i == other_i) {
          to_next = cluster.strokes[j].points[floor_other_i] -
                    cluster.strokes[j].points[floor_other_i - 1];
        }

        potential_ints.push_back(
          {{cluster.strokes[j].stroke_ind, cluster.strokes[j].cluster_ind, j,
            intersection.i, intersection.pt, to_next,
            tangent(cluster.strokes[j].points, intersection.i)},
           glm::dot(intersection.pt - origin, ortho)});
        if (glm::any(glm::isnan(potential_ints.back().point.tangent))) {
          potential_ints.pop_back();
        }
      }
    }
  }

  std::sort(potential_ints.begin(), potential_ints.end(),
            [](const PotentialInt &a, const PotentialInt &b) {
              return a.signed_dist < b.signed_dist;
            });

  if (potential_ints.size() > 1) {
    // Cut off intersections after any bad tangent angle
    int first_bad_tangent = -1;
    for (int n = potential_ints.size() - 1; n >= 0; --n) {
      if (potential_ints[n].signed_dist < 0 &&
          glm::dot(potential_ints[n].point.tangent, tan) < 0) {
        first_bad_tangent = n;
        break;
      }
    }
    int last_bad_tangent = potential_ints.size();
    for (int n = 0; n < potential_ints.size(); ++n) {
      if (potential_ints[n].signed_dist > 0 &&
          glm::dot(potential_ints[n].point.tangent, tan) < 0) {
        last_bad_tangent = n;
        break;
      }
    }

    potential_ints = std::vector<PotentialInt>(
      potential_ints.begin() + (1 + first_bad_tangent),
      potential_ints.begin() + last_bad_tangent);
  }

  Cluster::XSec xsec = {{}, {}, 0, 0.0, false};
  for (auto &potential_int : potential_ints) {
    xsec.points.push_back(potential_int.point);
    if (xsec.points.back().stroke_idx == stroke && xsec.points.back().i == i) {
      xsec.center_idx = xsec.points.size() - 1;
    }
  }

  return xsec;
}

std::vector<Cluster::XSec>
Parameterization::xsecs_from_params(const Cluster &cluster,
                                    bool ensure_sample_reference,
                                    double sample_rate, bool nonlinear) {
  double max_u = cluster.max_u();
  sample_rate = std::max(sample_rate, cluster.max_u() / 1e4);
  size_t num_xsecs = std::ceil(max_u / sample_rate);
  if (nonlinear)
    num_xsecs *= 3;

  std::vector<Cluster::XSec> result;
  std::unordered_set<double> seen_u;
  result.reserve(num_xsecs);

  std::vector<std::vector<bool>> sampled;
  std::vector<std::map<std::pair<size_t, size_t>, bool>> binary_sampled;
  for (auto &stroke : cluster.strokes) {
    sampled.emplace_back(stroke.points.size(), false);
    binary_sampled.emplace_back();
    for (size_t i = 0; i + 1 < stroke.points.size(); ++i) {
      binary_sampled.back()[std::make_pair(i, i + 1)] = false;
    }
  }

  for (size_t n = 0; n < num_xsecs; ++n) {
    double t = double(n) / double(num_xsecs - 1);
    double u;
    if (nonlinear) {
      u = poly_in_out(t, 2.25) * max_u;
    } else {
      u = t * max_u;
    }
    auto xsec = xsec_at_u(cluster, u);
    if (xsec.points.size() > 0) {
      result.push_back(xsec);

      for (auto &pt : xsec.points) {
        sampled[pt.stroke_idx][std::floor(pt.i)] = true;
        sampled[pt.stroke_idx][std::ceil(pt.i)] = true;

        int a = std::floor(pt.i);
        int b = std::ceil(pt.i);
        if (a == b) {
          if (b == sampled[pt.stroke_idx].size() - 1) {
            --a;
          } else {
            ++b;
          }
        }

        assert(a + 1 == b);
        auto binary_p = std::make_pair(a, b);
        assert(binary_sampled[pt.stroke_idx].count(binary_p));
        binary_sampled[pt.stroke_idx][binary_p] = true;
      }
    }
  }

  // for (const auto &xsec : result) {
  //   seen_u.emplace(xsec.u);
  // }

  // Make sure all sample points have a cross section referencing them
  for (size_t stroke_idx = 0;
       ensure_sample_reference && stroke_idx < cluster.strokes.size();
       ++stroke_idx) {
    for (size_t i = 0; i < cluster.strokes[stroke_idx].u.size(); ++i) {
      if (!sampled[stroke_idx][i]) {
        auto xsec = xsec_at_u(cluster, cluster.strokes[stroke_idx].u[i]);
        if (xsec.points.empty()) {
          glm::dvec2 to_next(0.0, 0.0);
          if (i < cluster.strokes[stroke_idx].points.size() - 1) {
            to_next = cluster.strokes[stroke_idx].points[i + 1] -
                      cluster.strokes[stroke_idx].points[i];
          } else {
            to_next = cluster.strokes[stroke_idx].points[i] -
                      cluster.strokes[stroke_idx].points[i - 1];
          }
          xsec = Cluster::XSec{
            {{cluster.strokes[stroke_idx].stroke_ind,
              cluster.strokes[stroke_idx].cluster_ind, stroke_idx, double(i),
              cluster.strokes[stroke_idx].points[i], to_next,
              tangent(cluster.strokes[stroke_idx].points, double(i))}},
            {},
            0,
            cluster.strokes[stroke_idx].u[i],
            false};
        }
        // if (!seen_u.count(xsec.u))
        result.push_back(xsec);
      }
      // Need to add a xsec between i and i + 1
      if (ensure_sample_reference &&
          i + 1 < cluster.strokes[stroke_idx].u.size() &&
          !binary_sampled[stroke_idx][std::make_pair(i, i + 1)]) {
        size_t itr = 0;
        double max_itr = 100;
        bool seen_s_xsec_pt = false;
        double itr_u = cluster.strokes[stroke_idx].u[i];
        Cluster::XSec xsec;
        while (itr < max_itr) {
          xsec = xsec_at_u(cluster, itr_u);

          for (auto &pt : xsec.points) {
            if (pt.stroke_ind == cluster.strokes[stroke_idx].stroke_ind &&
                pt.i >= i && pt.i < (i + 1)) {
              seen_s_xsec_pt = true;
            }
          }
          if (xsec.points.empty() || seen_s_xsec_pt)
            break;
          itr++;
          itr_u = itr / max_itr *
                    (cluster.strokes[stroke_idx].u[i + 1] -
                     cluster.strokes[stroke_idx].u[i]) +
                  cluster.strokes[stroke_idx].u[i];
        }
        assert(xsec.points.empty() || seen_s_xsec_pt);

        if (xsec.points.empty() || seen_s_xsec_pt) {
          glm::dvec2 to_next(0.0, 0.0);
          if (i < cluster.strokes[stroke_idx].points.size() - 1) {
            to_next = cluster.strokes[stroke_idx].points[i + 1] -
                      cluster.strokes[stroke_idx].points[i];
          } else {
            to_next = cluster.strokes[stroke_idx].points[i] -
                      cluster.strokes[stroke_idx].points[i - 1];
          }
          xsec = Cluster::XSec{
            {{cluster.strokes[stroke_idx].stroke_ind,
              cluster.strokes[stroke_idx].cluster_ind, stroke_idx, double(i),
              cluster.strokes[stroke_idx].points[i], to_next,
              tangent(cluster.strokes[stroke_idx].points, double(i))}},
            {},
            0,
            cluster.strokes[stroke_idx].u[i],
            false};
        }
        result.push_back(xsec);
      }
    }
  }

  // Sort based on U values
  std::sort(result.begin(), result.end(),
            [](const Cluster::XSec &xsec1, const Cluster::XSec &xsec2) -> bool {
              return xsec1.u < xsec2.u;
            });

  return result;
}

Cluster::XSec Parameterization::xsec_at_u(const Cluster &cluster, double u) {
  Cluster::XSec xsec = {{}, {}, 0, u, false};
  double period = cluster.max_u() + 1e-1;
  for (size_t stroke_idx = 0; stroke_idx < cluster.strokes.size();
       ++stroke_idx) {
    auto &stroke = cluster.strokes[stroke_idx];
    auto from = stroke.u.begin();

    // Periodic parameterizations are no longer monotinic; instead, they are
    // composed of monotonic chunks. We want to iterate over each chunk
    // individually.
    auto find_to = [](std::vector<double>::const_iterator it,
                      std::vector<double>::const_iterator end) {
      do {
        ++it;
      } while (it != end && *it > *(it - 1));
      return it;
    };

    for (auto to = from; from != stroke.u.end(); from = to) {
      to = find_to(from, stroke.u.end());

      auto it_gt = std::lower_bound(from, to, u);
      auto it_le = it_gt;
      if (it_gt == from) {
        // This will only be a problem if we're at the beginning of a stroke.
        // Periodic strokes will be special-cased later.
        if (it_gt == stroke.u.begin()) {
          continue;
        } else {
          --it_le;
        }
      } else if (it_gt != to && *it_gt - u < 1e-6) {
        ++it_gt;
      } else {
        --it_le;
      }
      if (it_gt == to) {
        // This is only a problem if we're at the end of a stroke. Periodic
        // strokes will be special-cased later.
        if (to == stroke.u.end()) {
          if (*it_le == u) {
            it_gt = it_le;
          } else {
            continue;
          }
        }
      }

      size_t a = (it_le - stroke.u.begin());
      size_t b = (it_gt - stroke.u.begin());

      double u_le = *it_le;
      double u_gt = *it_gt;
      if (u_gt < u_le) {
        if (u_le > u) {
          u_le -= period;
        } else {
          u_gt += period;
        }
      }

      double mix = 0;
      if (u_gt != u_le) {
        mix = (u - u_le) / (u_gt - u_le);
      }
      if (mix < 0 || mix > 1) {
        continue;
      }
      double i = a + mix;
      glm::dvec2 point = (1. - mix) * stroke.points[a] + mix * stroke.points[b];

      glm::dvec2 to_next(0.0, 0.0);
      if (it_gt != it_le) {
        to_next = stroke.points[b] - point;
      } else {
        to_next = stroke.points[b] - stroke.points[b - 1];
      }

      xsec.points.push_back({cluster.strokes[stroke_idx].stroke_ind,
                             cluster.strokes[stroke_idx].cluster_ind,
                             stroke_idx, i, point, to_next,
                             tangent(stroke.points, i)});
    }
  }

  for (size_t a = 0; a < xsec.points.size(); ++a) {
    for (size_t b = a + 1; b < xsec.points.size(); ++b) {
      xsec.connections.push_back({a, b, 1.});
    }
  }

  return xsec;
}

std::vector<Cluster::XSec>
Parameterization::xsecs_from_single_stroke_params(const Cluster &cluster,
                                                  bool ensure_sample_reference,
                                                  double sample_rate,
                                                  bool nonlinear) {
  double max_u = cluster.max_u();
  sample_rate = std::max(sample_rate, cluster.max_u() / 1e4);
  size_t num_xsecs = std::ceil(max_u / sample_rate);
  if (nonlinear)
    num_xsecs *= 3;

  std::vector<Cluster::XSec> result;
  std::unordered_set<double> seen_u;
  result.reserve(num_xsecs);

  std::vector<std::vector<bool>> sampled;
  std::vector<std::map<std::pair<size_t, size_t>, bool>> binary_sampled;
  for (auto &stroke : cluster.strokes) {
    sampled.emplace_back(stroke.points.size(), false);
    binary_sampled.emplace_back();
    for (size_t i = 0; i + 1 < stroke.points.size(); ++i) {
      binary_sampled.back()[std::make_pair(i, i + 1)] = false;
    }
  }

  for (size_t n = 0; n < num_xsecs; ++n) {
    double t = double(n) / double(num_xsecs - 1);
    double u;
    if (nonlinear) {
      u = poly_in_out(t, 2.25) * max_u;
    } else {
      u = t * max_u;
    }
    auto xsec = xsec_at_single_stroke_u(cluster, u);
    if (xsec.points.size() > 0) {
      result.push_back(xsec);

      for (auto &pt : xsec.points) {
        sampled[pt.stroke_idx][std::floor(pt.i)] = true;
        sampled[pt.stroke_idx][std::ceil(pt.i)] = true;

        int a = std::floor(pt.i);
        int b = std::ceil(pt.i);
        if (a == b) {
          if (b == sampled[pt.stroke_idx].size() - 1) {
            --a;
          } else {
            ++b;
          }
        }

        assert(a + 1 == b);
        auto binary_p = std::make_pair(a, b);
        assert(binary_sampled[pt.stroke_idx].count(binary_p));
        binary_sampled[pt.stroke_idx][binary_p] = true;
      }
    }
  }

  // for (const auto &xsec : result) {
  //   seen_u.emplace(xsec.u);
  // }

  // Make sure all sample points have a cross section referencing them
  for (size_t stroke_idx = 0;
       ensure_sample_reference && stroke_idx < cluster.strokes.size();
       ++stroke_idx) {
    for (size_t i = 0; i < cluster.strokes[stroke_idx].u.size(); ++i) {
      if (!sampled[stroke_idx][i]) {
        auto xsec =
          xsec_at_single_stroke_u(cluster, cluster.strokes[stroke_idx].u[i]);
        if (xsec.points.empty()) {
          glm::dvec2 to_next(0.0, 0.0);
          if (i < cluster.strokes[stroke_idx].points.size() - 1) {
            to_next = cluster.strokes[stroke_idx].points[i + 1] -
                      cluster.strokes[stroke_idx].points[i];
          } else {
            to_next = cluster.strokes[stroke_idx].points[i] -
                      cluster.strokes[stroke_idx].points[i - 1];
          }
          xsec = Cluster::XSec{
            {{cluster.strokes[stroke_idx].stroke_ind,
              cluster.strokes[stroke_idx].cluster_ind, stroke_idx, double(i),
              cluster.strokes[stroke_idx].points[i], to_next,
              tangent(cluster.strokes[stroke_idx].points, double(i))}},
            {},
            0,
            cluster.strokes[stroke_idx].u[i],
            false};
        }
        // if (!seen_u.count(xsec.u))
        result.push_back(xsec);
      }
      // Need to add a xsec between i and i + 1
      if (ensure_sample_reference &&
          i + 1 < cluster.strokes[stroke_idx].u.size() &&
          !binary_sampled[stroke_idx][std::make_pair(i, i + 1)]) {
        size_t itr = 0;
        double max_itr = 100;
        bool seen_s_xsec_pt = false;
        double itr_u = cluster.strokes[stroke_idx].u[i];
        Cluster::XSec xsec;
        while (itr < max_itr) {
          xsec = xsec_at_single_stroke_u(cluster, itr_u);

          for (auto &pt : xsec.points) {
            if (pt.stroke_ind == cluster.strokes[stroke_idx].stroke_ind &&
                pt.i >= i && pt.i < (i + 1)) {
              seen_s_xsec_pt = true;
            }
          }
          if (xsec.points.empty() || seen_s_xsec_pt)
            break;
          itr++;
          itr_u = itr / max_itr *
                    (cluster.strokes[stroke_idx].u[i + 1] -
                     cluster.strokes[stroke_idx].u[i]) +
                  cluster.strokes[stroke_idx].u[i];
        }
        assert(xsec.points.empty() || seen_s_xsec_pt);

        if (xsec.points.empty() || seen_s_xsec_pt) {
          glm::dvec2 to_next(0.0, 0.0);
          if (i < cluster.strokes[stroke_idx].points.size() - 1) {
            to_next = cluster.strokes[stroke_idx].points[i + 1] -
                      cluster.strokes[stroke_idx].points[i];
          } else {
            to_next = cluster.strokes[stroke_idx].points[i] -
                      cluster.strokes[stroke_idx].points[i - 1];
          }
          xsec = Cluster::XSec{
            {{cluster.strokes[stroke_idx].stroke_ind,
              cluster.strokes[stroke_idx].cluster_ind, stroke_idx, double(i),
              cluster.strokes[stroke_idx].points[i], to_next,
              tangent(cluster.strokes[stroke_idx].points, double(i))}},
            {},
            0,
            cluster.strokes[stroke_idx].u[i],
            false};
        }
        result.push_back(xsec);
      }
    }
  }

  // Sort based on U values
  std::sort(result.begin(), result.end(),
            [](const Cluster::XSec &xsec1, const Cluster::XSec &xsec2) -> bool {
              return xsec1.u < xsec2.u;
            });

  return result;
}

Cluster::XSec Parameterization::xsec_at_single_stroke_u(const Cluster &cluster,
                                                        double u) {
  Cluster::XSec xsec = {{}, {}, 0, u, false};
  double period = cluster.max_u() + 1e-1;
  for (size_t stroke_idx = 0; stroke_idx < cluster.strokes.size();
       ++stroke_idx) {
    auto &stroke = cluster.strokes[stroke_idx];
    auto from = stroke.u.begin();

    // Periodic parameterizations are no longer monotinic; instead, they are
    // composed of monotonic chunks. We want to iterate over each chunk
    // individually.
    auto find_to = [](std::vector<double>::const_iterator it,
                      std::vector<double>::const_iterator end) {
      do {
        ++it;
      } while (it != end && *it > *(it - 1));
      return it;
    };

    for (auto to = from; from != stroke.u.end(); from = to) {
      to = find_to(from, stroke.u.end());

      auto it_gt = std::lower_bound(from, to, u);
      auto it_le = it_gt;
      if (it_gt == from) {
        // This will only be a problem if we're at the beginning of a stroke.
        // Periodic strokes will be special-cased later.
        if (it_gt == stroke.u.begin()) {
          continue;
        } else {
          --it_le;
        }
      } else if (it_gt != to && *it_gt - u < 1e-6) {
        ++it_gt;
      } else {
        --it_le;
      }
      if (it_gt == to) {
        // This is only a problem if we're at the end of a stroke. Periodic
        // strokes will be special-cased later.
        if (to == stroke.u.end()) {
          if (*it_le == u) {
            it_gt = it_le;
          } else {
            continue;
          }
        }
      }

      size_t a = (it_le - stroke.u.begin());
      size_t b = (it_gt - stroke.u.begin());

      double u_le = *it_le;
      double u_gt = *it_gt;
      if (u_gt < u_le) {
        if (u_le > u) {
          u_le -= period;
        } else {
          u_gt += period;
        }
      }

      double mix = 0;
      if (u_gt != u_le) {
        mix = (u - u_le) / (u_gt - u_le);
      }
      if (mix < 0 || mix > 1) {
        continue;
      }
      double i = a + mix;
      glm::dvec2 point = (1. - mix) * stroke.points[a] + mix * stroke.points[b];

      glm::dvec2 to_next(0.0, 0.0);
      if (it_gt != it_le) {
        to_next = stroke.points[b] - point;
      } else {
        to_next = stroke.points[b] - stroke.points[b - 1];
      }

      double ref_u = average_u(stroke.u, i);
      if (std::abs(ref_u - u) > 2)
        continue;

      xsec.points.push_back({cluster.strokes[stroke_idx].stroke_ind,
                             cluster.strokes[stroke_idx].cluster_ind,
                             stroke_idx, i, point, to_next,
                             tangent(stroke.points, i)});
    }
  }

  for (size_t a = 0; a < xsec.points.size(); ++a) {
    for (size_t b = a + 1; b < xsec.points.size(); ++b) {
      xsec.connections.push_back({a, b, 1.});
    }
  }

  return xsec;
}

void Parameterization::ensure_connected(Cluster *cluster, Cluster::XSec *cut) {
  struct StrokeSegment {
    size_t stroke_idx;
    double begin;
    double end;
    bool small;
    bool cut;
  };
  std::unordered_map<size_t, std::vector<StrokeSegment>> segments;

  auto connection_crosses_cut = [&](const Cluster::XSec &xsec) {
    for (auto &conn : xsec.connections) {
      for (size_t idx : {conn.a_idx, conn.b_idx}) {
        size_t stroke = xsec.points[idx].stroke_idx;
        double i = xsec.points[idx].i;

        for (auto &pt : cut->points) {
          if (pt.stroke_idx == stroke && std::ceil(pt.i) >= std::ceil(i) &&
              std::floor(pt.i) <= std::floor(i)) {
            return true;
          }
        }
      }
    }
    return false;
  };

  for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
       ++stroke_idx) {
    std::vector<double> cuts;
    if (cut) {
      for (auto &pt : cut->points) {
        if (pt.stroke_idx == stroke_idx) {
          cuts.push_back(pt.i);
        }
      }
    }
    std::sort(cuts.begin(), cuts.end());
    StrokeSegment segment{stroke_idx, -1., 0.0, false, !cuts.empty()};
    for (double cut_i : cuts) {
      segment.end = cut_i;
      segment.small =
        std::min<int>(cluster->strokes[stroke_idx].points.size() - 1,
                      std::floor(segment.end)) -
          std::max<int>(0, std::ceil(segment.begin)) <=
        2;
      // if (std::min<int>(cluster->strokes[stroke_idx].points.size() - 1,
      // std::floor(segment.end)) - 	std::max<int>(0,
      // std::ceil(segment.begin)) >
      // 2) {
      segments[stroke_idx].push_back(segment);
      //}
      segment.begin = cut_i;
    }
    segment.end = cluster->strokes[stroke_idx].points.size();
    segment.small =
      std::min<int>(cluster->strokes[stroke_idx].points.size() - 1,
                    std::floor(segment.end)) -
        std::max<int>(0, std::ceil(segment.begin)) <=
      2;
    // if (std::min<int>(cluster->strokes[stroke_idx].points.size() - 1,
    // std::floor(segment.end)) - 	std::max<int>(0,
    // std::ceil(segment.begin))
    // >
    // 2) {
    segments[stroke_idx].push_back(segment);
    //}
  }

  StrokeSegment *first_non_cut_segment = &segments[0][0];
  for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
       ++stroke_idx) {
    if (segments[stroke_idx].size() == 1) {
      first_non_cut_segment = &segments[stroke_idx][0];
      break;
    }
  }

  bool ok = false;
  while (!ok) {
    std::unordered_map<StrokeSegment *, bool> connected;
    for (auto &kv : segments) {
      for (auto &seg : kv.second) {
        connected[&seg] = false;
      }
    }

    std::function<void(StrokeSegment *)> visit =
      [&](StrokeSegment *seg) -> void {
      if (connected[seg])
        return;
      connected[seg] = true;

      for (auto &xsec : cluster->xsecs) {
        for (auto &c : xsec.connections) {
          if (c.weight < 1e-4)
            continue;
          bool connected_a =
            xsec.points[c.a_idx].stroke_idx == seg->stroke_idx &&
            xsec.points[c.a_idx].i > seg->begin &&
            xsec.points[c.a_idx].i < seg->end;
          bool connected_b =
            xsec.points[c.b_idx].stroke_idx == seg->stroke_idx &&
            xsec.points[c.b_idx].i > seg->begin &&
            xsec.points[c.b_idx].i < seg->end;

          if (connected_a || connected_b) {
            size_t other = connected_a ? c.b_idx : c.a_idx;
            size_t other_stroke = xsec.points[other].stroke_idx;
            double other_i = xsec.points[other].i;

            auto it = std::find_if(
              segments[other_stroke].begin(), segments[other_stroke].end(),
              [&](const StrokeSegment &s) {
                return s.begin <= other_i && s.end > other_i;
              });
            if (it != segments[other_stroke].end()) {
              auto &other_seg = *it;
              visit(&other_seg);
            } else {
              std::cout << "Can't find it" << std::endl;
            }
          }
        }
      }
    };
    visit(first_non_cut_segment);

    double closest_dist = std::numeric_limits<double>::infinity();
    Cluster::XSec connector = {{}, {{0, 1, 1.0}}, 0, 0.0, true};
    ok = true;
    for (auto &disconnected : connected) {
      if (disconnected.second)
        continue;
      if (disconnected.first->cut && disconnected.first->small)
        continue;
      ok = false;
      size_t stroke_a = disconnected.first->stroke_idx;
      double min_i = disconnected.first->begin;
      double max_i = disconnected.first->end;

      for (auto &kv : connected) {
        if (!kv.second)
          continue;
        if (kv.first->cut && kv.first->small)
          continue;
        size_t stroke_b = kv.first->stroke_idx;
        double min_j = kv.first->begin;
        double max_j = kv.first->end;

        for (int ia :
             {std::max<int>(0, std::ceil(min_i)),
              std::min<int>(cluster->strokes[stroke_a].points.size() - 1,
                            std::floor(max_i))}) {
          for (int ib :
               {std::max<int>(0, std::ceil(min_j)),
                std::min<int>(cluster->strokes[stroke_b].points.size() - 1,
                              std::floor(max_j))}) {

            double dist = glm::distance2(cluster->strokes[stroke_a].points[ia],
                                         cluster->strokes[stroke_b].points[ib]);
            if (dist < closest_dist) {
              glm::dvec2 to_next_a(0.0, 0.0);
              if (ia == 0)
                to_next_a = cluster->strokes[stroke_a].points[ia + 1] -
                            cluster->strokes[stroke_a].points[ia];
              else
                to_next_a = cluster->strokes[stroke_a].points[ia] -
                            cluster->strokes[stroke_a].points[ia - 1];

              glm::dvec2 to_next_b(0.0, 0.0);
              if (ib == 0)
                to_next_b = cluster->strokes[stroke_b].points[ib + 1] -
                            cluster->strokes[stroke_b].points[ib];
              else
                to_next_b = cluster->strokes[stroke_b].points[ib] -
                            cluster->strokes[stroke_b].points[ib - 1];

              auto new_connector = connector;
              new_connector.points = {
                {cluster->strokes[stroke_a].stroke_ind,
                 cluster->strokes[stroke_a].cluster_ind, stroke_a, double(ia),
                 cluster->strokes[stroke_a].points[ia], to_next_a,
                 tangent(cluster->strokes[stroke_a].points, double(ia))},
                {cluster->strokes[stroke_b].stroke_ind,
                 cluster->strokes[stroke_b].cluster_ind, stroke_b, double(ib),
                 cluster->strokes[stroke_b].points[ib], to_next_b,
                 tangent(cluster->strokes[stroke_b].points, double(ib))},
              };
              new_connector.u = 0.5 * (cluster->strokes[stroke_a].u[ia] +
                                       cluster->strokes[stroke_b].u[ib]);

              if (!cut || !connection_crosses_cut(new_connector)) {
                closest_dist = dist;
                connector = new_connector;
              }
            }
          }
        }
      }
    }

    if (!ok) {
      if (connector.points.size() == 0) {
        break;
      }
      cluster->xsecs.push_back(connector);
    }
  }
}

void Parameterization::check_periodic(Cluster *cluster) {
  std::vector<Cluster::XSec> ortho_xsecs = orthogonal_xsecs(*cluster, true);

  struct Node;

  struct Edge {
    Node *to;
    int xsec_idx;
  };

  struct Node {
    size_t stroke_idx;
    int i;
    std::vector<Edge> edges;
  };

  std::vector<std::vector<Node>> nodes;
  nodes.reserve(cluster->strokes.size());
  for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size();
       ++stroke_idx) {
    nodes.emplace_back(cluster->strokes[stroke_idx].u.size(),
                       Node{stroke_idx, 0, {}});
    for (int i = 0; i < nodes[stroke_idx].size(); ++i) {
      nodes[stroke_idx][i].i = i;

      // Connections due to monotonicity
      if (i < nodes[stroke_idx].size() - 1) {
        nodes[stroke_idx][i].edges.push_back({&nodes[stroke_idx][i + 1], -1});
      }
    }
  }
  // Known isolines
  for (int xsec_idx = 0; xsec_idx < cluster->xsecs.size(); xsec_idx += 3) {
    auto &xsec = cluster->xsecs[xsec_idx];
    for (auto &conn : xsec.connections) {
      size_t sa = xsec.points[conn.a_idx].stroke_idx;
      size_t sb = xsec.points[conn.b_idx].stroke_idx;
      double raw_ia = xsec.points[conn.a_idx].i;
      double raw_ib = xsec.points[conn.b_idx].i;
      int ia = std::round(raw_ia);
      int ib = std::round(raw_ib);
      nodes[sa][ia].edges.push_back({&nodes[sb][ib], -1});
      nodes[sb][ib].edges.push_back({&nodes[sa][ia], -1});
    }
  }
  // Orthogonal jumps
  for (int xsec_idx = 0; xsec_idx < ortho_xsecs.size(); ++xsec_idx) {
    auto &xsec = ortho_xsecs[xsec_idx];
    for (auto &conn : xsec.connections) {
      size_t sa = xsec.points[conn.a_idx].stroke_idx;
      size_t sb = xsec.points[conn.b_idx].stroke_idx;
      double raw_ia = xsec.points[conn.a_idx].i;
      double raw_ib = xsec.points[conn.b_idx].i;
      int ia = std::round(raw_ia);
      int ib = std::round(raw_ib);
      nodes[sa][ia].edges.push_back({&nodes[sb][ib], xsec_idx});
      nodes[sb][ib].edges.push_back({&nodes[sa][ia], xsec_idx});
    }
  }

  auto pair = std::minmax_element(
    cluster->xsecs.begin(), cluster->xsecs.end(),
    [](const Cluster::XSec &a, const Cluster::XSec &b) { return a.u < b.u; });
  Cluster::XSec &from = *pair.second;
  Cluster::XSec &to = *pair.first;
  double xsec_step = 3;
  if (cluster->xsecs.size() > 2)
    xsec_step = (cluster->xsecs[1].u - cluster->xsecs[0].u) * 3;
  double close_theshold =
    std::ceil(std::max(std::min(10., (from.u - to.u) * 0.02), 4.) / xsec_step);
  auto is_close = [&](const Node &a, const Cluster::XSec &b) {
    for (auto &p_b : b.points) {
      // if (a.stroke_idx == p_b.stroke_idx &&
      //     std::abs(std::double_t(a.i) - p_b.i) < 5.0) {
      if (a.stroke_idx == p_b.stroke_idx &&
          std::abs(std::double_t(a.i) - p_b.i) < close_theshold) {
        return true;
      }
    }
    return false;
  };

  auto find_path = [&](const Node &from,
                       const Cluster::XSec &to) -> std::vector<int> {
    int allowed_jumps = 2;
    std::unordered_map<int, std::vector<std::vector<bool>>> visited;
    for (int jumped = 0; jumped <= allowed_jumps; ++jumped) {
      visited[jumped].reserve(nodes.size());
      for (auto &stroke_nodes : nodes) {
        visited[jumped].emplace_back(stroke_nodes.size(), false);
      }
    }

    std::function<std::vector<int>(const Node *, int)> visit =
      [&](const Node *n, int jumps) -> std::vector<int> {
      if (visited[jumps][n->stroke_idx][n->i])
        return {};
      visited[jumps][n->stroke_idx][n->i] = true;

      for (auto &edge : n->edges) {
        if (jumps == 0 && edge.xsec_idx != -1)
          continue;

        if (is_close(*edge.to, to)) {
          if (false && context.debug_viz) {
            add_debug_line(
              {cluster->strokes[n->stroke_idx].points[n->i],
               cluster->strokes[edge.to->stroke_idx].points[edge.to->i],
               edge.xsec_idx == -1 ? "#000" : "#F0F"});
          }
          return {edge.xsec_idx};
        } else {
          auto res = visit(edge.to, jumps - (edge.xsec_idx != -1 ? 1 : 0));
          if (!res.empty()) {
            if (false && context.debug_viz) {
              add_debug_line(
                {cluster->strokes[n->stroke_idx].points[n->i],
                 cluster->strokes[edge.to->stroke_idx].points[edge.to->i],
                 edge.xsec_idx == -1 ? "#000" : "#F0F"});
            }
            res.push_back(edge.xsec_idx);
            return res;
          }
        }
      }

      return {};
    };

    return visit(&from, allowed_jumps);
  };

  auto get_u = [&](const Cluster::XSecPoint &pt) {
    double mix = pt.i - std::floor(pt.i);
    return (1. - mix) * cluster->strokes[pt.stroke_idx].u[std::floor(pt.i)] +
           mix * cluster->strokes[pt.stroke_idx].u[std::ceil(pt.i)];
  };

  std::vector<double> period_samples;
  for (auto &stroke_nodes : nodes) {
    for (auto &node : stroke_nodes) {
      if (is_close(node, from)) {
        auto path = find_path(node, to);
        if (!path.empty()) {
          for (int xsec_idx : path) {
            if (xsec_idx == -1)
              continue;
            auto &step = ortho_xsecs[xsec_idx];
            for (auto &conn : step.connections) {
              double u_a = get_u(step.points[conn.a_idx]);
              double u_b = get_u(step.points[conn.b_idx]);
              double diff = std::abs(u_a - u_b);
              if (diff > 10.) {
                period_samples.push_back(diff);
              }
            }
          }
        }
      }
    }
  }

  // Not periodic
  if (period_samples.empty())
    return;

  cluster->periodic = true;

  int off = 0.2 * (period_samples.size() - 1);
  std::nth_element(period_samples.begin(), period_samples.begin() + off,
                   period_samples.end());
  double period = period_samples[off];

  // if (period < 0.8 * (from.u - to.u)) {
  //   cluster->periodic = false;
  //   return;
  // }

  // Make all u values be \in [0, period)
  for (auto &stroke : cluster->strokes) {
    for (auto &u : stroke.u) {
      while (u >= period) {
        u -= period;
      }
    }
  }
  bool ensure_sample_reference = true;
  cluster->xsecs =
    xsecs_from_params(*cluster, ensure_sample_reference,
                      stroke_sampling_size / stroke_sampling_size_ref);
}
