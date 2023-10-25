#include <cmath>
#include <future>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "Utils.h"

// For spiral cutting (assuming we have fewer strokes than spiral_offset)
size_t spiral_offset = 10000;

std::vector<Intersection> intersections(const std::vector<glm::dvec2> &polyline,
                                        glm::dvec2 from, glm::dvec2 to) {
  bool vertical = std::abs(from.x - to.x) < 1e-4;
  std::vector<double> idx_b;
  std::vector<double> x;
  std::vector<double> y;
  const int numSegments = polyline.size() - 1;

  double x0a = from.x;
  double x1a = to.x;
  double y0a = from.y;
  double y1a = to.y;

  std::vector<Intersection> intersections;
  for (int s = 0; s < numSegments; s++) {
    double x0b = polyline[s].x;
    double x1b = polyline[s + 1].x;
    double y0b = polyline[s].y;
    double y1b = polyline[s + 1].y;

    bool other_vertical = std::abs(x0b - x0a) < 1e-4;

    // Ignore if outside bbox
    if (std::min(x0b, x1b) > std::max(x0a, x1a) ||
        std::max(x0b, x1b) < std::min(x0a, x1a) ||
        std::min(y0b, y1b) > std::max(y0a, y1a) ||
        std::max(y0b, y1b) < std::min(y0a, y1a))
      continue;

    double tb = 0;
    if (vertical) {
      if (other_vertical) {
        bool intersection = std::max(std::min(y0a, y1a), std::min(y0b, y1b)) <
                            std::min(std::max(y0a, y1a), std::max(y0b, y1b));
        if (intersection) {
          tb = (std::max(std::min(y0a, y1a), std::min(y0b, y1b)) - y0b) /
               (y1b - y0b);
        } else {
          continue;
        }
      } else {
        tb = (x0a - x0b) / (x1b - x0b);
        if (tb < 0 || tb > 1)
          continue;
      }

    } else {
      double tb_num = y0a + (y1a - y0a) * (x0b - x0a) / (x1a - x0a) - y0b;
      double tb_denom = y1b - (y1a - y0a) * (x1b - x0b) / (x1a - x0a) - y0b;
      tb = tb_num / tb_denom;
      if (tb < 0 || tb > 1)
        continue;

      double ta = (x0b + (x1b - x0b) * tb - x0a) / (x1a - x0a);
      if (ta < 0 || ta > 1)
        continue;
    }

    intersections.emplace_back(
      Intersection{static_cast<double>(s) + tb,
                   glm::dvec2(x0b + tb * (x1b - x0b), y0b + tb * (y1b - y0b))});
  }

  return intersections;
}

double gaussian(double x, double sigma, double mu) {
  return std::exp(-0.5 * std::pow((x - mu) / sigma, 2));
}

double poly_in_out(double t, double k) {
  double t2 = 2 * t;
  if (t2 <= 1) {
    return std::pow(t2, k) / 2.;
  } else {
    return (2. - std::pow(2. - t2, k)) / 2.;
  }
}

double LinExpr::getValue(const std::vector<GRBVar> &param_vars) const {
  double value = 0;
  for (const auto &v_coef : var_idx_coeff) {
    value += v_coef.second * param_vars[v_coef.first].get(GRB_DoubleAttr_X);
  }
  value += const_coeff;
  return value;
}

GRBLinExpr l1_norm(GRBModel *model, const std::vector<GRBLinExpr> &x) {
  // min ||x||_1
  //
  // ...is equivalent to:
  //
  // min t
  // s.t.
  // x_i <= y_i,
  // -x_i <= y_i,
  // \sum_i y_i = t

  GRBLinExpr sum_y = 0.0;
  auto t =
    model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "L1");
  std::vector<GRBVar> y;
  y.reserve(x.size());
  size_t l1_i = 0;
  for (auto &term : x) {
    std::string name = std::string("l1_" + std::to_string(l1_i++));
    y.push_back(
      model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, name));
    sum_y += y.back();
    model->addConstr(term <= y.back());
    model->addConstr(-term <= y.back());
  }
  model->addConstr(sum_y == t);

  return t;
}

GRBLinExpr l1_norm(GRBModel *model, const std::vector<LinExpr> &x,
                   const std::vector<GRBVar> &param_vars, bool add_length) {
  GRBLinExpr sum_y = 0.0;
  // auto t =
  //   model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, "L1");
  std::vector<GRBVar> y;
  y.reserve(x.size());
  size_t l1_i = 0;
  for (auto &term : x) {
    GRBLinExpr grb_term = 0;
    for (const auto &v_coef : term.var_idx_coeff) {
      grb_term += v_coef.second * param_vars[v_coef.first];
    }
    if (add_length)
      grb_term += term.const_coeff;

    std::string name = std::string("l1_" + std::to_string(l1_i++));
    y.push_back(
      model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, name));
    sum_y += y.back();
    model->addConstr(grb_term <= y.back());
    model->addConstr(-grb_term <= y.back());
  }
  // model->addConstr(sum_y == t);

  // return t;
  return sum_y;
}

GRBLinExpr l1_norm_constr(GRBModel *model, const std::vector<LinExpr> &x,
                          const std::vector<GRBVar> &param_vars,
                          bool add_length) {
  GRBLinExpr sum_y = 0.0;
  std::vector<GRBVar> y;
  y.reserve(x.size());
  std::vector<GRBVar> diff;
  diff.reserve(x.size());
  size_t l1_i = 0;
  for (auto &term : x) {
    GRBLinExpr grb_term = 0;
    for (const auto &v_coef : term.var_idx_coeff) {
      grb_term += v_coef.second * param_vars[v_coef.first];
    }
    if (add_length)
      grb_term += term.const_coeff;

    std::string name = std::string("l1_" + std::to_string(l1_i++));
    y.push_back(
      model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, name));
    std::string diff_name = std::string("u_diff_" + std::to_string(l1_i++));
    diff.push_back(model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0,
                                 GRB_CONTINUOUS, diff_name));
    sum_y += y.back();
    model->addConstr(grb_term == diff.back());
    model->addGenConstrAbs(y.back(), diff.back());
  }

  return sum_y;
}

GRBQuadExpr l2_norm_sq(GRBModel *model, const std::vector<GRBLinExpr> &x) {
  GRBQuadExpr result = 0.0;
  for (auto &term : x) {
    result += term * term;
  }

  return result;
}

GRBQuadExpr l2_norm_sq(GRBModel *model, const std::vector<LinExpr> &x,
                       const std::vector<GRBVar> &param_vars) {
  // Sum terms
  std::map<size_t, double> quad_coeff;
  std::map<std::pair<size_t, size_t>, double> cross_coeff;
  std::map<size_t, double> lin_coeff;
  GRBQuadExpr result = 0.0;

  for (const auto &term : x) {
    for (const auto &v_coef : term.quad_coeff)
      quad_coeff[v_coef.first] += v_coef.second;
    for (const auto &v_coef : term.cross_coeff)
      cross_coeff[v_coef.first] += v_coef.second;
    for (const auto &v_coef : term.lin_coeff)
      lin_coeff[v_coef.first] += v_coef.second;

    // The constant doesn't matter to the optimization. Ignore.
    // result += term.const_coeff * term.const_coeff;
  }

  // Assemble equation
  std::vector<double> quad_vec, lin_vec;
  quad_vec.resize(param_vars.size(), 0);
  lin_vec.resize(param_vars.size(), 0);

  for (const auto &v_coef : quad_coeff)
    quad_vec[v_coef.first] = v_coef.second;
  for (const auto &v_coef : lin_coeff)
    lin_vec[v_coef.first] = v_coef.second;

  result.addTerms(quad_vec.data(), param_vars.data(), param_vars.data(),
                  param_vars.size());
  result.addTerms(lin_vec.data(), param_vars.data(), param_vars.size());

  for (const auto &v_coef : cross_coeff)
    result.addTerm(v_coef.second, param_vars[v_coef.first.first],
                   param_vars[v_coef.first.second]);

  return result;
}

void center_cluster(Cluster &cluster, glm::dvec2 &center) {
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &s : cluster.strokes) {
    for (const auto &point : s.points) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }
  }
  center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);

  for (auto &s : cluster.strokes) {
    for (auto &point : s.points) {
      point -= center;
    }
  }
  for (auto &s : cluster.xsecs) {
    for (auto &p : s.points) {
      p.point -= center;
    }
  }
}

void decenter_cluster(Cluster &cluster, const glm::dvec2 &center) {
  for (auto &s : cluster.strokes) {
    for (auto &point : s.points) {
      point += center;
    }
  }
  for (auto &s : cluster.xsecs) {
    for (auto &p : s.points) {
      p.point += center;
    }
  }
  for (auto &ss : cluster.orthogonal_xsecs_list) {
    for (auto &s : ss) {
      for (auto &p : s.points) {
        p.point += center;
      }
    }
  }
}

void center_input(const Input &input, double &w, double &h, double &x,
                  double &y, glm::dvec2 &center) {
  double padding = 5 * input.thickness;

  // Center
  center = glm::dvec2(0, 0);
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (auto &kv : input.clusters) {
    for (auto &stroke : kv.second.strokes) {
      for (const auto &point : stroke.points) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
      }
    }
  }
  center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
  w = (max_x - min_x) * input.thickness + 2 * padding;
  h = (max_y - min_y) * input.thickness + 2 * padding;
  x = -w / 2;
  y = -h / 2;
}

bool is_roughly_circular(const Cluster::Stroke &s) {
  // std::cout << "\t In is_roughly_circular" << std::endl;
  if (s.points.empty())
    return false;
  if (glm::distance(s.points.front(), s.points.back()) <
      std::numeric_limits<double>::epsilon())
    return true;

  double turning = 0;
  double min_angle_dist = 5 * M_PI;
  glm::dvec2 near_round_p = s.points.at(0);
  double closest_dist = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i + 1 < s.points.size(); i++) {
    glm::dvec2 t0 = s.points.at(i) - s.points.at(i - 1);
    glm::dvec2 t1 = s.points.at(i + 1) - s.points.at(i);
    t0 = glm::normalize(t0);
    t1 = glm::normalize(t1);
    glm::dvec2 up0(-t0.y, t0.x);
    double sign = (glm::dot(up0, t1) > 0) ? 1 : -1;
    turning += sign * std::acos(std::max(-1., std::min(glm::dot(t0, t1), 1.)));
    double angle_dist = std::abs(std::abs(turning) - 2 * M_PI);
    if (angle_dist < min_angle_dist) {
      near_round_p = s.points.at(i);
      min_angle_dist = angle_dist;
    }
    double c_dist = glm::distance(s.points.at(i), s.points.at(0));
    if (i > s.points.size() / 2 && c_dist < closest_dist) {
      closest_dist = c_dist;
    }
  }
  bool angle_check = (std::abs(turning) > 300. / 180 * M_PI);
  bool dist_check = glm::distance(near_round_p, s.points.at(0)) < 20;
  bool c_dist_check = closest_dist < 10;
  // std::cout << "\t Turn: " << s.stroke_ind << " - " << turning / M_PI * 180
  //           << "; " << glm::distance(near_round_p, s.points.at(0)) << "; "
  //           << closest_dist << std::endl;
  return angle_check || c_dist_check;
}

void compute_spiral_angles(Cluster::Stroke &s) {
  s.spiral_angles.clear();
  s.spiral_angles.emplace_back(0);
  double curv = 0;
  for (size_t i = 1; i + 1 < s.points.size(); ++i) {
    glm::dvec2 prev = glm::normalize(s.points[i] - s.points[i - 1]);
    glm::dvec2 cur = glm::normalize(s.points[i + 1] - s.points[i]);
    glm::dvec2 prev_norm(-prev.y, prev.x);

    double sign = (glm::dot(prev_norm, cur) > 0) ? 1 : -1;
    double dot = glm::dot(prev, cur);
    double angle = sign * std::acos(std::min(std::max(-1., dot), 1.));
    curv += angle;

    s.spiral_angles.emplace_back(curv);
  }
  s.spiral_angles.emplace_back(s.spiral_angles.back());
}

void cut_spiral_stroke(const Cluster::Stroke &spiral_s,
                       std::vector<Cluster::Stroke> &cut_ss,
                       const double cut_curv) {
  assert(spiral_s.spiral);
  if (spiral_s.points.empty()) {
    cut_ss.emplace_back(spiral_s);
    return;
  }

  double curv = 0;
  cut_ss.emplace_back();
  cut_ss.back().spiral = true;
  cut_ss.back().stroke_ind = spiral_s.stroke_ind;
  cut_ss.back().cluster_ind = spiral_s.cluster_ind;
  cut_ss.back().points.emplace_back(spiral_s.points[0]);
  cut_ss.back().u.emplace_back(spiral_s.u[0]);
  cut_ss.back().spiral_angles.emplace_back(spiral_s.spiral_angles[0]);
  for (size_t i = 1; i + 1 < spiral_s.points.size(); ++i) {
    glm::dvec2 prev =
      glm::normalize(spiral_s.points[i] - spiral_s.points[i - 1]);
    glm::dvec2 cur =
      glm::normalize(spiral_s.points[i + 1] - spiral_s.points[i]);
    glm::dvec2 prev_norm(-prev.y, prev.x);

    double sign = (glm::dot(prev_norm, cur) > 0) ? 1 : -1;
    double dot = glm::dot(prev, cur);
    double angle = sign * std::acos(std::min(std::max(-1., dot), 1.));
    curv += angle;

    // Add to current resulting stroke
    cut_ss.back().points.emplace_back(spiral_s.points[i]);
    cut_ss.back().u.emplace_back(spiral_s.u[i]);
    cut_ss.back().spiral_angles.emplace_back(spiral_s.spiral_angles[i]);

    // Cut
    if (std::abs(curv) > cut_curv) {
      cut_ss.emplace_back();
      cut_ss.back().spiral = true;
      cut_ss.back().stroke_ind =
        spiral_s.stroke_ind + spiral_offset * (cut_ss.size() - 1);
      cut_ss.back().cluster_ind = spiral_s.cluster_ind;
      cut_ss.back().points.emplace_back(spiral_s.points[i]);
      cut_ss.back().u.emplace_back(spiral_s.u[i]);
      cut_ss.back().spiral_angles.emplace_back(spiral_s.spiral_angles[i]);

      curv = 0;
    }
  }

  cut_ss.back().points.emplace_back(spiral_s.points.back());
  cut_ss.back().u.emplace_back(spiral_s.u.back());
  cut_ss.back().spiral_angles.emplace_back(spiral_s.spiral_angles.back());
}

void assemble_spiral_stroke(const std::vector<Cluster::Stroke> &cut_ss,
                            std::vector<Cluster::Stroke> &spiral_ss,
                            std::vector<Cluster::XSec> &xsecs) {
  std::map<size_t, std::vector<const Cluster::Stroke *>> to_assemble;

  // Order strokes
  for (const auto &s : cut_ss) {
    if (s.spiral) {
      size_t real_sid = s.stroke_ind % spiral_offset;
      to_assemble[real_sid].emplace_back(&s);
    } else {
      spiral_ss.emplace_back(s);
    }
  }

  if (to_assemble.empty())
    return;

  // Assemble
  std::map<size_t, std::vector<size_t>> strokes_cut_xsec_counts,
    strokes_cut_point_counts;
  for (const auto &xsec : xsecs) {
    for (const auto &p : xsec.points) {
      size_t real_sid = p.stroke_ind % spiral_offset;
      if (to_assemble.count(real_sid)) {
        if (!strokes_cut_xsec_counts.count(real_sid)) {
          strokes_cut_xsec_counts[real_sid].resize(to_assemble[real_sid].size(),
                                                   0);
        }

        size_t i = p.stroke_ind / spiral_offset;
        strokes_cut_xsec_counts[real_sid][i] = p.stroke_xsec_count;
      }
    }
  }
  for (auto &sid_ss : to_assemble) {
    std::sort(sid_ss.second.begin(), sid_ss.second.end(),
              [](const Cluster::Stroke *a, const Cluster::Stroke *b) {
                return a->stroke_ind < b->stroke_ind;
              });

    auto &cut_point_counts = strokes_cut_point_counts[sid_ss.first];
    spiral_ss.emplace_back();
    spiral_ss.back().spiral = true;
    spiral_ss.back().cluster_ind = sid_ss.second.front()->cluster_ind;
    spiral_ss.back().stroke_ind = sid_ss.first;
    for (const auto &s : sid_ss.second) {
      cut_point_counts.emplace_back(s->points.size());

      for (size_t i = 0; i + 1 < s->points.size(); ++i) {
        spiral_ss.back().points.emplace_back(s->points[i]);
        spiral_ss.back().u.emplace_back(s->u[i]);
        spiral_ss.back().spiral_angles.emplace_back(s->spiral_angles[i]);
      }
    }
    spiral_ss.back().points.emplace_back(sid_ss.second.back()->points.back());
    spiral_ss.back().u.emplace_back(sid_ss.second.back()->u.back());
    spiral_ss.back().spiral_angles.emplace_back(
      sid_ss.second.back()->spiral_angles.back());
  }

  // Map xsecs
  std::map<size_t, size_t> cut_xsec_totals;
  for (const auto &sid_xsecs : strokes_cut_xsec_counts) {
    for (const auto &xsec_c : sid_xsecs.second)
      cut_xsec_totals[sid_xsecs.first] += xsec_c - 1;
    cut_xsec_totals[sid_xsecs.first]++;
  }
  for (auto &xsec : xsecs) {
    for (auto &p : xsec.points) {
      size_t real_sid = p.stroke_ind % spiral_offset;
      size_t cut_index = p.stroke_ind / spiral_offset;
      if (to_assemble.count(real_sid)) {
        p.stroke_ind = real_sid;
        size_t cut_i_offset = 0, cut_xsec_offset = 0;
        for (size_t i = 0; i < cut_index; ++i) {
          assert(i < strokes_cut_point_counts[real_sid].size());
          assert(i < strokes_cut_xsec_counts[real_sid].size());
          cut_i_offset += strokes_cut_point_counts[real_sid][i] - 1;
          cut_xsec_offset += strokes_cut_xsec_counts[real_sid][i] - 1;
          if (i + 1 == cut_index) {
            cut_i_offset++;
            cut_xsec_offset++;
          }
        }
        p.i += cut_i_offset;
        p.stroke_within_idx += cut_xsec_offset;
      }
    }
  }

  std::sort(spiral_ss.begin(), spiral_ss.end(),
            [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
              return a.stroke_ind < b.stroke_ind;
            });
}

void connect_stroke(
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    &matching_pair,
  std::vector<Cluster::Stroke> &strokes,
  std::map<size_t, size_t> *connect_map_ptr) {
  // Merge strokes
  std::map<size_t, Cluster::Stroke *> s_map;
  for (auto &s : strokes) {
    s_map[s.stroke_ind] = &s;
  }

  std::set<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    matching_pair_indices;
  for (const auto &pair : matching_pair) {
    auto s_head = std::make_pair(pair.first.first, pair.first.second);
    auto s_tail = std::make_pair(pair.second.first, pair.second.second);
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
      matching_pair_indices.emplace(
        std::make_pair(std::make_pair(pair.first.first, pair.first.second),
                       std::make_pair(pair.second.first, pair.second.second)));
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
  std::vector<Cluster::Stroke> result_strokes;

  auto merge_strokes = [](Cluster::Stroke &to_s, const Cluster::Stroke &from_s,
                          bool is_to_tail, bool is_from_tail,
                          std::pair<size_t, bool> &s_head,
                          std::pair<size_t, bool> &s_tail) {
    // Move one point away from the endpoint to avoid segment with 0-length.
    if (from_s.points.size() <= 1)
      return;
    if (is_to_tail && !is_from_tail) {
      to_s.points.insert(to_s.points.end(), from_s.points.begin() + 1,
                         from_s.points.end());
      if (!to_s.u.empty())
        to_s.u.insert(to_s.u.end(), from_s.u.begin() + 1, from_s.u.end());
      s_tail.first = from_s.stroke_ind;
      s_tail.second = true;
    } else if (is_to_tail && is_from_tail) {
      to_s.points.insert(to_s.points.end(), from_s.points.rbegin() + 1,
                         from_s.points.rend());
      if (!to_s.u.empty())
        to_s.u.insert(to_s.u.end(), from_s.u.rbegin() + 1, from_s.u.rend());
      s_tail.first = from_s.stroke_ind;
      s_tail.second = false;
    } else if (!is_to_tail && is_from_tail) {
      to_s.points.erase(to_s.points.begin());
      to_s.points.insert(to_s.points.begin(), from_s.points.begin(),
                         from_s.points.end());
      if (!to_s.u.empty()) {
        to_s.u.erase(to_s.u.begin());
        to_s.u.insert(to_s.u.begin(), from_s.u.begin(), from_s.u.end());
      }
      s_head.first = from_s.stroke_ind;
      s_head.second = false;
    } else if (!is_to_tail && !is_from_tail) {
      to_s.points.erase(to_s.points.begin());
      to_s.points.insert(to_s.points.begin(), from_s.points.rbegin(),
                         from_s.points.rend());
      if (!to_s.u.empty()) {
        to_s.u.erase(to_s.u.begin());
        to_s.u.insert(to_s.u.begin(), from_s.u.rbegin(), from_s.u.rend());
      }
      s_head.first = from_s.stroke_ind;
      s_head.second = true;
    }

    // If any of the merged stroke is spiral then the resulting stroke is spiral
    to_s.spiral |= from_s.spiral;
  };

  std::map<size_t, size_t> cut_map;
  std::set<std::pair<size_t, bool>> valid_cut_points;
  for (const auto &g : merge_group) {
    if (g.empty())
      continue;
    auto s = *s_map[*g.begin()];
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
    result_strokes.emplace_back(s);
  }

  // Add not merged strokes
  for (const auto &s : strokes) {
    if (!merge_group_index.count(s.stroke_ind)) {
      result_strokes.emplace_back(s);
      valid_cut_points.emplace(s.stroke_ind, false);
      valid_cut_points.emplace(s.stroke_ind, true);
    }
  }
  std::sort(result_strokes.begin(), result_strokes.end(),
            [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
              return a.stroke_ind < b.stroke_ind;
            });
  if (connect_map_ptr) {
    std::map<size_t, size_t> s2i;
    for (size_t i = 0; i < strokes.size(); ++i)
      s2i[strokes[i].stroke_ind] = i;
    for (size_t i = 0; i < result_strokes.size(); ++i)
      (*connect_map_ptr)[result_strokes[i].stroke_ind] =
        s2i[result_strokes[i].stroke_ind];
  }
  strokes = result_strokes;
}

void sort_xsec(Cluster::XSec &xsec) {
  if (xsec.points.size() <= 1)
    return;
  struct PointPos {
    Cluster::XSecPoint p;
    size_t idx;
  };
  auto PointPos_compare = [](const std::pair<double, PointPos> &x,
                             const std::pair<double, PointPos> &y) {
    return x.first < y.first ||
           (!(y.first < x.first) && x.second.idx < y.second.idx);
  };
  std::set<std::pair<double, PointPos>, decltype(PointPos_compare)> sort_xsecs(
    PointPos_compare);
  size_t center_idx = xsec.center_idx;
  glm::dvec2 v_axis = xsec.points[center_idx].tangent;
  glm::dvec2 origin = xsec.points[center_idx].point;
  v_axis = glm::normalize(v_axis);
  v_axis = glm::dvec2(-v_axis.y, v_axis.x);

  for (size_t i = 0; i < xsec.points.size(); ++i) {
    double v = glm::dot(xsec.points[i].point - origin, v_axis);
    PointPos pp;
    pp.p = xsec.points[i];
    pp.idx = i;
    sort_xsecs.emplace(std::make_pair(v, pp));
  }

  std::map<size_t, size_t> reindexing;
  xsec.points.clear();
  for (const auto &pp : sort_xsecs) {
    reindexing[pp.second.idx] = xsec.points.size();
    xsec.points.emplace_back(pp.second.p);
  }

  for (auto conn : xsec.connections) {
    conn.a_idx = reindexing[conn.a_idx];
    conn.b_idx = reindexing[conn.b_idx];
  }
  xsec.center_idx = reindexing[xsec.center_idx];
}
