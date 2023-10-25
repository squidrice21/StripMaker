#include "JunctionMatching.h"

#include "../Util.h"

#include "Junction.h"
#include "JunctionProposal.h"
#include "SketchConnectivity/classifier.h"
#include "SketchConnectivity/closest.h"
#include "SketchConnectivity/intersect.h"
#include "SketchConnectivity/sketching.h"
#include "glm/detail/type_vec.hpp"

#include <map>
#include <set>
#include <vector>

double junc_threshold = 0.5;
double stroke_stroke_distance_ratio = 1.5;
double junc_close_width_ratio = 5;

bool to_distance_check = true;

double sample_u(Cluster::XSecPoint sample, const Cluster &merged_cluster) {
  std::map<size_t, const Cluster::Stroke *> sid2s;
  for (const auto &s : merged_cluster.strokes) {
    sid2s[s.stroke_ind] = &s;
  }

  double dist;
  size_t si = stroke_sample_projection_index(sample.point,
                                             *sid2s[sample.stroke_ind], dist);
  assert(dist < 2);
  return sid2s[sample.stroke_ind]->u[si];
}

// -1: Junction on the side with smaller U; 1: Larger U
int junction_side(const Cluster &c, const Cluster &merged_cluster,
                  const Junction &junc, double &end_u) {
  std::map<size_t, const Cluster::Stroke *> sid2s;
  for (const auto &s : merged_cluster.strokes) {
    sid2s[s.stroke_ind] = &s;
  }
  bool is_fit_head = (junc.from.first == c.fit.cluster_idx)
                       ? junc.from.second < 0.5
                       : junc.to.second < 0.5;
  Cluster::XSecPoint end_p = (is_fit_head)
                               ? c.fit.fit_sample_xsecs.front().points.front()
                               : c.fit.fit_sample_xsecs.back().points.front();
  double dist;
  size_t si =
    stroke_sample_projection_index(end_p.point, *sid2s[end_p.stroke_ind], dist);
  assert(dist < 2);
  end_u = sid2s[end_p.stroke_ind]->u[si];

  return (std::abs(end_u - sid2s[end_p.stroke_ind]->u.front()) <
          std::abs(end_u - sid2s[end_p.stroke_ind]->u.back()))
           ? -1
           : 1;
}

bool intersect(const sketching::StrokeBVH &bvh, size_t c1, size_t c2,
               std::vector<sketching::Vec2> &intersections) {
  sketching::PolylineBVH poly_bvh = bvh.polyline_bvh();
  sketching::PolylineBVHLeaf *cur_n_ptr = nullptr;
  sketching::PolylineBVHLeaf *other_n_ptr = nullptr;
  for (auto &n : poly_bvh.nodes) {
    if (n.geometry->index == c1)
      cur_n_ptr = &n;
    if (n.geometry->index == c2)
      other_n_ptr = &n;
  }

  cur_n_ptr->geometry->ensure_arclengths();
  other_n_ptr->geometry->ensure_arclengths();
  intersect_different(*cur_n_ptr, *other_n_ptr, intersections);
  return !intersections.empty();
}

void find_junctions(const std::vector<Junction> &junctions,
                    const FittedCurve &fit1, const FittedCurve &fit2,
                    sketching::JunctionType::Type type,
                    std::vector<Junction> &out_junctions) {
  sketching::Stroke s1;
  convert_stroke(fit1, s1);

  for (auto const &in_junc : junctions) {
    Junction cur_junc = in_junc;
    if (cur_junc.probability < junc_threshold ||
        (cur_junc.from.first != fit1.cluster_idx &&
         cur_junc.to.first != fit1.cluster_idx) ||
        (cur_junc.from.first != fit2.cluster_idx &&
         cur_junc.to.first != fit2.cluster_idx))
      continue;
    if (cur_junc.type != type)
      continue;
    if (cur_junc.from.first != fit1.cluster_idx) {
      std::swap(cur_junc.from, cur_junc.to);
      std::swap(cur_junc.from_pos, cur_junc.to_pos);
    }
    // Not a valid end
    if (!is_roughly_end(s1, cur_junc.from.second)) {
      std::swap(cur_junc.from, cur_junc.to);
      std::swap(cur_junc.from_pos, cur_junc.to_pos);
    }
    if (cur_junc.from.first != fit1.cluster_idx ||
        (cur_junc.type == sketching::JunctionType::X &&
         !is_roughly_end(s1, cur_junc.from.second)) ||
        (cur_junc.type == sketching::JunctionType::T &&
         cur_junc.from.second != 0 && cur_junc.from.second != 1)) {
      continue;
    }
    out_junctions.emplace_back(cur_junc);
  }
}

double closest_points(const sketching::Stroke &stroke1,
                      const sketching::Stroke &stroke2, double &out_arclen1,
                      double &out_arclen2) {
  assert(stroke1.size() > 0);
  assert(stroke2.size() > 0);
  assert(&stroke1 != &stroke2);

  stroke1.ensure_arclengths();
  stroke2.ensure_arclengths();

  const auto n1 = stroke1.size();
  const auto n2 = stroke2.size();
  auto best_dist = double(INFINITY);
  out_arclen1 = NAN;
  out_arclen2 = NAN;
  auto proj = sketching::Vec2(NAN, NAN);
  auto s = double(NAN);

  for (size_t i = 0; i < n1; ++i) {
    const auto dist = closest_point(stroke2, stroke1.xy(i), proj, s);
    if (dist < best_dist) {
      best_dist = dist;
      out_arclen1 = stroke1.arclength(i);
      out_arclen2 = s;
    }
  }
  for (size_t i = 0; i < n2; ++i) {
    const auto dist = closest_point(stroke1, stroke2.xy(i), proj, s);
    if (dist < best_dist) {
      best_dist = dist;
      out_arclen1 = s;
      out_arclen2 = stroke2.arclength(i);
    }
  }
  assert(std::isfinite(out_arclen1));
  assert(std::isfinite(out_arclen2));
  return (stroke1.pos(out_arclen1) - stroke2.pos(out_arclen2)).norm();
}

bool is_separate_forward_check(const Cluster &c1, const Cluster &c2,
                               const Cluster &merged_cluster,
                               const std::vector<FittedCurve> &fits,
                               const std::vector<Junction> &junctions,
                               std::vector<Junction> &decision_junctions) {
  std::vector<Junction> cur_decision_junctions;

  std::map<size_t, size_t> s2sep_c;
  for (auto const &s : c1.strokes) {
    s2sep_c[s.stroke_ind] = 0;
  }
  for (auto const &s : c2.strokes) {
    s2sep_c[s.stroke_ind] = 1;
  }

  // 1. Check if overlapping
  bool is_overlapping = false;
  for (const auto &xsec : merged_cluster.xsecs) {
    std::set<size_t> s_indices;
    for (const auto &p : xsec.points) {
      if (s2sep_c.count(p.stroke_ind))
        s_indices.emplace(s2sep_c[p.stroke_ind]);
      if (s_indices.size() == 2) {
        is_overlapping = true;
        break;
      }
    }
  }

  if (!is_overlapping)
    return false;

  sketching::Stroke s1;
  convert_stroke(c1.fit, s1);
  sketching::Stroke s2;
  convert_stroke(c2.fit, s2);

  // 2. Check the side of the junction
  struct JunctionHit {
    Junction junc;
    double end_u;
    int junc_side;
    double probability;
  };
  std::vector<JunctionHit> hits;

  for (auto const &in_junc : junctions) {
    Junction cur_junc = in_junc;
    if (cur_junc.probability < junc_threshold ||
        (cur_junc.from.first != c1.fit.cluster_idx &&
         cur_junc.to.first != c1.fit.cluster_idx))
      continue;
    if (cur_junc.from.first != c1.fit.cluster_idx) {
      std::swap(cur_junc.from, cur_junc.to);
      std::swap(cur_junc.from_pos, cur_junc.to_pos);
    }
    // Not a valid end
    if (!is_roughly_end(s1, cur_junc.from.second)) {
      std::swap(cur_junc.from, cur_junc.to);
      std::swap(cur_junc.from_pos, cur_junc.to_pos);
    }
    if (cur_junc.from.first != c1.fit.cluster_idx ||
        (cur_junc.type == sketching::JunctionType::X &&
         !is_roughly_end(s1, cur_junc.from.second)) ||
        (cur_junc.type == sketching::JunctionType::T &&
         cur_junc.from.second != 0 && cur_junc.from.second != 1)) {
      continue;
    }
    JunctionHit hit;
    hit.probability = cur_junc.probability;
    hit.junc_side = junction_side(c1, merged_cluster, cur_junc, hit.end_u);
    hit.junc = cur_junc;
    hits.emplace_back(hit);
  }

  if (hits.empty())
    return false;

  std::sort(hits.begin(), hits.end(),
            [](const JunctionHit &a, const JunctionHit &b) {
              return a.probability > b.probability;
            });

  std::map<size_t, FittedCurve const *> c2fit;
  for (auto const &fit2 : fits) {
    c2fit[fit2.cluster_idx] = &fit2;
  }

  double end_u2;
  int junc_side2;
  for (auto const &hit : hits) {
    Junction junc = hit.junc;
    double end_u1 = hit.end_u;
    int junc_side1 = hit.junc_side;

    for (auto const &in_junc : junctions) {
      Junction junc2 = in_junc;
      if (junc2.probability < junc_threshold ||
          (junc2.from.first != c2.fit.cluster_idx &&
           junc2.to.first != c2.fit.cluster_idx))
        continue;
      if (junc2.from.first != c2.fit.cluster_idx) {
        std::swap(junc2.from, junc2.to);
        std::swap(junc2.from_pos, junc2.to_pos);
      }
      // Not a valid end
      if (!is_roughly_end(s2, junc2.from.second)) {
        std::swap(junc2.from, junc2.to);
        std::swap(junc2.from_pos, junc2.to_pos);
      }
      if (junc2.from.first != c2.fit.cluster_idx ||
          (junc2.type == sketching::JunctionType::X &&
           !is_roughly_end(s2, junc2.from.second)) ||
          (junc2.type == sketching::JunctionType::T && junc2.from.second != 0 &&
           junc2.from.second != 1)) {
        continue;
      }

      // Not far away enough
      if (to_distance_check) {
        double ratio = stroke_stroke_distance_ratio;
        if (c1.strokes.size() == 1 && c2.strokes.size() == 1) {
          ratio = 1;
        }
        double to_length = total_length(c2fit[junc.to.first]->centerline);
        double to_width =
          *std::max_element(c2fit[junc.to.first]->widths.begin(),
                            c2fit[junc.to.first]->widths.end());
        if (junc.to.first == junc2.to.first &&
            std::abs(junc.to.second - junc2.to.second) * to_length <
              ratio * to_width)
          continue;
        if (glm::distance(junc.to_pos, junc2.to_pos) < ratio * to_width)
          continue;

        auto proj = sketching::Vec2::Empty();
        double s;
        auto dist1 = sketching::closest_point(
          s1, sketching::Vec2(junc2.from_pos.x, junc2.from_pos.y), proj, s);
        auto dist2 = sketching::closest_point(
          s2, sketching::Vec2(junc.from_pos.x, junc.from_pos.y), proj, s);
        double to_width1 =
          *std::max_element(c2fit[junc.from.first]->widths.begin(),
                            c2fit[junc.from.first]->widths.end());
        double to_width2 =
          *std::max_element(c2fit[junc2.from.first]->widths.begin(),
                            c2fit[junc2.from.first]->widths.end());
        if (std::min(dist1, dist2) < ratio * std::max(to_width1, to_width2))
          continue;
      }

      junc_side2 = junction_side(c2, merged_cluster, junc2, end_u2);
      if (junc_side2 * junc_side1 > 0 && (end_u2 - end_u1) * junc_side1 >= 0) {
        cur_decision_junctions.emplace_back(junc);
        cur_decision_junctions.emplace_back(junc2);
        break;
      }
    }
  }

  decision_junctions.insert(decision_junctions.end(),
                            cur_decision_junctions.begin(),
                            cur_decision_junctions.end());

  return !cur_decision_junctions.empty();
}
