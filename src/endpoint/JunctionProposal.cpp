#include "JunctionProposal.h"

#include "../solve/SolveUtil.h"
#include "Junction.h"
#include "SketchConnectivity/closest.h"
#include "SketchConnectivity/intersect.h"

#include "SketchConnectivity/3rdparty/span-lite/src/nonstd/span.hpp"
#include "SketchConnectivity/sketching.h"

double endpoint_distance_ratio = 3;

struct Hit {
  glm::dvec2 pos;
  sketching::StrokeTime st;
};

// Have to rewrite line_of_sight since span causes linking error...
bool line_of_sight_junc(
  sketching::Vec2 p, sketching::Vec2 q,
  nonstd::span<const sketching::Stroke> strokes,
  nonstd::span<const sketching::BoundingBox> centerline_bbs) {
  assert(strokes.size() == centerline_bbs.size());
  double s, t;
  const auto n = strokes.size();
  for (size_t i = 0; i < n; ++i) {
    if (intersect_segment_stroke_exclusive(
          p, q, {strokes[i], centerline_bbs[i]}, s, t)) {
      return false;
    }
  }
  return true;
}

bool is_roughly_end(const sketching::Stroke &s1, double t1) {
  double max_w = s1.pen_width();
  double arclen_i_mag = std::min(t1, 1 - t1);
  if (solve_context.stage == SolveContext::SolveStage::Endpoint)
    return (arclen_i_mag < 0.20) &&
           (arclen_i_mag * s1.length()) < endpoint_distance_ratio * max_w;
  return (arclen_i_mag < 0.15) &&
         (arclen_i_mag * s1.length()) < endpoint_distance_ratio * max_w;
}

void propose_end_end_intersections(
  size_t cur_c_idx, const sketching::StrokeBVH &bvh,
  std::vector<std::pair<sketching::StrokeTime, sketching::StrokeTime>>
    &out_intersections) {
  sketching::PolylineBVH poly_bvh = bvh.polyline_bvh();

  auto hits = std::multimap<double, Hit>(); // Sort by distance.

  auto to_glm = [](const sketching::Vec2 &vec) -> glm::dvec2 {
    return glm::dvec2(vec.x(), vec.y());
  };
  auto to_Vec = [](const glm::dvec2 &v) -> sketching::Vec2 {
    return sketching::Vec2(v.x, v.y);
  };

  sketching::PolylineBVHLeaf *cur_n_ptr = nullptr;
  for (auto &n : poly_bvh.nodes) {
    if (n.geometry->index == cur_c_idx)
      cur_n_ptr = &n;
  }
  assert(cur_n_ptr);
  cur_n_ptr->geometry->ensure_arclengths();

  const auto n = poly_bvh.nodes.size();
  for (size_t j = 0; j < n; ++j) {
    if (poly_bvh.nodes[j].geometry->index == cur_c_idx)
      continue;

    auto intersections = std::vector<sketching::Vec2>();
    intersect_different(*cur_n_ptr, poly_bvh.nodes[j], intersections);
    poly_bvh.nodes[j].geometry->ensure_arclengths();
    for (auto const &inter : intersections) {
      double arclen_i = inter.x() / cur_n_ptr->geometry->length();
      double arclen_j = inter.y() / poly_bvh.nodes[j].geometry->length();
      sketching::StrokeTime st1{cur_c_idx, arclen_i};
      sketching::StrokeTime st2{poly_bvh.nodes[j].geometry->index, arclen_j};
      if (is_roughly_end(*cur_n_ptr->geometry, arclen_i) ||
          is_roughly_end(*poly_bvh.nodes[j].geometry, arclen_j))
        out_intersections.emplace_back(
          (cur_c_idx < poly_bvh.nodes[j].geometry->index) ? st1 : st2,
          (cur_c_idx < poly_bvh.nodes[j].geometry->index) ? st2 : st1);
    }
  }
}

void propose_end_end_envelope_intersections(
  size_t cur_c_idx, const sketching::StrokeBVH &bvh,
  std::vector<std::pair<sketching::StrokeTime, sketching::StrokeTime>>
    &out_intersections) {
  sketching::PolylineBVH poly_bvh = bvh.polyline_bvh();

  auto hits = std::multimap<double, Hit>(); // Sort by distance.

  auto to_glm = [](const sketching::Vec2 &vec) -> glm::dvec2 {
    return glm::dvec2(vec.x(), vec.y());
  };
  auto to_Vec = [](const glm::dvec2 &v) -> sketching::Vec2 {
    return sketching::Vec2(v.x, v.y);
  };

  sketching::PolylineBVHLeaf *cur_n_ptr = nullptr;
  for (auto &n : poly_bvh.nodes) {
    if (n.geometry->index == cur_c_idx)
      cur_n_ptr = &n;
  }
  assert(cur_n_ptr);
  cur_n_ptr->geometry->ensure_arclengths();

  auto intersect_envelope = [&bvh](sketching::Vec2 xy, double w,
                                   const sketching::Stroke &stroke,
                                   double &arclen) {
    arclen = -1;
    stroke.ensure_arclengths();
    auto proj = sketching::Vec2::Empty();
    double s;
    auto dist = sketching::closest_point(stroke, xy, proj, s);
    if (dist < std::numeric_limits<double>::infinity()) {
      assert(s <= stroke.length());
      double time = s / stroke.length();

      // Invisible
      double vis_mag = (stroke.pos(s) - xy).norm() - 1e-5;
      if (!line_of_sight_junc(xy + (stroke.pos(s) - xy).normalized() * 1e-5,
                              xy + (stroke.pos(s) - xy).normalized() * vis_mag,
                              bvh.strokes(), bvh.centerline_bbs())) {
        return;
      }

      // Envelopes intersect
      if (dist < 0.5 * (w + stroke.width_at(s)))
        arclen = s;
    }
  };

  const auto n = poly_bvh.nodes.size();
  for (size_t j = 0; j < n; ++j) {
    if (poly_bvh.nodes[j].geometry->index == cur_c_idx)
      continue;

    poly_bvh.nodes[j].geometry->ensure_arclengths();
    auto intersections = std::vector<sketching::Vec2>();
    double arclen;
    intersect_envelope(cur_n_ptr->geometry->pos(0),
                       cur_n_ptr->geometry->width_at(0),
                       *poly_bvh.nodes[j].geometry, arclen);
    if (arclen >= 0) {
      intersections.emplace_back(sketching::Vec2(0, arclen));
    }

    intersect_envelope(
      cur_n_ptr->geometry->pos(cur_n_ptr->geometry->length()),
      cur_n_ptr->geometry->width_at(cur_n_ptr->geometry->length()),
      *poly_bvh.nodes[j].geometry, arclen);
    if (arclen >= 0) {
      intersections.emplace_back(
        sketching::Vec2(cur_n_ptr->geometry->length(), arclen));
    }

    for (auto const &inter : intersections) {
      double arclen_i = inter.x() / cur_n_ptr->geometry->length();
      double arclen_j = inter.y() / poly_bvh.nodes[j].geometry->length();
      sketching::StrokeTime st1{cur_c_idx, arclen_i};
      sketching::StrokeTime st2{poly_bvh.nodes[j].geometry->index, arclen_j};
      if (is_roughly_end(*cur_n_ptr->geometry, arclen_i) ||
          is_roughly_end(*poly_bvh.nodes[j].geometry, arclen_j))
        out_intersections.emplace_back(
          (cur_c_idx < poly_bvh.nodes[j].geometry->index) ? st1 : st2,
          (cur_c_idx < poly_bvh.nodes[j].geometry->index) ? st2 : st1);
    }
  }
}

void propose_end_end_candidates(const glm::dvec2 &v, size_t cur_c_idx,
                                size_t knn, sketching::StrokeBVH &bvh,
                                std::vector<sketching::StrokeTime> &candidates,
                                int to_c_idx) {
  const auto centerline_bbs = bvh.centerline_bbs();

  auto hits = std::multimap<double, Hit>(); // Sort by distance.

  auto to_glm = [](const sketching::Vec2 &vec) -> glm::dvec2 {
    return glm::dvec2(vec.x(), vec.y());
  };
  auto to_Vec = [](const glm::dvec2 &v) -> sketching::Vec2 {
    return sketching::Vec2(v.x, v.y);
  };
  sketching::Stroke *cur_s_ptr = nullptr;
  for (size_t i = 0; i < bvh.strokes().size(); ++i) {
    auto &s = bvh.stroke(i);
    if (s.index == cur_c_idx)
      cur_s_ptr = &s;
  }
  assert(cur_s_ptr);
  cur_s_ptr->ensure_arclengths();

  for (const auto &s : bvh.strokes()) {
    if (s.index == cur_c_idx || (to_c_idx >= 0 && s.index != to_c_idx))
      continue;

    auto cen_dist = glm::distance(v, to_glm(s.pos(0)));
    s.ensure_arclengths();

    // Close enough
    if (cen_dist <= std::min(std::min(s.length(), cur_s_ptr->length()), 15.0))
      hits.insert(
        {cen_dist, Hit{to_glm(s.pos(0)), sketching::StrokeTime{s.index, 0.0}}});

    cen_dist = glm::distance(v, to_glm(s.pos(s.length())));
    // Close enough
    if (cen_dist <= std::min(std::min(s.length(), cur_s_ptr->length()), 15.0))
      hits.insert({cen_dist, Hit{to_glm(s.pos(s.length())),
                                 sketching::StrokeTime{s.index, 1.0}}});
  }

  size_t i = 0;
  auto n_hits_looked_at = 0;
  for (const auto &hit : hits) {
    const auto other_vertex = hit.second.pos;
    auto v1 = to_Vec(v);
    auto v2 = to_Vec(other_vertex);

    if (line_of_sight_junc(v1, v2, bvh.strokes(), centerline_bbs)) {
      candidates.emplace_back(hit.second.st);
    }
    n_hits_looked_at++;
    if (n_hits_looked_at >= knn) {
      break;
    }
  }
}

void propose_end_stroke_candidates(
  const glm::dvec2 &v, size_t cur_c_idx, size_t knn,
  const sketching::StrokeBVH &bvh,
  std::vector<sketching::StrokeTime> &candidates, int to_c_idx) {
  const auto centerline_bbs = bvh.centerline_bbs();

  auto hits = std::multimap<double, Hit>(); // Sort by distance.

  auto to_glm = [](const sketching::Vec2 &vec) -> glm::dvec2 {
    return glm::dvec2(vec.x(), vec.y());
  };
  auto to_Vec = [](const glm::dvec2 &v) -> sketching::Vec2 {
    return sketching::Vec2(v.x, v.y);
  };
  auto xy = to_Vec(v);
  for (const auto &other_stroke : bvh.strokes()) {
    if (other_stroke.index == cur_c_idx ||
        (to_c_idx >= 0 && other_stroke.index != to_c_idx))
      continue;

    other_stroke.ensure_arclengths();
    auto proj = sketching::Vec2::Empty();
    double s;
    auto dist = sketching::closest_point(other_stroke, xy, proj, s);
    if (dist < std::numeric_limits<double>::infinity()) {
      assert(s <= other_stroke.length());
      double time = s / other_stroke.length();

      // Projects to an endpoint
      if (time == 0 || time == 1)
        continue;
      // Too close to an endpoint
      if (std::min(s, other_stroke.length() - s) <= other_stroke.width_at(s))
        continue;

      // Invisible
      double vis_mag = (other_stroke.pos(s) - xy).norm() - 1e-5;
      if (!line_of_sight_junc(
            xy + (other_stroke.pos(s) - xy).normalized() * 1e-5,
            xy + (other_stroke.pos(s) - xy).normalized() * vis_mag,
            bvh.strokes(), centerline_bbs)) {
        continue;
      }

      hits.insert({dist, Hit{to_glm(other_stroke.pos(time)),
                             sketching::StrokeTime{other_stroke.index, time}}});
    }
  }

  size_t i = 0;
  auto n_hits_looked_at = 0;
  for (const auto &hit : hits) {
    const auto other_vertex = hit.second.pos;
    auto v1 = to_Vec(v);
    auto v2 = to_Vec(other_vertex);
    candidates.emplace_back(hit.second.st);
    n_hits_looked_at++;
    if (n_hits_looked_at >= knn) {
      break;
    }
  }
}
