#include "junction_features_impl.h"

#include "../busyness.h"
#include "../bvh.h"
#include "../closest.h"
#include "../detail/suppress_warning.h"
#include "../detail/util.h"
#include "../force_assert.h"
#include "../intersect.h"
#include "../junction_type.h"
#include "../stroke_graph_extra.h"

namespace sketching {

void JunctionFeature::init(const StrokeGraph &graph) {
  force_assert(graph.orig_bvh_);
  init(*graph.orig_bvh_);
};

} // namespace sketching

namespace sketching::features {

namespace {

Float stepaway_projection_dist(const Stroke &s1, const bool s1_head,
                               const Float stepaway_amount, const Stroke &s2,
                               const Float connection_arclen2) {
  const auto stepaway_location =
    (s1_head ? std::min(stepaway_amount, s1.length())
             : std::max(s1.length() - stepaway_amount, 0.0));
  const auto stepaway = s1.pos(stepaway_location);
  constexpr auto near_factor = 2.0;
  // Start and end of search range.
  auto start =
    std::max(connection_arclen2 - near_factor * stepaway_amount, 0.0);
  auto stop =
    std::min(connection_arclen2 + near_factor * stepaway_amount, s2.length());

  auto proj = Vec2();
  auto s = Float(NAN);
  if (start == stop) {
    double small_stepaway_amount = 0.1;
    start =
      std::max(connection_arclen2 - near_factor * small_stepaway_amount, 0.0);
    stop = std::min(connection_arclen2 + near_factor * small_stepaway_amount,
                    s2.length());
  }
  force_assert(start < stop);
  const auto dist = closest_point_substroke(s2, start, stop, stepaway, proj, s);
  force_assert(std::isfinite(dist));
  return dist;
}

Float extra_s2_length(const StrokeGraph &graph, const Stroke &s2, Float arclen2,
                      size_t idx2, size_t he_idx) {
  s2.ensure_arclengths();
  StrokeTime interior_st((int)idx2, arclen2 / s2.length());
  const auto ok = convert_orig2strokes(graph, interior_st);
  force_assert(ok && "couldn't map from orig to strokes");
  // We can only reach this point if it's a corner (aka a high-valence vertex)
  if (!(interior_st.second == 0.0 || interior_st.second == 1.0)) {
    SPDLOG_DEBUG("Mapped position {} on stroke {} not at vertex",
                 interior_st.second, idx2);
  }

  // Round instead
  interior_st.second = (interior_st.second < 0.5) ? 0.0 : 1.0;

  // Get the sum of the part containing the shared vertex
  Endpoint end(interior_st.first, interior_st.second == 0.0);
  auto v = endpoint_to_vertex(graph, end);
  if (!v.is_active() || !v.is_valid())
    abort(); // Or assert here?

  Float dist_to_v = 0;
  const auto he = v.hedge();
  auto it = he;
  auto valence = 0;
  bool seen_he = false;
  do {
    if (it.index_ == he_idx) {
      seen_he = true;
      StrokeTime it_other_end(it.stroke_idx(), it.forward());
      const auto ok2 = convert_strokes2orig(graph, it_other_end);
      force_assert(ok2 && "couldn't map from strokes to orig");
      force_assert(it_other_end.first == idx2);

      dist_to_v = (it_other_end.second > (arclen2 / s2.length()))
                    ? arclen2
                    : s2.length() - arclen2;
      break;
    }

    it = it.twin().next();
    valence++;
    force_assert(valence < 1024 && "likely infinite loop found");
  } while (it != he);

  assert(seen_he);

  return dist_to_v;
}

} // namespace

Vec2 stepaway_tangent(const Stroke &s1, const bool s1_head,
                      const Float stepaway_amount) {
  s1.ensure_arclengths();
  const auto p = (s1_head ? s1.xy(0) : s1.xy(Back));
  auto furthest_sq_dist = -infinity;
  auto furthest_arclen = 0.0;
  if (s1_head) {
    for (Index i = 1; i < s1.size(); ++i) {
      const auto sq_dist = (s1.xy(i) - p).squaredNorm();
      if (sq_dist > furthest_sq_dist) {
        furthest_sq_dist = sq_dist;
        furthest_arclen = s1.arclength(i);
      }
    }
  } else {
    for (Index i = s1.size() - 2; i >= 0; --i) {
      const auto sq_dist = (s1.xy(i) - p).squaredNorm();
      if (sq_dist > furthest_sq_dist) {
        furthest_sq_dist = sq_dist;
        furthest_arclen = s1.arclength(i);
      }
    }
  }

  const auto stepaway =
    (s1_head
       ? s1.pos(std::min(stepaway_amount, furthest_arclen))
       : s1.pos(std::max(s1.length() - stepaway_amount, furthest_arclen)));
  return (p - stepaway).normalized();
}

std::string description(Normalization n) {
  switch (n) {
    case Normalization::BoundingBox:
      return "bounding box";
    case Normalization::PenWidth1:
      return "1st pen width";
    case Normalization::PenWidth2:
      return "2nd pen width";
    case Normalization::PenWidthPairwiseMax:
      return "max pen width of pair";
    case Normalization::PenWidthPairwiseMean:
      return "mean pen width of pair";
    case Normalization::StrokeLength1:
      return "1st stroke length";
    case Normalization::StrokeLength2:
      return "2nd stroke length";
    case Normalization::StrokeLengthPairwiseMin:
      return "min stroke length of pair";
    case Normalization::StrokeLengthPairwiseMax:
      return "max stroke length of pair";
    case Normalization::StrokeLengthPairwiseMean:
      return "mean stroke length of pair";
    default:
      std::abort();
  }
}

Float AbsEnvelopeDistance::operator()( //
  const Stroke &s1, const Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, size_t /*idx2*/) const {

  const auto centerline_dist = (s1.pos(arclen1) - s2.pos(arclen2)).norm();
  const auto iv1 = s1.fractional_index(arclen1);
  const auto iv2 = s2.fractional_index(arclen2);
  // TODO: This distance computation is very naive.
  const auto w1 =
    fast_lerp(s1.width(iv1.first), s1.width(iv1.first + 1), iv1.second);
  const auto w2 =
    fast_lerp(s2.width(iv2.first), s2.width(iv2.first + 1), iv2.second);
  return (centerline_dist - 0.5 * (w1 + w2));
}

Float AbsProjectionToClosestEndp1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float /*arclen2*/, const size_t /*idx2*/) const {

  Float s;
  smart_projection_dist(s1, arclen1, s2, s);
  return std::min(s, s2.length() - s);
}

Float AbsProjectionDist1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float /*arclen2*/, const size_t /*idx2*/) const {

  Float s;
  return smart_projection_dist(s1, arclen1, s2, s);
}

Float AbsStepawayDist1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();
  return stepaway_projection_dist(s1, /*s1_head=*/arclen1 < 0.5 * s1.length(),
                                  /*stepaway_amount=*/centerline_dist, s2,
                                  arclen2);
}

Float EndEndJunctionType::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t /*idx2*/) const {

  return end_end_junction_type(s1, arclen1, s2, arclen2);
}

Float EndStrokeJunctionType::operator()( //
  const Stroke & /*s1*/, const Float /*arclen1*/, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t /*idx2*/) const {

  return (Float)end_stroke_junction_type(s2, arclen2);
}

void EnvelopeDistance::init(const PolylineBVH &bvh) {
  if (m_normalization_scheme == Normalization::BoundingBox)
    m_inv_scale = 1.0 / bvh.bb.largest_axis_length();
  else
    m_inv_scale = std::numeric_limits<Float>::quiet_NaN();
}

void EnvelopeDistance::init(const StrokeGraph &graph) {
  graph_ptr_ = &graph;
  if (m_normalization_scheme == Normalization::BoundingBox)
    m_inv_scale = 1.0 / graph.orig_bvh_->bb.largest_axis_length();
  else
    m_inv_scale = std::numeric_limits<Float>::quiet_NaN();
}

Float EnvelopeDistance::operator()(const Stroke &s1, const Float arclen1,
                                   const size_t, const Stroke &s2,
                                   const Float arclen2, const size_t) const {
  const auto centerline_dist = (s1.pos(arclen1) - s2.pos(arclen2)).norm();
  const auto iv1 = s1.fractional_index(arclen1);
  const auto iv2 = s2.fractional_index(arclen2);
  // TODO: This distance computation is very naive.
  const auto w1 =
    fast_lerp(s1.width(iv1.first), s1.width(iv1.first + 1), iv1.second);
  const auto w2 =
    fast_lerp(s2.width(iv2.first), s2.width(iv2.first + 1), iv2.second);
  const auto abs_dist = (centerline_dist - 0.5 * (w1 + w2));
  switch (m_normalization_scheme) {
    case Normalization::PenWidth1: {
      return abs_dist / s1.pen_width();
    }
    case Normalization::PenWidth2: {
      return abs_dist / s2.pen_width();
    }
    case Normalization::PenWidthPairwiseMax: {
      return abs_dist / std::max(s1.pen_width(), s2.pen_width());
    }
    case Normalization::PenWidthPairwiseMean: {
      double mean_width = (0.5 * (s1.pen_width() + s2.pen_width()));
      if (graph_ptr_ && !graph_ptr_->orig_bvh_->combine_he_indices.empty()) {
        std::vector<std::pair<size_t, bool>> adjacent_stroke_indices;
        for (const auto he : graph_ptr_->orig_bvh_->combine_he_indices) {
          get_forward_chain_original_indices(*graph_ptr_, he,
                                             adjacent_stroke_indices);
        }
        double max_width = 0;
        for (const auto [sid, forward] : adjacent_stroke_indices) {
          max_width =
            std::max(max_width, graph_ptr_->orig_strokes_[sid].pen_width());
        }
        mean_width = (0.5 * (max_width + s1.pen_width()));
      }
      return abs_dist / mean_width;
    }
    case Normalization::StrokeLength1: {
      return abs_dist / s1.length();
    }
    case Normalization::StrokeLength2: {
      double length = s2.length();
      if (graph_ptr_ && !graph_ptr_->orig_bvh_->combine_he_indices.empty()) {
        length = 0;
        Float dist_to_v = 0;
        std::vector<std::pair<size_t, bool>> adjacent_stroke_indices;
        for (const auto he : graph_ptr_->orig_bvh_->combine_he_indices) {
          get_forward_chain_original_indices(*graph_ptr_, he,
                                             adjacent_stroke_indices);

          auto it = graph_ptr_->hedge(he);
          StrokeTime st((int)it.stroke_idx(), !it.forward());
          const auto ok = convert_strokes2orig(*graph_ptr_, st);
          force_assert(ok && "couldn't map from strokes to orig");

          // Subtract the section not adjacent to the dangling part
          auto &s = graph_ptr_->orig_strokes_[st.first];
          s.ensure_arclengths();
          dist_to_v += extra_s2_length(*graph_ptr_, s, s.length() * st.second,
                                       st.first, he);
        }

        std::unordered_set<size_t> seen_sid;
        for (const auto [sid, forward] : adjacent_stroke_indices) {
          if (seen_sid.count(sid))
            continue;
          seen_sid.emplace(sid);
          length += graph_ptr_->orig_bvh_->nodes[sid].geometry->length();
        }

        if (adjacent_stroke_indices.size() > 1)
          length -= dist_to_v;
      }
      return abs_dist / length;
    }
    case Normalization::StrokeLengthPairwiseMin: {
      return abs_dist / std::min(s1.length(), s2.length());
    }
    case Normalization::StrokeLengthPairwiseMax: {
      return abs_dist / std::max(s1.length(), s2.length());
    }
    case Normalization::StrokeLengthPairwiseMean: {
      return abs_dist / (0.5 * (s1.length() + s2.length()));
    }
    case Normalization::BoundingBox: {
      assert(!std::isnan(m_inv_scale));
      return abs_dist * m_inv_scale;
    }
    default:
      std::abort();
  }
}

Float ProjectionOverConnection1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  Float proj = smart_projection_dist(s1, arclen1, s2) / (q - p).norm();
  return proj;
}

Float ClosestAnyOverConnection1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();
  auto closest_dist = centerline_dist;
  Vec2 _proj = Vec2::Empty();
  Float _s;
  for (const auto &node : bvh_.nodes) {
    // TODO: Should we do something special in the node.geometry == &s1 case?
    if (node.bb.distanceTo(p) < closest_dist && node.geometry != &s1) {
      if (node.geometry == &s1) {
        closest_dist = std::min(
          closest_dist, smart_projection_dist(s1, arclen1, s1, _proj, _s));
      } else {
        closest_dist =
          std::min(closest_dist, closest_point(*node.geometry, p, _proj, _s));
      }
    }
  }
  return closest_dist / centerline_dist;
}

Float ClosestAnyOtherOverConnection1::operator()( //
  const Stroke &s1, Float arclen1, size_t idx1, //
  const Stroke &s2, Float arclen2, size_t idx2) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();

  const auto closest_p = closest(s1, arclen1, idx1, s2, arclen2, idx2);
  auto closest_dist = (closest_p - p).norm();
  if (!std::isfinite(closest_dist)) {
    closest_dist = std::max(s1.length(), s2.length());
  }
  return closest_dist / centerline_dist;
}

Vec2 ClosestAnyOtherOverConnection1::closest( //
  const Stroke &s1, Float arclen1, size_t idx1, //
  const Stroke &s2, Float arclen2, size_t idx2) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);

  std::vector<std::pair<size_t, bool>> adjacent_stroke_indices;

  // Graph based computation
  if (!graph_ptr_) {
    adjacent_stroke_indices.emplace_back(idx1, true);
    adjacent_stroke_indices.emplace_back(idx2, true);
  } else { // Stroke based computation
    // The strokes and indices here are about original strokes
    auto find_adjacent_strokes = [&](const Stroke &s, size_t idx,
                                     Float arclen) {
      s.ensure_arclengths();
      StrokeTime st((int)idx, arclen / s.length());
      const auto ok = convert_orig2strokes(*graph_ptr_, st);
      if (!ok)
        return;
      std::unordered_set<size_t> he_frontiers;

      if (st.second == 0.0 || st.second == 1.0) {
        // Get all adjacent half edges of this vertex
        Endpoint end(st.first, st.second == 0.0);
        auto v = endpoint_to_vertex(*graph_ptr_, end);
        if (!v.is_active() || !v.is_valid())
          return; // Or assert here?

        const auto he = v.hedge();
        auto it = he;
        auto valence = 0;
        do {
          if (!(it.flags() & StrokeGraph::HedgeRecord::Bridge)) {
            he_frontiers.emplace(it.index_);
          }
          it = it.twin().next();
          valence++;
          force_assert(valence < 1024 && "likely infinite loop found");
        } while (it != he);
      } else {
        // Get both directions if the end is in the interior
        const auto he = graph_ptr_->hedge(st.first * 2);
        he_frontiers.emplace(he.index_);
        he_frontiers.emplace(he.twin().index_);
      }

      for (const auto he : he_frontiers)
        get_forward_chain_original_indices(*graph_ptr_, he,
                                           adjacent_stroke_indices);
    };
    find_adjacent_strokes(s1, idx1, arclen1);
    find_adjacent_strokes(s2, idx2, arclen2);

    // Avoid the masked original strokes (the forward or not doesn't matter
    // here)
    for (const auto orig_sid : graph_ptr_->orig_bvh_->masked_nodes)
      adjacent_stroke_indices.emplace_back(orig_sid, true);
  }

  auto closest_dist = infinity;
  auto closest_p = Vec2(infinity, infinity);
  Vec2 proj = Vec2::Empty();
  Float _s = 0, _u = 0, _v = 0;
  force_assert((!graph_ptr_ || graph_ptr_->orig_bvh_) &&
               "Need to initialize orig_bvh_.");
  const auto &bvh = (graph_ptr_ ? *graph_ptr_->orig_bvh_ : bvh_);
  const auto blocking_bvh_node =
    ((idx2 == (size_t)-1 || idx2 >= bvh.nodes.size())
       ? PolylineBVHLeaf(s2, bounds(s2))
       : bvh.nodes[idx2]);
  for (size_t i = 0; i < bvh.nodes.size(); ++i) {
    const auto &node = bvh.nodes[i];
    // TODO: Should we do something special in the node.geometry == &s1 case?
    if (node.bb.distanceTo(p) < closest_dist && //
        std::find_if(adjacent_stroke_indices.begin(),
                     adjacent_stroke_indices.end(),
                     [=](const std::pair<size_t, bool> &p) {
                       return p.first == i;
                     }) == adjacent_stroke_indices.end()) {
      node.geometry->ensure_arclengths();
      const auto dist = closest_point(*node.geometry, p, proj, _s);
      if (dist < closest_dist &&

          // At a graph vertex, we want to avoid having the closest distance be
          // 0.
          (proj - p).squaredNorm() > 1e-10 &&
          (proj - q).squaredNorm() > 1e-10) {

        if (!(limit_to_visible_ && intersect_segment_stroke_exclusive(
                                     p, proj, blocking_bvh_node, _u, _v))) {
          closest_dist = dist;
          closest_p = proj;
        } else {
          const auto n_samples = 20;
          for (int k = 0; k < n_samples; ++k) {
            const auto other_arclen =
              ((Float)k / (n_samples - 1)) * node.geometry->length();
            MSVC_WARNING_SUPPRESS(
              4456); // declaration hides previous local declaration
            const auto proj = node.geometry->pos(other_arclen);
            MSVC_WARNING_SUPPRESS(
              4456); // declaration hides previous local declaration
            const auto dist = (proj - p).norm();
            if (dist < closest_dist && !intersect_segment_stroke_exclusive(
                                         p, proj, blocking_bvh_node, _u, _v)) {
              closest_dist = dist;
              closest_p = proj;
            }
          }
        }
      }
    }
  }
  return closest_p;
}

Float OtherEndpointClosestAnyEnvOverEnvConnection1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();
  const auto iv1 = s1.fractional_index(arclen1);
  const auto iv2 = s2.fractional_index(arclen2);
  // TODO: This distance computation is very naive.
  const auto w1 =
    fast_lerp(s1.width(iv1.first), s1.width(iv1.first + 1), iv1.second);
  const auto w2 =
    fast_lerp(s2.width(iv2.first), s2.width(iv2.first + 1), iv2.second);
  const auto conn_env_dist = (centerline_dist - 0.5 * (w1 + w2));
  // TODO: Re-enable this one we re-snap after hook removal.
  // assert(conn_env_dist > 0.0);

  const auto other_end_index =
    (arclen1 < 0.5 * s1.length() ? s1.size() - 1 : 0);
  const auto other_end = s1.xy(other_end_index);
  const auto other_end_radius = 0.5 * s1.width(other_end_index);
  // Find the closest point on any stroke to other_end.
  auto closest_dist = infinity;
  for (const auto &node : bvh_.nodes) {
    // TODO: Should we do something special in the node.geometry == &s1 case?
    if (node.bb.distanceTo(other_end) < closest_dist && node.geometry != &s1) {
      closest_dist =
        std::min(closest_dist, signed_dist_stroke(other_end, *node.geometry) -
                                 other_end_radius);
      if (closest_dist <= 0.0) {
        return 0.0;
      }
    }
  }

  return closest_dist / conn_env_dist;
}

Float ClosestEndpointOverConnection1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t idx1, //
  const Stroke &s2, const Float arclen2, const size_t idx2) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();

  const auto closest_p = closest(s1, arclen1, idx1, s2, arclen2, idx2);
  auto closest_dist = (closest_p - p).norm();
  if (!std::isfinite(closest_dist)) {
    closest_dist = std::max(s1.length(), s2.length());
  }
  return closest_dist / centerline_dist;
}

Vec2 ClosestEndpointOverConnection1::closest( //
  const Stroke &s1, const Float arclen1, const size_t idx1, //
  const Stroke & /*s2*/, const Float /*arclen2*/, const size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  auto closest_dist = Float(INFINITY);
  auto closest_p = Vec2(INFINITY, INFINITY);
  for (size_t si = 0; si < bvh_.nodes.size(); ++si) {
    if (bvh_.nodes[si].geometry->size() <= 1) {
      continue;
    }

    for (const auto head : {true, false}) {
      const auto &other_stroke = *bvh_.nodes[si].geometry;
      const auto other_p = other_stroke.xy(head ? 0 : other_stroke.size() - 1);
      const auto dist = (p - other_p).norm();
      if (dist <= 1e-10) {
        continue; // Ignore own endpoint.
      }
      if (dist < closest_dist) {
        Float _u = 0, _v = 0;
        auto has_line_of_sight = true;
        for (size_t other_si = 0; other_si < bvh_.nodes.size(); ++other_si) {
          if (other_si != idx1 && intersect_segment_stroke_exclusive(
                                    p, other_p, bvh_.nodes[other_si], _u, _v)) {
            has_line_of_sight = false;
            break;
          }
        }
        if (has_line_of_sight) {
          closest_dist = dist;
          closest_p = other_p;
        }
      }
    }
  }
  return closest_p;
}

Float StepawayOverConnection1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {
  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();
  auto stepaway_proj_dist =
    stepaway_projection_dist(s1, /*s1_head=*/arclen1 < 0.5 * s1.length(),
                             /*stepaway_amount=*/centerline_dist, s2, arclen2);
  std::unordered_set<size_t> seen_sid;
  if (graph_ptr_) {
    for (const auto he : graph_ptr_->orig_bvh_->combine_he_indices) {
      StrokeTime orig_st(graph_ptr_->hedge(he).stroke_idx(),
                         (graph_ptr_->hedge(he).forward()) ? 0 : 1);
      const auto ok = convert_strokes2orig(*graph_ptr_, orig_st);
      force_assert(ok && "couldn't map from strokes to orig");

      size_t sid = orig_st.first;
      if (seen_sid.count(sid))
        continue;
      seen_sid.emplace(sid);
      auto proj_dist = stepaway_projection_dist(
        s1, /*s1_head=*/arclen1 < 0.5 * s1.length(),
        /*stepaway_amount=*/centerline_dist, graph_ptr_->orig_strokes_[sid],
        orig_st.second);
      stepaway_proj_dist = std::min(stepaway_proj_dist, proj_dist);
    }
  }

  return (stepaway_proj_dist / centerline_dist);
}

Float NearestEndpointOverStepaway1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();
  const auto stepaway_proj_dist =
    stepaway_projection_dist(s1, /*s1_head=*/arclen1 < 0.5 * s1.length(),
                             /*stepaway_amount=*/centerline_dist, s2, arclen2);
  assert(stepaway_proj_dist != 0.0 && "will result in division by 0");
  Float s;
  smart_projection_dist(s1, arclen1, s2, s);
  const auto dist_nearest_endp = std::min(s, s2.length() - s);
  return dist_nearest_endp / stepaway_proj_dist;
}

Float StepawayOverProjection1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {

  const auto p = s1.pos(arclen1);
  const auto q = s2.pos(arclen2);
  const auto centerline_dist = (p - q).norm();
  const auto stepaway_proj_dist =
    stepaway_projection_dist(s1, /*s1_head=*/arclen1 < 0.5 * s1.length(),
                             /*stepaway_amount=*/centerline_dist, s2, arclen2);
  return stepaway_proj_dist / smart_projection_dist(s1, arclen1, s2);
}

void TangentAngle::init(const PolylineBVH &bvh) {
  const auto n = bvh.nodes.size();
  m_tangents.resize(n, Eigen::NoChange);
  // TODO: This is really inefficient; we only need to compute tangents for
  // strokes that
  //       could possibly need it.
  for (auto i = decltype(n){0}; i < n; ++i) {
    const auto &stroke = *bvh.nodes[i].geometry;
    std::vector<Bezier> bezier;
    if (stroke.size() <= 1) {
      m_tangents.row(i).fill(0.0);
    } else if (stroke.size() < 21) {
      fit_bezier_spline(stroke, 3.0, bezier);
      m_tangents.block(i, 0, 1, 2) =
        to_eigen(bezier[0].normalized_head_tangent()).transpose();
      m_tangents.block(i, 2, 1, 2) =
        to_eigen(bezier.back().normalized_tail_tangent()).transpose();
    } else {
      fit_bezier_spline(ConstStrokeView(stroke, 0, 10), 3.0, bezier);
      m_tangents.block(i, 0, 1, 2) =
        to_eigen(bezier[0].normalized_head_tangent()).transpose();
      bezier.clear();
      fit_bezier_spline(
        ConstStrokeView(stroke, stroke.size() - 10, stroke.size()), 3.0,
        bezier);
      m_tangents.block(i, 2, 1, 2) =
        to_eigen(bezier.back().normalized_tail_tangent()).transpose();
    }
  }
}

Float InteriorTangentAngle::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t idx2) const {

  const auto p1 = s1.pos(arclen1);
  const auto p2 = s2.pos(arclen2);
  if ((p1 - p2).squaredNorm() < 1e-10) {
    // I assume that a value of 90deg is the "best" value for a T-junction.
    return M_PI_2;
  } else {
    const auto connection = (p2 - p1).normalized();
    const auto &fit = smooth_fit(s2, idx2);
    // TODO: Use corner_tangents always.
    auto avg_tangent = Vec2::Empty();
    if (idx2 == (size_t)-1) {
      const auto [tangent1, tangent2] = fit.corner_tangents(arclen2);
      // Use average of tangents if we are at a corner.
      avg_tangent = (0.5 * (tangent1 + tangent2)).normalized();
    } else {
      avg_tangent = fit.tangent(arclen2);
    }
    const auto angle =
      std::acos(std::clamp(connection.dot(avg_tangent), -1.0, 1.0));
    assert(!std::isnan(angle));
    // Unlike endpoint tangents, the interior tangent angle doesn't have a
    // specific orientation.
    return std::min(angle, M_PI - angle);
  }
}

Vec2 InteriorTangentAngle::head_tangent(const size_t stroke_idx) const {
  const auto &fit = fits_[stroke_idx + 1];
  assert(fit.n_segments_ > 0);
  return -fit.tangent(0.0);
}

Vec2 InteriorTangentAngle::tail_tangent(const size_t stroke_idx) const {
  const auto &fit = fits_[stroke_idx + 1];
  assert(fit.n_segments_ > 0);
  return fit.tangent(fit.length());
}

const BezierSpline &
InteriorTangentAngle::smooth_fit(const Stroke &stroke,
                                 const size_t stroke_idx) const {
  if (stroke_idx == (size_t)-1) {
    fits_[0] = fit_bezier_spline_with_corners(stroke, 0.5 * stroke.pen_width());
    return fits_[0];
  }
  if (fits_[stroke_idx + 1].n_segments_ == 0) {
    fits_[stroke_idx + 1] =
      fit_bezier_spline_with_corners(stroke, 0.5 * stroke.pen_width());
  }
  return fits_[stroke_idx + 1];
}

Float TangentAngle1::operator()(const Stroke &s1, Float arclen1, size_t idx1,
                                const Stroke &s2, Float arclen2, size_t) const {
  const Vec2 p1 = s1.pos(arclen1);
  const Vec2 p2 = s2.pos(arclen2);
  if ((p1 - p2).squaredNorm() < 1e-10) {
    return 0.0;
  } else {
    const Vec2 v = (p2 - p1).normalized();
    const auto tangent = from_eigen(
      m_tangents.block(idx1, (arclen1 < 0.5 * s1.length() ? 0 : 2), 1, 2)
        .transpose());
    const auto angle = std::acos(std::clamp(v.dot(tangent), -1.0, 1.0));
    assert(!std::isnan(angle));
    return angle;
  }
}

Float TangentAngle2::operator()(const Stroke &s1, Float arclen1,
                                size_t /*idx1*/, const Stroke &s2,
                                Float arclen2, size_t idx2) const {
  const Vec2 p1 = s1.pos(arclen1);
  const Vec2 p2 = s2.pos(arclen2);
  if ((p1 - p2).squaredNorm() < 1e-10) {
    return 0.0;
  } else {
    const Vec2 v = (p1 - p2).normalized();
    const auto tangent = from_eigen(
      m_tangents.block(idx2, (arclen2 < 0.5 * s2.length() ? 0 : 2), 1, 2)
        .transpose());
    const auto angle = std::acos(std::clamp(v.dot(tangent), -1.0, 1.0));
    assert(!std::isnan(angle));
    return angle;
  }
}

Vec2 TangentAngle::head_tangent(const Index stroke_idx) const {
  return from_eigen(m_tangents.block(stroke_idx, 0, 1, 2).transpose());
}

Vec2 TangentAngle::tail_tangent(const Index stroke_idx) const {
  return from_eigen(m_tangents.block(stroke_idx, 2, 1, 2).transpose());
}

Float StepawayTangentAngle1::operator()( //
  const Stroke &s1, Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, Float arclen2, size_t /*idx2*/) const {

  const Vec2 p = s1.pos(arclen1);
  const Vec2 q = s2.pos(arclen2);
  if ((p - q).squaredNorm() < 1e-10) {
    return 0.0;
  } else {
    auto v = q - p;
    const auto centerline_dist = v.norm();
    v *= 1.0 / centerline_dist;
    const auto tangent =
      stepaway_tangent(s1, arclen1 < 0.5 * s1.length(), centerline_dist);
    const auto angle = std::acos(std::clamp(v.dot(tangent), -1.0, 1.0));
    assert(!std::isnan(angle));
    return angle;
  }
}

Float ClosestDistanceOnExtension1::operator()( //
  const Stroke &s1, const Float arclen1, size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t idx2) const {

  const Vec2 p = s1.pos(arclen1);
  const Vec2 q = s2.pos(arclen2);
  if ((p - q).squaredNorm() < 1e-10) {
    return 0.0;
  } else {
    // Compute tangent angle via stepaway.
    auto v = q - p;
    const auto centerline_dist = v.norm();
    v *= 1.0 / centerline_dist;
    const auto tangent =
      stepaway_tangent(s1, arclen1 < 0.5 * s1.length(), centerline_dist);

    // Find intersection if it exists.
    constexpr auto max_step = 1.5; // Or 1.0.
    Float s, t;
    if (intersect_segment_stroke_exclusive(
          p, p + max_step * centerline_dist * tangent, bvh_.nodes[idx2], s,
          t)) {
      return 0.0;
    }

    // Sample at a few points.
    const auto pw = s1.pen_width();
    auto closest_dist = smart_projection_dist(s1, arclen1, s2);
    auto proj = Vec2::Empty();
    for (const auto step : {0.5, 1.0, max_step}) {
      const auto new_dist =
        closest_point(bvh_.nodes[idx2], p + step * centerline_dist * tangent,
                      closest_dist, proj, s);
      if (&s1 != &s2 || std::abs(s - arclen1) > pw)
        closest_dist = std::min(closest_dist, new_dist);
    }
    return closest_dist / centerline_dist;
  }
}

void Busyness1::init(const PolylineBVH &bvh) {
  bvh_ = bvh;
  const auto n = bvh.nodes.size();
  const auto bufn = 3 * n;
  memory_.reset(new Float[bufn]);
  for (size_t i = 0; i < bufn; ++i) {
    memory_[i] = -1.0;
  }
  head_busyness_ = &memory_[0];
  tail_busyness_ = &memory_[n];
  avg_widths_ = &memory_[2 * n];
  average_widths(bvh, {avg_widths_, n});
}

Float Busyness1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t idx2) const {
  // Bypass cache.
  const auto p = s1.pos(arclen1);
  auto busyness =
    busyness_at(bvh_, p, {avg_widths_, bvh_.nodes.size()}, busyness_falloff_);
  if (idx2 == (size_t)-1) {
    // We want to pretend there is no endpoint at the dissolved stroke, so we
    // reverse the effects of that.
    busyness -= gaussian_factor(
      (p - s2.pos(arclen2)).norm() / average_width(s2), busyness_falloff_);
  }
  return busyness;
}

void PenWidth1::init(const PolylineBVH &bvh) {
  const auto n = bvh.nodes.size();
  auto acc = 0.0;
  auto denom = 0.0;
  for (auto i = decltype(n){0}; i < n; ++i) {
    const auto &stroke = *bvh.nodes[i].geometry;
    if (stroke.size() > 0) {
      const auto pw = stroke.pen_width();
      acc += pw * stroke.length();
      denom += stroke.length();
    }
  }
  drawing_inv_average_width_ = denom / acc;
  force_assert(acc > 0);
  force_assert(denom > 0);
  force_assert(!std::isnan(drawing_inv_average_width_));
}

Float ConnectedDistanceToEndpoint::operator()(const Stroke & /*s1*/,
                                              Float /*arclen1*/,
                                              size_t /*idx1*/, const Stroke &s2,
                                              Float arclen2,
                                              size_t idx2) const {
  if (!graph_ptr_ || graph_ptr_->orig_bvh_->combine_he_indices.empty()) {
    return std::min(arclen2, s2.length() - arclen2) / s2.length();
  }
  s2.ensure_arclengths();

  int he_idx = -1;

  Float length = 0;
  for (const auto he_i : graph_ptr_->orig_bvh_->combine_he_indices) {
    const auto it = graph_ptr_->hedge(he_i);
    StrokeTime end(it.stroke_idx(), it.forward());
    const auto ok2 = convert_strokes2orig(*graph_ptr_, end);
    if (end.first == idx2) {
      he_idx = he_i;
      break;
    }
  }
  assert(he_idx >= 0);
  {
    std::vector<std::pair<size_t, bool>> adjacent_stroke_indices;
    get_forward_chain_original_indices(*graph_ptr_, he_idx,
                                       adjacent_stroke_indices);
    std::unordered_set<size_t> seen_sid;
    for (const auto [sid, forward] : adjacent_stroke_indices) {
      if (seen_sid.count(sid))
        continue;
      seen_sid.emplace(sid);
      length += graph_ptr_->orig_bvh_->nodes[sid].geometry->length();
    }
  }

  Float other_length = 0;
  {
    std::vector<std::pair<size_t, bool>> adjacent_stroke_indices;
    Float dist_to_v = 0;
    for (const auto he : graph_ptr_->orig_bvh_->combine_he_indices) {
      if (he != he_idx) {
        get_forward_chain_original_indices(*graph_ptr_, he,
                                           adjacent_stroke_indices);

        // Subtract the section not adjacent to the dangling part
        auto it = graph_ptr_->hedge(he);
        StrokeTime st((int)it.stroke_idx(), !it.forward());
        const auto ok = convert_strokes2orig(*graph_ptr_, st);
        force_assert(ok && "couldn't map from strokes to orig");

        auto &s = graph_ptr_->orig_strokes_[st.first];
        s.ensure_arclengths();
        dist_to_v +=
          extra_s2_length(*graph_ptr_, s, s.length() * st.second, st.first, he);
      }
    }
    std::unordered_set<size_t> seen_sid;
    for (const auto [sid, forward] : adjacent_stroke_indices) {
      if (seen_sid.count(sid))
        continue;
      seen_sid.emplace(sid);
      other_length += graph_ptr_->orig_bvh_->nodes[sid].geometry->length();
    }

    other_length -= dist_to_v;
  }

  Float dist_to_v = extra_s2_length(*graph_ptr_, s2, arclen2, idx2, he_idx);
  Float ratio = (length - dist_to_v) / ((length - dist_to_v) + other_length);

  return std::min(ratio, 1 - ratio);
}

Float ConnectedLocationOverConnection::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t /*idx2*/) const {

  const auto centerline_dist = (s1.pos(arclen1) - s2.pos(arclen2)).norm();
  return std::min(arclen2, s2.length() - arclen2) / centerline_dist;
}

Float ProjectionToEndpointRatio1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float /*arclen2*/, const size_t /*idx2*/) const {

  Float proj_arclen;
  smart_projection_dist(s1, arclen1, s2, proj_arclen);
  return std::min(proj_arclen, s2.length() - proj_arclen) / s2.length();
}

Float ProjectionToEndpointOverConnection1::operator()( //
  const Stroke &s1, const Float arclen1, const size_t /*idx1*/, //
  const Stroke &s2, const Float arclen2, const size_t /*idx2*/) const {

  Float proj_arclen;
  smart_projection_dist(s1, arclen1, s2, proj_arclen);
  const auto centerline_dist = (s1.pos(arclen1) - s2.pos(arclen2)).norm();
  return std::min(proj_arclen, s2.length() - proj_arclen) / centerline_dist;
}

// FIXME: These drawing order features are broken if we use edge indices.

Float DrawingOrder::operator()( //
  const Stroke & /*s1*/, const Float arclen1, const size_t idx1, //
  const Stroke & /*s2*/, const Float /*arclen2*/, const size_t idx2) const {

  if (idx1 == idx2) {
    return arclen1 ==
           0.0; // Head trying to connect to tail. Head is drawn first.
  }
  return (idx1 < idx2 ? 0.0 : 1.0);
}

Float DrawingOrderUnsignedDiff::operator()( //
  const Stroke & /*s1*/, const Float /*arclen1*/, const size_t idx1, //
  const Stroke & /*s2*/, const Float /*arclen2*/, const size_t idx2) const {

  return (Float)std::abs((int64_t)idx2 - (int64_t)idx1);
}

Float DrawingOrderBucketedUnsignedDiff::operator()( //
  const Stroke & /*s1*/, const Float /*arclen1*/, const size_t idx1, //
  const Stroke & /*s2*/, const Float /*arclen2*/, const size_t idx2) const {

  const auto dt = std::abs((int64_t)idx2 - (int64_t)idx1);
  if (dt == 0) {
    return 0.0;
  } else if (dt == 1) {
    return 0.1;
  } else if (dt <= 5) {
    return 0.5;
  } else {
    return 1.0;
  }
}

} // namespace sketching::features
