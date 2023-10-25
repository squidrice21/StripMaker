#include "deform.h"

#include "bounding_box.h"
#include "closest.h"
#include "degrees.h"
#include "detail/util.h"
#include "incremental_util.h"
#include "intersect.h"

namespace sketching {

namespace {
int max_deformation_trials = 5;

bool similar_gaps(const Stroke &stroke1, const Float s_start, const Float s_end,
                  const Stroke &stroke2, const Float t_start, const Float t_end,
                  const Float tol) {
  // Let's not merge if the shared amount is too small and there is little
  // overlap.
  if (std::abs(s_end - s_start) < stroke1.pen_width() ||
      std::abs(t_end - t_start) < stroke2.pen_width()) {
    auto width_start =
      0.5 * (stroke1.width_at(s_start) + stroke2.width_at(t_start));
    auto width_end = 0.5 * (stroke1.width_at(s_end) + stroke2.width_at(t_end));
    if ((stroke1.pos(s_start) - stroke2.pos(t_start)).norm() > width_start ||
        (stroke1.pos(s_end) - stroke2.pos(t_end)).norm() > width_end) {
      return false;
    }
  }

  // Check lengths are similar.
  const auto ratio = std::abs(s_end - s_start) / std::abs(t_end - t_start);
  if (!(1 / 1.2 < ratio && ratio < 1.2)) {
    return false;
  }

  // Check angle.
  const auto v1 = (stroke1.pos(s_start) - stroke1.pos(s_end)).normalized();
  const auto v2 = (stroke2.pos(t_start) - stroke2.pos(t_end)).normalized();
  const auto angle = std::abs(std::acos(std::clamp(v1.dot(v2), -1.0, 1.0)));
  if (angle > 20_deg) {
    return false;
  }

  // Check distance between corresponding points isn't too large.
  const auto n_samples =
    std::min(10, (int)std::ceil(0.5 * (std::abs(s_end - s_start) +
                                       std::abs(t_end - t_start) / tol)));
  auto width1 = 0.0;
  auto width2 = 0.0;
  for (int i = 0; i <= n_samples; ++i) {
    const auto u = i / Float(n_samples);
    const auto s = clamped_lerp(s_start, s_end, u);
    const auto t = clamped_lerp(t_start, t_end, u);
    width1 = std::max(width1, stroke1.width_at(s));
    width2 = std::max(width2, stroke2.width_at(t));
  }
  const auto new_tol = 0.17 * std::max(width1, width2);
  for (int i = 1; i < n_samples; ++i) {
    const auto u = i / Float(n_samples);
    const auto s = clamped_lerp(s_start, s_end, u);
    const auto t = clamped_lerp(t_start, t_end, u);
    const auto w1 = stroke1.width_at(s);
    const auto w2 = stroke2.width_at(t);
    // Check if envelope distance is greater than a threshold.
    if ((stroke1.pos(s) - stroke2.pos(t)).norm() - 0.5 * (w1 + w2) >= new_tol) {
      return false;
    }
  }

  return true;
}

enum Easing : uint8_t {
  FadeIn = 1, // Start with 0 contribution from correction stroke.
  FadeOut = 2, // End with 0 contribution from correction stroke.
};

/// Merge the segment defined from vertices [src_arclen_start, src_arclen_end]
/// or [src_arclen_end, src_arclen_start] (depending on if `src_arclen_start <
/// src_arclen_end`) from the src stroke into the segment defined from vertices
/// [dst_start, dst_end] of the dst stroke.
void merge_runs(Stroke &dst, const Index dst_start, const Index dst_end,
                const Float dst_arclen_start, const Float dst_arclen_end,
                const Index offset, const Stroke &src,
                const Float src_arclen_start, const Float src_arclen_end,
                const uint8_t easing) {
  assert(src.has_arclengths());
  assert(dst.has_arclengths());
  assert(dst_arclen_start != dst_arclen_end);

  for (Index i = dst_start; i <= dst_end; ++i) {
    const auto sample_u = (dst.arclength(i - offset) - dst_arclen_start) /
                          (dst_arclen_end - dst_arclen_start);
    const auto sample_location =
      clamped_lerp(src_arclen_start, src_arclen_end, sample_u);
    const auto [j, u] = src.fractional_index(sample_location);
    const auto x2 = fast_lerp(src.x(j), src.x(j + 1), u);
    const auto y2 = fast_lerp(src.y(j), src.y(j + 1), u);
    const auto w2 = fast_lerp(src.width(j), src.width(j + 1), u);
    const auto d = (dst.xy(i) - Vec2(x2, y2)).norm();
    if (d < 1e-8) {
      dst.width(i) = std::max(dst.width(i), w2);
    } else if (d + 0.5 * dst.width(i) <= 0.5 * w2) {
      // src vertex completely covers dst.
      dst.width(i) = w2;
      dst.x(i) = x2;
      dst.y(i) = y2;
    } else if (d + 0.5 * w2 > 0.5 * dst.width(i)) {
      // Neither completely covers the other.
      const auto envelope_span = 0.5 * (dst.width(i) + w2) + d;
      if (envelope_span <= dst.width(i) + w2) {
        // We place a vertex in the middle of the "inked" area.
        dst.width(i) = envelope_span;
        const auto influence = 0.5 + (w2 - dst.width(i)) / (4 * d);
        dst.x(i) = fast_lerp(dst.x(i), x2, influence);
        dst.y(i) = fast_lerp(dst.y(i), y2, influence);
      } else {
        // There is a gap in-between the two strokes.
        Float blend_factor;
        if ((easing & FadeIn) && (easing & FadeOut)) {
          // Reflect function around u = 0.5.
          const auto x = (sample_u < 0.5 ? sample_u : 1 - sample_u);
          blend_factor = 3 * square(2 * x) - 2 * cube(2 * x);
        } else if (easing & FadeIn) {
          blend_factor = 3 * square(sample_u) - 2 * cube(sample_u);
        } else if (easing & FadeOut) {
          blend_factor = 3 * square(1 - sample_u) - 2 * cube(1 - sample_u);
        } else {
          std::abort(); // Bad argument.
        }
        // We don't allow the centerline to move to a point that would cause the
        // width to expand outside the envelope span.
        const auto influence =
          fast_lerp(0.5 * w2 / d, 1.0 - 0.5 * dst.width(i) / d, blend_factor);
        dst.width(i) = dst.width(i) + w2;
        dst.x(i) = fast_lerp(dst.x(i), x2, influence);
        dst.y(i) = fast_lerp(dst.y(i), y2, influence);
      }
    }
  }
}

} // namespace

void straight_line_connect(Stroke &stroke, const Vec2 &head_pos,
                           const Vec2 &tail_pos) {
  if ((head_pos - stroke.xy(0)).norm() > 1e-6) {
    stroke.insert(0, head_pos.x_, head_pos.y_, stroke.width(0));
  } else {
    stroke.x(0) = head_pos.x_;
    stroke.y(0) = head_pos.y_;
  }
  if ((tail_pos - stroke.xy(Back)).norm() > 1e-6) {
    stroke.push_back(tail_pos.x_, tail_pos.y_, stroke.width(Back));
  } else {
    stroke.x(Back) = tail_pos.x_;
    stroke.y(Back) = tail_pos.y_;
  }
  stroke.invalidate_arclengths();
}

void smart_deform(Stroke &stroke, const Vec2 &head_pos, const Vec2 &tail_pos) {
  auto succeeded = false;
#ifdef HAS_GUROBI
  succeeded = snap_endpoints_divided(stroke, head_pos, tail_pos, false);
#else
  static auto warned = false;
  if (!warned) {
    SPDLOG_WARN(
      "Gurobi not linked; falling back to smooth_deform for snapping");
    warned = true;
  }
#endif
  if (!succeeded) {
    smooth_deform(stroke, head_pos, tail_pos);
  }
}

bool non_intersecting_smart_deform(const StrokeGraph::VertexView v,
                                   Stroke &stroke, const Vec2 &head_pos,
                                   const Vec2 &tail_pos) {
  auto succeeded = false;
  std::vector<Vec2> out;

#ifdef HAS_GUROBI
  assert(v.is_valid());

  // 1. Try deform and test if any new intersection exists.
  Stroke varying_stroke;
  int deformation_size = -1;
  int decrease_step = -1;
  do {
    succeeded = false;
    varying_stroke = stroke.clone();
    bool forward;
    succeeded = snap_endpoints_adaptive(varying_stroke, head_pos, tail_pos,
                                        false, deformation_size, forward);
    assert(deformation_size < varying_stroke.size());
    // Determine the decrease step size
    if (decrease_step < 0) {
      decrease_step = std::max(1, deformation_size / max_deformation_trials);
    }
    if (succeeded) {
      varying_stroke.ensure_arclengths();
      const auto node1 =
        PolylineBVHLeaf(varying_stroke, bounds(varying_stroke));

      // Check for self-intersection
      {
        out.clear();
        intersect_self(node1, out);
        Float check_length = node1.geometry->length();
        for (const auto &out_v : out) {
          if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
            goto failure;
          }
        }
      }

      // Check against all strokes (except for this stroke itself)
      for (size_t i = 0; i < v.graph_->strokes_.size(); ++i) {
        if (&stroke == &v.graph_->strokes_[i] ||
            v.graph_->strokes_[i].size() == 0)
          continue;
        out.clear();
        const Stroke &s = v.graph_->strokes_[i];
        s.ensure_arclengths();
        Float check_length = s.length();
        const auto node2 = PolylineBVHLeaf(s, bounds(s));
        intersect_different(node1, node2, out);
        for (const auto &out_v : out) {
          if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
            succeeded = false;
            goto failure;
          }
        }
      }
    }

failure:
    if (!succeeded) {
      assert(decrease_step > 0);
      deformation_size = deformation_size - decrease_step;
    }
  } while (!succeeded && deformation_size > 0);

  if (deformation_size > 0) {
    stroke = std::move(varying_stroke);
    return true;
  }

  // 2. If fails, resolve to do the simple connection
  SPDLOG_WARN("Failed to adaptively deform a stroke");
  straight_line_connect(stroke, head_pos, tail_pos);
#else
  static auto warned = false;
  if (!warned) {
    SPDLOG_WARN(
      "Gurobi not linked; falling back to smooth_deform for snapping");
    warned = true;
  }
  if (!succeeded) {
    smooth_deform(stroke, head_pos, tail_pos);
  }
#endif

  // Check intersections again (necessary for high-valence vertex)
  auto node1 = PolylineBVHLeaf(stroke, bounds(stroke));
  node1.geometry->ensure_arclengths();

  for (size_t i = 0; i < v.graph_->strokes_.size(); ++i) {
    if (node1.geometry == &v.graph_->strokes_[i] ||
        v.graph_->strokes_[i].size() == 0)
      continue;
    out.clear();
    const Stroke &s = v.graph_->strokes_[i];
    s.ensure_arclengths();
    Float check_length = s.length();
    const auto node2 = PolylineBVHLeaf(s, bounds(s));
    intersect_different(node1, node2, out);
    for (const auto &out_v : out) {
      if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
        SPDLOG_WARN(
          "Failed to straight_line_connect a stroke without intersections");
        return false;
      }
    }
  }

  return true;
}

void smooth_deform(Stroke &s, const Vec2 &start_endpoint,
                   const Vec2 &end_endpoint, Float radius) {
  const auto n = s.size();
  assert(n > 1);
  if (n <= 2 || radius == 0.0 ||
      (s.xy(0).isApprox(start_endpoint) &&
       s.xy(n - 1).isApprox(end_endpoint))) {
    s.x(0) = start_endpoint.x();
    s.y(0) = start_endpoint.y();
    s.x(n - 1) = end_endpoint.x();
    s.y(n - 1) = end_endpoint.y();
  } else {
    const Vec2 disp_start = start_endpoint - s.xy(0);
    const Vec2 disp_end = end_endpoint - s.xy(n - 1);
    const auto l = s.length();
    if (radius < 0.0) {
      radius = std::min(10 * s.pen_width(), 0.5 * l);
    }
    radius = std::min(radius, l);
    s.ensure_arclengths();
    s.x(0) = start_endpoint.x();
    s.y(0) = start_endpoint.y();
    for (Index i = 1; i < n - 1; ++i) {
      const auto t_s = std::max(1.0 - s.arclength(i) / radius, 0.0);
      const auto t_e = std::max(1.0 - (l - s.arclength(i)) / radius, 0.0);
      const auto influence_s = 3 * square(t_s) - 2 * cube(t_s);
      const auto influence_e = 3 * square(t_e) - 2 * cube(t_e);
      const Vec2 new_p = influence_s * disp_start + influence_e * disp_end;
      s.x(i) += new_p.x();
      s.y(i) += new_p.y();
    }
    s.x(n - 1) = end_endpoint.x();
    s.y(n - 1) = end_endpoint.y();
  }
  s.invalidate_arclengths();
}

bool non_intersecting_smooth_deform(const StrokeGraph::VertexView v,
                                    Stroke &stroke, const Vec2 &start_endpoint,
                                    const Vec2 &end_endpoint,
                                    Float max_radius) {
  assert((stroke.xy(0) == start_endpoint || stroke.xy(Back) == end_endpoint) &&
         "two-endpoint simultaneous deform not implemented");

  const size_t si = &stroke - &v.graph_->strokes_[0];
  const auto forward = stroke.xy(0) != start_endpoint;
  if (max_radius < 0) {
    // Note this is a little different from the default radius in smooth_deform.
    const size_t limit_index = 0;
    max_radius = (forward ? stroke.arclength(limit_index)
                          : stroke.length() - stroke.arclength(limit_index));
    // Keep this line since it's faster to have a smaller range.
    max_radius = std::min(max_radius, 4 * stroke.pen_width());
  }
  const auto radius_step = max_radius / max_deformation_trials;

  auto intersections = std::vector<Vec2>();
  for (Index trial_idx = max_deformation_trials; trial_idx > 0; --trial_idx) {
    auto varying_stroke = stroke.clone();
    const auto radius = trial_idx * radius_step;
    smooth_deform(varying_stroke, start_endpoint, end_endpoint, radius);
    varying_stroke.compute_arclengths();

    // Check intersection against strokes adjacent to the vertex which this
    // substroke is snapped to.
    auto node1 = PolylineBVHLeaf(varying_stroke, bounds(varying_stroke));
    node1.geometry->ensure_arclengths();
    std::unordered_set<size_t> seen_strokes;

    // Check for self-intersection
    {
      intersections.clear();
      intersect_self(node1, intersections);
      Float check_length = node1.geometry->length();
      for (const auto &out_v : intersections) {
        if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
          goto failure;
        }
      }
    }

    // Check against all strokes (except for this stroke itself)
    for (size_t i = 0; i < v.graph_->strokes_.size(); ++i) {
      if (&stroke == &v.graph_->strokes_[i] ||
          v.graph_->strokes_[i].size() == 0)
        continue;
      intersections.clear();
      const Stroke &s = v.graph_->strokes_[i];
      s.ensure_arclengths();
      Float check_length = s.length();
      const auto node2 = PolylineBVHLeaf(s, bounds(s));
      intersect_different(node1, node2, intersections);
      for (const auto &out_v : intersections) {
        if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
          goto failure;
        }
      }
    }

    // Success.
    stroke = std::move(varying_stroke);
    return true;

failure:;
  }

  // SPDLOG_WARN("Failed to smooth_deform a stroke without intersections");
  straight_line_connect(stroke, start_endpoint, end_endpoint);

  // Check intersections again (necessary for high-valence vertex)
  auto node1 = PolylineBVHLeaf(stroke, bounds(stroke));
  node1.geometry->ensure_arclengths();

  for (size_t i = 0; i < v.graph_->strokes_.size(); ++i) {
    if (node1.geometry == &v.graph_->strokes_[i] ||
        v.graph_->strokes_[i].size() == 0)
      continue;
    intersections.clear();
    const Stroke &s = v.graph_->strokes_[i];
    s.ensure_arclengths();
    Float check_length = s.length();
    const auto node2 = PolylineBVHLeaf(s, bounds(s));
    intersect_different(node1, node2, intersections);
    for (const auto &out_v : intersections) {
      if (out_v.y() > 1e-5 && std::abs(check_length - out_v.y()) > 1e-5) {
        SPDLOG_WARN(
          "Failed to straight_line_connect a stroke without intersections");
        return false;
      }
    }
  }

  return true;
}

bool oversketch_deform(Stroke &stroke, const Stroke &correction) {
  if (stroke.size() <= 1 || correction.size() <= 1) {
    return false;
  }
  assert(&stroke != &correction && "not supported");
  assert(stroke.has_time() == correction.has_time());

  // There are two major cases, distinguishable by the endpoint projections:
  //
  // 1.    ────────────
  //              ─────────────
  //
  // 2.    ────────────────────
  //            ──────────

  // TODO: Avoid creating strokes where one vertex is very close to another.

  stroke.ensure_arclengths();
  correction.ensure_arclengths();
  const auto tol = std::max(stroke.pen_width(), correction.pen_width());
  auto _proj = Vec2::Empty();
  Float s1, s2, t1, t2;
  const auto d_s1 = closest_point(stroke, correction.xy(0), _proj, s1) - tol;
  const auto d_s2 = closest_point(stroke, correction.xy(Back), _proj, s2) - tol;
  const auto d_t1 = closest_point(correction, stroke.xy(0), _proj, t1) - tol;
  const auto d_t2 = closest_point(correction, stroke.xy(Back), _proj, t2) - tol;
  const auto n = stroke.size();
  if (d_s1 <= 0 && d_t2 <= 0 && s1 > 0 && t2 > 0 && //
      s1 < stroke.length() && t2 < correction.length() &&
      similar_gaps(stroke, s1, stroke.length(), correction, 0, t2, tol)) {
    //                   s1
    //     stroke  ───────────►
    // correction         ────────────►
    //                        t2
    const auto smax = stroke.length();
    const auto [i, u] = stroke.fractional_index(s1);
    auto [j, v] = correction.fractional_index(t2);
    j++;
    const auto new_size = n + (correction.size() - j);
    // Reserve instead of resize, because resize destroys the arclengths array.
    stroke.reserve(new_size);
    stroke.size_ = new_size;
    merge_runs(stroke, (Index)std::ceil(i + u), n - 1, s1, smax, 0, correction,
               0.0, t2, FadeIn);
    copy_run(stroke, n, new_size - 1, correction, j, correction.size() - 1);
    stroke.invalidate_arclengths();
    return true;
  } else if (d_s2 <= 0 && d_t2 <= 0 && s2 > 0 && t2 > 0 && //
             s2 < stroke.length() && t2 < correction.length() &&
             similar_gaps(stroke, s2, stroke.length(), correction,
                          correction.length(), t2, tol)) {
    //                   s2
    //     stroke  ───────────►
    // correction         ◄────────────
    //                        t2
    //
    //     result  ───────────────────►
    const auto smax = stroke.length();
    const auto [i, u] = stroke.fractional_index(s2);
    const auto [j, v] = correction.fractional_index(t2);
    const auto margin = (Index)std::ceil(j + v);
    const auto new_size = n + margin;
    stroke.reserve(new_size);
    stroke.size_ = new_size;
    merge_runs(stroke, (Index)std::ceil(i + u), n - 1, s2, smax, 0, correction,
               correction.length(), t2, FadeIn);
    copy_run(stroke, n, new_size - 1, correction, margin - 1, 0);
    stroke.invalidate_arclengths();
    return true;
  } else if (d_s1 <= 0 && d_t1 <= 0 && s1 > 0 && t1 > 0 && //
             s1 < stroke.length() && t1 < correction.length() &&
             similar_gaps(stroke, s1, 0, correction, 0, t1, tol)) {
    //                   s1
    //     stroke  ◄───────────
    // correction         ────────────►
    //                        t1
    //
    //     result  ◄───────────────────
    const auto [i, u] = stroke.fractional_index(s1);
    auto [j, v] = correction.fractional_index(t1);
    j++;
    const auto offset = correction.size() - j;
    const auto new_size = n + offset;
    stroke.reserve(new_size);
    stroke.size_ = new_size;
    const auto dst_start = 0 + offset;
    const auto dst_end = i + offset;
    copy_run(stroke, new_size - 1, new_size - n, stroke, n - 1, 0);
    merge_runs(stroke, dst_start, dst_end, 0.0, s1, offset, correction, t1, 0.0,
               FadeOut);
    copy_run(stroke, 0, dst_start - 1, correction, correction.size() - 1, j);
    stroke.invalidate_arclengths();
    return true;
  } else if (d_s2 <= 0 && d_t1 <= 0 && s2 > 0 && t1 > 0 && //
             s2 < stroke.length() && t1 < correction.length() &&
             similar_gaps(stroke, s2, 0, correction, correction.length(), t1,
                          tol)) {
    //                   s2
    //     stroke  ◄───────────
    // correction         ◄────────────
    //                        t1
    //
    //     result  ◄───────────────────
    const auto [i, u] = stroke.fractional_index(s2);
    const auto [j, v] = correction.fractional_index(t1);
    const auto offset = (Index)std::ceil(j + v);
    const auto new_size = n + offset;
    stroke.reserve(new_size);
    stroke.size_ = new_size;
    const auto dst_start = 0 + offset;
    const auto dst_end = i + offset;
    copy_run(stroke, new_size - 1, new_size - n, stroke, n - 1, 0);
    merge_runs(stroke, dst_start, dst_end, 0.0, s2, offset, correction, t1,
               correction.length(), FadeOut);
    copy_run(stroke, 0, dst_start - 1, correction, 0, dst_start - 1);
    stroke.invalidate_arclengths();
    return true;
  } else if (d_s1 <= 0 && d_s2 <= 0 && s1 > 0 && s2 < stroke.length() &&
             // If we cut an original stroke at a sharp turn, we want to be able
             // to consolidate the two pieces if applicable.
             (s1 < stroke.length() || t2 < 1e-6) &&
             (s2 > 0 || t1 >= correction.length() - 1e-6) &&
             similar_gaps(stroke, s1, s2, correction, 0, correction.length(),
                          tol)) {
    //                s1         s2
    //     stroke  ───────────────────►
    // correction      ──────────►
    //
    // OR
    //
    //                s2         s1
    //     stroke  ───────────────────►
    // correction      ◄──────────
    //                t1         t2
    if (s2 < s1) {
      std::swap(s1, s2);
    }
    const auto [start_i, u] = stroke.fractional_index(s1);
    const auto [stop_i, v] = stroke.fractional_index(s2);
    const auto dst_start = (Index)std::ceil(start_i + u);
    const auto dst_end = (Index)std::floor(stop_i + v);
    Float src_start_arclen = 0.0, src_end_arclen = correction.length();
    if (t2 < t1) {
      std::swap(src_start_arclen, src_end_arclen);
    }
    merge_runs(stroke, dst_start, dst_end, s1, s2, 0, correction,
               src_start_arclen, src_end_arclen, FadeIn | FadeOut);
    stroke.invalidate_arclengths();
    return true;
  } else if (d_t1 <= 0 && d_t2 <= 0 && t1 > 0 && t2 < correction.length() &&
             // If we cut an original stroke at a sharp turn, we want to be able
             // to consolidate the two pieces if applicable.
             (t1 < correction.length() || s2 < 1e-6) &&
             (t2 > 0 || s1 >= stroke.length() - 1e-6) &&
             similar_gaps(correction, t1, t2, stroke, 0, stroke.length(),
                          tol)) {
    //     stroke      ──────────►
    // correction  ───────────────────►
    //                t1         t2
    //
    // OR
    //
    //                s2         s1
    //     stroke      ──────────►
    // correction  ◄───────────────────
    //                t1         t2
    const auto reversed = s2 < s1;
    // @optimize
    const auto temp = std::move(stroke);
    stroke = correction.clone();
    const auto rtn = oversketch_deform(stroke, temp);
    if (rtn && reversed) {
      stroke.reverse();
    }
    return rtn;
  }
  return false;
}

} // namespace sketching
