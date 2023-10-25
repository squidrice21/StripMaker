#include "render.h"

#include "bounding_box.h"
#include "detail/util.h"
#include "eigen_compat.h"
#include "force_assert.h"
#include "intersect.h"
#include "resample.h"
#include "stroke_graph_extra.h"

#include <spdlog/spdlog.h>

namespace sketching {

namespace {

#define ARRAY_SIZE(arr) (sizeof((arr)) / sizeof((arr)[0]))

struct Circle {
  Float x;
  Float y;
  Float r;
};

/// Find the points on circles c0 and c1 that lie on lines tangent to both.
void tangent_points(const Circle& c0, const Circle& c1, Vec2& p0_0, Vec2& p0_1,
                    Vec2& p1_0, Vec2& p1_1) {
  // See http://www.ambrsoft.com/TrigoCalc/Circles2/Circles2Tangent_.htm.
  const auto r0 = c0.r, r1 = c1.r;
  const auto r0r1 = std::abs(r1 - r0);
  if (r0r1 < 1e-4) {
    // Parallel tangent lines.
    const Vec2 centre0 = Vec2(c0.x, c0.y);
    const Vec2 centre1 = Vec2(c1.x, c1.y);
    const Vec2 tangent = centre1 - centre0;
    const Vec2 normal = Vec2(-tangent.y(), tangent.x()).normalized();
    p0_0 = centre0 + r0 * normal;
    p0_1 = centre0 - r0 * normal;
    p1_0 = centre1 + r1 * normal;
    p1_1 = centre1 - r1 * normal;
  } else {
    // Converging tangent lines.
    const auto a = c0.x;
    const auto b = c0.y;
    const auto c = c1.x;
    const auto d = c1.y;
    const auto x_p = (c * r0 - a * r1) / (r0 - r1);
    const auto y_p = (d * r0 - b * r1) / (r0 - r1);
    if (r0 > 1e-6) {
      const auto sqrt1 = std::sqrt(square(x_p - a) + square(y_p - b) - r0 * r0);
      p0_0.x() = (r0 * r0 * (x_p - a) + r0 * (y_p - b) * sqrt1) /
                   (square(x_p - a) + square(y_p - b)) +
                 a;
      p0_0.y() = (r0 * r0 * (y_p - b) - r0 * (x_p - a) * sqrt1) /
                   (square(x_p - a) + square(y_p - b)) +
                 b;
      p0_1.x() = (r0 * r0 * (x_p - a) - r0 * (y_p - b) * sqrt1) /
                   (square(x_p - a) + square(y_p - b)) +
                 a;
      p0_1.y() = (r0 * r0 * (y_p - b) + r0 * (x_p - a) * sqrt1) /
                   (square(x_p - a) + square(y_p - b)) +
                 b;
    } else {
      p0_0 = p0_1 = {c0.x, c0.y};
    }
    if (r1 > 1e-6) {
      const auto sqrt2 = std::sqrt(square(x_p - c) + square(y_p - d) - r1 * r1);
      p1_0.x() = (r1 * r1 * (x_p - c) + r1 * (y_p - d) * sqrt2) /
                   (square(x_p - c) + square(y_p - d)) +
                 c;
      p1_0.y() = (r1 * r1 * (y_p - d) - r1 * (x_p - c) * sqrt2) /
                   (square(x_p - c) + square(y_p - d)) +
                 d;
      p1_1.x() = (r1 * r1 * (x_p - c) - r1 * (y_p - d) * sqrt2) /
                   (square(x_p - c) + square(y_p - d)) +
                 c;
      p1_1.y() = (r1 * r1 * (y_p - d) + r1 * (x_p - c) * sqrt2) /
                   (square(x_p - c) + square(y_p - d)) +
                 d;
    } else {
      p1_0 = p1_1 = {c1.x, c1.y};
    }
    if (r1 < r0) {
      // Try to keep points stable across all configurations of r0, r1.
      std::swap(p0_0, p0_1);
      std::swap(p1_0, p1_1);
    }
  }
}

auto should_decimate = true;
auto cap_vertices = 7;

void partial_cap(const Vec2& c, const Float r, const Float start_theta,
                 const Float end_theta, std::vector<Vec2>& out_vertices) {
  const auto num_steps =
    (int)std::round(std::abs(end_theta - start_theta) * (cap_vertices / M_PI));
  const auto step = (end_theta - start_theta) / num_steps;
  for (auto i = 1; i < num_steps; ++i) {
    const auto t = start_theta + i * step;
    out_vertices.emplace_back(c.x() + r * std::cos(t), c.y() + r * std::sin(t));
  }
}
void partial_cap_ccw(const Vec2& c, const Float r, const Float start_theta,
                     Float end_theta, std::vector<Vec2>& out_vertices) {
  if (end_theta < start_theta) {
    end_theta += 2 * M_PI;
  }
  // For a noisy almost-straight line, we will sometimes turn "backwards" angle-wise.
  // However, we do not need a full ~360deg join here; a straight line connection will
  // suffice---hence the 1.99 * pi threshold.
  if (end_theta - start_theta <= 1.99 * M_PI) {
    partial_cap(c, r, start_theta, end_theta, out_vertices);
  }
}
void partial_cap_cw(const Vec2& c, const Float r, Float start_theta,
                    const Float end_theta, std::vector<Vec2>& out_vertices) {
  if (end_theta > start_theta) {
    start_theta += 2 * M_PI;
  }
  // For a noisy almost-straight line, we will sometimes turn "backwards" angle-wise.
  // However, we do not need a full ~360deg join here; a straight line connection will
  // suffice---hence the 1.99 * pi threshold.
  if (start_theta - end_theta <= 1.99 * M_PI) {
    partial_cap(c, r, start_theta, end_theta, out_vertices);
  }
}
void cap_full_circle(const Index v, const ConstStrokeView& s,
                     std::vector<CoordMat>& out_vertices) {
  const Vec2 p = s.xy(v);
  const Float w = 0.5 * s.width(v);
  auto& coords = out_vertices.emplace_back(2 * cap_vertices, 2);
  for (auto i = 0; i < 2 * cap_vertices; ++i) {
    const auto t = i * M_PI / cap_vertices;
    coords(i, 0) = p.x() + w * std::cos(t);
    coords(i, 1) = p.y() + w * std::sin(t);
  }
}

void finalize_polygon(std::vector<Vec2>& forward_vertices,
                      std::vector<Vec2>& reverse_vertices, const Float pen_width,
                      std::vector<CoordMat>& out_vertices, const Float dec_acc) {
  if (forward_vertices.empty() && reverse_vertices.empty())
    return;
  auto& polygon =
    out_vertices.emplace_back(forward_vertices.size() + reverse_vertices.size(), 2);
  auto row = Index(0);
  for (size_t i = 0; i < forward_vertices.size(); ++i, ++row) {
    polygon.row(row) = to_eigen(forward_vertices[i]);
  }
  for (size_t i = reverse_vertices.size(); i-- > 0; ++row) {
    polygon.row(row) = to_eigen(reverse_vertices[i]);
  }
  if (should_decimate) {
    if (dec_acc < 0)
      decimate(polygon, 0.02 * pen_width);
    else
      decimate(polygon, dec_acc);
  }
  forward_vertices.clear();
  reverse_vertices.clear();
}

Float find_theta(const Vec2 center, const Vec2 p) {
  Vec2 q = p - center;
  return std::atan2(q.y(), q.x());
}

#ifdef PASTEL_PALETTE
const Col3 palette[] = {Col3::from_hex_rgb(0x8dd3c7), Col3::from_hex_rgb(0xffffb3),
                        Col3::from_hex_rgb(0xbebada), Col3::from_hex_rgb(0xfb8072),
                        Col3::from_hex_rgb(0x80b1d3), Col3::from_hex_rgb(0xfdb462),
                        Col3::from_hex_rgb(0xb3de69), Col3::from_hex_rgb(0xfccde5),
                        Col3::from_hex_rgb(0xd9d9d9), Col3::from_hex_rgb(0xbc80bd),
                        Col3::from_hex_rgb(0xccebc5), Col3::from_hex_rgb(0xe5c494)};
#else
// const Col3 palette[] = {
//  Col3::from_hex_rgb(0xbcbd22), Col3::from_hex_rgb(0xc5b0d5),
//  Col3::from_hex_rgb(0x17becf), Col3::from_hex_rgb(0xffbb78),
//  Col3::from_hex_rgb(0xc49c94), Col3::from_hex_rgb(0xaec7e8),
//  Col3::from_hex_rgb(0x9be4ff), Col3::from_hex_rgb(0x800000),
//  Col3::from_hex_rgb(0xff9896), Col3::from_hex_rgb(0x2ca02c),
//  Col3::from_hex_rgb(0xd62728), Col3::from_hex_rgb(0xe377c2),
//  Col3::from_hex_rgb(0x54eeee), Col3::from_hex_rgb(0x9467bd),
//  Col3::from_hex_rgb(0xdbdb8d), Col3::from_hex_rgb(0xf7b6d2),
//  Col3::from_hex_rgb(0x98df8a), Col3::from_hex_rgb(0x1f77b4),
//  Col3::from_hex_rgb(0xf2f200), Col3::from_hex_rgb(0xff7f0e),
//  Col3::from_hex_rgb(0xee00ee), Col3::from_hex_rgb(0x8c564b),
//  Col3::from_hex_rgb(0x00ee00),
//};
const Col3 palette[] = {
  Col3::from_hex_rgb(0x00a391), Col3::from_hex_rgb(0xffd700),
  Col3::from_hex_rgb(0xa52a2a), Col3::from_hex_rgb(0x45769e),
  Col3::from_hex_rgb(0xf77ff2), Col3::from_hex_rgb(0xff0000),
  Col3::from_hex_rgb(0xdb7093), Col3::from_hex_rgb(0x1e90ff),
  Col3::from_hex_rgb(0x00ee00), Col3::from_hex_rgb(0xff8c00),
  Col3::from_hex_rgb(0x9046cc), Col3::from_hex_rgb(0x6a7b53),
  Col3::from_hex_rgb(0xff1493), Col3::from_hex_rgb(0x7fffd4),
  Col3::from_hex_rgb(0x9acd32), Col3::from_hex_rgb(0xba55d3),
  Col3::from_hex_rgb(0xff00ff), Col3::from_hex_rgb(0xd2691e),
  Col3::from_hex_rgb(0x458c45), Col3::from_hex_rgb(0x954595),
  Col3::from_hex_rgb(0xf0e68c), Col3::from_hex_rgb(0x00bfff),
  Col3::from_hex_rgb(0xffa07a),
};
#endif

} // namespace

bool get_render_decimation() { return should_decimate; }

void set_render_decimation(bool enabled) { should_decimate = enabled; }

int get_render_n_cap_vertices() { return cap_vertices; }

void set_render_n_cap_vertices(int n) { cap_vertices = n; }

void outline_to_polygons(const ConstStrokeView& s, std::vector<CoordMat>& out_vertices,
                         const Float dec_acc) {
  const auto n = s.size();
  if (n == 0) {
    return;
  } else if (n == 1) {
    cap_full_circle(0, s, out_vertices);
    return;
  }
  auto forward_vertices = std::vector<Vec2>();
  auto reverse_vertices = std::vector<Vec2>();
  auto need_start_cap = true;

  auto pen_width = Float(0.0);
  for (Index i = 0; i < n; ++i) {
    pen_width = std::max(pen_width, s.width(i));
  }

  const auto handle_overlapping_circles = [&](const Index i) {
    // In this case, we just split the polygon here.
    // Use the bigger cap (the smaller will be completely hidden).
    if (s.width(i - 1) > s.width(i)) {
      if (need_start_cap) {
        finalize_polygon(forward_vertices, reverse_vertices, pen_width, out_vertices,
                         dec_acc);
        cap_full_circle(i - 1, s, out_vertices);
      } else if (!reverse_vertices.empty()) {
        const Vec2 c = s.xy(i - 1);
        const Vec2 start = forward_vertices.back();
        const Vec2 end = reverse_vertices.back();
        partial_cap_ccw(c, 0.5 * s.width(i - 1), find_theta(c, start), find_theta(c, end),
                        forward_vertices);
        finalize_polygon(forward_vertices, reverse_vertices, pen_width, out_vertices,
                         dec_acc);
      }
      need_start_cap = false;
    } else {
      finalize_polygon(forward_vertices, reverse_vertices, pen_width, out_vertices,
                       dec_acc);
      need_start_cap = true;
    }
  };

  for (Index i = 1; i < n; ++i) {
    const Vec2 tangent = s.xy(i) - s.xy(i - 1);
    if (tangent.squaredNorm() > 1e-9) {
      const auto r0r1 = 0.5 * std::abs(s.width(i) - s.width(i - 1));
      // TODO: This check should be against the last *non-eclipsed* point, not the last
      //       point.
      if ((s.xy(i) - s.xy(i - 1)).squaredNorm() <= square(r0r1 * (1 + 1e-6))) {
        // Overlapping circles case.
        handle_overlapping_circles(i);
      } else {
        const auto c0 = Circle{s.x(i - 1), s.y(i - 1), 0.5 * s.width(i - 1)};
        const auto c1 = Circle{s.x(i), s.y(i), 0.5 * s.width(i)};
        auto for_p = Vec2::Empty(), rev_p = Vec2::Empty(), next_for_p = Vec2::Empty(),
             next_rev_p = Vec2::Empty();
        tangent_points(c0, c1, rev_p, for_p, next_rev_p, next_for_p);
        if (for_p.hasNaN() || rev_p.hasNaN() || next_for_p.hasNaN() ||
            next_rev_p.hasNaN()) {
          // Overlapping circles case (due to numerical imprecision).
          SPDLOG_WARN("NaN in tangent computation");
          handle_overlapping_circles(i);
        } else {
          if (need_start_cap) {
            // Start cap.
            partial_cap_ccw(s.xy(i - 1), 0.5 * s.width(i - 1),
                            find_theta(s.xy(i - 1), rev_p),
                            find_theta(s.xy(i - 1), for_p), forward_vertices);
            need_start_cap = false;
          }
          if (!reverse_vertices.empty()) {
            const auto a1 =
              signed_area(reverse_vertices.back(), forward_vertices.back(), for_p);
            const auto a2 =
              signed_area(reverse_vertices.back(), forward_vertices.back(), rev_p);
            if (a1 < 0 && a2 >= 0) {
              // Outer join.
              partial_cap_cw(s.xy(i - 1), 0.5 * s.width(i - 1),
                             find_theta(s.xy(i - 1), reverse_vertices.back()),
                             find_theta(s.xy(i - 1), rev_p), reverse_vertices);
            } else if (a2 < 0 && a1 >= 0) {
              partial_cap_ccw(s.xy(i - 1), 0.5 * s.width(i - 1),
                              find_theta(s.xy(i - 1), forward_vertices.back()),
                              find_theta(s.xy(i - 1), for_p), forward_vertices);
            } else {
              // Outer joins (a1, a2 > 0) or inner joins (a1, a2 < 0).  Same as above.
              partial_cap_cw(s.xy(i - 1), 0.5 * s.width(i - 1),
                             find_theta(s.xy(i - 1), reverse_vertices.back()),
                             find_theta(s.xy(i - 1), rev_p), reverse_vertices);
              partial_cap_ccw(s.xy(i - 1), 0.5 * s.width(i - 1),
                              find_theta(s.xy(i - 1), forward_vertices.back()),
                              find_theta(s.xy(i - 1), for_p), forward_vertices);
            }
          }
          // Start of isosceles trapezoid.
          forward_vertices.emplace_back(for_p);
          reverse_vertices.emplace_back(rev_p);
          // End of isosceles trapezoid.
          forward_vertices.emplace_back(next_for_p);
          reverse_vertices.emplace_back(next_rev_p);
        }
      }
    } else {
      // Skip near-duplicate vertex.
      // TODO: Should use maximum of widths of duplicate vertices.
    }
  }

  // End cap.
  if (need_start_cap) {
    finalize_polygon(forward_vertices, reverse_vertices, pen_width, out_vertices,
                     dec_acc);
    cap_full_circle(n - 1, s, out_vertices);
    need_start_cap = false; // noop
  } else if (!reverse_vertices.empty()) {
    const Vec2 c = s.xy(n - 1);
    const Vec2 start = forward_vertices.back();
    const Vec2 end = reverse_vertices.back();
    partial_cap_ccw(c, 0.5 * s.width(n - 1), find_theta(c, start), find_theta(c, end),
                    forward_vertices);
  }
  finalize_polygon(forward_vertices, reverse_vertices, pen_width, out_vertices, dec_acc);
}

void outline_to_triangulation(const ConstStrokeView& s,
                              std::vector<Float>& out_coordinates,
                              std::vector<unsigned>& out_indices) {

  // Helpers.
  const auto cap_full_circle = [&](Index v) {
    const auto start_index = unsigned(out_coordinates.size() / 2);
    const auto r = 0.5 * s.width(v);
    const auto x = s.x(v);
    const auto y = s.y(v);
    for (auto i = 0; i < 2 * cap_vertices; ++i) {
      const auto t = i * M_PI / cap_vertices;
      out_coordinates.push_back(x + r * std::cos(t));
      out_coordinates.push_back(y + r * std::sin(t));
    }
    // Triangle fan.
    for (unsigned i = 2; i < 2 * cap_vertices; ++i) {
      out_indices.push_back(start_index);
      out_indices.push_back(start_index + i - 1);
      out_indices.push_back(start_index + i);
    }
  };
  const auto cap_partial = [&](const Index v, const Float start_theta,
                               const Float end_theta, const unsigned start_idx,
                               const unsigned end_idx, const unsigned fanout_idx) {
    const auto cap_start_idx = unsigned(out_coordinates.size() / 2);
    const auto x = s.x(v);
    const auto y = s.y(v);
    const auto r = 0.5 * s.width(v);
    const auto num_steps =
      (unsigned)std::round(std::abs(end_theta - start_theta) * (cap_vertices / M_PI));
    const auto step = (end_theta - start_theta) / num_steps;
    for (unsigned i = 1; i < num_steps; ++i) {
      const auto t = start_theta + i * step;
      out_coordinates.push_back(x + r * std::cos(t));
      out_coordinates.push_back(y + r * std::sin(t));
    }
    // Triangle fan.
    if (num_steps <= 1) {
      out_indices.push_back(fanout_idx);
      out_indices.push_back(start_idx);
      out_indices.push_back(end_idx);
      if (end_theta < start_theta)
        std::swap(out_indices.back(), out_indices[out_indices.size() - 2]);
    } else {
      // Need to push for_p and rev_p after this function.
      out_indices.push_back(fanout_idx);
      out_indices.push_back(start_idx);
      out_indices.push_back(cap_start_idx);
      if (end_theta < start_theta)
        std::swap(out_indices.back(), out_indices[out_indices.size() - 2]);
      for (unsigned i = 2; i < num_steps; ++i) {
        out_indices.push_back(fanout_idx);
        out_indices.push_back(cap_start_idx + i - 2);
        out_indices.push_back(cap_start_idx + i - 1);
        if (end_theta < start_theta)
          std::swap(out_indices.back(), out_indices[out_indices.size() - 2]);
      }
      if (fanout_idx != end_idx) {
        out_indices.push_back(fanout_idx);
        out_indices.push_back(cap_start_idx + num_steps - 2);
        out_indices.push_back(end_idx);
        if (end_theta < start_theta)
          std::swap(out_indices.back(), out_indices[out_indices.size() - 2]);
      }
    }
  };
  const auto cap_partial_ccw = [&](Index v, unsigned start_idx, unsigned end_idx,
                                   unsigned fanout_idx) {
    const auto start_pt =
      Vec2(out_coordinates[2 * start_idx], out_coordinates[2 * start_idx + 1]);
    const auto end_pt =
      Vec2(out_coordinates[2 * end_idx], out_coordinates[2 * end_idx + 1]);
    auto start_theta = find_theta(s.xy(v), start_pt);
    auto end_theta = find_theta(s.xy(v), end_pt);
    if (end_theta < start_theta) {
      end_theta += 2 * M_PI;
    }
    cap_partial(v, start_theta, end_theta, start_idx, end_idx, fanout_idx);
  };
  const auto cap_partial_cw = [&](Index v, unsigned start_idx, unsigned end_idx,
                                  unsigned fanout_idx) {
    const auto start_pt =
      Vec2(out_coordinates[2 * start_idx], out_coordinates[2 * start_idx + 1]);
    const auto end_pt =
      Vec2(out_coordinates[2 * end_idx], out_coordinates[2 * end_idx + 1]);
    auto start_theta = find_theta(s.xy(v), start_pt);
    auto end_theta = find_theta(s.xy(v), end_pt);
    if (end_theta > start_theta) {
      start_theta += 2 * M_PI;
    }
    cap_partial(v, start_theta, end_theta, start_idx, end_idx, fanout_idx);
  };

  // Base cases.
  if (s.size() == 0) {
    return;
  } else if (s.size() == 1) {
    cap_full_circle(0);
    return;
  }

  auto need_start_cap = true;
  auto last_for_p_idx = (unsigned)~0;

  // More helpers.
  const auto handle_overlapping_circles = [&](const Index i) {
    // Use the bigger cap (the smaller will be completely hidden).
    if (s.width(i - 1) > s.width(i)) {
      if (need_start_cap) {
        cap_full_circle(i - 1);
      } else if (last_for_p_idx != (unsigned)~0) {
        cap_partial_cw(i - 1, last_for_p_idx + 3, last_for_p_idx + 2, last_for_p_idx + 2);
      }
      need_start_cap = false;
    } else {
      need_start_cap = true;
    }
    last_for_p_idx = unsigned(~0);
  };

  for (Index i = 1; i < s.size(); ++i) {
    const Vec2 tangent = s.xy(i) - s.xy(i - 1);
    if (tangent.squaredNorm() > 1e-9) {
      const auto r0r1 = 0.5 * std::abs(s.width(i) - s.width(i - 1));
      if ((s.xy(i) - s.xy(i - 1)).norm() <= r0r1 + 1e-9) {
        // Overlapping circles case.
        handle_overlapping_circles(i);
      } else {
        const auto c0 = Circle{s.x(i - 1), s.y(i - 1), 0.5 * s.width(i - 1)};
        const auto c1 = Circle{s.x(i), s.y(i), 0.5 * s.width(i)};
        auto for_p = Vec2::Empty(), rev_p = Vec2::Empty(), next_for_p = Vec2::Empty(),
             next_rev_p = Vec2::Empty();
        tangent_points(c0, c1, rev_p, for_p, next_rev_p, next_for_p);
        if (for_p.hasNaN() || rev_p.hasNaN() || next_for_p.hasNaN() ||
            next_rev_p.hasNaN()) {
          // Overlapping circles case (due to numerical imprecision).
          SPDLOG_WARN("NaN in tangent computation");
          handle_overlapping_circles(i);
        } else {
          const auto for_p_idx = unsigned(out_coordinates.size() / 2);
          // Start of isosceles trapezoid.
          out_coordinates.push_back(for_p.x()); // for_p_idx
          out_coordinates.push_back(for_p.y());
          out_coordinates.push_back(rev_p.x()); // for_p_idx + 1
          out_coordinates.push_back(rev_p.y());
          // End of isosceles trapezoid.
          out_coordinates.push_back(next_for_p.x()); // for_p_idx + 2
          out_coordinates.push_back(next_for_p.y());
          out_coordinates.push_back(next_rev_p.x()); // for_p_idx + 3
          out_coordinates.push_back(next_rev_p.y());
          // Triangle 1.
          out_indices.push_back(for_p_idx);
          out_indices.push_back(for_p_idx + 2);
          out_indices.push_back(for_p_idx + 1);
          // Triangle 2.
          out_indices.push_back(for_p_idx + 1);
          out_indices.push_back(for_p_idx + 2);
          out_indices.push_back(for_p_idx + 3);
          if (need_start_cap) {
            // Start cap.
            cap_partial_ccw(i - 1, for_p_idx + 1, for_p_idx, for_p_idx);
            need_start_cap = false;
          } else if (last_for_p_idx != (unsigned)~0) {
            const Vec2 for_back(out_coordinates[2 * (size_t)last_for_p_idx + 4],
                                out_coordinates[2 * (size_t)last_for_p_idx + 5]);
            const Vec2 rev_back(out_coordinates[2 * (size_t)last_for_p_idx + 6],
                                out_coordinates[2 * (size_t)last_for_p_idx + 7]);
            const auto a1 = signed_area(rev_back, for_back, for_p);
            const auto a2 = signed_area(rev_back, for_back, rev_p);
            if (a1 < 0 && a2 >= 0) {
              cap_partial_cw(i - 1, last_for_p_idx + 3, for_p_idx + 1,
                             last_for_p_idx + 2);
            } else if (a2 < 0 && a1 >= 0) {
              cap_partial_ccw(i - 1, last_for_p_idx + 2, for_p_idx, last_for_p_idx + 3);
            } else {
              // Outer joins (a1, a2 > 0) or inner joins (a1, a2 < 0).  Same as above.
              cap_partial_cw(i - 1, last_for_p_idx + 3, for_p_idx + 1,
                             last_for_p_idx + 2);
              cap_partial_ccw(i - 1, last_for_p_idx + 2, for_p_idx, for_p_idx + 1);
            }
          }
          last_for_p_idx = for_p_idx;
        }
      }
    } else {
      // Skip near-duplicate vertex.
      // TODO: Should use maximum of widths of duplicate vertices.
    }
  }

  // End cap.
  if (need_start_cap) {
    cap_full_circle(s.size() - 1);
    need_start_cap = false; // noop
  } else if (last_for_p_idx != (unsigned)~0) {
    cap_partial_cw(s.size() - 1, last_for_p_idx + 3, last_for_p_idx + 2,
                   last_for_p_idx + 2);
  }
}

Float rasterize_scale(span<const Stroke> strokes, const int width, const int height,
                      BoundingBox& bb) {
  if (bb.isEmpty()) {
    bb = visual_bounds(strokes);
  }
  return 0.9 * std::min(width / bb.width(), height / bb.height());
}

span<const Col3> get_color_palette() {
  return {palette, sizeof(palette) / sizeof(palette[0])};
}

} // namespace sketching
