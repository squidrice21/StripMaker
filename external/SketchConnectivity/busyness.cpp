#include "busyness.h"

#include "bvh.h"
#include "detail/util.h"
#include "force_assert.h"
#include "sketching.h"

namespace sketching {

double gaussian_factor(const double x, const double sigma) {
  return std::exp(-0.5 * square(x / sigma));
}

Float average_width(const Stroke& s) {
  if (s.size() == 0) {
    return std::numeric_limits<Float>::quiet_NaN();
  } else if (s.size() == 1) {
    return s.width(0);
  } else if (s.length() < 1e-8) {
    const auto width = s.width();
    return *std::max_element(width.begin(), width.end());
  }
  s.ensure_arclengths();
  double acc = 0.0;
  for (auto i = 1; i < s.size(); ++i) {
    const auto mid = 0.5 * (s.width(i) + s.width(i - 1));
    const auto seglen = s.arclength(i) - s.arclength(i - 1);
    acc += mid * seglen;
  }
  const auto res = acc / s.length();
  force_assert(res >= 0.0);
  return res;
}

void average_widths(const PolylineBVH& bvh, span<Float> out_widths) {
  const auto nstrokes = bvh.nodes.size();
  assert(out_widths.size() == nstrokes);
  for (auto i = decltype(nstrokes){0}; i < nstrokes; ++i) {
    out_widths[i] = average_width(*bvh.nodes[i].geometry);
  }
}

Float busyness_at(const PolylineBVH& bvh, const Vec2& p, span<const Float> avg_widths,
                  const Float busyness_falloff) {
  double busyness = 0.0;
  auto at_endpoint = false;
  const auto nstrokes = bvh.nodes.size();
  for (auto i = decltype(nstrokes){0}; i < nstrokes; ++i) {
    const auto& stroke = *bvh.nodes[i].geometry;
    const auto npoints = stroke.size();
    // Avoid the endpoints on the masked strokes
    if (npoints > 1 && std::find(bvh.masked_nodes.begin(), bvh.masked_nodes.end(), i) ==
                         bvh.masked_nodes.end()) {
      const auto w = avg_widths[i];
      if (w > 1e-6) {
        for (auto j : {Index(0), npoints - 1}) {
          if (stroke.xy(j).isApprox(p)) {
            // We want to avoid counting endpoints at appear at a graph vertex multiple
            // times.
            at_endpoint = true;
          } else {
            const auto d = Float{(stroke.xy(j) - p).norm()};
            busyness += gaussian_factor(d / w, busyness_falloff);
          }
        }
      }
    }
  }
  return (at_endpoint ? busyness + 1 : busyness);
}

Eigen::Matrix<Float, Eigen::Dynamic, 2>
busyness_at_endpoints(const PolylineBVH& bvh, const Float busyness_falloff) {
  const auto nstrokes = bvh.nodes.size();
  Eigen::Matrix<Float, Eigen::Dynamic, 2> busyness;
  busyness.resize(nstrokes, 2);
  auto avg_widths_mem = std::make_unique<Float[]>(nstrokes);
  auto avg_widths = span<Float>{&avg_widths_mem[0], nstrokes};
  average_widths(bvh, avg_widths);
  for (auto i = decltype(nstrokes){0}; i < nstrokes; ++i) {
    const auto& s = *bvh.nodes[i].geometry;
    if (s.size() <= 1) {
      busyness.row(i).fill(std::numeric_limits<Float>::quiet_NaN());
    } else {
      busyness(i, 0) = busyness_at(bvh, s.xy(0), avg_widths, busyness_falloff);
      busyness(i, 1) = busyness_at(bvh, s.xy(s.size() - 1), avg_widths, busyness_falloff);
    }
  }
  force_assert(busyness.minCoeff() > 0.99);
  return busyness;
}

} // namespace sketching
