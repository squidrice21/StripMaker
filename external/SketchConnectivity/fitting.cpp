#include "fitting.h"

#include "detail/util.h"
#include "eigen_compat.h"
#include "force_assert.h"
#include "resample.h"

#ifdef HAS_GUROBI
#include <StrokeStrip/FittingEigenSparse.h>
#include <StrokeStrip/Parameterization.h>
#include <StrokeStrip/SketchInfo.h>
#include <StrokeStrip/StrokeOrientation.h>
#endif

#include <GraphicsGems/FitCurves.h>
#include <spdlog/spdlog.h>

namespace sketching {

namespace {

/// Cubic Bézier basis coefficients.
Float B0(Float u) {
  Float tmp = 1.0 - u;
  return (tmp * tmp * tmp);
}
Float B1(Float u) {
  Float tmp = 1.0 - u;
  return (3 * u * (tmp * tmp));
}
Float B2(Float u) {
  Float tmp = 1.0 - u;
  return (3 * u * u * tmp);
}
Float B3(Float u) { return (u * u * u); }

Vec2 evaluate_bezier(const int degree, const Vec2 control_pts[], Float t) {
  assert(degree <= 3);

  // Local copy of control points.
  Vec2 Vtemp[4];
  for (int i = 0; i < 4; ++i)
    Vtemp[i] = control_pts[i];

  // Triangle computation.
  for (int i = 1; i <= degree; i++) {
    for (int j = 0; j <= degree - i; j++) {
      Vtemp[j] = (1.0 - t) * Vtemp[j] + t * Vtemp[j + 1];
    }
  }

  return Vtemp[0]; // Point on curve at parameter t.
}

/// @param bez Current fitted curve.
/// @param p Digitized point.
/// @param u Parameter value for p.
Float newton_raphson_root_find(const Bezier& bez, const Vec2 p, const Float u) {
  // Compute bez(u).
  const Vec2 Q_u = evaluate_bezier(3, bez.pts, u);

  Vec2 Q1[3]; // Control vertices for bez'.
  for (int i = 0; i < 3; i++) {
    Q1[i] = 3.0 * (bez.pts[i + 1] - bez.pts[i]);
  }

  Vec2 Q2[2]; // Control vertices for bez''.
  for (int i = 0; i < 2; i++) {
    Q2[i] = 2.0 * (Q1[i + 1] - Q1[i]);
  }

  // Compute Q'(u) and Q''(u).
  const Vec2 Q1_u = evaluate_bezier(2, Q1, u);
  const Vec2 Q2_u = evaluate_bezier(1, Q2, u);

  // Compute f(u) / f'(u).
  const auto denominator = (Q1_u.x_) * (Q1_u.x_) + (Q1_u.y_) * (Q1_u.y_) +
                           (Q_u.x_ - p.x_) * (Q2_u.x_) + (Q_u.y_ - p.y_) * (Q2_u.y_);
  if (denominator < 1e-7)
    return u;
  const auto numerator = (Q_u.x_ - p.x_) * (Q1_u.x_) + (Q_u.y_ - p.y_) * (Q1_u.y_);

  // u <- u - f(u) / f'(u)
  return u - (numerator / denominator);
}

void reparameterize(span<const Float> x, span<const Float> y,
                    span<Float> parameterization, const Bezier& bez) {
  const auto n = parameterization.size();
  for (size_t i = 0; i < n; ++i) {
    parameterization[i] =
      newton_raphson_root_find(bez, Vec2(x[i], y[i]), parameterization[i]);
  }
}

Bezier bezier_iteration(span<const Float> x, span<const Float> y,
                        span<const Float> parameterization, const Vec2 start_tangent,
                        const Vec2 end_tangent) {
  const auto n = (int)x.size();
  const auto p_first = Vec2(x[0], y[0]);
  const auto p_last = Vec2(x[n - 1], y[n - 1]);

  // Create the C and X matrices.
  double C[2][2];
  C[0][0] = 0.0;
  C[0][1] = 0.0;
  C[1][0] = 0.0;
  C[1][1] = 0.0;
  double X[2];
  X[0] = 0.0;
  X[1] = 0.0;
  for (int i = 0; i < n; ++i) {
    // RHS
    const auto A0 = start_tangent * B1(parameterization[i]);
    const auto A1 = -end_tangent * B2(parameterization[i]);

    C[0][0] += A0.dot(A0);
    C[0][1] += A0.dot(A1);
    // C[1][0] += A0.dot(A1);
    C[1][0] = C[0][1];
    C[1][1] += A1.dot(A1);

    const auto tmp =
      Vec2(x[i], y[i]) -
      (B0(parameterization[i]) * p_first + B1(parameterization[i]) * p_first +
       B2(parameterization[i]) * p_last + B3(parameterization[i]) * p_last);
    X[0] += A0.dot(tmp);
    X[1] += A1.dot(tmp);
  }

  // Compute the determinants of C and X.
  const auto det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
  const auto det_C0_X = C[0][0] * X[1] - C[1][0] * X[0];
  const auto det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1];

  // Finally, derive left and right alpha values.
  auto alpha_l = (det_C0_C1 == 0.0 ? 0.0 : det_X_C1 / det_C0_C1);
  auto alpha_r = (det_C0_C1 == 0.0 ? 0.0 : det_C0_X / det_C0_C1);

  // If alpha negative, use the Wu/Barsky heuristic (see text).
  // (If alpha is 0, you get coincident control points that lead to divide by zero in any
  // subsequent newton_raphson_root_find call.)
  const auto seglen = (p_last - p_first).norm();
  const auto epsilon = 1e-6 * seglen;
  if (alpha_l < epsilon || alpha_r < epsilon) {
    // Fall back to standard (probably inaccurate) formula.
    alpha_l = alpha_r = seglen / 3.0;
  }

  // First and last control points of the Bézier curve are positioned exactly at the first
  // and last data points.
  // Control points 1 and 2 are positioned an alpha distance out on the tangent vectors,
  // left and right, respectively.
  auto bez = Bezier();
  bez.pts[0] = p_first;
  bez.pts[3] = p_last;
  bez.pts[1] = bez.pts[0] + alpha_l * start_tangent;
  bez.pts[2] = bez.pts[3] - alpha_r * end_tangent;
  return bez;
}

// NOTE: Not thread-safe.
std::vector<Bezier> g_bezier_out;
std::vector<int> g_start_indices;
std::vector<double> g_parameterization;
void bezier_out_append(const int n, Point2* control_pts, const int first, int last,
                       double* parameterization) {
  auto& bezier = g_bezier_out.emplace_back();
  assert(n == 3);
  for (int i = 0; i <= n; ++i) {
    bezier.pts[i] = Vec2(control_pts[i].x, control_pts[i].y);
  }
  if (parameterization) {
    for (int i = 0; i < last - first; ++i) {
      g_parameterization.push_back(parameterization[i]);
    }
  } else {
    assert(last - first == 1);
    g_parameterization.push_back(0.0);
  }
  g_start_indices.push_back(first);
}

std::pair<Index, Float> fractional_index(const span<const Float> arclengths,
                                         const Float s) {
  if (s == 0.0) {
    return {0, 0.0};
  } else if (s == arclengths.back()) {
    return {arclengths.size() - 2, 1.0};
  }

  // Essentially copied from VPaint.
  auto i = Index(0);
  auto j = Index(arclengths.size() - 1);
  auto si = arclengths[i];
  auto sj = arclengths[j];

  while (j - i > 1) {
    // compute an index hopefully close to s
    const auto u = (s - si) / (sj - si);
    auto k = Index(std::floor((1 - u) * i + u * j));
    // make sure i < k < j
    k = std::min(j - 1, std::max(i + 1, k));
    // recurse
    const auto sk = arclengths[k];
    if (sk > s) {
      j = k;
      sj = sk;
    } else {
      i = k;
      si = sk;
    }
  }
  assert(j == i + 1);
  assert(arclengths[i] - 1e-6 <= s);
  assert(arclengths[j] + 1e-6 >= s);
  const auto u = (s - si) / (sj - si);
  assert(u >= -1e-6);
  return {i, u};
}

} // namespace

Vec2 head_tangent_lagrange(const ConstStrokeView& s, const Float step_factor) {
  assert(s.size() > 1);
  assert(s.length() > 1e-6);
  if (s.size() == 2) {
    return (s.xy(0) - s.xy(1)).normalized();
  }
  const auto double_step = std::min(s.length(), 2 * step_factor * s.stroke().pen_width());
  const auto step = 0.5 * double_step;
  if (step <= s.arclength(1)) {
    return (s.xy(0) - s.xy(1)).normalized();
  }
  const Vec2 p0 = s.xy(0);
  const Vec2 p1 = s.pos(step);
  const Vec2 p2 = s.pos(double_step);
  // Compute tangent of quadratic Lagrange polynomial.
  constexpr double w0 = 3.0 / 2.0;
  constexpr double w1 = -2.0;
  constexpr double w2 = 0.5;
  return (w0 * p0 + w1 * p1 + w2 * p2).normalized();
}

Vec2 tail_tangent_lagrange(const ConstStrokeView& s, const Float step_factor) {
  assert(s.size() > 1);
  assert(s.length() > 1e-6);
  const auto n = s.size();
  if (n == 2) {
    return (s.xy(n - 1) - s.xy(n - 2)).normalized();
  }
  const auto double_step = std::min(s.length(), 2 * step_factor * s.stroke().pen_width());
  const auto step = 0.5 * double_step;
  if (step <= s.length() - s.arclength(n - 2)) {
    return (s.xy(n - 1) - s.xy(n - 2)).normalized();
  }
  const Vec2 p0 = s.xy(Back);
  const Vec2 p1 = s.pos(s.length() - step);
  const Vec2 p2 = s.pos(s.length() - double_step);
  // Compute tangent of quadratic Lagrange polynomial.
  constexpr double w0 = 3.0 / 2.0;
  constexpr double w1 = -2.0;
  constexpr double w2 = 0.5;
  return (w0 * p0 + w1 * p1 + w2 * p2).normalized();
}

Vec2 BezierSpline::pos(const Float s) const {
  const auto [i, u] = fractional_index({&arclength_[0], (size_t)n_segments_ + 1}, s);
  return segments_[i].pos(u);
}

Vec2 BezierSpline::tangent(const Float s) const {
  const auto [i, u] = fractional_index({&arclength_[0], (size_t)n_segments_ + 1}, s);
  return segments_[i].tangent(u);
}

std::pair<Vec2, Vec2> BezierSpline::corner_tangents(const Float s) const {
  const auto [i, u] = fractional_index({&arclength_[0], (size_t)n_segments_ + 1}, s);
  const auto tangent2 = segments_[i].tangent(u);
  if (i > 0 && u < 1e-8) {
    const auto tangent1 = segments_[i - 1].tangent(1.0);
    return {tangent1, tangent2};
  } else {
    return {tangent2, tangent2};
  }
}

void fit_bezier_spline(const ConstStrokeView& s, const Float error,
                       std::vector<Bezier>& out_spline,
                       std::vector<int>* out_start_indices,
                       std::vector<double>* out_parameterization) {
  const auto n = s.size();
  const auto points = std::make_unique<Point2[]>(n);
  for (auto i = 0; i < n; ++i) {
    points[i] = {s.x(i), s.y(i)};
  }
  g_bezier_out.clear();
  g_start_indices.clear();
  g_parameterization.clear();
  FitCurve(points.get(), (int)n, error, bezier_out_append);
  if (out_parameterization) {
    std::swap(*out_parameterization, g_parameterization);
    out_parameterization->push_back(1.0); // Missing last point.
  }
  if (out_start_indices) {
    std::swap(*out_start_indices, g_start_indices);
    out_start_indices->push_back((int)s.size());
  }
  std::swap(out_spline, g_bezier_out);
}

BezierSpline fit_bezier_spline_with_corners(const Stroke& s, const Float tolerance) {
  const auto n = s.size();
  if (n < 2) {
    return BezierSpline{0, nullptr, nullptr};
  } else if (s.size() == 2) {
    auto spline = BezierSpline();
    spline.n_segments_ = 1;
    spline.arclength_ = std::make_unique<Float[]>(2);
    spline.arclength_[0] = 0.0;
    spline.arclength_[1] = s.length();
    spline.segments_ = std::make_unique<Bezier[]>(1);
    auto& bez = spline.segments_[0];
    bez.pts[0] = s.xy(0);
    bez.pts[1] = lerp(s.xy(0), s.xy(1), 1 / 3.0);
    bez.pts[2] = lerp(s.xy(0), s.xy(1), 2 / 3.0);
    bez.pts[3] = s.xy(1);
    return spline;
  }
  auto corner = std::make_unique<bool[]>(n);
  ramer_douglas_peucker(s, tolerance, corner.get());
  assert(corner[0] && corner[n - 1]);

  auto spline = BezierSpline();
  spline.n_segments_ = 1;
  for (Index i = 1; i < n - 1; ++i) {
    if (corner[i])
      spline.n_segments_++;
  }
  spline.arclength_ = std::make_unique<Float[]>(spline.n_segments_ + 1);
  spline.arclength_[0] = 0.0;
  spline.segments_ = std::make_unique<Bezier[]>(spline.n_segments_);

  s.ensure_arclengths();

  auto segment_i = Index(0);
  auto start = Index(0);
  auto parameterization = std::vector<Float>();
  for (Index i = start + 1; i < n; ++i) {
    if (corner[i]) {
      const auto end = i;
      const auto seg_size = end - start + 1;
      const auto view = ConstStrokeView(s, start, end + 1);
      const auto start_tangent = -head_tangent_lagrange(view);
      const auto end_tangent = tail_tangent_lagrange(view);
      parameterization.resize(seg_size);
      parameterization[0] = 0.0;
      const auto inv_len = 1.0 / (s.arclength(end) - s.arclength(start));
      for (Index k = 1; k < seg_size - 1; ++k) {
        parameterization[k] = (s.arclength(start + k) - s.arclength(start)) * inv_len;
      }
      parameterization[seg_size - 1] = 1.0;
      auto bez =
        fit_bezier({&s.x_[start], (size_t)seg_size}, {&s.y_[start], (size_t)seg_size},
                   parameterization, start_tangent, end_tangent);
      spline.segments_[segment_i] = bez;
      spline.arclength_[segment_i + 1] = s.arclength(end);
      segment_i++;
      // Now continue.
      start = end;
    }
  }
  assert(segment_i == spline.n_segments_);

  return spline;
}

Bezier fit_bezier(span<const Float> x, span<const Float> y, span<Float> parameterization,
                  const Vec2 start_tangent, const Vec2 end_tangent) {
  // Based on "An Algorithm for Automatically Fitting Digitized Curves" (Schneider 1990).

  assert(x.size() == y.size());
  assert(x.size() == parameterization.size());
  assert(parameterization[0] == 0.0);
  assert(parameterization.back() == 1.0);
  const auto last = int(x.size()) - 1;
  if (x.size() == 2) {
    auto bez = Bezier();
    bez.pts[0] = Vec2(x[0], y[0]);
    bez.pts[3] = Vec2(x[last], y[last]);
    const auto dist = (bez.pts[3] - bez.pts[0]).norm() / 3.0;
    bez.pts[1] = bez.pts[0] + dist * start_tangent;
    bez.pts[2] = bez.pts[3] - dist * end_tangent;
    return bez;
  }

  auto bez = bezier_iteration(x, y, parameterization, start_tangent, end_tangent);
  // Key insight from Schneider (1990): we can iteratively update the parameterization
  // after each fit to get a better fit.
  constexpr auto num_iter = 4;
  for (int i = 0; i < num_iter; ++i) {
    reparameterize(x, y, parameterization, bez);
    bez = bezier_iteration(x, y, parameterization, start_tangent, end_tangent);
  }
  return bez;
}

Clustering::Clustering(std::unique_ptr<int[]>&& stroke2cluster,
                       std::map<int, std::vector<int>>&& cluster2stroke,
                       const int n_strokes)
  : cluster2stroke_(std::move(cluster2stroke))
  , stroke2cluster_(std::move(stroke2cluster))
  , n_strokes_(n_strokes) {

  // For strokes without a specified cluster, we must assign them a cluster index
  // without colliding with any existing clusters.
  // It is for this reason we need to "seal" (i.e. make immutable) this object, because
  // we do not want `max_cluster_index` changing underneath us.
  const auto max_cluster_index = cluster2stroke_.rbegin()->first;
  for (int stroke_index = 0; stroke_index < n_strokes; ++stroke_index) {
    if (stroke2cluster_[stroke_index] == -1) {
      const auto cluster_index = max_cluster_index + 1 + stroke_index;
      cluster2stroke_[cluster_index].push_back(stroke_index);
      stroke2cluster_[stroke_index] = cluster_index;
    }
  }

  // Enforce a consistent ordering.
  for (auto& [cluster_index, stroke_indices] : cluster2stroke_) {
    std::sort(stroke_indices.begin(), stroke_indices.end());
  }
}

std::string Clustering::repr() const {
  auto ss = std::stringstream();
  ss << "Clustering([";
  for (auto i = 0; i < n_strokes_ - 1; ++i) {
    ss << i << ", ";
  }
  ss << stroke2cluster_[n_strokes_ - 1] << "])";
  return ss.str();
}

std::vector<int> Clustering::clusters() const {
  auto out = std::vector<int>();
  out.reserve(cluster2stroke_.size());
  for (const auto& [cluster_index, _] : cluster2stroke_) {
    out.push_back(cluster_index);
  }
  return out;
}

ClusteringBuilder::ClusteringBuilder(const int n_strokes)
  : stroke2cluster_(new int[n_strokes])
  , n_strokes_(n_strokes) {
  std::fill_n(stroke2cluster_.get(), n_strokes, -1);
}

void ClusteringBuilder::add(const int stroke_index, const int cluster_index) {
  if (stroke2cluster_[stroke_index] != -1) {
    std::stringstream ss;
    ss << "cut strokes not supported, found multiple clusters associated with stroke "
       << stroke_index;
    throw std::runtime_error(ss.str());
  }
  stroke2cluster_[stroke_index] = cluster_index;
  cluster2stroke_[cluster_index].push_back(stroke_index);
}

#ifdef HAS_GUROBI

namespace {

std::vector<double> compute_arclengths(const SketchUI::Polyline2D& s) {
  std::vector<double> out;
  out.reserve(s.points.size());
  auto sofar = Float(0.0);
  out.push_back(0.0);
  for (size_t i = 1; i < s.points.size(); ++i) {
    const auto p_i = s.points[i].first;
    const auto p_i1 = s.points[i - 1].first;
    const auto dist = (p_i - p_i1).Length();
    sofar += dist;
    out.push_back(sofar);
  }
  return out;
}

} // namespace

void strokes_to_capture(const std::vector<const Stroke*>& strokes, const Float accuracy,
                        const bool cut_strokes, std::vector<StrokeMapping>& out_mappings,
                        Capture& capture) {
  auto widths = 0.0;
  auto num_vertices = 0;
  auto max_length = 0.0;
  auto polyline_index = 0u;
  for (auto i = size_t(0); i < strokes.size(); ++i) {
    const auto& stroke = *strokes[i];
    const auto n = stroke.size();
    auto& mapping = out_mappings.emplace_back();
    if (n > 1 && (strokes.size() == 1 || stroke.length() > 1.5 * stroke.pen_width())) {
      stroke.ensure_arclengths();
      const Vec dx = as_eigen(stroke.x()).tail(n - 1) - as_eigen(stroke.x()).head(n - 1);
      const Vec dy = as_eigen(stroke.y()).tail(n - 1) - as_eigen(stroke.y()).head(n - 1);

      auto* polyline = &capture.sketchedPolylines.emplace_back();
      polyline->width = as_eigen(stroke.width()).mean();
      polyline->stroke_ind = (int)polyline_index++;
      polyline->additional_ind = (int)i;
      polyline->group_ind = 0;
      polyline->points.emplace_back(SketchUI::Point2D(stroke.x(0), stroke.y(0)), 0);
      mapping.domain_arclens_.push_back(0.0);
      for (Index j = 1; j < n - 1; ++j) {
        polyline->points.emplace_back(SketchUI::Point2D(stroke.x(j), stroke.y(j)), 0);
        polyline->stroke_ind = (int)polyline_index;
        polyline->additional_ind = (int)i;
        if (cut_strokes) {
          const Vec2 a(dx(j - 1), dy(j - 1));
          const Vec2 b(dx(j), dy(j));
          const auto turn_angle = std::abs(
            std::atan2(a.x() * b.y() - a.y() * b.x(), a.x() * b.x() + a.y() * b.y()));
          if (turn_angle > 0.5 * M_PI) {
            // Split stroke here.
            polyline = &capture.sketchedPolylines.emplace_back();
            polyline->width = as_eigen(stroke.width()).mean();
            polyline->stroke_ind = (int)polyline_index++;
            polyline->additional_ind = (int)i;
            polyline->group_ind = 0;
            polyline->points.emplace_back(SketchUI::Point2D(stroke.x(j), stroke.y(j)), 0);
            mapping.domain_arclens_.push_back(stroke.arclength(j));
          }
        }
      }
      polyline->points.emplace_back(SketchUI::Point2D(stroke.x(n - 1), stroke.y(n - 1)),
                                    0);
      mapping.domain_arclens_.push_back(stroke.length());

      widths += as_eigen(stroke.width()).sum();
      num_vertices += (int)stroke.size();
      max_length = std::max(max_length, stroke.length());
    } else {
      if (n <= 1)
        SPDLOG_WARN("singular stroke");

      mapping.domain_arclens_.push_back(0.0);
      mapping.domain_arclens_.push_back(stroke.length());
    }
  }
  // Clamp width to prevent very long computations for very thin strokes.
  capture.thickness =
    std::max(widths / ((double)num_vertices * accuracy), 0.001 * max_length);
}

Input from_capture(Capture& capture, glm::dvec2& out_center) {
  Input input;
  auto& clusters = input.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input.thickness = capture.thickness;

  double max_width = -1;
  for (auto& polyline : capture.sketchedPolylines) {
    max_width = std::max(max_width, polyline.width);
  }

  double rate = 2.75;
  double reparam_rate = rate * capture.thickness;

  double width_ratio = 5;
  if (max_width > width_ratio * reparam_rate)
    reparam_rate = max_width / width_ratio;

  for (auto& polyline : capture.sketchedPolylines) {
    polyline.reparameterize(std::min(reparam_rate, polyline.totalLen() / 3));

    for (const auto& point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  out_center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
  input.width = (max_x - min_x) / capture.thickness;
  input.height = (max_y - min_y) / capture.thickness;

  for (auto& polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;
    clusters[polyline.group_ind].strokes.emplace_back();
    auto& stroke = clusters[polyline.group_ind].strokes.back();
    stroke.points.reserve(polyline.points.size());
    stroke.u.reserve(polyline.points.size());
    for (size_t i = 0; i < polyline.points.size(); ++i) {
      const auto& point = polyline.points[i];

      // Normalize to stroke width
      stroke.points.push_back((glm::dvec2(point.first.x, point.first.y) - out_center) /
                              capture.thickness);
      stroke.u.push_back(0);
      stroke.max_width = polyline.width / capture.thickness;
    }
  }

  return input;
}

Stroke fit_stroke_to_cluster(const std::vector<const Stroke*>& strokes,
                             const Float accuracy, const bool cut_strokes,
                             std::vector<StrokeMapping>& out_mappings) {
  // TODO: Should just convert from Strokes to Input directly.  But for now, we rely on
  //       their resampling...
  auto skipped_stroke = std::vector<bool>(strokes.size(), false);
  Capture capture;
  for (auto i = size_t(0); i < strokes.size(); ++i) {
    const auto& stroke = *strokes[i];
    const auto n = stroke.size();
    if (!(n > 1 && stroke.length() > 1.5 * stroke.pen_width())) {
      skipped_stroke[i] = true;
    }
  }
  strokes_to_capture(strokes, accuracy, cut_strokes, out_mappings, capture);

  if (capture.sketchedPolylines.empty()) {
    throw std::runtime_error("no strokes in cluster were long enough");
  } else if (capture.sketchedPolylines.size() == 1) {
    const Stroke* unskipped = nullptr;
    for (auto i = size_t(0); i < strokes.size(); ++i) {
      if (skipped_stroke[i]) {
        for (const auto _ : out_mappings[i].domain_arclens_)
          out_mappings[i].range_arclens_.push_back(0.0);
      } else {
        const auto& stroke = *strokes[i];
        out_mappings[i].range_arclens_.push_back(0.0);
        out_mappings[i].range_arclens_.push_back(stroke.length());
        unskipped = &stroke;
      }
    }
    return unskipped->clone();
  }

  glm::dvec2 center;
  Input input = from_capture(capture, center);
  const Context context(true);
  auto orientations = std::vector<int>();
  {
    StrokeOrientation orientation(context);
    orientation.orient_strokes(&input);
    orientation.flip_strokes(&input);

    orientations = std::move(orientation.orientations.begin()->second);
  }

  {
    Parameterization param(context);
    param.parameterize(&input);
  }

  {
    FittingEigenSparse fitting(context);
    auto out = fitting.fit(&input, true);
    force_assert(out.size() == 1);
    const auto& arrays = out.begin()->second;
    const auto& fitted = arrays[capture.sketchedPolylines.size()];
    auto strokefit = Stroke(fitted.size(), false);

    // Copy points.
    for (size_t i = 0; i < fitted.size(); ++i) {
      const auto& p = fitted[i];
      strokefit.x(i) = p.x * capture.thickness + center.x;
      strokefit.y(i) = p.y * capture.thickness + center.y;
      strokefit.width(i) =
        accuracy * capture.thickness; // TODO: Transfer thickness better.
    }

    // Compute mapping.
    const auto fit_length = strokefit.length();
    auto cut_stroke_index = 0u;
    for (size_t i = 0; i < strokes.size(); ++i) {
      auto& our_mapping = out_mappings[i];
      if (!skipped_stroke[i]) {
        for (size_t j = 0; j < our_mapping.domain_arclens_.size() - 1; ++j) {
          const auto& ss_mapping = arrays[cut_stroke_index];
          auto range_start = ss_mapping[0].t * capture.thickness;
          auto range_stop = ss_mapping.back().t * capture.thickness;
          if (orientations[cut_stroke_index] == -1)
            std::swap(range_start, range_stop);
          range_start = std::min(range_start, fit_length);
          range_stop = std::min(range_stop, fit_length);
          if (our_mapping.range_arclens_.empty()) {
            our_mapping.range_arclens_.push_back(range_start);
          } else {
            auto& last_pushed = our_mapping.range_arclens_.back();
            last_pushed = 0.5 * (last_pushed + range_start); // Average.
          }
          our_mapping.range_arclens_.push_back(range_stop);

          cut_stroke_index++;
        }
      } else {
        for (const auto _ : our_mapping.domain_arclens_)
          our_mapping.range_arclens_.push_back(0.0);
      }
      force_assert(our_mapping.domain_arclens_.size() ==
                   our_mapping.range_arclens_.size());
    }
    return decimated(strokefit, 0.1 * accuracy * capture.thickness);
  }
}

#endif

} // namespace sketching
