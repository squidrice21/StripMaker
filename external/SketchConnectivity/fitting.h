#pragma once

#include "mapping.h"
#include "stroke_view.h"

#include <map>

#ifdef HAS_GUROBI
#include <StrokeStrip/SketchInfo.h>
struct Input;
#endif

namespace sketching {

/**
 * Samples control points at step_factor * pen_width and 2 * step_factor * pen_width away
 * from the endpoint.
 *
 * Result is normalized and points away from the stroke.
 */
Vec2 head_tangent_lagrange(const ConstStrokeView& s, Float step_factor = 2.5);

/**
 * Samples control points at step_factor * pen_width and 2 * step_factor * pen_width away
 * from the endpoint.
 *
 * Result is normalized and points away from the stroke.
 */
Vec2 tail_tangent_lagrange(const ConstStrokeView& s, Float step_factor = 2.5);

/// Cubic Bézier curve.
struct Bezier {
  static constexpr auto degree = 3;
  Vec2 pts[degree + 1];

  /// Return the point on the curve at parameter t.
  Vec2 pos(const Float t) const {
    // Local copy of control points.
    Vec2 tmp[degree + 1];
    for (int i = 0; i < degree + 1; ++i)
      tmp[i] = pts[i];

    // Triangle computation.
    for (int i = 1; i <= degree; i++) {
      for (int j = 0; j <= degree - i; j++) {
        tmp[j] = (1.0 - t) * tmp[j] + t * tmp[j + 1];
      }
    }

    return tmp[0];
  }

  Vec2 tangent(const Float t) const {
    return (-3 * (1 - t) * (1 - t) * pts[0] //
            + 3 * (1 - t) * (1 - t) * pts[1] //
            - 6 * t * (1 - t) * pts[1] //
            - 3 * t * t * pts[2] //
            + 6 * t * (1 - t) * pts[2] //
            + 3 * t * t * pts[3])
      .normalized();
  }

  /// Return the negative normalized head tangent.  This points away from the stroke.
  Vec2 normalized_head_tangent() const { return -tangent(0.0).normalized(); }

  Vec2 normalized_tail_tangent() const { return tangent(1.0).normalized(); }
};

struct BezierSpline {
  /// @param s Arc length value.  Ranges from 0 to spline length.
  Vec2 pos(Float s) const;

  /// @param s Arc length value.  Ranges from 0 to spline length.
  Vec2 tangent(Float s) const;

  /// Get the two tangents at a corner.
  /// If s corresponds to a non-corner, the same tangent will be returned twice.
  std::pair<Vec2, Vec2> corner_tangents(Float s) const;

  /// This value is approximate.
  Float length() const { return arclength_[n_segments_]; }

  /// Number of Bézier segments.
  Index n_segments_ = 0;
  /// Size is n_segments_ + 1 (first value is 0, last is total length).
  std::unique_ptr<Float[]> arclength_;
  std::unique_ptr<Bezier[]> segments_;
};

/**
 * Returns (bezier_segments, start_indices).
 *
 * @deprecated @ref fit_bezier_spline_with_corners should probably be used instead.
 */
void fit_bezier_spline(const ConstStrokeView& s, Float error,
                       std::vector<Bezier>& out_spline,
                       std::vector<int>* out_start_indices = nullptr,
                       std::vector<double>* out_parameterization = nullptr);

/**
 * Fit a Bézier spline, aligning individual segment boundaries with sharp turns (corners).
 */
BezierSpline fit_bezier_spline_with_corners(const Stroke& s, Float tolerance);

/**
 * Fit a single Bézier segment to an ordered collection of points described by the
 * parallel arrays x, y, and inout_parameterization.
 *
 * @param x X-coordinates of original curve.
 * @param y Y-coordinates of original curve.  Must be the same size as x.
 * @param[in, out] inout_parameterization The initial value of this buffer should contain
 *     the initial parameterization for the fitted curve.  This buffer will be updated
 *     with the final parameterization when the function returns.  Must be the same size
 *     as x.
 * @param start_tangent Estimated tangent of the beginning of the fit.  Should point
 *     towards the curve.
 * @param end_tangent Estimated tangent of the end of the fit.  Should point away from the
 *     curve.
 */
Bezier fit_bezier(span<const Float> x, span<const Float> y,
                  span<Float> inout_parameterization, Vec2 start_tangent,
                  Vec2 end_tangent);

struct Clustering {
  Clustering()
    : n_strokes_(0) {}

  Clustering(std::unique_ptr<int[]>&& stroke2cluster,
             std::map<int, std::vector<int>>&& cluster2stroke, int n_strokes);

  int get_cluster_index(const int stroke_index) const {
    if (stroke_index < 0 || stroke_index >= n_strokes_)
      throw std::out_of_range("");
    return stroke2cluster_[stroke_index];
  }

  size_t get_cluster_index(const size_t stroke_index) const {
    return stroke2cluster_[stroke_index];
  }

  const std::vector<int>* get_stroke_indices(const int cluster_index) const {
    const auto it = cluster2stroke_.find(cluster_index);
    if (it != cluster2stroke_.end()) {
      return &it->second;
    }
    throw std::out_of_range("");
  }

  explicit operator bool() const { return !cluster2stroke_.empty(); }

  std::string repr() const;

  Eigen::Map<const Eigen::VectorXi> array() const {
    return Eigen::Map<const Eigen::VectorXi>(stroke2cluster_.get(), n_strokes_);
  }

  std::vector<int> clusters() const;

  const std::map<int, std::vector<int>>& cluster2stroke() const {
    return cluster2stroke_;
  }

  int n_strokes() const { return n_strokes_; }

private:
  std::map<int, std::vector<int>> cluster2stroke_;
  std::unique_ptr<int[]> stroke2cluster_;
  const int n_strokes_;
};

struct ClusteringBuilder {
  explicit ClusteringBuilder(int n_strokes);

  void add(int stroke_index, int cluster_index);

  /// The ClusteringBuilder instance becomes invalid after this.
  Clustering finalize() {
    return Clustering(std::move(stroke2cluster_), std::move(cluster2stroke_), n_strokes_);
  }

  explicit operator bool() const { return !cluster2stroke_.empty(); }

private:
  std::map<int, std::vector<int>> cluster2stroke_;
  std::unique_ptr<int[]> stroke2cluster_;
  const int n_strokes_;
};

#ifdef HAS_GUROBI

// Helper functions to convert between data structs used by different libraries
void strokes_to_capture(const std::vector<const Stroke*>& strokes, const Float accuracy,
                        const bool cut_strokes, std::vector<StrokeMapping>& out_mappings,
                        Capture& capture);
Input from_capture(Capture& capture, glm::dvec2& out_center);
Input from_capture(Capture& capture);

/**
 * Consolidate a stroke cluster into a single representative stroke.
 *
 * @param strokes Strokes in the stroke cluster.
 * @param accuracy How closely the fitted curve should follow the original strokes. Trades
 *                 fidelity for smoothness.  1.0 is a good starting point.
 * @param cut_strokes Whether to cut strokes at sharp turns.
 * @param[out] out_mappings Mappings from each stroke to the fitted curve.
 */
Stroke fit_stroke_to_cluster(const std::vector<const Stroke*>& strokes, Float accuracy,
                             bool cut_strokes, std::vector<StrokeMapping>& out_mappings);

#endif

} // namespace sketching
