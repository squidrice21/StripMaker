#pragma once

#include "../closest.h"
#include "../eigen_compat.h"
#include "../fitting.h"
#include "../junction_features.h"
#include "../sketching.h"
#include <string>

namespace sketching::features {

enum class Normalization {
  BoundingBox,
  PenWidth1,
  PenWidth2,
  PenWidthPairwiseMean,
  PenWidthPairwiseMax,
  StrokeLength1,
  StrokeLength2,
  StrokeLengthPairwiseMean,
  StrokeLengthPairwiseMin,
  StrokeLengthPairwiseMax,
};

std::string description(Normalization n);

Vec2 stepaway_tangent(const Stroke &s1, const bool s1_head,
                      const Float stepaway_amount);

///////////////////
// Meta-features //
///////////////////

/// Feature which reverses the roles of candidate 1 and candidate 2.
template <typename Feature>
struct ReverseFeature : JunctionFeature {
  template <typename... Args>
  explicit ReverseFeature(Args &&...args)
    : feature_(std::forward<Args>(args)...) {}

  void init(const PolylineBVH &bvh) final {
    static_cast<JunctionFeature *>(&feature_)->init(bvh);
  }
  void init(const StrokeGraph &graph) final {
    static_cast<JunctionFeature *>(&feature_)->init(graph);
  }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2,
                   size_t idx2) const override {
    return feature_(s2, arclen2, idx2, s1, arclen1, idx1);
  }

  void human_readable(span<Float> v) const override {
    feature_.human_readable(v);
  }

private:
  Feature feature_;
};

template <typename Feature>
struct MinFeature : JunctionFeature {
  template <typename... Args>
  explicit MinFeature(Args &&...args)
    : feature_(std::forward<Args>(args)...) {}

  void init(const PolylineBVH &bvh) final {
    static_cast<JunctionFeature *>(&feature_)->init(bvh);
  }
  void init(const StrokeGraph &graph) final {
    static_cast<JunctionFeature *>(&feature_)->init(graph);
  }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2,
                   size_t idx2) const override {
    return std::min(feature_(s1, arclen1, idx1, s2, arclen2, idx2),
                    feature_(s2, arclen2, idx2, s1, arclen1, idx1));
  }

  void human_readable(span<Float> v) const override {
    feature_.human_readable(v);
  }

private:
  Feature feature_;
};

template <typename Feature>
struct MaxFeature : JunctionFeature {
  template <typename... Args>
  explicit MaxFeature(Args &&...args)
    : feature_(std::forward<Args>(args)...) {}

  void init(const PolylineBVH &bvh) final {
    static_cast<JunctionFeature *>(&feature_)->init(bvh);
  }
  void init(const StrokeGraph &graph) final {
    static_cast<JunctionFeature *>(&feature_)->init(graph);
  }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2,
                   size_t idx2) const override {
    return std::max(feature_(s1, arclen1, idx1, s2, arclen2, idx2),
                    feature_(s2, arclen2, idx2, s1, arclen1, idx1));
  }

  void human_readable(span<Float> v) const override {
    feature_.human_readable(v);
  }

private:
  Feature feature_;
};

template <typename Feature>
struct MeanFeature : JunctionFeature {
  void init(const PolylineBVH &bvh) final {
    static_cast<JunctionFeature *>(&feature_)->init(bvh);
  }
  void init(const StrokeGraph &graph) final {
    static_cast<JunctionFeature *>(&feature_)->init(graph);
  }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2,
                   size_t idx2) const override {
    return 0.5 * (feature_(s1, arclen1, idx1, s2, arclen2, idx2) +
                  feature_(s2, arclen2, idx2, s1, arclen1, idx1));
  }

  void human_readable(span<Float> v) const override {
    feature_.human_readable(v);
  }

private:
  Feature feature_;
};

///////////////////////
// Absolute features //
///////////////////////

struct AbsCenterlineDistance final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t /*idx1*/, //
                   const Stroke &s2, Float arclen2,
                   size_t /*idx2*/) const final {
    return (s1.pos(arclen1) - s2.pos(arclen2)).norm();
  }

  std::string description() const final {
    return "Absolute centerline connection dist";
  }
};

struct AbsEnvelopeDistance final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Absolute envelope connection dist";
  }
};

struct AbsStrokeLength1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float /*arclen1*/, size_t /*idx1*/, //
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t /*idx2*/) const final {
    return s1.length();
  }

  std::string description() const final { return "Absolute stroke length 1"; }
};

struct AbsStrokeLength2 final : ReverseFeature<AbsStrokeLength1> {
  std::string description() const final { return "Absolute stroke length 2"; }
};

struct AbsPenWidth1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float /*arclen1*/, size_t /*idx1*/, //
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t /*idx2*/) const final {
    return s1.pen_width();
  }

  std::string description() const final {
    return "Absolute pen width of first stroke";
  }
};

struct AbsPenWidth2 final : ReverseFeature<AbsPenWidth1> {
  std::string description() const final {
    return "Absolute pen width of second stroke";
  }
};

struct AbsWidth1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t /*idx1*/, //
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t /*idx2*/) const final {
    return s1.width_at(arclen1);
  }

  std::string description() const final {
    return "Absolute width at 1st point";
  }
};

struct AbsWidth2 final : ReverseFeature<AbsWidth1> {
  std::string description() const final {
    return "Absolute width at 2nd point";
  }
};

struct AbsProjectionDist1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Abs projection distance of 1st endp";
  }
};

struct AbsProjectionDist2 final : ReverseFeature<AbsProjectionDist1> {
  std::string description() const final {
    return "Abs projection distance of 2nd endp";
  }
};

struct AbsProjectionToClosestEndp1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Abs dist along stroke from proj of 1st endp to nearest endp";
  }
};

struct AbsProjectionToClosestEndp2 final
  : ReverseFeature<AbsProjectionToClosestEndp1> {
  std::string description() const final {
    return "Abs dist along stroke from proj of 2nd endp to nearest endp";
  }
};

struct AbsStepawayDist1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Absolute stepaway projection dist (1)";
  }
};

struct AbsStepawayDist2 final : ReverseFeature<AbsStepawayDist1> {
  std::string description() const final {
    return "Absolute stepaway projection dist (2)";
  }
};

//////////////////////
// Regular features //
//////////////////////

struct EndEndJunctionType final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final { return "End-end junction type"; }
};

struct EndStrokeJunctionType final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final { return "End-stroke junction type"; }
};

struct EnvelopeDistance final : JunctionFeature {
  explicit EnvelopeDistance(const Normalization scheme)
    : m_normalization_scheme(scheme) {}

  void init(const PolylineBVH &bvh) final;
  void init(const StrokeGraph &graph) final;

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1,
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Envelope dist divided by " +
           ::sketching::features::description(m_normalization_scheme);
  }

private:
  Float m_inv_scale = std::numeric_limits<Float>::quiet_NaN();
  Normalization m_normalization_scheme;
  const StrokeGraph *graph_ptr_ = nullptr;
};

struct ProjectionOverConnection1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1,
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Dist to projection divided by conn dist (1)";
  }
};

struct ProjectionOverConnection2 final
  : ReverseFeature<ProjectionOverConnection1> {
  std::string description() const final {
    return "Dist to projection divided by conn dist (2)";
  }
};

struct ProjectionOverConnectionMin final
  : MinFeature<ProjectionOverConnection1> {
  std::string description() const final {
    return "Dist to projection divided by conn dist (min)";
  }
};

struct ProjectionOverConnectionMax final
  : MaxFeature<ProjectionOverConnection1> {
  std::string description() const final {
    return "Dist to projection divided by conn dist (max)";
  }
};

/**
 * Return the shortest distance from the connecting point to any other stroke,
 * divided by the connection distance.
 */
struct ClosestAnyOverConnection1 final : JunctionFeature {
  void init(const PolylineBVH &bvh) final { bvh_ = bvh; }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1,
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Dist to closest stroke divided by conn dist (1)";
  }

  PolylineBVH bvh_;
};

/**
 * Return the shortest distance from the connected-to point to any other stroke,
 * divided by the connection distance.
 */
struct ClosestAnyOverConnection2 final
  : ReverseFeature<ClosestAnyOverConnection1> {
  std::string description() const final {
    return "Dist to closest stroke divided by conn dist (2)";
  }
};

struct ClosestAnyOverConnectionMin final
  : MinFeature<ClosestAnyOverConnection1> {
  std::string description() const final {
    return "Dist to closest stroke divided by conn dist (min)";
  }
};

struct ClosestAnyOverConnectionMax final
  : MaxFeature<ClosestAnyOverConnection1> {
  std::string description() const final {
    return "Dist to closest stroke divided by conn dist (max)";
  }
};

struct ClosestAnyOtherOverConnection1 final : JunctionFeature {
  explicit ClosestAnyOtherOverConnection1(bool limit_to_visible = false)
    : limit_to_visible_(limit_to_visible) {}

  void init(const PolylineBVH &bvh) final { bvh_ = bvh; }
  void init(const StrokeGraph &graph) final { graph_ptr_ = &graph; }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  Vec2 closest(const Stroke &s1, Float arclen1, size_t idx1, //
               const Stroke &s2, Float arclen2, size_t idx2) const;

  std::string description() const final {
    return "Dist to closest any other " +
           std::string(limit_to_visible_ ? "visible " : "") +
           "stroke divided by conn dist (1)";
  }

  bool limit_to_visible_;
  PolylineBVH bvh_;
  const StrokeGraph *graph_ptr_ = nullptr;
};

struct ClosestAnyOtherOverConnection2 final
  : ReverseFeature<ClosestAnyOtherOverConnection1> {
  explicit ClosestAnyOtherOverConnection2(bool limit_to_visible = false)
    : ReverseFeature(limit_to_visible)
    , limit_to_visible_(limit_to_visible) {}

  std::string description() const final {
    return "Dist to closest any other " +
           std::string(limit_to_visible_ ? "visible " : "") +
           "stroke divided by conn dist (2)";
  }

  bool limit_to_visible_;
};

struct ClosestAnyOtherOverConnectionMin final
  : MinFeature<ClosestAnyOtherOverConnection1> {
  explicit ClosestAnyOtherOverConnectionMin(bool limit_to_visible = false)
    : MinFeature(limit_to_visible)
    , limit_to_visible_(limit_to_visible) {}

  std::string description() const final {
    return "Dist to closest any other " +
           std::string(limit_to_visible_ ? "visible " : "") +
           "stroke divided by conn dist (min)";
  }

  bool limit_to_visible_;
};

struct ClosestAnyOtherOverConnectionMax final
  : MaxFeature<ClosestAnyOtherOverConnection1> {
  explicit ClosestAnyOtherOverConnectionMax(bool limit_to_visible = false)
    : MaxFeature(limit_to_visible)
    , limit_to_visible_(limit_to_visible) {}

  std::string description() const final {
    return "Dist to closest any other " +
           std::string(limit_to_visible_ ? "visible " : "") +
           "stroke divided by conn dist (max)";
  }

  bool limit_to_visible_;
};

struct OtherEndpointClosestAnyEnvOverEnvConnection1 final : JunctionFeature {
  void init(const PolylineBVH &bvh) final { bvh_ = bvh; }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Other endp env dist to closest stroke divided by env conn dist (1)";
  }

  PolylineBVH bvh_;
};

struct OtherEndpointClosestAnyEnvOverEnvConnection2 final
  : ReverseFeature<OtherEndpointClosestAnyEnvOverEnvConnection1> {
  std::string description() const final {
    return "Other endp env dist to closest stroke divided by env conn dist (2)";
  }
};

struct ClosestEndpointOverConnection1 final : JunctionFeature {
  void init(const PolylineBVH &bvh) final { bvh_ = bvh; }

  Vec2 closest(const Stroke &s1, Float arclen1, size_t idx1, //
               const Stroke &s2, Float arclen2, size_t idx2) const;

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Distance to closest visible endpoint divided by conn dist (1)";
  }

  PolylineBVH bvh_;
};

struct StepawayOverConnection1 final : JunctionFeature {
  explicit StepawayOverConnection1(Float factor = 1.0)
    : factor_(factor) {}

  void init(const PolylineBVH &) final {}
  void init(const StrokeGraph &graph) final { graph_ptr_ = &graph; }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1,
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    if (factor_ == 1.0) {
      return "Dist when stepping away divided by conn dist (1)";
    } else if (std::ceil(factor_) == factor_) { // Integer.
      return "Dist when stepping away (" + std::to_string(int(factor_)) +
             ") divided by conn dist (1)";
    }
    return "Dist when stepping away (" + std::to_string(factor_) +
           ") divided by conn dist (1)";
  }

private:
  Float factor_;
  const StrokeGraph *graph_ptr_ = nullptr;
};

struct StepawayOverConnection2 final : ReverseFeature<StepawayOverConnection1> {
  explicit StepawayOverConnection2(Float factor = 1.0)
    : ReverseFeature(factor) {}

  std::string description() const final {
    return "Dist when stepping away divided by conn dist (2)";
  }
};

struct StepawayOverConnectionMin final : MinFeature<StepawayOverConnection1> {
  explicit StepawayOverConnectionMin(Float factor = 1.0)
    : MinFeature(factor) {}

  std::string description() const final {
    return "Dist when stepping away divided by conn dist (min)";
  }
};

struct StepawayOverConnectionMax final : MaxFeature<StepawayOverConnection1> {
  explicit StepawayOverConnectionMax(Float factor = 1.0)
    : MaxFeature(factor) {}

  std::string description() const final {
    return "Dist when stepping away divided by conn dist (max)";
  }
};

struct NearestEndpointOverStepaway1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Dist from stepaway proj to nearest endp divided by dist when "
           "stepping away "
           "(1)";
  }
};

struct NearestEndpointOverStepaway2 final
  : ReverseFeature<NearestEndpointOverStepaway1> {
  std::string description() const final {
    return "Dist from stepaway proj to nearest endp divided by dist when "
           "stepping away "
           "(2)";
  }
};

struct NearestEndpointOverStepawayMin final
  : MinFeature<NearestEndpointOverStepaway1> {
  std::string description() const final {
    return "Dist from stepaway proj to nearest endp divided by dist when "
           "stepping away "
           "(min)";
  }
};

struct NearestEndpointOverStepawayMax final
  : MaxFeature<NearestEndpointOverStepaway1> {
  std::string description() const final {
    return "Dist from stepaway proj to nearest endp divided by dist when "
           "stepping away "
           "(max)";
  }
};

struct StepawayOverProjection1 final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Dist when stepping away divided by endp projection dist (1)";
  }
};

struct StepawayOverProjection2 final : ReverseFeature<StepawayOverProjection1> {
  std::string description() const final {
    return "Dist when stepping away divided by endp projection dist (2)";
  }
};

struct StepawayOverProjectionMin final : MinFeature<StepawayOverProjection1> {
  std::string description() const final {
    return "Dist when stepping away divided by endp projection dist (min)";
  }
};

struct StepawayOverProjectionMax final : MaxFeature<StepawayOverProjection1> {
  std::string description() const final {
    return "Dist when stepping away divided by endp projection dist (max)";
  }
};

struct TangentAngle : JunctionFeature {
  void init(const PolylineBVH &) final;

  const Eigen::Matrix<Float, Eigen::Dynamic, 4> &tangents() const {
    return m_tangents;
  }

  void human_readable(span<Float> v) const final {
    for (auto &x : v) {
      x *= 180 / M_PI;
    }
  }

  Vec2 head_tangent(Index stroke_idx) const;

  Vec2 tail_tangent(Index stroke_idx) const;

protected:
  Eigen::Matrix<Float, Eigen::Dynamic, 4> m_tangents;
};

/**
 * Return the angle between the endpoint tangent of the first stroke and the
 * point on the other stroke.
 */
struct TangentAngle1 final : TangentAngle {
  Float operator()(const Stroke &s1, Float arclen1, size_t idx1,
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Angle between 1st endpoint tangent, connection";
  }
};

/**
 * Return the angle between the endpoint tangent of the second stroke and the
 * point on the other stroke.
 */
struct TangentAngle2 final : TangentAngle {
  Float operator()(const Stroke &s1, Float arclen1, size_t idx1,
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Angle between 2nd endpoint tangent, connection";
  }
};

struct StepawayTangentAngle1 : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  void human_readable(span<Float> v) const final {
    for (auto &x : v) {
      x *= 180 / M_PI;
    }
  }

  std::string description() const final {
    return "Angle between 1st endp tangent (stepaway), conn.";
  }
};

struct StepawayTangentAngle2 final : ReverseFeature<StepawayTangentAngle1> {
  std::string description() const final {
    return "Angle between 2nd endp tangent (stepaway), conn.";
  }
};

struct StepawayTangentAngleMin final : MinFeature<StepawayTangentAngle1> {
  std::string description() const final {
    return "Min angle between endpoint tangents (stepaway), conn.";
  }
};

struct StepawayTangentAngleMax final : MaxFeature<StepawayTangentAngle1> {
  std::string description() const final {
    return "Max angle between endpoint tangents (stepaway), conn.";
  }
};

struct InteriorTangentAngle final : JunctionFeature {
  void init(const PolylineBVH &bvh) final {
    fits_.clear();
    fits_.resize(bvh.nodes.size() + 1);
  }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  void human_readable(span<Float> v) const final {
    for (auto &x : v) {
      x *= 180 / M_PI;
    }
  }

  /// For testing purposes only.  They are only valid for fits already stored in
  /// fits_.
  Vec2 head_tangent(size_t stroke_idx) const;

  /// For testing purposes only.  They are only valid for fits already stored in
  /// fits_.
  Vec2 tail_tangent(size_t stroke_idx) const;

  std::string description() const final {
    return "Smallest angle between interior tangent, connection";
  }

private:
  /// The fit for stroke i is at index i + 1.
  /// Index 0 holds a temporary, so that we can always return a reference in
  /// smooth_fit.
  mutable std::vector<BezierSpline> fits_;

  const BezierSpline &smooth_fit(const Stroke &, size_t stroke_idx) const;
};

struct ClosestDistanceOnExtension1 final : JunctionFeature {
  void init(const PolylineBVH &bvh) final { bvh_ = bvh; }

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Closest dist on extension divided by conn dist (1)";
  }

private:
  PolylineBVH bvh_;
};

struct ClosestDistanceOnExtension2 final
  : ReverseFeature<ClosestDistanceOnExtension1> {
  std::string description() const final {
    return "Closest dist on extension divided by conn dist (2)";
  }
};

struct ClosestDistanceOnExtensionMin final
  : MinFeature<ClosestDistanceOnExtension1> {
  std::string description() const final {
    return "Closest dist on extension divided by conn dist (min)";
  }
};

struct ClosestDistanceOnExtensionMax final
  : MaxFeature<ClosestDistanceOnExtension1> {
  std::string description() const final {
    return "Closest dist on extension divided by conn dist (max)";
  }
};

/**
 * Return the value of the busyness metric at the endpoint of the first stroke.
 */
struct Busyness1 final : JunctionFeature {
  explicit Busyness1(const Float busyness_falloff)
    : busyness_falloff_(busyness_falloff) {}

  void init(const PolylineBVH &) final;

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const final;

  std::string description() const final {
    return "Busyness at connecting point";
  }

private:
  PolylineBVH bvh_;
  Float busyness_falloff_;
  std::unique_ptr<Float[]> memory_;
  Float *head_busyness_ = nullptr;
  Float *tail_busyness_ = nullptr;
  Float *avg_widths_ = nullptr;
};

/**
 * Return the value of the busyness metric at the endpoint of the second stroke.
 */
struct Busyness2 final : ReverseFeature<Busyness1> {
  explicit Busyness2(const Float busyness_falloff)
    : ReverseFeature(busyness_falloff) {}

  std::string description() const final {
    return "Busyness at connected-to point";
  }
};

struct BusynessMin final : MinFeature<Busyness1> {
  explicit BusynessMin(const Float busyness_falloff)
    : MinFeature(busyness_falloff) {}

  std::string description() const final { return "Min busyness at two points"; }
};

struct BusynessMax final : MaxFeature<Busyness1> {
  explicit BusynessMax(const Float busyness_falloff)
    : MaxFeature(busyness_falloff) {}

  std::string description() const final { return "Max busyness at two points"; }
};

/**
 * Return the width of the first stroke, relative to the average.
 */
struct PenWidth1 final : JunctionFeature {
  void init(const PolylineBVH &) final;

  Float operator()(const Stroke &s1, Float /*arclen1*/, size_t /*idx1*/, //
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t /*idx2*/) const final {
    const auto pen_width = s1.pen_width() * drawing_inv_average_width_;
    assert(pen_width <= 10000 && pen_width > 0);
    return pen_width;
  }

  std::string description() const final {
    return "Width of stroke 1 relative to avg";
  }

private:
  Float drawing_inv_average_width_ = 0.0;
};

/**
 * Return the width of the second stroke, relative to the average.
 */
struct PenWidth2 final : ReverseFeature<PenWidth1> {
  std::string description() const final {
    return "Width of stroke 2 relative to avg";
  }
};

struct PenWidthMax final : MaxFeature<PenWidth1> {
  std::string description() const final {
    return "Max width of both strokes relative to avg";
  }
};

struct PenWidthMean final : MeanFeature<PenWidth1> {
  std::string description() const final {
    return "Mean width of both strokes relative to avg";
  }
};

/**
 * Return the distance from the connecting point to the closest endpoint along
 * the stroke being connected to, as a proportion of the total length of the
 * stroke.
 */
struct ConnectedDistanceToEndpoint final : JunctionFeature {
  void init(const PolylineBVH &) override {}
  void init(const StrokeGraph &graph) final { graph_ptr_ = &graph; }

  Float operator()(const Stroke & /*s1*/, Float /*arclen1*/, size_t /*idx1*/, //
                   const Stroke &s2, Float arclen2,
                   size_t /*idx2*/) const override;

  std::string description() const override {
    return "Connected-to location as ratio";
  }

  const StrokeGraph *graph_ptr_ = nullptr;
};

struct ConnectedLocationOverConnection final : JunctionFeature {
  void init(const PolylineBVH &) override {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const override;

  std::string description() const override {
    return "Connected-to location over conn. dist.";
  }
};

/**
 * Return the distance from the closest point to the closest endpoint along the
 * stroke being connected to, as a proportion of the total length of the stroke.
 */
struct ProjectionToEndpointRatio1 final : JunctionFeature {
  void init(const PolylineBVH &) override {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const override;

  std::string description() const override {
    return "Projection location of 1st endpoint as ratio";
  }
};

struct ProjectionToEndpointRatio2 final
  : ReverseFeature<ProjectionToEndpointRatio1> {
  std::string description() const final {
    return "Projection location of 2nd endpoint as ratio";
  }
};

struct ProjectionToEndpointRatioMin final
  : MinFeature<ProjectionToEndpointRatio1> {
  std::string description() const final {
    return "Projection location of endpoints as ratio (min)";
  }
};

struct ProjectionToEndpointRatioMax final
  : MaxFeature<ProjectionToEndpointRatio1> {
  std::string description() const final {
    return "Projection location of endpoints as ratio (max)";
  }
};

struct ProjectionToEndpointOverConnection1 final : JunctionFeature {
  void init(const PolylineBVH &) override {}

  Float operator()(const Stroke &s1, Float arclen1, size_t idx1, //
                   const Stroke &s2, Float arclen2, size_t idx2) const override;

  std::string description() const override {
    return "Projection location of endpoint over conn. dist (1)";
  }
};

struct ProjectionToEndpointOverConnection2 final
  : ReverseFeature<ProjectionToEndpointOverConnection1> {
  std::string description() const override {
    return "Projection location of endpoint over conn. dist (2)";
  }
};

struct ProjectionToEndpointOverConnectionMax final
  : MaxFeature<ProjectionToEndpointOverConnection1> {
  std::string description() const override {
    return "Projection location of endpoint over conn. dist (max)";
  }
};

struct ProjectionToEndpointOverConnectionMin final
  : MinFeature<ProjectionToEndpointOverConnection1> {
  std::string description() const override {
    return "Projection location of endpoint over conn. dist (min)";
  }
};

struct DrawingOrder final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke & /*s1*/, Float /*arclen1*/, size_t idx1,
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t idx2) const final;

  std::string description() const final { return "Binary drawing order"; }
};

struct DrawingOrderUnsignedDiff final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke & /*s1*/, Float /*arclen1*/, size_t idx1,
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t idx2) const final;

  std::string description() const final {
    return "Drawing order unsigned difference";
  }
};

struct DrawingOrderBucketedUnsignedDiff final : JunctionFeature {
  void init(const PolylineBVH &) final {}

  Float operator()(const Stroke & /*s1*/, Float /*arclen1*/, size_t idx1,
                   const Stroke & /*s2*/, Float /*arclen2*/,
                   size_t idx2) const final;

  std::string description() const final {
    return "Drawing order bucketed unsigned difference";
  }
};

Vec2 stepaway_tangent(const Stroke &s1, const bool s1_head,
                      const Float stepaway_amount);

} // namespace sketching::features
