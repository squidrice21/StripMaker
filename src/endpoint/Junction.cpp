#include "Junction.h"

#include "../solve/SolveUtil.h"
#include "JunctionProposal.h"
#include "SketchConnectivity/bvh.h"
#include "SketchConnectivity/classifier.h"
#include "SketchConnectivity/junction_features.h"
#include "SketchConnectivity/types.h"

#include "glm/detail/type_vec.hpp"

#include <limits>
#include <map>
#include <vector>

void convert_stroke(const FittedCurve &fit, sketching::Stroke &junc_stroke) {
  junc_stroke = sketching::Stroke(fit.centerline.size(), false);
  for (size_t i = 0; i < fit.centerline.size(); ++i) {
    assert(fit.widths.size() == fit.centerline.size());
    junc_stroke.x(i) = fit.centerline[i].x;
    junc_stroke.y(i) = fit.centerline[i].y;
    junc_stroke.width(i) = fit.widths[i];

    junc_stroke.index = fit.cluster_idx;
  }
}

void construct_bvh(const std::vector<FittedCurve> &fits,
                   sketching::StrokeBVH &bvh) {
  for (auto const &f : fits) {
    sketching::Stroke s;
    convert_stroke(f, s);
    bvh.add(std::move(s));
  }
}

void prepare_fits(const std::map<size_t, Cluster> &clusters,
                  std::vector<FittedCurve> &fits) {
  for (auto const &cc : clusters) {
    fits.emplace_back(cc.second.fit);
    fits.back().cluster_idx = cc.first;
  }
}

void predict_candidate(const sketching::Stroke &stroke, bool is_head,
                       const FittedCurve &fit, sketching::StrokeBVH &bvh,
                       sketching::StrokeTime cand, Junction &junc) {
  double arclen1 = (is_head) ? 0 : 1;
  arclen1 *= stroke.length();

  std::map<size_t, size_t> sid2bvh;
  for (size_t i = 0; i < bvh.strokes().size(); ++i) {
    sid2bvh[bvh.stroke(i).index] = i;
  }
  static const auto end_end_features = sketching::get_end_end_features();
  static const auto end_stroke_features = sketching::get_end_stroke_features();
  for (auto &feat : end_end_features) {
    feat->init(bvh.polyline_bvh());
  }
  for (auto &feat : end_stroke_features) {
    feat->init(bvh.polyline_bvh());
  }

  sketching::FeatureVector fea_vec;
  auto type = (cand.second == 0 || cand.second == 1)
                ? sketching::JunctionType::Type::R
                : sketching::JunctionType::Type::T;
  if (type == sketching::JunctionType::Type::R) {
    for (size_t fi = 0; fi < sketching::EndEndFeatures::n_features_; ++fi) {
      const auto &feat = *end_end_features[fi];
      const auto &s2 = bvh.stroke(sid2bvh[cand.first]);
      fea_vec.ee_fea_.data_[fi] =
        feat(stroke, arclen1, sid2bvh[fit.cluster_idx], s2,
             cand.second * s2.length(), sid2bvh[cand.first]);
    }
    fea_vec.type_ = sketching::FeatureVector::EndEnd;
  } else if (type == sketching::JunctionType::Type::T) {
    for (size_t fi = 0; fi < sketching::EndStrokeFeatures::n_features_; ++fi) {
      const auto &feat = *end_stroke_features[fi];
      const auto &s2 = bvh.stroke(sid2bvh[cand.first]);
      fea_vec.es_fea_.data_[fi] =
        feat(stroke, arclen1, sid2bvh[fit.cluster_idx], s2,
             cand.second * s2.length(), sid2bvh[cand.first]);
    }
    fea_vec.type_ = sketching::FeatureVector::EndStroke;
  }

  junc.from = sketching::StrokeTime{fit.cluster_idx, (is_head) ? 0 : 1};
  junc.from_pos = (is_head) ? fit.centerline.front() : fit.centerline.back();
  junc.to = cand;
  const auto &s2 = bvh.stroke(sid2bvh[cand.first]);
  const auto v2 = s2.pos(cand.second * s2.length());
  junc.to_pos = glm::dvec2(v2.x(), v2.y());
  junc.type = type;

  if (type == sketching::JunctionType::Type::R) {
    junc.probability = sketching::clf_endpoint(fea_vec.ee_fea_);
  } else if (type == sketching::JunctionType::Type::T) {
    junc.probability = sketching::clf_tjunction(fea_vec.es_fea_);
  }
}

void predict_junctions(const FittedCurve &fit, bool is_head, size_t knn,
                       sketching::StrokeBVH &bvh,
                       std::vector<Junction> &junctions) {
  assert(!fit.widths.empty());
  sketching::Stroke stroke;
  convert_stroke(fit, stroke);
  double arclen1 = (is_head) ? 0 : 1;
  stroke.ensure_arclengths();
  arclen1 *= stroke.length();

  // 1. Propose junction candidates
  std::vector<sketching::StrokeTime> candidates;
  glm::dvec2 v = (is_head) ? fit.centerline.front() : fit.centerline.back();
  propose_end_end_candidates(v, fit.cluster_idx, knn, bvh, candidates);
  propose_end_stroke_candidates(v, fit.cluster_idx, knn, bvh, candidates);

  std::vector<std::pair<sketching::StrokeTime, sketching::StrokeTime>>
    out_intersections;
  propose_end_end_intersections(fit.cluster_idx, bvh, out_intersections);

  if (solve_context.stage == SolveContext::SolveStage::Endpoint) {
    propose_end_end_envelope_intersections(fit.cluster_idx, bvh,
                                           out_intersections);
  }

  // Close to intersection

  // 2. Predict
  std::map<size_t, size_t> sid2bvh;
  for (size_t i = 0; i < bvh.strokes().size(); ++i) {
    sid2bvh[bvh.stroke(i).index] = i;
  }
  for (auto const &cand : candidates) {
    junctions.emplace_back();
    predict_candidate(stroke, is_head, fit, bvh, cand, junctions.back());
  }

  // 3. Add intersections
  for (auto const &intersection : out_intersections) {
    auto compute_pos = [&](const sketching::StrokeTime &st) -> glm::dvec2 {
      const auto &s = bvh.stroke(sid2bvh[st.first]);
      const auto v = s.pos(st.second * s.length());
      return glm::dvec2(v.x(), v.y());
    };
    junctions.emplace_back();
    junctions.back().from = intersection.first;
    junctions.back().from_pos = compute_pos(intersection.first);
    junctions.back().to = intersection.second;
    junctions.back().to_pos = compute_pos(intersection.second);
    junctions.back().type = sketching::JunctionType::Type::X;
    junctions.back().probability = 1;
  }
}
