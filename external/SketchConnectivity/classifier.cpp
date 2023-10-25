#include "classifier.h"

#include "detail/alloca.h"
#include "features/junction_features_impl.h"
#include "force_assert.h"
#include "stroke_graph_extra.h"

#include <spdlog/spdlog.h>

namespace sketching {

namespace {

constexpr auto invalid = StrokeGraph::invalid;

} // namespace

Index n_classifier_evaluations_vv(const StrokeGraph& graph, const Endpoint e1,
                                  const Endpoint e2) {
  auto n_1st_endp_entries = Index(0);
  {
    const auto he = endpoint_to_hedge(graph, e1);
    auto it = he;
    do {
      assert(it.continuity() == invalid);
      if (!(it.flags() & StrokeGraph::HedgeRecord::Bridge)) {
        n_1st_endp_entries++;
      }
      it = it.twin().next();
      assert(n_1st_endp_entries < 1024 && "likely infinite loop found");
    } while (it != he);
  }
  auto n_2nd_endp_entries = Index(0);
  {
    const auto he = endpoint_to_hedge(graph, e2);
    auto it = he;
    do {
      if (!(it.flags() & StrokeGraph::HedgeRecord::Bridge)) {
        // Deduplicate.  Note the condition also works for invalid continuity.
        if (it.continuity() >= it.index_) {
          n_2nd_endp_entries++;
        }
      }
      it = it.twin().next();
      assert(n_2nd_endp_entries < 1024 && "likely infinite loop found");
    } while (it != he);
  }
  return n_1st_endp_entries * n_2nd_endp_entries;
}

Index n_classifier_evaluations_vs(const StrokeGraph& graph, const Endpoint e1) {
  auto count = Index(0);
  const auto he = endpoint_to_hedge(graph, e1);
  auto it = he;
  do {
    assert(it.continuity() == invalid);
    count++;
    it = it.twin().next();
    assert(count < 1024 && "likely infinite loop found");
  } while (it != he);
  return count;
}

void expand_candidates_vv(const StrokeGraph& graph, const Endpoint e1, const Endpoint e2,
                          const span<Endpoint> out_cand1,
                          const span<StrokeTime> out_cand2) {
  assert(out_cand1.size() == out_cand2.size());
  assert((Index)out_cand1.size() == n_classifier_evaluations_vv(graph, e1, e2));

  auto out = size_t(0);
  const auto he1 = endpoint_to_hedge(graph, e1);
  const auto he2 = endpoint_to_hedge(graph, e2);
  auto it1 = he1;
  do {
    auto it2 = he2;
    do {
      assert(it1.continuity() == invalid);
      if (!(it1.flags() & StrokeGraph::HedgeRecord::Bridge) &&
          !(it2.flags() & StrokeGraph::HedgeRecord::Bridge)) {
        // Deduplicate.  Note the condition also works for invalid continuity.
        if (it2.continuity() >= it2.index_) {
          out_cand1[out] = hedge_to_endpoint(it1);
          const auto endp2 = hedge_to_endpoint(it2);
          out_cand2[out] =
            StrokeTime((int)endp2.stroke_idx(), endp2.is_head() ? 0.0 : 1.0);
          // Note that if we have a continuous edge, it will always be in cand2.
          out++;
        }
      }
      it2 = it2.twin().next();
    } while (it2 != he2);
    it1 = it1.twin().next();
  } while (it1 != he1);
  force_assert(out == out_cand1.size() &&
               "expand_candidates_vv out of sync with n_classifier_evaluations_vv");
}

void expand_candidates_vs(const StrokeGraph& graph, const Endpoint e1,
                          const span<Endpoint> out_cand) {
  assert((Index)out_cand.size() == n_classifier_evaluations_vs(graph, e1));

  auto out = size_t(0);
  const auto he = endpoint_to_hedge(graph, e1);
  auto it = he;
  do {
    assert(it.continuity() == invalid);
    out_cand[out++] = hedge_to_endpoint(it);
    it = it.twin().next();
  } while (it != he);
  force_assert(out == out_cand.size() &&
               "expand_candidates_vs out of sync with n_classifier_evaluations_vs");
}

void compute_junction_probabilities(
  const StrokeGraph& graph, const span<const Endpoint> cand1,
  const span<const StrokeTime> cand2, FeatureType feature_type,
  const span<Float> out_prob, std::vector<ClassifierPrediction>* const out_predictions,
  bool expand, const std::vector<JunctionType::Type>& junction_types) {
  assert(cand1.size() == cand2.size());
  assert(cand1.size() == out_prob.size());

  std::vector<Endpoint> expanded_cand1;
  std::vector<StrokeTime> expanded_cand2;
  std::vector<Index> n_evaluations;
  auto n_evaluations_total = Index(0);
  if (expand) {
    // Compute the number of classifier evaluations we need to do.
    n_evaluations = std::vector<Index>(cand2.size());
    for (size_t i = 0; i < cand1.size(); ++i) {
      const auto [si2, narclen2] = cand2[i];
      force_assert(0.0 <= narclen2 && narclen2 <= 1.0 && "invalid normalized arc length");
      if (narclen2 == 0.0 || narclen2 == 1.0) {
        const auto cand2_endp = Endpoint(si2, narclen2 == 0.0);
        n_evaluations[i] = n_classifier_evaluations_vv(graph, cand1[i], cand2_endp);
      } else {
        n_evaluations[i] = n_classifier_evaluations_vs(graph, cand1[i]);
      }
      n_evaluations_total += n_evaluations[i];
    }

    expanded_cand1 = std::vector<Endpoint>(n_evaluations_total, Endpoint(0));
    expanded_cand2 = std::vector<StrokeTime>(n_evaluations_total);
    for (size_t i = 0, k = 0; i < cand1.size(); ++i) {
      const auto [si2, narclen2] = cand2[i];
      const auto n = (size_t)n_evaluations[i];
      if (narclen2 == 0.0 || narclen2 == 1.0) {
        const auto cand2_endp = Endpoint(si2, narclen2 == 0.0);
        expand_candidates_vv(graph, cand1[i], cand2_endp, {&expanded_cand1[k], n},
                             {&expanded_cand2[k], n});
      } else {
        expand_candidates_vs(graph, cand1[i], {&expanded_cand1[k], n});
        for (size_t j = 0; j < n; ++j) {
          expanded_cand2[k + j] = cand2[i];
        }
      }
      k += n;
    }
  } else {
    n_evaluations_total = Index(cand2.size());
    n_evaluations.reserve(cand2.size());
    expanded_cand1.reserve(cand2.size());
    for (auto const& c1 : cand1)
      expanded_cand1.emplace_back(c1);
    expanded_cand2.reserve(cand2.size());
    for (auto const& c2 : cand2) {
      expanded_cand2.emplace_back(c2);
      n_evaluations.emplace_back(1);
    }
  }

  auto feature_mat = std::vector<FeatureVector>(n_evaluations_total);
  if (feature_type == FeatureType::OrigStroke)
    compute_features(graph, expanded_cand1, expanded_cand2, feature_mat, junction_types);
  else if (feature_type == FeatureType::Graph)
    compute_features(graph, graph.bvh_.polyline_bvh(), expanded_cand1, expanded_cand2,
                     feature_mat);
  else
    abort();

  auto expanded_prob = std::vector<Float>(n_evaluations_total);
  for (Index i = 0; i < n_evaluations_total; ++i) {
    if (feature_mat[i].type_ == FeatureVector::EndEnd) {
      expanded_prob[i] = clf_endpoint(feature_mat[i].ee_fea_);
    } else if (feature_mat[i].type_ == FeatureVector::EndStroke) {
      expanded_prob[i] = clf_tjunction(feature_mat[i].es_fea_);
    } else {
      std::abort();
    }
  }

  if (out_predictions) {
    out_predictions->reserve(n_evaluations_total);
    for (Index i = 0; i < n_evaluations_total; ++i) {
      auto& pred_info = out_predictions->emplace_back();
      const auto v1 = endpoint_to_vertex(graph, expanded_cand1[i]);
      // We assume end-stroke predictions have the stroke as the second member.
      const auto [si2, narclen2] = expanded_cand2[i];
      if (narclen2 == 0.0 || narclen2 == 1.0) {
        const auto he2 = endpoint_to_hedge(graph, Endpoint(si2, narclen2 == 0.0));
        const auto v2 = he2.origin();
        // Note that the type will be wrong for end-stroke evaluations across continuous
        // edges, but we want the visualization to be able to group them together. Can't
        // be helped for now...
        pred_info.key.type = JunctionType::R;
        pred_info.key.cand1 = (int)v1.index_;
        pred_info.key.cand2 = (int)v2.index_;
        pred_info.p_b = v2.pos();
      } else {
        pred_info.key.type = JunctionType::T;
        pred_info.key.cand1 = (int)v1.index_;
        pred_info.key.cand2 = (int)si2;
        pred_info.p_b = graph.strokes_[si2].pos_norm(narclen2);
      }
      pred_info.p_a = v1.pos();
      pred_info.prob = expanded_prob[i];
      pred_info.fea = feature_mat[i];

      auto orig_a = expanded_cand1[i].as_pair();
      auto ok = convert_strokes2orig(graph, orig_a);
      force_assert(ok && "failed to map from strokes to orig");
      pred_info.orig_a = orig_a;
      auto orig_b =
        std::make_pair((int)expanded_cand2[i].first, expanded_cand2[i].second);
      ok = convert_strokes2orig(graph, orig_b);
      force_assert(ok && "failed to map from strokes to orig");
      pred_info.orig_b = orig_b;
    }
  }

  // Aggregate probabilities (necessary for high-valence junctions).
  for (size_t i = 0, k = 0; i < cand2.size(); ++i) {
    const auto n = (size_t)n_evaluations[i];
    out_prob[i] = 0.0;
    for (size_t j = 0; j < n; ++j) {
      out_prob[i] = std::max(out_prob[i], expanded_prob[k]);
      k++;
    }
  }
}

std::vector<std::unique_ptr<JunctionFeature>> get_end_end_features() {
  using namespace ::sketching::features;

  auto features = std::vector<std::unique_ptr<JunctionFeature>>();
  features.push_back(std::make_unique<EndEndJunctionType>());
  features.push_back(
    std::make_unique<EnvelopeDistance>(Normalization::PenWidthPairwiseMean));
  features.push_back(
    std::make_unique<EnvelopeDistance>(Normalization::StrokeLengthPairwiseMin));
  features.push_back(
    std::make_unique<EnvelopeDistance>(Normalization::StrokeLengthPairwiseMax));
  features.push_back(std::make_unique<StepawayTangentAngleMin>());
  features.push_back(std::make_unique<StepawayTangentAngleMax>());
  features.push_back(std::make_unique<ClosestAnyOtherOverConnectionMin>());
  features.push_back(std::make_unique<ClosestAnyOtherOverConnectionMax>());
  features.push_back(std::make_unique<ProjectionToEndpointRatioMin>());
  features.push_back(std::make_unique<ProjectionToEndpointRatioMax>());
  features.push_back(std::make_unique<StepawayOverConnectionMin>());
  features.push_back(std::make_unique<StepawayOverConnectionMax>());
  features.push_back(std::make_unique<ProjectionOverConnectionMin>());
  features.push_back(std::make_unique<ProjectionOverConnectionMax>());
  assert(features.size() == EndEndFeatures::n_features_);
  return features;
}

std::vector<std::unique_ptr<JunctionFeature>> get_end_stroke_features() {
  using namespace ::sketching::features;

  auto features = std::vector<std::unique_ptr<JunctionFeature>>();
  features.push_back(
    std::make_unique<EnvelopeDistance>(Normalization::PenWidthPairwiseMean));
  features.push_back(std::make_unique<EnvelopeDistance>(Normalization::StrokeLength1));
  features.push_back(std::make_unique<EnvelopeDistance>(Normalization::StrokeLength2));
  features.push_back(std::make_unique<StepawayTangentAngle1>());
  features.push_back(std::make_unique<Busyness1>(1.0));
  features.push_back(
    std::make_unique<ClosestAnyOtherOverConnection1>(/*limit_to_visible=*/true));
  features.push_back(std::make_unique<ConnectedDistanceToEndpoint>());
  features.push_back(std::make_unique<StepawayOverConnection1>());
  assert(features.size() == EndStrokeFeatures::n_features_);
  return features;
}

std::vector<std::string> get_end_end_feature_descriptions() {
  static const auto features = get_end_end_features();
  auto descriptions = std::vector<std::string>();
  for (const auto& feat : features) {
    descriptions.emplace_back(feat->description());
  }
  return descriptions;
}

std::vector<std::string> get_end_stroke_feature_descriptions() {
  static const auto features = get_end_stroke_features();
  auto descriptions = std::vector<std::string>();
  for (const auto& feat : features) {
    descriptions.emplace_back(feat->description());
  }
  return descriptions;
}

void human_readable_end_end_features(const span<Float> feature_vec) {
  static const auto features = get_end_end_features();
  for (size_t fi = 0; fi < EndEndFeatures::n_features_; ++fi) {
    const auto& feat = features[fi];
    feat->human_readable({&feature_vec[fi], 1});
  }
}

void human_readable_end_stroke_features(const span<Float> feature_vec) {
  static const auto features = get_end_stroke_features();
  for (size_t fi = 0; fi < EndStrokeFeatures::n_features_; ++fi) {
    const auto& feat = features[fi];
    feat->human_readable({&feature_vec[fi], 1});
  }
}

void compute_features(const StrokeGraph& graph, const span<const Endpoint> cand1,
                      const span<const StrokeTime> cand2,
                      const span<FeatureVector> out_features,
                      const std::vector<JunctionType::Type>& junction_types) {
  assert(cand1.size() == cand2.size());
  assert(cand1.size() == out_features.size());
  static const auto end_end_features = get_end_end_features();
  static const auto end_stroke_features = get_end_stroke_features();

  for (auto& feat : end_end_features) {
    feat->init(graph);
  }
  for (auto& feat : end_stroke_features) {
    feat->init(graph);
  }

  for (size_t i = 0; i < cand1.size(); ++i) {
    force_assert(!(graph.hedges_[2 * cand1[i].stroke_idx()].flags_ &
                   StrokeGraph::HedgeRecord::Bridge));
    force_assert(
      !(graph.hedges_[2 * cand2[i].first].flags_ & StrokeGraph::HedgeRecord::Bridge));
    const auto is_stroke_interior = (cand2[i].second != 0.0 && cand2[i].second != 1.0);
    const auto orig_pos1 = as_orig_position(graph, cand1[i].as_pair());
    auto orig_pos2 = StrokeTime();
    if (is_stroke_interior) {
      orig_pos2 = reproject_cand2_on_original(graph, orig_pos1, cand2[i]);
    } else {
      orig_pos2 = as_orig_position(graph, cand2[i]);
    }

    const auto idx1 = orig_pos1.first;
    const auto idx2 = orig_pos2.first;
    const auto& s1 = graph.orig_strokes_[idx1];
    const auto& s2 = graph.orig_strokes_[idx2];
    const auto arclen1 = orig_pos1.second * s1.length();
    const auto arclen2 = orig_pos2.second * s2.length();
    const auto endp2 =
      Endpoint(cand2[i].first, cand2[i].second < 0.5); // cand2 might not be an endpoint.

    bool is_end_stroke =
      is_stroke_interior || endpoint_to_hedge(graph, endp2).continuity() != invalid;
    if (junction_types.size() == cand1.size())
      is_end_stroke = (junction_types[i] == JunctionType::T);
    if (is_end_stroke) {
      for (size_t fi = 0; fi < EndStrokeFeatures::n_features_; ++fi) {
        const auto& feat = *end_stroke_features[fi];
        out_features[i].es_fea_.data_[fi] = feat(s1, arclen1, idx1, s2, arclen2, idx2);
      }
      out_features[i].type_ = FeatureVector::EndStroke;
    } else {
      for (size_t fi = 0; fi < EndEndFeatures::n_features_; ++fi) {
        const auto& feat = *end_end_features[fi];
        out_features[i].ee_fea_.data_[fi] = feat(s1, arclen1, idx1, s2, arclen2, idx2);
      }
      out_features[i].type_ = FeatureVector::EndEnd;
    }
  }
}

void compute_features(const StrokeGraph& graph, const PolylineBVH& bvh,
                      const span<const Endpoint> cand1,
                      const span<const StrokeTime> cand2,
                      const span<FeatureVector> out_features) {
  assert(cand1.size() == cand2.size());
  assert(cand1.size() == out_features.size());
  static const auto end_end_features = get_end_end_features();
  static const auto end_stroke_features = get_end_stroke_features();

  for (auto& feat : end_end_features) {
    feat->init(bvh);
  }
  for (auto& feat : end_stroke_features) {
    feat->init(bvh);
  }

  for (size_t i = 0; i < cand1.size(); ++i) {
    const auto idx1 = cand1[i].stroke_idx();
    const auto idx2 = cand2[i].first;
    const auto& s1 = *bvh.nodes[idx1].geometry;
    const auto& s2 = *bvh.nodes[idx2].geometry;
    const auto arclen1 = (cand1[i].is_head() ? 0.0 : s1.length());
    const auto arclen2 = cand2[i].second * s2.length();
    const auto endp2 = Endpoint(idx2, arclen2 == 0.0); // cand2 might not be an endpoint.

    if (arclen2 != 0.0 && arclen2 != s2.length()) {
      for (size_t fi = 0; fi < EndStrokeFeatures::n_features_; ++fi) {
        const auto& feat = *end_stroke_features[fi];
        out_features[i].es_fea_.data_[fi] = feat(s1, arclen1, idx1, s2, arclen2, idx2);
      }
      out_features[i].type_ = FeatureVector::EndStroke;
    } else if (endpoint_to_hedge(graph, endp2).continuity() != invalid) {
      const auto he = endpoint_to_hedge(graph, endp2);
      Float join_arclen;
      const auto dissolved = concatenate_edges(graph, he.continuity_edge().twin().index_,
                                               he.index_, join_arclen);
      force_assert(dissolved.size() > 0);
      for (size_t fi = 0; fi < EndStrokeFeatures::n_features_; ++fi) {
        const auto& feat = *end_stroke_features[fi];
        out_features[i].es_fea_.data_[fi] =
          feat(s1, arclen1, idx1, dissolved, join_arclen, (size_t)-1);
      }
      out_features[i].type_ = FeatureVector::EndStroke;
    } else {
      for (size_t fi = 0; fi < EndEndFeatures::n_features_; ++fi) {
        const auto& feat = *end_end_features[fi];
        out_features[i].ee_fea_.data_[fi] = feat(s1, arclen1, idx1, s2, arclen2, idx2);
      }
      out_features[i].type_ = FeatureVector::EndEnd;
    }
  }
}

} // namespace sketching
