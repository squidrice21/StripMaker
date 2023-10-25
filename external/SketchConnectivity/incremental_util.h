#pragma once

#include "incremental.h"
#include "stroke_graph.h"
#include "types.h"

#include <vector>

namespace sketching {
struct Junction;
struct Stroke;

struct GraphState {
  int prev_state_idx_ = -1;
  size_t state_idx_ = 0;
  size_t curr_stroke_step_ = 0;
  Float obj_ = 0;
  StrokeGraph graph_;
  std::vector<int> face_colors_;
  IncrementalCache cache_;
  std::vector<Junction> candidates_;
  std::vector<bool> connectivity_;
  std::unordered_map<std::string, Float> junc_distance_map_;
  std::unordered_map<std::string, Float> region_size_cache_;

  GraphState()
    : graph_(StrokeGraph::SnappingType::Deformation) {}
  GraphState(size_t curr_stroke_step, const StrokeGraph &graph,
             const IncrementalCache &cache,
             const std::vector<Junction> &candidates,
             const std::vector<bool> &connectivities,
             const std::unordered_map<std::string, Float> &junc_distance_map,
             const std::unordered_map<std::string, Float> &region_size_cache)
    : curr_stroke_step_(curr_stroke_step)
    , graph_(graph.clone())
    , cache_(cache)
    , candidates_(candidates)
    , connectivity_(connectivities)
    , junc_distance_map_(junc_distance_map)
    , region_size_cache_(region_size_cache) {
    for (auto &s : graph_.strokes_)
      s.ensure_arclengths();
  }
  GraphState(const GraphState &other) { *this = other; }

  GraphState &operator=(const GraphState &other) {
    prev_state_idx_ = other.prev_state_idx_;
    state_idx_ = other.state_idx_;
    obj_ = other.obj_;
    curr_stroke_step_ = other.curr_stroke_step_;
    graph_ = other.graph_.clone();
    for (auto &s : graph_.strokes_)
      s.ensure_arclengths();
    face_colors_ = other.face_colors_;
    cache_ = other.cache_;
    candidates_ = other.candidates_;
    connectivity_ = other.connectivity_;
    junc_distance_map_ = other.junc_distance_map_;
    region_size_cache_ = other.region_size_cache_;

    return *this;
  }
};

struct RegionSolution {
  StrokeGraph graph_;
  IncrementalCache cache_;
  std::vector<Junction> candidates_;
  std::vector<bool> connectivity_;
  std::unordered_map<std::string, Float> junc_distance_map_;
  std::unordered_map<std::string, Float> region_size_cache_;
};

int to_v_idx(const StrokeGraph &stroke_graph, const StrokeTime &end);

void print_junctions(const StrokeGraph &stroke_graph,
                     const std::vector<Junction> &candidates);
void print_junctions(const std::vector<Junction> &candidates);

int end_degree(const StrokeGraph &stroke_graph, const StrokeTime &end);

bool has_bridge_vertex(const StrokeGraph &stroke_graph, const Junction &junc);

bool is_snap_valid(const StrokeGraph &stroke_graph, const Junction &junc);

bool is_high_valence_junction_on_boundary_cycle(
  const StrokeGraph &graph, const span<const Junction> candidates, size_t fi,
  const StrokeGraph::VertexView v, const std::string &vid_str);
bool is_high_valence_junction_in_region_cycle(
  const StrokeGraph &graph, const span<const Junction> candidates, size_t fi,
  const StrokeGraph::VertexView v, const std::string &vid_str);

bool is_candidate_valid(const StrokeGraph &stroke_graph,
                        span<const Stroke> strokes,
                        const StrokeGraph &final_stroke_graph,
                        size_t cur_stroke_step, const Junction &junc,
                        const std::vector<Junction> &final_predictions);
bool is_weak_candidate_valid(const StrokeGraph &stroke_graph,
                             span<const Stroke> strokes,
                             const StrokeGraph &final_stroke_graph,
                             size_t cur_stroke_step, const Junction &junc);

void expand_junction(const StrokeGraph &graph, const Junction &in_junc,
                     std::vector<Junction> &bind_junc);
void expand_adjacent_junction(const StrokeGraph &graph, const Junction &in_junc,
                              std::vector<Junction> &bind_junc);

bool update_disconnected_junction_predictions(
  const StrokeGraph &stroke_graph, std::vector<Junction> &varying_candidates,
  FeatureType feature_type, bool to_expand = false,
  bool use_input_junc_type = false);
bool update_high_valence_junction_predictions(
  const StrokeGraph &stroke_graph, std::vector<Junction> &varying_candidates,
  FeatureType feature_type, bool to_dedup = false);

void complete_candidate_graph(
  const StrokeGraph &graph, const std::vector<Junction> &candidates,
  std::vector<Junction> &out_candidates,
  const std::vector<bool> &connectivity = std::vector<bool>());

bool reproject_t_junc(const StrokeGraph &graph, Junction &junc);
bool reproject_t_junc_orig(const StrokeGraph &graph, Junction &junc);

bool replay_deformations(span<const Stroke> strokes,
                         span<const Junction> candidates,
                         StrokeGraph &stroke_graph);
bool replay_connection(const std::vector<Junction> &candidates,
                       const std::vector<bool> &junction_connected,
                       StrokeGraph &stroke_graph,
                       std::vector<Junction> &out_junctions);

void color_endpoint_graph(const StrokeGraph &stroke_graph,
                          const span<Junction> candidates,
                          std::vector<int> &junc_components,
                          const Float min_prob = 0.0);

void trim_overshoots(StrokeGraph &stroke_graph);
void build_plane_graph(span<const Stroke> strokes, StrokeGraph &stroke_graph);
StrokeSnapInfo increment_strokes(StrokeGraph &graph, span<const Stroke> strokes,
                                 size_t start_idx, size_t end_idx);

void vanilla_candidates(const StrokeGraph &stroke_graph,
                        std::vector<Junction> &dangling_predictions,
                        bool train_time, int num_cand = -1);
void vanilla_candidates(const StrokeGraph &end_stroke_graph,
                        const StrokeGraph &stroke_graph,
                        std::vector<Junction> &dangling_predictions,
                        std::vector<Junction> &corner_predictions,
                        bool train_time);
void predicted_corner_candidates(const StrokeGraph &end_stroke_graph,
                                 const StrokeGraph &stroke_graph,
                                 const FeatureType feature_type,
                                 std::vector<Junction> &corner_predictions,
                                 bool to_dedup = false,
                                 bool to_include_prev_connections = false);
void complete_predicted_corner_candidates(
  const StrokeGraph &end_stroke_graph, const StrokeGraph &stroke_graph,
  const FeatureType feature_type, std::vector<Junction> &corner_predictions,
  Float min_probability = -1, bool to_include_prev_connections = false);

bool is_corner_candidate(const StrokeGraph &stroke_graph, const Junction &junc);
void expand_corner_candidates(
  const StrokeGraph &stroke_graph,
  const std::vector<Junction> &corner_predictions,
  std::vector<Junction> &expanded_corner_predictions);

void connect_graph(StrokeGraph &stroke_graph,
                   std::vector<Junction> &candidates);

///////////////////////
//   Visualization.  //
///////////////////////
struct VizProbabilities {
  Float prob_ = -1, future_prob_ = -1, env_dist_ = -1;
  std::string repr = "";
};
void junc2pred(
  const StrokeGraph &graph, const std::vector<Junction> &junctions,
  std::vector<std::pair<ClassifierPrediction, VizProbabilities>> &predictions,
  const std::vector<Junction> &future_junctions = std::vector<Junction>());
Junction pred2junc(const ClassifierPrediction &pred);
void viz_connections(const GraphState &state,
                     const std::vector<Junction> &future_junctions,
                     std::vector<std::string> &connected_junc_strs);
void viz_connected_junctions(
  const StrokeGraph &graph,
  std::vector<std::pair<std::string, Vec2>> &junctions);

} // namespace sketching
