#pragma once

#include "incremental.h"
#include "stroke_graph.h"
#include "types.h"

namespace sketching {
struct Stroke;
struct StrokeGraph;
struct Junction;
struct GraphState;

void group_connected_junctions(const StrokeGraph &stroke_graph,
                               span<const Junction> candidates,
                               const std::vector<bool> &junction_connected,
                               std::vector<size_t> &connected_indices,
                               std::vector<size_t> &component_labels);

void propose_candidates_incremental(
  span<const Stroke> strokes, const StrokeGraph &final_graph,
  const std::vector<Junction> &final_predictions, const GraphState &in_state,
  size_t cur_stroke_step, std::vector<Junction> &new_candidates,
  bool train_time);
void propose_candidates_final(span<const Stroke> strokes,
                              const StrokeGraph &final_graph,
                              const std::vector<Junction> &final_predictions,
                              const GraphState &in_state,
                              size_t cur_stroke_step,
                              std::vector<Junction> &new_candidates);
void complete_graph_candidates(
  const StrokeGraph &stroke_graph, span<const Stroke> strokes,
  const StrokeGraph &final_graph, size_t cur_stroke_step,
  const std::vector<Junction> &to_complete_candidates,
  std::vector<Junction> &new_candidates);
void update_candidate_record(const StrokeGraph &stroke_graph,
                             const std::vector<Junction> &in_candidates,
                             const std::vector<bool> &in_connectivity,
                             const std::vector<Junction> &new_candidates,
                             std::vector<Junction> &varying_candidates,
                             std::vector<Junction> &out_candidates,
                             std::vector<bool> &out_connectivity);

void expand_trial_candidate(const StrokeGraph &graph,
                            span<const Junction> candidates, const size_t idx,
                            std::set<size_t> &binding);

bool is_junction_on_cycle(const StrokeGraph::VertexView v,
                          const Junction &candidate);
bool is_assignment_valid(span<const Junction> candidates,
                         const std::vector<bool> &junction_connected,
                         const StrokeGraph::VertexView v);

Float junction_distance_init(const StrokeGraph &stroke_graph,
                             const Junction &junc);
Float junction_distance(const StrokeGraph &stroke_graph, const Junction &junc);

std::string face_id(const StrokeGraph::FaceView f);
Float face_max_stroke_width(const StrokeGraph::FaceView f);

/**
 * Return a new graph with the requested modifications (candidates[i] where
 * junction_connected[i] is true).
 *
 * Returns nullptr if the modification could not be made.
 *
 * TODO: Rename to modified_graph to reflect that this function is not in-place.
 */
std::unique_ptr<StrokeGraph>
modify_graph(const StrokeGraph &stroke_graph, span<Junction> candidates,
             const std::vector<bool> &junction_connected,
             std::vector<std::pair<size_t, size_t>> &adj_faces,
             std::vector<Float> &junc_distances,
             std::vector<StrokeGraph::VertexID> &junc_vertices,
             StrokeGraph::SnappingType snapping_type);

// TODO: Rename to modified_graph to reflect that this function is not in-place.
std::unique_ptr<StrokeGraph>
modify_graph(const StrokeGraph &stroke_graph, span<Junction> candidates,
             const std::vector<bool> &junction_connected,
             std::vector<std::pair<size_t, size_t>> &adj_faces,
             std::vector<Float> &junc_distances,
             std::vector<StrokeGraph::VertexID> &junc_vertices);

void junction_distance_sort(const StrokeGraph &stroke_graph,
                            std::vector<Junction> &candidates);
void junction_drawing_order_sort(std::vector<Junction> &candidates);
void junction_probability_sort(std::vector<Junction> &candidates);

Float face_maximum_inscribing_circle_radius(const StrokeGraph &stroke_graph,
                                            size_t face_idx,
                                            Eigen::Vector2d &center);
Float face_maximum_inscribing_circle_radius_clipping(
  const StrokeGraph &stroke_graph, size_t face_idx, Eigen::Vector2d &center);
Float face_perimeter(const StrokeGraph &stroke_graph, size_t face_idx,
                     bool include_interior_strokes = false);
Float face_area(const StrokeGraph &stroke_graph, size_t face_idx);
inline Float
face_maximum_inscribing_circle_radius(const StrokeGraph &stroke_graph,
                                      size_t face_idx) {
  Eigen::Vector2d center;
  return face_maximum_inscribing_circle_radius(stroke_graph, face_idx, center);
}
Float face_stroke_width_min(const StrokeGraph &graph, size_t fi);

} // namespace sketching
