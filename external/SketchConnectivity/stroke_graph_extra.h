#pragma once

#include "endpoint.h"
#include "stroke_graph.h"
#include "types.h"

namespace sketching {

/** The returned endpoint indexes into graph.strokes_ (unlike graph.vertex2endpoint_). */
inline Endpoint hedge_to_endpoint(StrokeGraph::HedgeView he) {
  return {he.stroke_idx(), he.forward()};
}

/** The returned endpoint indexes into graph.strokes_ (unlike graph.vertex2endpoint_). */
Endpoint vertex_to_endpoint(StrokeGraph::VertexView v);

/** `endp` indexes into graph.strokes_ (unlike graph.endpoint2vertex_). */
inline StrokeGraph::HedgeView endpoint_to_hedge(const StrokeGraph& graph,
                                                const Endpoint& edge_endp) {
  const auto hi = 2 * edge_endp.stroke_idx() + (edge_endp.is_head() ? 0 : 1);
  return graph.hedge(hi);
}

/** `endp` indexes into graph.strokes_ (unlike graph.endpoint2vertex_). */
inline StrokeGraph::VertexView endpoint_to_vertex(const StrokeGraph& graph,
                                                  const Endpoint& edge_endp) {
  return endpoint_to_hedge(graph, edge_endp).origin();
}

[[nodiscard]] bool is_bridge_vertex(const StrokeGraph& stroke_graph,
                                    const StrokeGraph::VertexView v);

/**
 * Return true iff snapping vertex vi1 to vi2 would result in a spurious self-connection
 * for some edge.
 */
[[nodiscard]] bool is_snap_spurious(const StrokeGraph& graph, const size_t vi1,
                                    const size_t vi2);

/**
 * Return true iff vertex vi1 can be snapped to vertex vi2 without forming bad topology.
 *
 * Doesn't do line-of-sight checks (because they are expensive), so make sure that the two
 * vertices have a line of sight to each other.
 */
[[nodiscard]] bool is_snap_valid(const StrokeGraph& graph, size_t vi1, size_t vi2);

/**
 * Return true iff vertex vi1 can be snapped to the specified stroke interior without
 * forming bad topology.
 *
 * Doesn't do line-of-sight checks (because they are expensive), so make sure that the two
 * points have a line of sight to each other.
 */
[[nodiscard]] bool is_snap_valid(const StrokeGraph& graph, size_t vi1,
                                 StrokeTime stroke_time);

/**
 * Return true iff snapping vertex vi1 to the specified stroke interior would result in a
 * sliver face
 *
 * `stroke_time` is specified as (stroke index, normalized arc length).
 */
[[nodiscard]] bool would_form_sliver_vert_interior(const StrokeGraph& graph, size_t vi,
                                                   StrokeTime stroke_time);

/**
 * Return true iff there is a continuous half-edge emanating from this vertex.
 */
bool has_continuous_edge(StrokeGraph::VertexView v);

/**
 * Requires `out_component_ids.size() == graph.strokes_.size()`.
 *
 * Returns the number of connected_components.
 */
int connected_components_edge(const StrokeGraph& graph, span<int> out_component_ids);

/**
 * Remove tiny hooks from dangling edges.
 */
void dehook_dangling_edges(StrokeGraph& graph);

/**
 * Remove tiny hooks from a dangling vertex.
 */
void dehook_dangling_vertex(StrokeGraph& graph, StrokeGraph::VertexView v);

/**
 * Append the positions on a half-edge cycle to out_vertices.
 */
void cycle_positions(StrokeGraph::HedgeView he, std::vector<Vec2>& out_vertices);

/**
 * Delete a stroke and its associated half-edges in constant time.  Indices will be
 * rearranged.
 *
 * The deleted stroke's associated half-edges must *not* have any references to them when
 * this function is called.  This is not an edge collapse function.
 *
 * @param graph
 * @param si Index of the stroke to delete.
 */
void delete_stroke(StrokeGraph& graph, size_t si);

/**
 * Delete a vertex in constant time.  Indices will be rearranged.
 *
 * The deleted vertex must *not* have any references to it when this function is called.
 *
 * @param graph
 * @param vi Index of the vertex to delete.
 */
void delete_vertex(StrokeGraph& graph, size_t vi);

/**
 * Dissolve a vertex by merging the edges on both sides.  For this operation to be valid,
 * the vertex must have valence 2.
 *
 * Returns the index of the merged stroke if the operation was successful.  If the
 * return value is `StrokeGraph::invalid`, the graph is left unmodified.
 */
size_t dissolve_vertex(StrokeGraph& graph, size_t vi);

/**
 * Return the stroke that would result from dissolving vertex vi, without actually
 * dissolving the vertex.  The arc length position on the dissolved stroke
 * corresponding to the join point (the vertex) will be output to `join_arclen`.
 *
 * Returns an empty stroke if the dissolve operation is invalid.
 */
Stroke dissolved_stroke(const StrokeGraph& graph, size_t vi, Float& join_arclen);

/**
 * Return the stroke that would result from joining the half-edges with indices hi1 and
 * hi2.  The destination vertex of hedge hi1 must be coincident with the origin vertex of
 * hi2.  For a more general chaining operation, see `consolidate_with_chaining`.
 *
 * The arc length position on the resulting stroke corresponding to the join point will be
 * output to `join_arclen`.
 *
 * The direction of the result is unspecified.
 */
Stroke concatenate_edges(const StrokeGraph& graph, size_t hi1, size_t hi2,
                         Float& join_arclen);

/**
 * Move a vertex to a new location. move_to_v is the reference to the vertex vi is moved
 * to.
 */
bool move_vertex(StrokeGraph& graph, size_t vi, Vec2 new_pos,
                 StrokeGraph::SnappingType type, StrokeGraph::VertexView move_to_v);

/**
 * Connect vertices vi1 and vi2 with an edge.
 *
 * The vertices must have a line of sight to each other, otherwise the planarity of the
 * graph will be violated.
 *
 * Returns the stroke index of the newly added bridge edge.
 */
size_t bridge_vertices(StrokeGraph& graph, size_t vi1, size_t vi2);

/**
 * Add a new mapping (orig stroke domain -> stroke range). Domain and range are normalized
 * arc length. While domain has the two interval boundaries: first < second, range can
 * have any normalized arc length values. This function can have several effects:
 *
 * 1- When no mapping info is saved for the given original stroke. It simply saves the
 * input mapping.
 * 2- When the input orig stroke domain partially overlaps the orig stroke domain of an
 * existing mapping. It split the existing one into two (note we don't allow the
 * overlapping that splits the domain into three).
 * 3- When the input orig stroke domain fully overlaps the orig stroke domain of an
 * existing mapping. It modifies the saved mapping info.
 *
 * Note case 2 should only be used by split_stroke_mapping. Other parts should only call
 * this function for case 1 and 3 (with si_remain and s_remain_range as default).
 *
 * @param graph Stroke graph.
 * @param orig_si Index of the original stroke.
 * @param orig_domain The input orig stroke domain. After calling this function, there
 * should be a mapping record with this exact orig stroke domain.
 * @param si Index of the stroke mapped to.
 * @param s_range The range of the stroke mapped to.
 * @param si_remain When orig_domain partially overlaps an existing mapping domain. The
 * stroke index of the remaining mapping after split. -1 means keeping it as the mapping
 * before the change.
 * @param s_remain_range When orig_domain partially overlaps an existing mapping domain.
 * The range of the remaining mapping after split. (-1, -1) means keeping it as the
 * mapping before the change.
 */
void add_stroke_mapping(
  StrokeGraph& graph, size_t orig_si, std::pair<Float, Float> orig_domain, size_t si,
  std::pair<Float, Float> s_range, int si_remain = -1,
  std::pair<Float, Float> s_remain_range = std::pair<Float, Float>(-1, -1));

/**
 * Higher-level mapping modification function. It handles all mapping changes when a
 * stroke (indexed as si) is split into a stroke si1 and a stroke si2. si1 is the part <
 * the split position and si2 is the part > the split position.
 *
 * @param graph Stroke graph.
 * @param si Index of the stroke that gets split.
 * @param split The split position indicated as a normalized arc length on si.
 * @param si1 Index of the substroke si1 after the split.
 * @param si2 Index of the substroke si2 after the split. si2 = -1 makes si2 = si.
 */
void split_stroke_mapping(StrokeGraph& graph, size_t si, Float split, size_t si1,
                          int si2 = -1);

inline void split_stroke_mapping(StrokeGraph& graph, size_t si, Float split, size_t si1,
                                 size_t si2) {
  split_stroke_mapping(graph, si, split, si1, (int)si2);
}

/**
 * Higher-level mapping modification function. It handles all mapping changes when two
 * strokes si1, si2, are merged into a stroke si, due to valence-2 vertex dissolving.
 *
 * Note it also keeps record of dissolved vertex IDs (save_vid) on a stroke. This is used
 * for the objective value computation in the integrated solve.
 *
 * @param graph Stroke graph.
 * @param si1 Index of stroke si1.
 * @param s1_to_s The range of the mapping s1 -> s in normalized arc length unit. The
 * domain is implied as [0, 1]. The first and second can have any lesser and greater
 * relationship.
 * @param si2 Index of stroke si2.
 * @param s2_to_s The range of the mapping s2 -> s in normalized arc length unit. The
 * domain is implied as [0, 1]. The first and second can have any lesser and greater
 * relationship.
 * @param si Index of the stroke that gets split.
 * @param save_vid The ID of the dissolved vertex.
 */
void merge_stroke_mapping(StrokeGraph& graph, size_t si1, std::pair<Float, Float> s1_to_s,
                          size_t si2, std::pair<Float, Float> s2_to_s, size_t si,
                          const std::vector<StrokeGraph::VertexID>& save_vids);

/**
 * Similar to trim_original_mapping, but doesn't trim the original stroke.
 * The result is similar to if the first/last section of the stroke got shortened via
 * deformation.
 */
void clamp_edge_mapping(StrokeGraph& graph, size_t edge_idx,
                        std::pair<Float, Float> new_domain);

/**
 * Check if mapping info is consistent. This is done both on original strokes and strokes:
 * 1- Check if the domain/range intervals are without overlapping and gaps.
 * 2- Check if the interval union is [0, 1].
 * Note we ignore empty mapping info since the StrokeGraph default init may cause it when
 * very short strokes exist.
 *
 * @param graph Stroke graph.
 */
bool check_mapping(const StrokeGraph& graph);

/**
 * Check if the stroke interiors are intersecting.
 *
 * @param graph Stroke graph.
 */
void check_interior_intersection(const StrokeGraph& graph);

/**
 * Convert a index-normalized-arc-length position between strokes and original strokes.
 *
 * Returns true on success.  Otherwise `end` is left unmodified.
 */
[[nodiscard]] bool convert_strokes2orig(const StrokeGraph& graph, StrokeTime& end);
[[nodiscard]] bool convert_orig2strokes(const StrokeGraph& graph, StrokeTime& end);

/**
 * Like convert_strokes2orig but returns a new StrokeTime.
 *
 * The program will abort if the mapping fails.
 */
[[nodiscard]] StrokeTime as_orig_position(const StrokeGraph& graph,
                                          const StrokeTime& edge_pos);
/**
 * Like convert_orig2strokes but returns a new StrokeTime.
 *
 * The program will abort if the mapping fails.
 */
[[nodiscard]] StrokeTime as_edge_position(const StrokeGraph& graph,
                                          const StrokeTime& orig_pos);

/** Returns true on success.  Otherwise the input junc is left unmodified. */
[[nodiscard]] bool stroke_to_original_stroke_indexing(const StrokeGraph& stroke_graph,
                                                      Junction& junc);
[[nodiscard]] bool original_stroke_to_stroke_indexing(const StrokeGraph& stroke_graph,
                                                      Junction& junc);

StrokeGraph::VertexView orig2endpoint(const StrokeGraph& graph, const StrokeTime& end);

/** Return value uses original indexing. */
StrokeTime reproject_cand2_on_original(const StrokeGraph& graph,
                                       const StrokeTime& orig_pos1,
                                       const StrokeTime& cand2);

/** Disconnect an edge from the stroke graph and then deactivate it. */
void drop_edge(StrokeGraph& graph, size_t si);

/**
 * Drop a vertex and all its neighbouring edges from the graph.
 *
 * The vertices and edges will merely be deactivated.  To remove their records from the
 * graph entirely, call `delete_vertex` and `delete_stroke` afterwards.
 */
void drop_vertex_and_neighbors(StrokeGraph& graph, size_t vi,
                               span<Junction> junctions_to_update);

/**
 * Return the apparent radius of this vertex, based on all emanating edges.
 */
Float vertex_radius(StrokeGraph::VertexView v);

/**
 * Print a nice string representation of the orig2strokes_ mapping.
 */
std::string orig2strokes_mapping_repr(const StrokeGraph& graph);

/**
 * Print a nice string representation of the strokes2orig_ mapping.
 */
std::string strokes2orig_mapping_repr(const StrokeGraph& graph);

bool similar_gaps_antisliver(const Stroke& stroke1, Float s_start, Float s_end,
                             const Stroke& stroke2, Float t_start, Float t_end);

/**
 * Returns true iff the vertex is a corner.
 */
bool is_corner(const StrokeGraph& graph, size_t v_idx, bool include_t = false);

/**
 * Get all corresponding original stroke end positions for a corner vertex.
 */
void get_corner_original_positions(const StrokeGraph& graph, size_t v_idx,
                                   std::vector<Vec2>& positions, bool include_t = false);

/**
 * Returns true iff the chain is cyclic.
 */
bool get_forward_chain_original_indices(
  const StrokeGraph& graph, size_t he_idx,
  std::vector<std::pair<size_t, bool>>& orig_indices, bool to_dissolve = false);

/**
 * Direction of the output stroke matches the half-edge direction.
 * The chain extends in both directions, unlike `get_forward_chain_original_indices`.
 */
Stroke get_chain(StrokeGraph::HedgeView he_index,
                 std::vector<std::pair<size_t, bool>>& orig_indices);

/**
 * Format of `out_orig_indices` is list of list of (original stroke index, direction).
 * Outer list has one element per stroke in the chained drawing.  Inner list has one
 * element per chain segment.
 */
std::vector<Stroke>
chained_drawing(const StrokeGraph& graph,
                std::vector<std::vector<std::pair<size_t, bool>>>& out_mapping);

/**
 * Query connectivity information about endpoints.
 */
struct Connections {
  using Connection = std::pair<StrokeTime, StrokeTime>;

  Connections()
    : max_orig_stroke_idx_(INT_MAX) {}

  explicit Connections(const StrokeGraph& graph);

  explicit Connections(const EnvelopeBVH& strokes);

  /**
   * Get all the connections associated with the input endpoint.
   */
  [[nodiscard]] span<const Connection> associated_connections(StrokeTime endp) const;

  /** Get all the connections. */
  [[nodiscard]] span<const Connection> connections() const { return connections_; }

  /**
   * Return the largest original stroke index with a connection to this endpoint.
   *
   * Unlike in `associated_connections`, `endp` is included as a connection for thus
   * function, making it straightforward to determine if you should delay a decision (by
   * comparing the return value to the largest original stroke index added so far).
   */
  [[nodiscard]] int largest_connected_stroke_idx(StrokeTime endp) const;

  /**
   * Return true iff the specified connection is present.
   *
   * Points must be specified in terms of original stroke indices, using normalized arc
   * lengths.  `a` must be an endpoint.  `b` can be an endpoint or stroke interior.
   */
  [[nodiscard]] bool contains(StrokeTime a, StrokeTime b) const;

  [[nodiscard]] bool empty() const { return connections_.empty(); }

  void clear() {
    connections_.clear();
    max_orig_stroke_idx_ = INT_MAX;
  }

  /**
   * For each pair (p_i, q_i) in (a, b), add (p_i, q_i) (and (q_i, p_i) for end-end
   * connections) to `connections_`.
   */
  void insert(span<const StrokeTime> a, span<const StrokeTime> b);

  /**
   * Create a subset of the current stored connections so only the ones involving the
   * input stroke (indexed as si) are in this returned set.
   */
  Connections subset_involving(size_t si) const;

private:
  /// Sorted and expressed in terms of original stroke indices.
  std::vector<Connection> connections_;

  /// Used for validation.
  int max_orig_stroke_idx_;
};

} // namespace sketching
