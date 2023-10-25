#pragma once

#include "bvh.h"
#include "classifier.h"
#include "closest.h"
#include "endpoint.h"
#include "junction.h"
#include "mapping.h"
#include "types.h"

#include <unordered_map>
#include <vector>

namespace sketching {

struct ExtendedJunction;
struct PolylineBVH;
struct Stroke;
struct StrokeCoverage;

/**
 * Small half-edge data structure implementation of a stroke graph with faces.
 */
struct StrokeGraph {
  static constexpr auto invalid = (size_t)-1;

  struct HedgeRecord;

  struct FaceRecord {
    /// Indices of half-edges, *one* per cycle.  In other words, for faces without holes,
    /// `cycles_` will have length 1.
    std::vector<size_t> cycles_;
  };

  struct VertexID {
    std::uint32_t connection_type_ = 0;
    enum Type : decltype(connection_type_) { //
      Initialization = 1,
      Junction = 2
    };

    size_t connection_index_ = 0;

    bool is_valid() const { return connection_type_ > 0; }

    std::string repr() const {
      return ((connection_type_ == Type::Initialization) ? "init_" : "junc_") +
             std::to_string(connection_index_);
    }

    bool operator==(const VertexID& id2) const {
      return (connection_type_ == id2.connection_type_) &&
             (connection_index_ == id2.connection_index_);
    }
  };

  struct VertexRecord {
    /**
     * @param p Location of vertex.
     */
    explicit VertexRecord(const Vec2& p)
      : p_(p) {}

    /// Index of an arbitrary outbound half-edge.
    /// When the vertex is a T-junction, hedge_ is pointed to the horizontal bar of the T
    /// (aka the part attached to the end point of the other stroke).
    size_t hedge_ = invalid;
    /// Location of vertex.
    Vec2 p_;
    /// Additional properties.
    std::uint32_t flags_ = 0;
    enum Flags : decltype(flags_) {
      Inactive = 1,
      Grazing = 1 << 1,
      OriginallyDangling = 1 << 2,
      Overshoot = 1 << 3,
      CenterlineIntersection = 1 << 4,
      // Topologically dangling but overlapping some other stroke.
      // Works around some snapping issues.
      Overlapping = 1 << 5,
    };
    /// Junction type.
    JunctionType::Type junc_type_;

    /// Unique IDs depending on whether the vertex is created during initialization or by
    /// connecting multiple specific junctions.
    std::vector<VertexID> ids_;

    bool is_active() const { return !(flags_ & Inactive); }

    void deactivate() { flags_ |= Inactive; }

    bool is_grazing() const { return flags_ & Grazing; }

    bool is_originally_dangling() const { return flags_ & OriginallyDangling; }
  };

  struct HedgeRecord {
    /// Index of next half-edge on the same face.
    size_t next_ = invalid;
    /// Index of previous half-edge on the same face.
    size_t prev_ = invalid;
    /// Index of the continuity half-edge (i.e. from the same original stroke).  The two
    /// continuity half-edges reference each other and share the same origin vertex.
    size_t continuity_ = invalid;
    /// Index of vertex which this half-edge emanates from.
    size_t origin_ = invalid;
    /// Index of face.
    size_t face_ = invalid;
    /// Additional properties.
    std::uint32_t flags_ = 0;
    enum Flags : decltype(flags_) {
      Inactive = 1,
      Grazing = 1 << 1,
      Bridge = 1 << 2,
    };

    bool is_active() const { return !(flags_ & Inactive); }

    void deactivate() { flags_ |= Inactive; }
  };

  struct HedgeView;

  struct VertexView {
    VertexView()
      : index_(invalid)
      , graph_(nullptr) {}

    VertexView(const StrokeGraph* graph, const size_t index)
      : index_(index)
      , graph_(graph) {}

    HedgeView hedge() const {
      return HedgeView(graph_, graph_->vertices_[index_].hedge_);
    }

    const Vec2 pos() const { return graph_->vertices_[index_].p_; }

    bool is_valid() const { return graph_ && index_ != invalid; }

    bool is_active() const { return graph_->vertices_[index_].is_active(); }

    bool is_originally_dangling() const {
      return graph_->vertices_[index_].is_originally_dangling();
    }

    auto flags() const { return graph_->vertices_[index_].flags_; }

    explicit operator bool() const { return is_valid() && is_active(); }

    bool operator==(VertexView other) const {
      return graph_ == other.graph_ && index_ == other.index_;
    }

    bool operator!=(VertexView other) const { return !(*this == other); }

    size_t valence() const;

    bool is_dangling() const;

    JunctionType::Type junction_type() const {
      return graph_->vertices_[index_].junc_type_;
    }

    const std::vector<VertexID>& vertex_ids() const {
      return graph_->vertices_[index_].ids_;
    }

    size_t index_;

    const StrokeGraph* graph_;
  };

  struct FaceView;

  struct HedgeView {
    HedgeView()
      : index_(invalid)
      , graph_(nullptr) {}

    HedgeView(const StrokeGraph* graph, const size_t index)
      : index_(index)
      , graph_(graph) {}

    HedgeView next() const { return HedgeView(graph_, graph_->hedges_[index_].next_); }

    HedgeView prev() const { return HedgeView(graph_, graph_->hedges_[index_].prev_); }

    VertexView origin() const {
      return VertexView(graph_, graph_->hedges_[index_].origin_);
    }

    VertexView dest() const { return twin().origin(); }

    HedgeView twin() const {
      return HedgeView(graph_, (forward() ? index_ + 1 : index_ - 1));
    }

    const Stroke& stroke() const { return graph_->strokes_[stroke_idx()]; }

    const Stroke& orig_stroke() const { return graph_->orig_strokes_[orig_stroke_idx()]; }

    size_t orig_stroke_idx() const;

    size_t stroke_idx() const { return index_ / 2; }

    FaceView face() const { return FaceView(graph_, face_idx()); }

    size_t face_idx() const { return graph_->hedges_[index_].face_; }

    bool forward() const { return index_ % 2 == 0; }

    bool is_valid() const { return graph_ && index_ != invalid; }

    auto flags() const { return graph_->hedges_[index_].flags_; }

    explicit operator bool() const { return is_valid() && is_active(); }

    bool operator==(HedgeView other) const {
      return graph_ == other.graph_ && index_ == other.index_;
    }

    bool operator!=(HedgeView other) const { return !(*this == other); }

    bool is_active() const { return graph_->hedges_[index_].is_active(); }

    size_t continuity() const { return graph_->hedges_[index_].continuity_; };

    // TODO: For consistency with the other methods, `continuity` should be renamed
    //       `continuity_idx` and `continuity_edge` should be renamed `continuity`.
    HedgeView continuity_edge() const {
      const auto cont = graph_->hedges_[index_].continuity_;
      return (cont != invalid ? HedgeView(graph_, cont) : HedgeView());
    };

    size_t index_;

    const StrokeGraph* graph_;
  };

  struct FaceView {
    FaceView()
      : index_(invalid)
      , graph_(nullptr) {}

    FaceView(const StrokeGraph* graph, const size_t index)
      : index_(index)
      , graph_(graph) {}

    const std::vector<size_t>& cycles() const { return graph_->faces_[index_].cycles_; }

    bool is_valid() const { return graph_ && index_ != invalid; }

    bool operator==(FaceView other) const {
      return graph_ == other.graph_ && index_ == other.index_;
    }

    bool operator!=(FaceView other) const { return !(*this == other); }

    /// Return the number of neighbouring faces.
    size_t n_neighbors() const;

    size_t n_edges() const;

    size_t index_;

    const StrokeGraph* graph_;
  };

  ///////////////////////
  // Low-level access. //
  ///////////////////////

  /// Acceleration structure containing the edge geometries.
  StrokeBVH bvh_;

  /// Stopgap to avoid a massive refactor.
  struct IndirectStrokes {
    explicit IndirectStrokes(StrokeGraph* graph)
      : graph_(graph) {}

    size_t size() const { return graph_->bvh_.strokes().size(); }
    const Stroke& back() const { return graph_->bvh_.strokes().back(); }
    Stroke& operator[](size_t i) { return graph_->bvh_.stroke(i); }
    const Stroke& operator[](size_t i) const { return graph_->bvh_.stroke(i); }
    const Stroke* begin() const { return graph_->bvh_.strokes().begin(); }
    const Stroke* end() const { return graph_->bvh_.strokes().end(); }
    operator span<const Stroke>() const { return graph_->bvh_.strokes(); }

  private:
    StrokeGraph* graph_;
  } strokes_;

  /// Container of vertices.
  std::vector<VertexRecord> vertices_;
  /// Container of half-edges.  They are arranged in the same order as the strokes they
  /// correspond to, such that twins are kept together, with the forward twin before the
  /// backward twin.  Thus the layout looks like [forward hedge of stroke 0, backward
  /// hedge of stroke 0, forward hedge of stroke 1, backward hedge of stroke 1, ...].
  std::vector<HedgeRecord> hedges_;
  /// Container of faces.
  std::vector<FaceRecord> faces_;
  /// Index of boundary face.  Usually the first one.
  size_t boundary_face_ = 0;

  /// Choose snapping method implementations
  enum SnappingType : std::uint32_t { //
    Connection = 0,
    Deformation = 1,
    LinearSolve = 2,
  };
  SnappingType snapping_method_type_ = Connection;

  ////////////////////////
  // High-level access. //
  ////////////////////////

  VertexView vertex(const size_t i) const {
    assert(i < vertices_.size() && "out of bounds");
    return {this, i};
  }
  HedgeView hedge(const size_t i) const {
    assert(i < hedges_.size() && "out of bounds");
    return {this, i};
  }
  FaceView face(const size_t i) const {
    assert(i < faces_.size() && "out of bounds");
    return {this, i};
  }

  /// Create an empty stroke graph.
  StrokeGraph(SnappingType type = Connection)
    : strokes_(this)
    , snapping_method_type_(type) {}

  /// Create a stroke graph, creating junctions conservatively.
  explicit StrokeGraph(span<const Stroke> strokes, SnappingType type = Connection);

  std::string repr() const;

  StrokeGraph clone() const;

  // Move constructors.
  StrokeGraph(StrokeGraph&& other) noexcept;
  StrokeGraph& operator=(StrokeGraph&& other) noexcept;

  // Use `clone` to make an explicit copy.
  StrokeGraph(const StrokeGraph& other) = delete;
  StrokeGraph& operator=(const StrokeGraph& other) = delete;

  /**
   * Build a face polygon.
   */
  std::vector<Vec2> face_positions(size_t face_index) const;

  /////////////////////////////
  // Topological operations. //
  /////////////////////////////

  /**
   * Snap the vertices with indices vi1 and vi2 together.
   *
   * High-valence snaps are supported.
   *
   * The probability `prob` will be recorded if specified.
   */
  bool snap_vertices(size_t vi1, size_t vi2, SnappingType type, Float prob = -1.0,
                     const Vec2* snap_pos = nullptr);

  bool snap_vertices(size_t vi1, size_t vi2) {
    return snap_vertices(vi1, vi2, snapping_method_type_);
  }

  /// Snap the input vertex to the given position of a stroke. Depending on the distance
  /// between the snapping position and existing vertices, this operation may or may not
  /// add a new vertex.
  VertexView snap_endpoint_to_edge(VertexView v, size_t stroke_idx, Float arclen,
                                   SnappingType type);

private:
  VertexView snap_endpoint_to_edge_legacy(VertexView v, size_t stroke_idx, Float arclen,
                                          SnappingType type);

public:
  ///////////////////////////////////
  // Book-keeping data structures. //
  ///////////////////////////////////
  /// Mapping from endpoints on the original stroke to vertex indices.
  std::unordered_map<Endpoint, size_t> endpoint2vertex_;
  std::unordered_map<size_t, Endpoint> vertex2endpoint_;

  // TODO: Rethink these.  These don't belong here.
  /// Copy of original strokes.
  std::vector<Stroke> orig_strokes_;
  /// BVH for original strokes.
  std::unique_ptr<PolylineBVH> orig_bvh_;

  /// Mapping from original stroke index to graph stroke indices.
  /// Original stroke index -> stroke index, mapping (can be multiple if an original
  /// stroke is split). The StrokeMappings are sorted and all domains should have first <
  /// second, while the ranges don't have this assumption. See check_mapping in
  /// stroke_graph_extra.h for more detailed explanations.
  std::vector<std::vector<std::pair<size_t, StrokeMapping>>> orig2strokes_;

  /// Stroke index -> original stroke index (can be multiple if valence-2 vertex is
  /// dissolved so original strokes are merged). The actual mapping info is saved in
  /// orig2strokes_, which can be refered using the original stroke index. See
  /// check_mapping in stroke_graph_extra.h for more detailed explanations.
  std::vector<std::vector<size_t>> strokes2orig_;

  /// Stroke index -> dissolved vertex position on this stroke (in normalized arc length),
  /// VertexID (as index for saved information used by integrated solve).
  struct VidComparator {
    bool operator()(const std::pair<Float, VertexID>& a,
                    const std::pair<Float, VertexID>& b) const {
      return a.first < b.first || ((a.first == b.first) && (a.second.connection_index_ <
                                                            b.second.connection_index_));
    }
  };
  std::vector<std::set<std::pair<Float, VertexID>, VidComparator>> strokes2vid_;

  std::vector<ClassifierPrediction> snap_history_;

  mutable std::unordered_map<std::uint64_t, bool> parallel_endpoints_cache_;

  void snap_endpoints();

  void snap_endpoints_legacy();

private:
  void init(span<const Stroke>, const StrokeCoverage&);
};

/**
 * Return true if the stroke will not be included in the stroke graph.
 */
bool skip_stroke(const Stroke& s);

bool point_in_face(const StrokeGraph::FaceView f, const Vec2 p);

void check_continuity(const StrokeGraph& graph);
void check_consistency(const StrokeGraph& graph);
void check_intersection(const StrokeGraph& graph);

/** Starting from a graph without faces, build the faces from scratch. */
void construct_faces(StrokeGraph& graph);

void insert_into_star(StrokeGraph& graph, size_t hi, size_t vi);

/** Re-uses existing vertices if close enough. */
bool add_vertex(StrokeGraph& graph, size_t stroke_idx, Float arclen,
                StrokeGraph::VertexView& new_vertexview);

/** Guaranteed to create a new vertex if arclen is not equal to 0.0 or stroke length. */
bool add_vertex_precise(StrokeGraph& graph, size_t stroke_idx, Float arclen,
                        StrokeGraph::VertexView& new_vertexview);

} // namespace sketching
