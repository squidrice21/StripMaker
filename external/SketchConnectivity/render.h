#pragma once

#include "stroke_view.h"

#include <vector>

namespace sketching {

struct BoundingBox;
struct Junction;
struct Stroke;

struct PolygonizeOptions {
  using type = std::uint8_t;
  enum Flags : type {
    None = 0,
    StartCapRound = 1,
    EndCapRound = 1 << 1,
    JoinsRound = 1 << 2,
    JoinsBevel = 1 << 3,
    RoundCaps = StartCapRound | EndCapRound,
  };
};

bool get_render_decimation();

/** Enable or disable decimation of output polygons. */
void set_render_decimation(bool enabled);

int get_render_n_cap_vertices();

/** Set number of vertices to use to approximate caps when outputting render polygons. */
void set_render_n_cap_vertices(int n);

/**
 * Convert the outline of a stroke into a collection of polygons.
 *
 * This method works by simulating sweeping a circle along the polyline, with the diameter
 * of the circle varying piecewise-linearly across the segments of the polyline.
 */
void outline_to_polygons(const ConstStrokeView& s, std::vector<CoordMat>& out_vertices,
                         const Float dec_acc);
inline void outline_to_polygons(const ConstStrokeView& s,
                                std::vector<CoordMat>& out_vertices) {
  outline_to_polygons(s, out_vertices, -1);
}

/**
 * Convert the outline of a stroke into a vertex buffer and index buffer, suitable for
 * rendering with e.g. OpenGL.
 *
 * This method works by simulating sweeping a circle along the polyline, with the diameter
 * of the circle varying piecewise-linearly across the segments of the polyline.
 */
void outline_to_triangulation(const ConstStrokeView&, std::vector<Float>& out_coordinates,
                              std::vector<unsigned int>& out_indices);

Float rasterize_scale(span<const Stroke> strokes, int width, int height, BoundingBox& bb);

using Bitmap = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ColorBitmap =
  Eigen::Matrix<std::uint32_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct Col3 {
  static constexpr Col3 from_hex_rgb(uint32_t hex) noexcept {
    auto rgb = Col3();
    rgb.r = Float((hex & 0xff0000) >> 16) / 0xff;
    rgb.g = Float((hex & 0xff00) >> 8) / 0xff;
    rgb.b = Float(hex & 0xff) / 0xff;
    return rgb;
  }

  Float r = 0.0;
  Float g = 0.0;
  Float b = 0.0;
};

span<const Col3> get_color_palette();

struct StrokeGraph;

} // namespace sketching
