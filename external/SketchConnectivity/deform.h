#pragma once

#include "sketching.h"
#include "stroke_graph.h"

namespace sketching {

void smart_deform(Stroke& stroke, const Vec2& head_pos, const Vec2& tail_pos);

bool non_intersecting_smart_deform(StrokeGraph::VertexView v, Stroke& stroke,
                                   const Vec2& head_pos, const Vec2& tail_pos);

void smooth_deform(Stroke& s, const Vec2& start_endpoint, const Vec2& end_endpoint,
                   Float radius = -1.0);

bool non_intersecting_smooth_deform(StrokeGraph::VertexView v, Stroke& stroke,
                                    const Vec2& start_endpoint, const Vec2& end_endpoint,
                                    Float radius = -1.0);

void straight_line_connect(Stroke& stroke, const Vec2& head_pos, const Vec2& tail_pos);

/**
 * Deform a stroke by applying the correction specified by the oversketched stroke
 * `correction`.
 *
 * @return True if the deformation was applied successfully.  If false, the input stroke
 *         will be left unmodified.
 */
bool oversketch_deform(Stroke& stroke, const Stroke& correction);

} // namespace sketching
