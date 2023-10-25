#pragma once

#include "types.h"

namespace sketching {

struct PolylineBVH;
struct Stroke;

double gaussian_factor(double x, double sigma);

/**
 * Compute the average width of a stroke.
 */
Float average_width(const Stroke&);

void average_widths(const PolylineBVH&, span<Float> out_widths);

/**
 * Compute the busyness metric at p.
 */
Float busyness_at(const PolylineBVH&, const Vec2& p, span<const Float> avg_widths,
                  Float busyness_falloff);

Eigen::Matrix<Float, Eigen::Dynamic, 2> busyness_at_endpoints(const PolylineBVH&,
                                                              Float busyness_falloff);

} // namespace sketching
