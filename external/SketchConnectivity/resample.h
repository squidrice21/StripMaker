#pragma once

#include "types.h"

namespace sketching {

struct Stroke;

/**
 * Returns true if and only if at least one vertex was removed.
 */
bool remove_duplicate_vertices(Stroke& s);

/**
 * Create a copy of `s` with duplicate, near-duplicate, eclipsed (vertices completely
 * covered by nearby vertices), and nearly-eclipsed vertices removed.
 *
 * Note that the returned stroke may look different than `s` because removing eclipsed
 * vertices can still change the transitions between vertices.  Only use this function for
 * evenly-sampled strokes if you wish to approximately preserve the envelope shape.
 */
Stroke duplicate_vertices_removed(const Stroke& s);

/**
 * @param s Stroke with n vertices.
 * @param tolerance
 * @param out_include Array of size n, representing whether each vertex should be included
 *                    in a decimation of the given tolerance.
 */
void ramer_douglas_peucker(const Stroke& stroke, double tolerance, bool out_include[]);

void ramer_douglas_peucker(const Float x[], const Float y[], Index n, double tolerance,
                           bool out_include[]);

/**
 * Similar to @ref ramer_douglas_peucker, but tries to keep potential cut points for hooks
 * in the output.
 *
 * @param s Stroke with n vertices.
 * @param tolerance
 * @param out_include Array of size n, representing whether each vertex should be included
 *                    in a decimation of the given tolerance.
 */
void corners_and_hooks(const Stroke& stroke, double tolerance, bool out_include[]);

/**
 * Decimate (simplify) a stroke using the Ramer-Douglas-Peucker algorithm by removing
 * low-contribution points within a certain tolerance.
 *
 * This function does not (yet) take into account the width.
 */
void decimate(Stroke& s, Float tolerance);

void decimate(CoordMat& xy, Float tolerance);

/**
 * Return the input stroke decimated (simplified) using the Ramer-Douglas-Peucker
 * algorithm by removing low-contribution points within a certain tolerance.
 *
 * This function does not (yet) take into account the width.
 */
Stroke decimated(const Stroke& s, Float tolerance);

} // namespace sketching
