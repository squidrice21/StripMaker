#pragma once

#include "GraphicsGems.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Fit a Bezier curve to a set of digitized points
 *
 * Args:
 *     d: Array of digitized points
 *     nPts: Number of digitized points
 *     error: User-defined error squared
 *     callback: Called for each Bezier segment.  Arguments are (degree, curve,
 *         start index, end index, segment parameterization).  Segment parameterization
 *         is nullptr for a segment with 2 points.
 */
void FitCurve(Point2* d, int nPts, double error, void(*callback)(int, Point2*, int, int, double*));

void FitSingleCubic(Point2* d, int nPts, double error, void(*callback)(int, Point2*, int, int, double*));

#ifdef __cplusplus
}
#endif
