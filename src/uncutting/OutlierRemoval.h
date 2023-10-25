#pragma once

#include "../endpoint/Junction.h"
#include "../stroke_strip_src/Cluster.h"

#include <map>

/**
 * Removes outliers in the post-processing.
 * @param input The input sketch data.
 * @param in_clusters The input clusters to remove outliers from.
 * @param in_fits The input fitted curves to use for outlier removal.
 * @param clusters The resulting clusters after outlier removal.
 * @param junctions The junctions in the sketch.
 * @param skip_cut Whether to skip cutting the sketch.
 * @return True if outliers were successfully removed, false otherwise.
 */
bool remove_outliers(const Input &input,
                     const std::map<size_t, Cluster> &in_clusters,
                     const std::map<int, FittedCurve> &in_fits,
                     std::map<size_t, Cluster> &clusters,
                     const std::vector<Junction> &junctions, bool skip_cut);
