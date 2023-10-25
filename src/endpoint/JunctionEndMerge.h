#pragma once

#include "../stroke_strip_src/Cluster.h"
#include "Junction.h"

#include <vector>

/**
 * Find junctions that are likely to be connected between input strips/clusters.
 *
 * @param input The input sketch data.
 * @param clusters A map of cluster IDs to Cluster objects.
 * @param junctions A vector of Junction objects.
 */
void junction_end_merge(const Input &input, std::map<size_t, Cluster> &clusters,
                        std::vector<Junction> &junctions);
/**
 * Checks if two fitted curves are continuations of each other at a junction.
 * @param fit1 The first fitted curve.
 * @param fit2 The second fitted curve.
 * @param junc The junction at which the curves meet.
 * @param angle The angle between the two curves at the junction.
 * @return True if the two curves are continuations of each other at the
 * junction, false otherwise.
 */
bool is_continuation(const FittedCurve &fit1, const FittedCurve &fit2,
                     const Junction &junc, double &angle);
/**
 * Resolves continuation of junctions based on the given junction angles.
 *
 * @param junc_angles The vector of junction angles to use for resolving
 * continuation.
 * @param junctions The vector of Junction objects to modify based on the
 * resolution.
 */
void resolve_continuation(const std::vector<double> &junc_angles,
                          std::vector<Junction> &junctions);
/**
 * Merges the ends of the input clusters based on the given junctions.
 *
 * @param input The input data.
 * @param clusters The clusters to merge.
 * @param merged_clusters The resulting merged clusters.
 * @param junctions The junctions to use for merging.
 */
void merge_ends(const Input &input, std::map<size_t, Cluster> &clusters,
                std::map<size_t, Cluster> &merged_clusters,
                const std::vector<Junction> &junctions);
