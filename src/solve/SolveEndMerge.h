#pragma once

#include "SolveLazy.h"

/**
 * Merges the sub-strips within the same initial local strip after the splitting
 * during the refinement.
 *
 * @param input The input data.
 * @param clusters The map of clusters.
 * @param vis_folder The path to the folder where visualization files will be
 * saved.
 * @param csv_str The CSV string to append the cluster information to.
 * @param csv_cluster_offset The offset to add to the cluster IDs when writing
 * to the CSV string.
 * @param merged_cluster_ptr A pointer to the merged cluster, if any. Note that
 * merged_cluster_ptr is only valid when clusters are initially from the same
 * cluster. In this case, merged_cluster_ptr points to the parameterization of
 * the initial cluster.
 */
void solve_end_merge(const Input &input, std::map<size_t, Cluster> &clusters,
                     std::string vis_folder, std::string &csv_str,
                     size_t csv_cluster_offset,
                     Cluster *merged_cluster_ptr = nullptr);
