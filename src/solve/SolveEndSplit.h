#pragma once

#include "SolveLazy.h"

/**
 * The splitting step of the refinement. It creates sub-strips from the initial
 * strip if the probability of merging them is relatively low. This step is
 * followed by a merging step. Note that this splitting uses a given
 * parameterization of the initial strip to avoid computation overhead.
 *
 * @param input The input to solve the end split problem for.
 * @param clusters A map of cluster IDs to Cluster objects to update with the
 * final clusters.
 * @param intermediate_clusters A map of cluster IDs to Cluster objects to
 * update with intermediate clusters.
 * @param intermediate_folder The folder to save intermediate visualizations to.
 * @param vis_folder The folder to save final visualizations to.
 */
void solve_end_split(const Input &input, std::map<size_t, Cluster> &clusters,
                     std::map<size_t, Cluster> &intermediate_clusters,
                     std::string intermediate_folder, std::string vis_folder);
