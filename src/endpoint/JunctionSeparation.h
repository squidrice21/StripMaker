#pragma once

#include "../stroke_strip_src/Cluster.h"

#include <map>
#include <vector>

/**
 * @brief Checks if substrips within a single local strip are separable and
 * returns the separable strips.
 *
 * @param input The input sketch.
 * @param clusters The initial strips.
 * @param sep_clusters The separable strips.
 * @param merged_cluster The merged strip (containing the parameterization of
 * the local strip).
 * @param vis_folder The visualization folder.
 * @return true if strips are separable, false otherwise.
 */
bool are_separable_clusters(const Input &input,
                            std::map<size_t, Cluster> clusters,
                            std::map<size_t, Cluster> &sep_clusters,
                            Cluster merged_cluster,
                            const std::string &vis_folder);

/**
 * @brief Separates the substrips within a single local strip and returns the
 * separated strips.
 *
 * @param input The input sketch.
 * @param init_clusters The initial strips.
 * @param sec_clusters The separated strips.
 * @param sep_cid The separated cluster IDs.
 * @param vis_folder The visualization folder.
 */
void separable_clusters(const Input &input,
                        std::map<size_t, Cluster> &init_clusters,
                        std::map<size_t, Cluster> &sec_clusters,
                        std::vector<std::pair<size_t, size_t>> &sep_cid,
                        const std::string &vis_folder = "");
