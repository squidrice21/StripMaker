#include "SolveEndMerge.h"

#include "../Logger.h"

void solve_end_merge(const Input &input, std::map<size_t, Cluster> &clusters,
                     std::string vis_folder, std::string &csv_str,
                     size_t csv_cluster_offset, Cluster *merged_cluster_ptr) {
  // Offset the cluster indices for visualization
  std::map<size_t, Cluster> offset_clusters;
  for (auto &cc : clusters) {
    offset_clusters[cc.first + csv_cluster_offset] = cc.second;
  }
  clusters = offset_clusters;

  int start_merge_cid = clusters.begin()->first;
  while (start_merge_cid <= clusters.rbegin()->first) {
    // Recursively merge this current cluster to other clusters
    std::vector<size_t> merge_history;
    bool in_place = true;
    compare_cluster_cluster_recursive(input, start_merge_cid, clusters,
                                      merge_history, vis_folder, csv_str,
                                      in_place, merged_cluster_ptr);

    // Move to the next merging seed cluster
    int new_start_merge_cid = -1;
    for (auto const &cc : clusters) {
      if (cc.first > start_merge_cid) {
        new_start_merge_cid = cc.first;
        break;
      }
    }

    if (new_start_merge_cid < 0)
      break;

    start_merge_cid = new_start_merge_cid;
  }
}
