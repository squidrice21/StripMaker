#include "../stroke_strip_src/Serialization.h"
#include "Logger.h"
#include "Util.h"
#include "endpoint/Junction.h"
#include "feature/FeatureIO.h"
#include "solve/SolveEndMerge.h"
#include "solve/SolveEndSplit.h"
#include "solve/SolveUtil.h"

#include <CLI/CLI.hpp>
#include <iomanip>
#include <string>
#include <vector>

std::string output_folder;
std::string intermediate_folder = "";

int width, height;
double input_thickness;
Input input;
Capture capture;

void try_solve_end_split(std::map<size_t, Cluster> clusters,
                         std::map<size_t, Cluster> &out_clusters,
                         std::map<size_t, Cluster> &split_clusters,
                         std::vector<double> thresholds,
                         std::string base_filename, std::string prefix,
                         std::string vis_folder) {
  std::vector<double> tmp_thresholds = solve_context.split_merge_cutoffs;
  std::map<std::string, Cluster> tmp_merged_cluster_cache =
    solve_context.merged_cluster_cache;
  solve_context.split_merge_cutoffs = thresholds;
  solve_context.separation_sid.clear();
  split_clusters.clear();

  logger().info("try_solve_end_split: {}", prefix);
  solve_end_split(input, clusters, split_clusters, intermediate_folder,
                  vis_folder);
  solve_context.split_merge_cutoffs = tmp_thresholds;
  out_clusters = clusters;

  dump_clusters(capture, out_clusters, input, width, height, input_thickness,
                output_folder, base_filename + "_" + prefix,
                prefix + "_prediction_fit");
  dump_clusters(capture, split_clusters, input, width, height, input_thickness,
                output_folder, base_filename + "_split",
                "split_prediction_fit");
}

void try_solve_end_merge(std::map<size_t, Cluster> clusters,
                         std::vector<double> thresholds,
                         std::string base_filename, std::string prefix,
                         std::string param_prefix, std::string vis_folder) {
  std::vector<double> tmp_thresholds{solve_context.single_final_cutoff,
                                     solve_context.cutoff,
                                     solve_context.single_single_final_cutoff};
  std::map<std::string, Cluster> tmp_merged_cluster_cache =
    solve_context.merged_cluster_cache;
  solve_context.single_final_cutoff = thresholds[0];
  solve_context.cutoff = thresholds[1];
  solve_context.single_single_final_cutoff = thresholds[2];

  Cluster *merged_cluster_ptr = nullptr;
  std::string merge_csv_str = "@";
  solve_end_merge(input, clusters, vis_folder, merge_csv_str, 0,
                  merged_cluster_ptr);
  dump_clusters(capture, clusters, input, width, height, input_thickness,
                output_folder, base_filename + "_" + prefix,
                prefix + "_prediction_fit");

  solve_context.single_final_cutoff = tmp_thresholds[0];
  solve_context.cutoff = tmp_thresholds[1];
  solve_context.single_single_final_cutoff = tmp_thresholds[2];
  solve_context.merged_cluster_cache = tmp_merged_cluster_cache;

  // Write the result parameterizations
  std::string cache_folder = intermediate_folder + "/cache/";
  for (auto &cc : clusters) {
    std::string hash_str = get_hash_str(cc.second);
    std::string json_filename_all =
      cache_folder + "/" + hash_str + "_" + param_prefix + ".json";

    if (!file_exists(json_filename_all)) {
      write_json(json_filename_all, cc.second);
    }
  }
}

/**
 * The entrance function for solving. It runs the initial local temporary step
 * and the refinement step (turned on by "-s,--split").
 */
int main(int argc, char *argv[]) {
  CLI::App app{"Solve."};

  std::string scap_filename;

  std::string fea_file;
  std::string sec_fea_file;
  std::string vis_folder = "";
  std::string cut_file = "";

  // Use lazy comparison by default
  bool lazy_cluster_comparison = true;
  bool use_gt_label = false;
  bool use_cached_label = false;

  // Extra stages
  bool to_end_merge = false;
  bool to_end_split = false;
  bool to_junc_sep = false;

  bool to_count_comparisons = false;

  app.add_option("-i,--input", scap_filename, "The input scap.")->required();
  app.add_option("-o,--output", output_folder, "The output folder.")
    ->required();
  app.add_option("--interm", intermediate_folder,
                 "The folder for intermediate results.");
  app.add_option("-f,--fea", fea_file, "The feature set.");
  app.add_option("--sec", sec_fea_file, "The secondary feature set.");
  app.add_option("-u,--uncut", cut_file,
                 "The records of cut points between strokes.");
  app.add_flag("-c,--cached", use_cached_label,
               "Use labels saved as the first stage result.");
  app.add_flag(
    "-s,--split", to_end_split,
    "Apply in-cluster spliting after the initial local temporary step.");
  app.add_flag(
    "-m,--merge", to_end_merge,
    "Apply cluster-cluster merging after the initial local temporary step.");
  app.add_flag("-j,--junc", to_junc_sep, "Apply junction separation.");

  CLI11_PARSE(app, argc, argv);

  solve_context.count_skip = to_count_comparisons;

  if (use_gt_label) {
    logger().warn(
      "The solve is running in the GT mode. It uses the GT labels as if they "
      "are predictions. This should only be used for debugging!");
  }

  if (intermediate_folder.empty() && !vis_folder.empty())
    intermediate_folder = vis_folder;

  // 1. Read. This also sorts the input strokes based on indices.
  if (!fea_file.empty())
    read_feature_definitions(fea_file, solve_context.feature_set);
  if (!sec_fea_file.empty())
    read_feature_definitions(sec_fea_file, solve_context.sec_feature_set);

  std::string base_filename = scap_filename;
  if (scap_filename.find_last_of("/\\") != std::string::npos) {
    base_filename = scap_filename.substr(scap_filename.find_last_of("/\\") + 1);
  }
  // Remove the "_out" from the name
  if (use_cached_label) {
    base_filename = base_filename.erase(base_filename.find_last_of("_"));
  }
  if (base_filename.find(".") != std::string::npos) {
    base_filename = base_filename.erase(base_filename.find("."));
  }

  bool to_preprocess = false;
  read_input(scap_filename, input, width, height, input_thickness,
             to_preprocess);
  read_input(scap_filename, capture, width, height, to_preprocess);

  // Sum up the sample count
  {
    double sum_points = 0;
    size_t s_count = 0;
    size_t max_s = 0;
    for (const auto &c : input.clusters) {
      for (const auto &s : c.second.original_input_strokes) {
        sum_points += s.points.size();
        max_s = std::max(max_s, s.points.size());
        s_count++;
      }
    }
    stroke_sampling_size =
      sum_points / stroke_sampling_total_ref * stroke_sampling_size_ref;
    stroke_sampling_size =
      std::min(input_thickness * 10,
               std::max(stroke_sampling_size_ref, stroke_sampling_size));

    if (stroke_sampling_size != stroke_sampling_size_ref) {
      std::cout << base_filename << " reparam: " << stroke_sampling_size
                << std::endl;
      input.clusters.clear();
      read_input(scap_filename, input, width, height, input_thickness,
                 to_preprocess);
      capture.sketchedPolylines.clear();
      read_input(scap_filename, capture, width, height, to_preprocess);
    }
  }

  // Read cut file (if enabled and exists)
  std::string cut_filename = cut_file;
  if (!cut_file.empty()) {
    std::ifstream file;
    file.open(cut_filename, std::ios_base::in);
    if (file.is_open()) {
      std::string line;
      while (std::getline(file, line)) {
        std::vector<std::string> cut_strings = SketchUI::split_str(line, ",");
        assert(cut_strings.size() == 4);

        solve_context.matching_pair.emplace_back();
        solve_context.matching_pair.back().first.first =
          std::stol(cut_strings[0]);
        solve_context.matching_pair.back().first.second =
          std::stol(cut_strings[1]);
        solve_context.matching_pair.back().second.first =
          std::stol(cut_strings[2]);
        solve_context.matching_pair.back().second.second =
          std::stol(cut_strings[3]);
      }
    }
  }

  // Compute per cluster features
  bool has_global = false;
  for (auto const &fea : solve_context.sec_feature_set) {
    if (fea.find("global_") != std::string::npos) {
      has_global = true;
      break;
    }
  }
  logger().info("has_global: {}", has_global);
  if (use_cached_label && has_global) {
    for (auto &c_cluster : input.clusters) {
      logger().info("Fit cluster {}", c_cluster.first);
      bool to_orient = true;
      bool test_time = false;
      fit_cluster(input.thickness, input.width, input.height, c_cluster.second,
                  to_orient, test_time, &solve_context.matching_pair);

      // Handling single stroke cluster
      if (c_cluster.second.strokes.size() == 1) {
        Input single_input;
        single_input.thickness = input.thickness;
        single_input.width = input.width;
        single_input.height = input.height;
        single_input.clusters[0] = c_cluster.second;
        single_input.clusters[0].xsecs.clear();
        auto obj_term_values = fea_param.parameterize(&single_input);
        c_cluster.second = single_input.clusters[0];
        c_cluster.second.obj_term_values = obj_term_values;
      }

      if (!c_cluster.second.strokes.empty() &&
          c_cluster.second.obj_term_values.size() >= 2) {
        bool reorient = false;
        solve_context.cluster_features.emplace(
          c_cluster.first, FeatureVector(input.width, input.height));
        solve_context.cluster_features[c_cluster.first].compute_typical(
          c_cluster.second, input.thickness, &solve_context.matching_pair,
          reorient, "");
      }
    }
  }

  // 2. Solve
  if (!vis_folder.empty())
    make_folder(vis_folder, false);
  if (!intermediate_folder.empty())
    make_cache_folder(intermediate_folder, false);

  auto begin = std::chrono::high_resolution_clock::now();

  if (use_gt_label) {
    for (const auto &c : input.clusters) {
      for (const auto &s : c.second.original_input_strokes) {
        solve_context.gt_sid2cid[s.stroke_ind] = c.first;
      }
    }
  }

  solve_context.stage = SolveContext::SolveStage::Initial;
  solve_context.width = input.width;
  solve_context.height = input.height;
  solve_context.norm_thickness = 1;
  solve_context.input_thickness = input_thickness;
  solve_context.to_junc_sep = to_junc_sep;
  std::map<size_t, Cluster> clusters, merge_clusters, split_clusters;

  if (!use_cached_label) {
    solve(input, clusters, lazy_cluster_comparison, vis_folder);

    {
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << base_filename << " stage1: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                         begin)
                       .count() /
                     1000.0
                << " s" << std::endl;
    }

    input.clusters.clear();
    for (auto &cc : clusters) {
      input.clusters[cc.first] = cc.second;
      input.clusters[cc.first].fit.cluster_idx = cc.first;
    }
    dump_clusters(capture, clusters, input, width, height, input_thickness,
                  output_folder, base_filename, "final_prediction_fit");

    // Write to cache folder
    // Parameterize saved clusters
    std::string cache_folder = vis_folder + "/cache/";
    if (!intermediate_folder.empty())
      cache_folder = intermediate_folder + "/cache/";
    auto file_exists = [](const std::string &file_name) {
      std::ifstream f(file_name.c_str());
      return f.good();
    };

    for (auto &cc : clusters) {
      std::string hash_str = get_hash_str(cc.second);
      std::string json_filename_all =
        cache_folder + "/" + hash_str + "_init.json";

      if (!file_exists(json_filename_all)) {
        write_json(json_filename_all, cc.second);
      }
    }

    // Precompute for global features
    if (has_global) {
      for (auto &c_cluster : clusters) {
        // Handling single stroke cluster
        if (c_cluster.second.strokes.size() == 1) {
          Input single_input;
          single_input.thickness = input.thickness;
          single_input.width = input.width;
          single_input.height = input.height;
          single_input.clusters[0] = c_cluster.second;
          single_input.clusters[0].xsecs.clear();
          auto obj_term_values = fea_param.parameterize(&single_input);
          c_cluster.second = single_input.clusters[0];
          c_cluster.second.obj_term_values = obj_term_values;
        }

        if (!c_cluster.second.strokes.empty() &&
            c_cluster.second.obj_term_values.size() >= 2) {
          bool reorient = false;
          solve_context.cluster_features.emplace(
            c_cluster.first, FeatureVector(input.width, input.height));
          solve_context.cluster_features[c_cluster.first].compute_typical(
            c_cluster.second, input.thickness, &solve_context.matching_pair,
            reorient, "");
        }
      }
    }
  } else {
    // Parameterize saved clusters
    std::string cache_folder = vis_folder + "/cache/";
    if (!intermediate_folder.empty())
      cache_folder = intermediate_folder + "/cache/";
    auto file_exists = [](const std::string &file_name) {
      std::ifstream f(file_name.c_str());
      return f.good();
    };

    for (auto const &cc : input.clusters) {
      clusters[cc.first] = cc.second;
      clusters[cc.first].fit.cluster_idx = cc.first;
    }
    for (auto &cc : clusters) {
      std::string hash_str = get_hash_str(cc.second);
      std::string json_filename_all =
        cache_folder + "/" + hash_str + "_init.json";

      if (file_exists(json_filename_all)) {
        read_json(json_filename_all, cc.second);
      } else {
        bool to_orient = true;
        bool test_time = true;
        fit_cluster(solve_context.norm_thickness, input.width, input.height,
                    cc.second, to_orient, test_time,
                    &solve_context.matching_pair);

        write_json(json_filename_all, cc.second);
      }

      // Write to the solve cache
      solve_context.merged_cluster_cache[hash_str] = cc.second;
    }
  }

  auto begin2 = std::chrono::high_resolution_clock::now();
  if (to_end_split) {
    {
      std::vector<std::vector<double>> threshold_sets{
        std::vector<double>{0.3, 0.5, 0.3}};
      for (size_t i = 0; i < threshold_sets.size(); ++i) {
        std::string prefix = (i == 0) ? "smerge" : "smerge" + std::to_string(i);
        // prefix = "smerge" + std::to_string(i + 8);
        std::map<size_t, Cluster> out_clusters;
        try_solve_end_split(clusters, out_clusters, split_clusters,
                            threshold_sets[i], base_filename, prefix,
                            vis_folder);
        if (i == 0)
          merge_clusters = out_clusters;
      }
      clusters = merge_clusters;
    }

    if (!to_end_merge) {
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << base_filename << " stage2: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                         begin2)
                       .count() /
                     1000.0
                << " s" << std::endl;
      std::cout << base_filename << " time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                         begin)
                       .count() /
                     1000.0
                << " s" << std::endl;

      dump_clusters(capture, split_clusters, input, width, height,
                    input_thickness, output_folder, base_filename + "_split",
                    "split_prediction_fit");
    }

    if (to_end_merge) {
      solve_context.stage = SolveContext::SolveStage::Merging;

      // Register single stroke clusters
      {
        solve_context.final_single_sid.clear();
        for (auto &c_cluster : clusters) {
          if (c_cluster.second.strokes.size() == 1)
            solve_context.final_single_sid.emplace(
              c_cluster.second.strokes.front().stroke_ind);
        }
      }

      // Recompute for global features
      if (has_global) {
        solve_context.cluster_features.clear();
        for (auto &c_cluster : clusters) {
          if (c_cluster.second.fit.widths.empty()) {
            Input input_fit;
            input_fit.clusters[0] = c_cluster.second;
            context.widths = true;
            input_fit.clusters[0].fit =
              fea_fitting.fit(&input_fit).begin()->second;
            c_cluster.second.fit = input_fit.clusters[0].fit;
            context.widths = false;
          }

          // Write to the solve cache
          std::string hash_str = get_hash_str(c_cluster.second);
          solve_context.merged_cluster_cache[hash_str] = c_cluster.second;

          // Handling single stroke cluster
          if (c_cluster.second.strokes.size() == 1) {
            Input single_input;
            single_input.thickness = input.thickness;
            single_input.width = input.width;
            single_input.height = input.height;
            single_input.clusters[0] = c_cluster.second;
            single_input.clusters[0].xsecs.clear();
            auto obj_term_values = fea_param.parameterize(&single_input);
            c_cluster.second = single_input.clusters[0];
            c_cluster.second.obj_term_values = obj_term_values;
          }

          if (!c_cluster.second.strokes.empty() &&
              c_cluster.second.obj_term_values.size() >= 2) {
            bool reorient = false;
            solve_context.cluster_features.emplace(
              c_cluster.first, FeatureVector(input.width, input.height));
            solve_context.cluster_features[c_cluster.first].compute_typical(
              c_cluster.second, input.thickness, &solve_context.matching_pair,
              reorient, "");
          }
        }
      }

      dump_clusters(capture, split_clusters, input, width, height,
                    input_thickness, output_folder, base_filename + "_split",
                    "split_prediction_fit");

      std::vector<std::vector<double>> threshold_sets{
        std::vector<double>{0.5, 0.55, 0.7}};
      for (size_t i = 0; i < threshold_sets.size(); ++i) {
        std::string prefix = (i == 0) ? "merge" : "merge" + std::to_string(i);
        std::string end_prefix = (i == 0) ? "end" : "end" + std::to_string(i);
        try_solve_end_merge(clusters, threshold_sets[i], base_filename, prefix,
                            end_prefix, vis_folder);
      }
      auto end = std::chrono::high_resolution_clock::now();
      std::cout << base_filename << " stage2: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                         begin2)
                       .count() /
                     1000.0
                << " s" << std::endl;
      std::cout << base_filename << " time: "
                << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                         begin)
                       .count() /
                     1000.0
                << " s" << std::endl;
    }
  } else {
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << base_filename << " stage2: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       begin2)
                     .count() /
                   1000.0
              << " s" << std::endl;
    std::cout << base_filename << " time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       begin)
                     .count() /
                   1000.0
              << " s" << std::endl;
  }

  {
    std::cout << base_filename << " orientt: " << std::fixed
              << std::setprecision(6) << context.orientation_timer << " s"
              << std::endl;
    std::cout << base_filename << " paramt: " << std::fixed
              << std::setprecision(6) << context.parameterization_timer << " s"
              << std::endl;
    std::cout << base_filename << " fitt: " << std::fixed
              << std::setprecision(6) << context.fitting_timer << " s"
              << std::endl;
    std::cout << base_filename << " featuret: " << std::fixed
              << std::setprecision(6) << solve_context.feature_timer << " s"
              << std::endl;
    std::cout << base_filename << " classifiert: " << std::fixed
              << std::setprecision(6) << solve_context.classifier_timer << " s"
              << std::endl;
    std::cout << base_filename << " skipped: " << solve_context.skip_count
              << "," << solve_context.single_count << ","
              << solve_context.ours_count << ";" << solve_context.call_count
              << "," << solve_context.single_call_count << ","
              << solve_context.ours_call_count << std::endl;
  }

  return 0;
}
