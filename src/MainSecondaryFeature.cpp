#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Serialization.h"
#include "../stroke_strip_src/StrokeCutting.h"
#include "Capture2Input.h"
#include "Logger.h"
#include "Util.h"
#include "feature/FeatureIO.h"
#include "solve/SolveUtil.h"

#include <CLI/CLI.hpp>

#include <fstream>
#include <string>

/**
 * The entrance function of global feature computation given generated training
 * examples.
 */
int main(int argc, char *argv[]) {
  CLI::App app{"Training feature generation."};

  std::string scap_filename;
  std::string example_csv_filename;
  std::string output_csv_filename;
  std::string output_vis_folder = "";
  std::string cache_folder = "";
  std::string fea_file;
  bool has_cut_file = false;
  bool compute_prob = false;

  app.add_option("-i,--input", scap_filename, "The ground truth scap input.");
  app.add_option("-e,--example", example_csv_filename,
                 "The csv file containing the training example.");
  app.add_option("-o,--output", output_csv_filename, "The output csv file.");
  app.add_option("-v,--vis", output_vis_folder,
                 "The output folder for visualizations.");
  app.add_option("-c,--cache", cache_folder,
                 "The output folder for visualizations.");
  app.add_option("-f,--fea", fea_file, "The feature set.");
  app.add_flag("-u,--uncut", has_cut_file,
               "To generate stroke-stroke pairs for measurement.");

  CLI11_PARSE(app, argc, argv);

  // 1. Preprocess
  Input input;
  int width, height;
  double input_thickness;
  // bool to_preprocess = true;
  bool to_preprocess = false;
  read_input(scap_filename, input, width, height, input_thickness,
             to_preprocess);

  // Read cut file (if enabled and exists)
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    matching_pair;
  if (has_cut_file) {
    std::string cut_filename =
      std::regex_replace(scap_filename, std::regex(".scap"), "_cut.csv");
    std::ifstream file;
    file.open(cut_filename, std::ios_base::in);
    if (file.is_open()) {
      std::string line;
      while (std::getline(file, line)) {
        std::vector<std::string> cut_strings = SketchUI::split_str(line, ",");
        assert(cut_strings.size() == 4);

        matching_pair.emplace_back();
        matching_pair.back().first.first = std::stol(cut_strings[0]);
        matching_pair.back().first.second = std::stol(cut_strings[1]);
        matching_pair.back().second.first = std::stol(cut_strings[2]);
        matching_pair.back().second.second = std::stol(cut_strings[3]);
      }
    }
  }

  if (compute_prob) {
    read_feature_definitions(fea_file, solve_context.sec_feature_set);
    solve_context.stage = SolveContext::SolveStage::Splitting;
  }

  // 2. Read examples
  std::vector<std::pair<Cluster, Cluster>> clusters;
  read_samples(input, example_csv_filename, clusters);

  // 3. Compute per cluster features
  std::map<size_t, FeatureVector> cluster_features;
  if (!compute_prob) {
    for (auto &c_cluster : input.clusters) {
      logger().info("Fit cluster {}", c_cluster.first);
      bool to_orient = true;
      bool test_time = false;
      fit_cluster(input.thickness, input.width, input.height, c_cluster.second,
                  to_orient, test_time, &matching_pair);

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
        cluster_features.emplace(c_cluster.first,
                                 FeatureVector(input.width, input.height));
        cluster_features[c_cluster.first].compute_typical(
          c_cluster.second, input.thickness, &matching_pair, reorient,
          output_vis_folder);
      }
    }
  }

  // 4. Compute features
  double thickness = 1;
  std::vector<FeatureVector> feas;
  bool cache_exists = false;
  for (auto &cc : clusters) {
    std::string stroke_str;
    for (auto const &s : cc.first.strokes) {
      if (!stroke_str.empty())
        stroke_str += ' ';
      stroke_str += std::to_string(s.stroke_ind);
    }
    stroke_str += ',';
    for (auto const &s : cc.second.strokes) {
      if (!stroke_str.empty())
        stroke_str += ' ';
      stroke_str += std::to_string(s.stroke_ind);
    }
    std::cout << "Compute on " + stroke_str << std::endl;

    // Look up cached parameterization
    std::vector<size_t> stroke_indices1, stroke_indices2;
    for (auto const &s : cc.first.strokes) {
      stroke_indices1.emplace_back(s.stroke_ind);
    }
    for (auto const &s : cc.second.strokes) {
      stroke_indices2.emplace_back(s.stroke_ind);
    }
    std::string svg_name = (cc.first.strokes.front().cluster_ind ==
                            cc.second.strokes.front().cluster_ind)
                             ? "pos_"
                             : "neg_";
    std::string hash_str = example_hash(stroke_indices1, stroke_indices2);
    std::string json_filename1, json_filename2, json_filename_all;

    Cluster merged_cluster;
    json_filename1 = cache_folder + "/" + svg_name + hash_str + "_c1.json";
    json_filename2 = cache_folder + "/" + svg_name + hash_str + "_c2.json";
    json_filename_all = cache_folder + "/" + svg_name + hash_str + "_comb.json";

    auto file_exists = [](const std::string &file_name) {
      std::ifstream f(file_name.c_str());
      return f.good();
    };

    bool to_write = !cache_folder.empty() && !(file_exists(json_filename1) &&
                                               file_exists(json_filename2) &&
                                               file_exists(json_filename_all));
    if (!cache_folder.empty())
      cache_exists |= file_exists(json_filename1) &&
                      file_exists(json_filename2) &&
                      file_exists(json_filename_all);

    if (!cache_folder.empty() && file_exists(json_filename1) &&
        file_exists(json_filename2) && file_exists(json_filename_all)) {
      read_json(json_filename1, cc.first);
      read_json(json_filename2, cc.second);
      read_json(json_filename_all, merged_cluster);

      // Make it compatible to use old cached file that doesn't save
      // original_input_strokes
      if (merged_cluster.original_input_strokes.empty()) {
        merged_cluster.original_input_strokes.insert(
          merged_cluster.original_input_strokes.end(),
          cc.first.original_input_strokes.begin(),
          cc.first.original_input_strokes.end());
        merged_cluster.original_input_strokes.insert(
          merged_cluster.original_input_strokes.end(),
          cc.second.original_input_strokes.begin(),
          cc.second.original_input_strokes.end());
      }

      // Combined parameterization failed
      if (!cc.first.fit.centerline.empty() &&
          !cc.second.fit.centerline.empty() && merged_cluster.strokes.empty())
        continue;
    } else {
      // Fit the two individual clusters
      bool to_orient = true;
      bool test_time = false;
      fit_cluster(input.thickness, input.width, input.height, cc.first,
                  to_orient, test_time, &matching_pair);
      fit_cluster(input.thickness, input.width, input.height, cc.second,
                  to_orient, test_time, &matching_pair);
    }

    // Either parameterizations fail
    if (cc.first.strokes.empty() || cc.second.strokes.empty()) {
      continue;
    }

    bool reorient = false;
    FeatureVector fea(input.width, input.height);
    if (compute_prob) {
      fea.compute_probability(cc.first, cc.second, thickness, merged_cluster,
                              &matching_pair, reorient, output_vis_folder);
    } else {
      fea.compute_secondary(cc.first, cc.second, thickness, merged_cluster,
                            cluster_features, &matching_pair, reorient,
                            output_vis_folder);
    }

    if (!merged_cluster.strokes.empty())
      feas.emplace_back(fea);

    if (to_write) {
      write_json(json_filename1, cc.first);
      write_json(json_filename2, cc.second);
      write_json(json_filename_all, merged_cluster);
    }
  }

  // 4. Write
  if (compute_prob) {
    write_probability_feature_csv(output_csv_filename, feas);
  } else {
    write_secondary_feature_csv(output_csv_filename, feas);
  }

  return 0;
}
