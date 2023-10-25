#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Serialization.h"
#include "../stroke_strip_src/StrokeCutting.h"
#include "Capture2Input.h"
#include "Util.h"
#include "feature/FeatureIO.h"

#include <CLI/CLI.hpp>

#include <fstream>
#include <string>

/**
 * The entrance function of feature computation given generated training
 * examples.
 */
int main(int argc, char *argv[]) {
  CLI::App app{"Training feature generation."};

  std::string scap_filename;
  std::string example_csv_filename;
  std::string output_csv_filename;
  std::string output_vis_folder = "";
  std::string cache_folder = "";
  bool has_cut_file = false;

  app.add_option("-i,--input", scap_filename, "The ground truth scap input.");
  app.add_option("-e,--example", example_csv_filename,
                 "The csv file containing the training example.");
  app.add_option("-o,--output", output_csv_filename, "The output csv file.");
  app.add_option("-v,--vis", output_vis_folder,
                 "The output folder for visualizations.");
  app.add_option("-c,--cache", cache_folder,
                 "The output folder for parameterization cache.");
  app.add_flag("-u,--uncut", has_cut_file,
               "TO generate stroke-stroke pairs for measurement.");

  CLI11_PARSE(app, argc, argv);

  // 1. Preprocess
  Input input;
  int width, height;
  double input_thickness;
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

  // 2. Read examples
  std::vector<std::pair<Cluster, Cluster>> clusters;
  read_samples(input, example_csv_filename, clusters);

  auto begin = std::chrono::high_resolution_clock::now();

  // 3. Compute features
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
    fea.compute(cc.first, cc.second, thickness, merged_cluster, false,
                &matching_pair, reorient, output_vis_folder);
    if (!merged_cluster.strokes.empty())
      feas.emplace_back(fea);

    if (to_write) {
      write_json(json_filename1, cc.first);
      write_json(json_filename2, cc.second);
      write_json(json_filename_all, merged_cluster);
    }
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::cout << scap_filename << " feature time: " << std::fixed
            << std::setprecision(6)
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count() /
                 1000.0
            << " s" << std::endl;

  // 4. Write
  write_feature_csv(output_csv_filename, feas);

  return 0;
}
