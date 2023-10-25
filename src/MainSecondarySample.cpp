#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/Parameterization.h"
#include "../stroke_strip_src/StrokeOrientation.h"
#include "../stroke_strip_src/Utils.h"
#include "Logger.h"
#include "Util.h"
#include "feature/FeatureIO.h"
#include "measure/Measurement.h"
#include "sample/SampleFilters.h"
#include "sample/SampleGeneration.h"
#include "sample/SampleSecondaryGeneration.h"

#include <CLI/CLI.hpp>

#include <regex>
#include <string>
#include <utility>

/**
 * The entrance function of global sample generation given ground truth scap
 * files.
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
  app.add_option("-o,--output", output_csv_filename, "The output csv file.");
  app.add_option("-c,--cache", cache_folder,
                 "The output folder for visualizations.");
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

  // Try fitting, the result fit curve is empty if the parameterization fails
  for (auto &c_cluster : input.clusters) {
    logger().info("Fit cluster {}", c_cluster.first);
    bool to_orient = true;
    bool test_time = false;
    fit_cluster(input.thickness, input.width, input.height, c_cluster.second,
                to_orient, test_time, &matching_pair);
  }

  // Reindex to have stroke indices going from 0
  std::map<size_t, size_t> inv_reindexing;
  {
    std::map<size_t, size_t> reindexing;

    for (const auto &c : input.clusters) {
      for (const auto &s : c.second.strokes) {
        reindexing[s.stroke_ind] = 0;
      }
    }
    size_t s_idx = 0;
    for (auto &ss : reindexing) {
      ss.second = s_idx++;
      inv_reindexing[ss.second] = ss.first;
    }

    for (auto &c : input.clusters) {
      for (auto &s : c.second.strokes) {
        s.stroke_ind = reindexing[s.stroke_ind];
      }
      for (auto &s : c.second.original_input_strokes) {
        s.stroke_ind = reindexing[s.stroke_ind];
      }
      for (auto &xsec : c.second.xsecs) {
        for (auto &p : xsec.points) {
          p.stroke_idx = reindexing[p.stroke_idx];
          p.stroke_ind = reindexing[p.stroke_ind];
        }
      }
    }
  }

  // 2. Generate positive and negative training samples
  std::vector<Sample> pos_stroke_samples, neg_stroke_samples,
    pos_cluster_samples, neg_cluster_samples;

  logger().info("generate samples");
  generate_secondary_stroke_positive_samples(input, pos_stroke_samples,
                                             &matching_pair);
  generate_secondary_cluster_positive_samples(input, pos_cluster_samples,
                                              &matching_pair);
  generate_stroke_negative_samples(input, neg_stroke_samples, &matching_pair);
  generate_cluster_negative_samples(input, neg_cluster_samples, &matching_pair);

  // Map the index back
  auto inv_indexing = [&inv_reindexing](std::vector<int> &cluster_indices) {
    for (auto &sid : cluster_indices) {
      sid = inv_reindexing[sid];
    }
  };
  auto inv_samples = [&inv_indexing](std::vector<Sample> &samples) {
    for (auto &sample : samples) {
      inv_indexing(sample.first);
      inv_indexing(sample.second);
    }
  };
  inv_samples(pos_stroke_samples);
  inv_samples(neg_stroke_samples);
  inv_samples(pos_cluster_samples);
  inv_samples(neg_cluster_samples);

  // 2.5 Filter samples
  for (auto &c_cluster : input.clusters) {
    for (auto &s : c_cluster.second.strokes) {
      s.stroke_ind = inv_reindexing[s.stroke_ind];
    }
    for (auto &s : c_cluster.second.original_input_strokes) {
      s.stroke_ind = inv_reindexing[s.stroke_ind];
    }
    for (auto &xsec : c_cluster.second.xsecs) {
      sort_xsec(xsec);
      for (auto &p : xsec.points)
        p.stroke_ind = inv_reindexing[p.stroke_ind];
    }
  }

  bool disable_short_overlapping = true;
  std::vector<FilteredSample> filtered_samples;
  logger().info("filter_sample_extra: stroke");

  filter_sample_extra(input, pos_stroke_samples, neg_stroke_samples,
                      filtered_samples, &matching_pair, cache_folder,
                      output_vis_folder, disable_short_overlapping);
  logger().info("filter_sample_extra: cluster");
  filter_sample_extra(input, pos_cluster_samples, neg_cluster_samples,
                      filtered_samples, &matching_pair, cache_folder,
                      output_vis_folder, disable_short_overlapping);

  // 3. Write
  write_samples(output_csv_filename + "_sec_pos_stroke.csv",
                pos_stroke_samples);
  write_samples(output_csv_filename + "_sec_neg_stroke.csv",
                neg_stroke_samples);
  write_samples(output_csv_filename + "_sec_pos_cluster.csv",
                pos_cluster_samples);
  write_samples(output_csv_filename + "_sec_neg_cluster.csv",
                neg_cluster_samples);

  write_filtered_samples(output_csv_filename + "_sec_filtered.csv",
                         filtered_samples);

  return 0;
}
