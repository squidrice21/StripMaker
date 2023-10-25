#include "Util.h"
#include "feature/FeatureIO.h"
#include "uncutting/Uncutting.h"

#include <CLI/CLI.hpp>

#include <numeric>
#include <string>
#include <type_traits>

int main(int argc, char *argv[]) {
  CLI::App app{"Input preprocessing."};

  std::string scap_filename;
  std::string orig_scap_filename;
  std::string output_scap_filename;
  std::string csv_filename;
  std::string vis_csv_filename;

  bool use_gt_label = false;
  bool use_angle = false;
  bool raw_limited = false;

  app.add_option("-i,--input", scap_filename, "The ground truth scap input.")
    ->required();
  app.add_option("-o,--output", output_scap_filename, "The output scap file.")
    ->required();
  app.add_option("-c,--csv", csv_filename, "The output cut csv file.");
  app.add_option("-v,--vis", vis_csv_filename,
                 "The output visualization csv file for uncut positions.");
  app.add_option("-r,--raw", orig_scap_filename,
                 "The original scap input (with no preprocessing).");
  app.add_flag("-g,--gt", use_gt_label,
               "Consider the GT label in the input file during uncutting.");
  app.add_flag("-l,--limit", raw_limited,
               "Limit the uncut to only cuts based on the input raw file.");
  app.add_flag("-a,--angle", use_angle,
               "Consider the angle threshold at the uncutting points.");

  CLI11_PARSE(app, argc, argv);

  // 1. Preprocess
  int width, height;
  bool to_preprocess = false;
  Capture capture;
  read_input(scap_filename, capture, width, height, to_preprocess);
  Input input;
  double input_thickness;
  read_input(scap_filename, input, width, height, input_thickness,
             to_preprocess);
  std::vector<std::pair<std::pair<Cluster::Stroke, bool>,
                        std::pair<Cluster::Stroke, bool>>>
    matching_pair, cut_pair;

  matching_point(input, matching_pair);
  std::cout << "finish matching point" << std::endl;

  cut_pair = matching_pair;

  if (use_gt_label) {
    check_GT_cluster(input, matching_pair);
    std::cout << "finish check_GT_cluster" << std::endl;
  }

  if (use_angle) {
    check_tangent(input, matching_pair);
    std::cout << "finish check_tangent" << std::endl;
  }

  if (!orig_scap_filename.empty()) {
    Capture capture_orig;
    int width_ref, height_ref;
    read_input(orig_scap_filename, capture_orig, width_ref, height_ref,
               to_preprocess);
    double w, h, w_ref, h_ref;
    glm::dvec2 in_center = capture_center(capture, w, h);
    glm::dvec2 ref_center = capture_center(capture_orig, w_ref, h_ref);

    for (auto &s : capture_orig.sketchedPolylines) {
      for (auto &point : s.points) {
        point.first.x -= ref_center.x;
        point.first.y -= ref_center.y;
        point.first.x *= w / w_ref;
        point.first.y *= w / w_ref;
        point.first.x += in_center.x;
        point.first.y += in_center.y;
      }
    }

    check_original_input(capture, capture_orig, cut_pair);

    if (raw_limited) {
      check_original_input(capture, capture_orig, matching_pair);
      std::cout << "finish check_original_input" << std::endl;
    }
  }

  connect_stroke(capture, matching_pair, cut_pair, vis_csv_filename,
                 csv_filename);

  // 2. Write preprocessed results
  std::ofstream scap_ofs(output_scap_filename);
  std::string out_buffer = capture.to_string();
  scap_ofs << "#" << width << "\t" << height << std::endl;
  scap_ofs.write(out_buffer.c_str(), out_buffer.size());
  scap_ofs.close();

  return 0;
}
