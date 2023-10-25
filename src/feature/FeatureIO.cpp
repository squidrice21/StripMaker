#include "FeatureIO.h"

#include "../Capture2Input.h"
#include "../Logger.h"
#include "../Util.h"
#include "../measure/Measurement.h"
#include "../sample/SampleFilters.h"
#include "../stroke_strip_src/Polyline2D.h"
#include "Filters.h"

#include <fstream>
#include <regex>
#include <set>
#include <string>
#include <vector>

static void write_single_sample(std::ofstream &file,
                                const std::vector<int> &cluster_indices) {
  for (size_t i = 0; i < cluster_indices.size(); ++i) {
    if (i != 0)
      file << " ";
    file << cluster_indices[i];
  }
}

void read_samples(const Input &input, const std::string &example_filename,
                  std::vector<std::pair<Cluster, Cluster>> &clusters) {
  std::ifstream file;
  file.open(example_filename, std::ios_base::in);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open example file {}", example_filename);
    return;
  }

  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : input.clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      sid2s[s.stroke_ind] = s;
    }
  }

  std::string line;
  while (std::getline(file, line)) {
    std::vector<std::string> cluster_strings = SketchUI::split_str(line, ",");
    assert(cluster_strings.size() == 2);

    auto read_cluster = [&sid2s](const std::string &c_str, Cluster &cluster) {
      std::set<size_t> gt_cid;
      std::vector<std::string> s_strings = SketchUI::split_str(c_str, " ");
      std::vector<size_t> s_indices;
      for (auto sstr : s_strings) {
        size_t s_idx = std::stol(sstr);
        s_indices.emplace_back(s_idx);
      }

      for (auto s_idx : s_indices) {
        assert(sid2s.count(s_idx));
        cluster.strokes.emplace_back(sid2s[s_idx]);
        gt_cid.emplace(sid2s[s_idx].cluster_ind);
      }
      cluster.original_input_strokes = cluster.strokes;
      assert(gt_cid.size() == 1);
    };

    clusters.emplace_back();
    read_cluster(cluster_strings[0], clusters.back().first);
    read_cluster(cluster_strings[1], clusters.back().second);
  }
}

void read_samples(const Input &input, const std::string &example_filename,
                  std::vector<Cluster> &clusters) {
  std::ifstream file;
  file.open(example_filename, std::ios_base::in);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open example file {}", example_filename);
    return;
  }

  std::map<size_t, Cluster::Stroke> sid2s;
  for (const auto &c : input.clusters) {
    for (const auto &s : c.second.original_input_strokes) {
      sid2s[s.stroke_ind] = s;
    }
  }

  std::string line;
  while (std::getline(file, line)) {
    // Turn binary example to single
    line = std::regex_replace(line, std::regex(","), " ");
    std::vector<std::string> cluster_strings = SketchUI::split_str(line, ",");
    assert(cluster_strings.size() == 1);

    auto read_cluster = [&sid2s](const std::string &c_str, Cluster &cluster) {
      std::set<size_t> gt_cid;
      std::vector<std::string> s_strings = SketchUI::split_str(c_str, " ");
      std::vector<size_t> s_indices;
      for (auto sstr : s_strings) {
        size_t s_idx = std::stol(sstr);
        s_indices.emplace_back(s_idx);
      }

      for (auto s_idx : s_indices) {
        assert(sid2s.count(s_idx));
        cluster.original_input_strokes.emplace_back(sid2s[s_idx]);
        cluster.strokes.emplace_back(sid2s[s_idx]);
        gt_cid.emplace(sid2s[s_idx].cluster_ind);
      }

      std::sort(cluster.strokes.begin(), cluster.strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });
      std::sort(cluster.original_input_strokes.begin(),
                cluster.original_input_strokes.end(),
                [](const Cluster::Stroke &a, const Cluster::Stroke &b) {
                  return a.stroke_ind < b.stroke_ind;
                });
    };

    clusters.emplace_back();
    read_cluster(cluster_strings[0], clusters.back());
  }
}

void read_feature_definitions(const std::string &fea_filename,
                              std::vector<std::string> &fea_csv_names) {
  fea_csv_names.clear();

  std::ifstream file;
  file.open(fea_filename, std::ios_base::in);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open feature definition file {}",
                   fea_filename);
    return;
  }

  std::string line;
  while (std::getline(file, line)) {
    // Lines would go like "- [feature name]"
    if (line.find("-") == std::string::npos)
      continue;
    fea_csv_names.emplace_back(line.substr(2));
  }
}

void write_filtered_samples(
  const std::string &filtered_example_filename,
  const std::vector<FilteredSample> &filtered_samples) {
  std::ofstream file(filtered_example_filename);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open training filtered sample file {}",
                   filtered_example_filename);
    return;
  }

  for (const auto &sample_row : filtered_samples) {
    write_single_sample(file, sample_row.sample.first);
    file << ",";
    write_single_sample(file, sample_row.sample.second);

    // Additional printouts for debugging
    file << ",";
    file << int(sample_row.is_pos);
    file << ",";
    file << int(sample_row.cause);

    file << "\n";
  }
}

void write_samples(const std::string &example_filename,
                   const std::vector<Sample> &pos_samples) {
  std::ofstream file(example_filename);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open training sample file {}",
                   example_filename);
    return;
  }

  for (const auto &sample_row : pos_samples) {
    write_single_sample(file, sample_row.first);
    file << ",";
    write_single_sample(file, sample_row.second);
    file << "\n";
  }
}

void write_feature_csv(const std::string &feature_filename,
                       const std::vector<FeatureVector> &fea) {
  std::ofstream file(feature_filename);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open feature file {}", feature_filename);
    return;
  }

  // Header
  file << "cluster_idx1,stroke_indices1,cluster_idx2,stroke_indices2,ans";
  const auto feature_csv = get_stroke_cluster_features();
  for (size_t i = 0; i < feature_csv.size(); ++i) {
    file << ",";
    if (feature_csv[i]->num_binned_ == 1)
      file << feature_csv[i]->csv();
    else {
      for (size_t j = 0; j < feature_csv[i]->num_binned_; ++j) {
        file << feature_csv[i]->csv() << "_b" << j;
        if (j + 1 < feature_csv[i]->num_binned_)
          file << ",";
      }
    }
  }
  std::vector<std::string> obj_term_values{
    "velocity",
    "alignment",
  };
  for (size_t i = 0; i < obj_term_values.size(); ++i) {
    file << ",";
    file << obj_term_values[i];
  }
  std::vector<std::string> vis_files{
    "input", "input_cluster", "isoline", "fit", "parameterization",
  };
  for (size_t i = 0; i < vis_files.size(); ++i) {
    file << ",";
    file << vis_files[i];
  }
  file << "\n";

  // Actual feature values
  for (const auto &f : fea) {
    // Filters: Filter out pairs without side-by-side section
    if (filter_feature(f))
      continue;

    assert(f.cluster_info_vec_.size() == 2);
    file << f.cluster_info_vec_[0].cluster_idx << ",";
    for (size_t i = 0; i < f.cluster_info_vec_[0].stroke_indices.size(); ++i) {
      if (i != 0)
        file << " ";
      file << f.cluster_info_vec_[0].stroke_indices[i];
    }
    file << ",";
    file << f.cluster_info_vec_[1].cluster_idx << ",";
    for (size_t i = 0; i < f.cluster_info_vec_[1].stroke_indices.size(); ++i) {
      if (i != 0)
        file << " ";
      file << f.cluster_info_vec_[1].stroke_indices[i];
    }
    file << ","
         << (f.cluster_info_vec_[0].cluster_idx ==
             f.cluster_info_vec_[1].cluster_idx);
    for (size_t i = 0; i < f.features_.size(); ++i) {
      file << ",";
      file << f.features_[i];
    }

    // Parameterization objective
    for (size_t i = 0; i < obj_term_values.size(); ++i) {
      file << ",";
      if (f.obj_term_values_.count(obj_term_values[i]))
        file << f.obj_term_values_.at(obj_term_values[i]);
    }

    // Visualization
    for (size_t i = 0; i < vis_files.size(); ++i) {
      file << ",";
      if (f.vis_files_.count(vis_files[i]))
        file << f.vis_files_.at(vis_files[i]);
    }

    file << "\n";
  }
}

void write_probability_feature_csv(const std::string &feature_filename,
                                   const std::vector<FeatureVector> &fea) {
  std::ofstream file(feature_filename);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open feature file {}", feature_filename);
    return;
  }

  // Header
  file << "cluster_idx1,stroke_indices1,cluster_idx2,stroke_indices2,ans";
  const auto feature_csv = get_stroke_cluster_probability_features();
  for (size_t i = 0; i < feature_csv.size(); ++i) {
    file << ",";
    if (feature_csv[i]->num_binned_ == 1)
      file << feature_csv[i]->csv();
    else {
      for (size_t j = 0; j < feature_csv[i]->num_binned_; ++j) {
        file << feature_csv[i]->csv() << "_b" << j;
        if (j + 1 < feature_csv[i]->num_binned_)
          file << ",";
      }
    }
  }
  std::vector<std::string> obj_term_values{
    "velocity",
    "alignment",
  };
  for (size_t i = 0; i < obj_term_values.size(); ++i) {
    file << ",";
    file << obj_term_values[i];
  }
  std::vector<std::string> vis_files{
    "input", "input_cluster", "isoline", "fit", "parameterization",
  };
  for (size_t i = 0; i < vis_files.size(); ++i) {
    file << ",";
    file << vis_files[i];
  }
  file << "\n";

  // Actual feature values
  for (const auto &f : fea) {
    // Filters: Filter out pairs without side-by-side section
    if (filter_feature(f))
      continue;

    assert(f.cluster_info_vec_.size() == 2);
    file << f.cluster_info_vec_[0].cluster_idx << ",";
    for (size_t i = 0; i < f.cluster_info_vec_[0].stroke_indices.size(); ++i) {
      if (i != 0)
        file << " ";
      file << f.cluster_info_vec_[0].stroke_indices[i];
    }
    file << ",";
    file << f.cluster_info_vec_[1].cluster_idx << ",";
    for (size_t i = 0; i < f.cluster_info_vec_[1].stroke_indices.size(); ++i) {
      if (i != 0)
        file << " ";
      file << f.cluster_info_vec_[1].stroke_indices[i];
    }
    file << ","
         << (f.cluster_info_vec_[0].cluster_idx ==
             f.cluster_info_vec_[1].cluster_idx);
    for (size_t i = 0; i < f.features_.size(); ++i) {
      file << ",";
      file << f.features_[i];
    }

    // Parameterization objective
    for (size_t i = 0; i < obj_term_values.size(); ++i) {
      file << ",";
      if (f.obj_term_values_.count(obj_term_values[i]))
        file << f.obj_term_values_.at(obj_term_values[i]);
    }

    // Visualization
    for (size_t i = 0; i < vis_files.size(); ++i) {
      file << ",";
      if (f.vis_files_.count(vis_files[i]))
        file << f.vis_files_.at(vis_files[i]);
    }

    file << "\n";
  }
}

void write_secondary_feature_csv(const std::string &feature_filename,
                                 const std::vector<FeatureVector> &fea) {
  std::ofstream file(feature_filename);
  if (!file.is_open()) {
    logger().error("ERROR: Unable to open feature file {}", feature_filename);
    return;
  }

  // Header
  file << "cluster_idx1,stroke_indices1,cluster_idx2,stroke_indices2,ans";
  const auto feature_csv = get_stroke_cluster_features();
  size_t first_count = 0;
  for (size_t i = 0; i < feature_csv.size(); ++i) {
    file << ",";
    if (feature_csv[i]->num_binned_ == 1) {
      file << feature_csv[i]->csv();
      first_count++;
    } else {
      for (size_t j = 0; j < feature_csv[i]->num_binned_; ++j) {
        file << feature_csv[i]->csv() << "_b" << j;
        if (j + 1 < feature_csv[i]->num_binned_)
          file << ",";
        first_count++;
      }
    }
  }
  std::vector<std::string> obj_term_values{
    "velocity",
    "alignment",
  };
  for (size_t i = 0; i < obj_term_values.size(); ++i) {
    file << ",";
    file << obj_term_values[i];
  }

  const auto global_feature_calculators = get_secondary_cluster_features();
  file << ",";
  for (size_t i = 0; i < global_feature_calculators.size(); ++i) {
    auto fea_csv = global_feature_calculators[i]->csvs();
    for (size_t j = 0; j < fea_csv.size(); ++j) {
      const auto &csv = fea_csv[j];
      file << csv;
      if (i + 1 < global_feature_calculators.size() || j + 1 < fea_csv.size())
        file << ",";
    }
  }
  std::vector<std::string> vis_files{
    "input", "input_cluster", "isoline", "fit", "parameterization",
  };
  for (size_t i = 0; i < vis_files.size(); ++i) {
    file << ",";
    file << vis_files[i];
  }
  file << "\n";

  // Actual feature values
  for (const auto &f : fea) {
    // Filters: Filter out pairs without side-by-side section
    if (filter_feature(f))
      continue;

    assert(f.cluster_info_vec_.size() == 2);
    file << f.cluster_info_vec_[0].cluster_idx << ",";
    for (size_t i = 0; i < f.cluster_info_vec_[0].stroke_indices.size(); ++i) {
      if (i != 0)
        file << " ";
      file << f.cluster_info_vec_[0].stroke_indices[i];
    }
    file << ",";
    file << f.cluster_info_vec_[1].cluster_idx << ",";
    for (size_t i = 0; i < f.cluster_info_vec_[1].stroke_indices.size(); ++i) {
      if (i != 0)
        file << " ";
      file << f.cluster_info_vec_[1].stroke_indices[i];
    }
    file << ","
         << (f.cluster_info_vec_[0].cluster_idx ==
             f.cluster_info_vec_[1].cluster_idx);
    for (size_t i = 0; i < f.features_.size(); ++i) {
      file << ",";
      file << f.features_[i];

      // Parameterization objective
      if (i == first_count) {
        for (size_t j = 0; j < obj_term_values.size(); ++j) {
          file << ",";
          if (f.obj_term_values_.count(obj_term_values[j]))
            file << f.obj_term_values_.at(obj_term_values[j]);
        }
      }
    }

    // Visualization
    for (size_t i = 0; i < vis_files.size(); ++i) {
      file << ",";
      if (f.vis_files_.count(vis_files[i]))
        file << f.vis_files_.at(vis_files[i]);
    }

    file << "\n";
  }
}

void save_fitting(const std::string &scap_filename,
                  const std::string &output_filename, const bool fit_width,
                  const std::string &vis_folder,
                  const std::string &cut_filename, const bool fit_single,
                  const bool to_spiral, const bool disable_cut_orientation) {
  context.tighter_fit = true;
  context.to_spiral = to_spiral;

  Input input_in;
  int width, height;
  double input_thickness;
  bool to_preprocess = false;
  read_input(scap_filename, input_in, width, height, input_thickness,
             to_preprocess);
  Input input;
  if (fit_single) {
    input.orig_center = input_in.orig_center;
    input.thickness = input_in.thickness;
    for (auto &c_cluster : input_in.clusters) {
      for (auto &s : c_cluster.second.strokes) {
        // s.spiral = true;
        s.spiral = false;
        set_spiral_cut_angle(s);
        size_t new_c = input.clusters.size();
        input.clusters[new_c].strokes.emplace_back(s);
        input.clusters[new_c].original_input_strokes.emplace_back(s);
      }
    }
  } else {
    read_input(scap_filename, input, width, height, input_thickness,
               to_preprocess);
  }

  // Get the cut points and use them to merge strokes (if all strokes are
  // consistently oriented)
  std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    matching_pair;
  if (!cut_filename.empty()) {
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

  if (!vis_folder.empty()) {
    context.debug_viz = true;
  }

  bool to_orient = true;
  bool test_time = true;
  if (disable_cut_orientation)
    context.ignore_endpoint = true;
  reorient_param_input(input, to_orient, test_time, matching_pair);
  if (disable_cut_orientation)
    context.ignore_endpoint = false;

  if (!vis_folder.empty()) {
    // Input cluster
    std::string final_output_name = vis_folder + "/input_cluster.svg";
    {
      std::ofstream input_svg(final_output_name);
      fea_param.save_non_isoline_svg(input_svg, input);
    }

    // Orientation
    final_output_name = vis_folder + "/orientations.svg";
    {
      std::ofstream orientations_svg(final_output_name);
      fea_orientation.orientation_debug(orientations_svg, input);
    }

    // Isoline
    final_output_name = vis_folder + "/isolines.svg";
    {
      std::ofstream isolines_svg(final_output_name);
      fea_param.isolines_svg(isolines_svg, -1, input);
    }

    // Orthogonal isoline
    size_t max_c = 0;
    for (const auto &c : input.clusters)
      for (size_t i = 0; i < c.second.orthogonal_xsecs_list.size(); ++i)
        max_c = std::max(max_c, i);
    for (size_t i = 0; i <= max_c; ++i) {
      if (i == 0)
        final_output_name = vis_folder + "/relaxed_isolines.svg";
      else
        final_output_name =
          vis_folder + "/itr" + std::to_string(i) + "_isolines.svg";
      {
        std::ofstream isolines_svg(final_output_name);
        fea_param.orthogonal_isolines_svg(isolines_svg, -1, input, i);
      }
    }

    // Parameterization
    final_output_name = vis_folder + "/parameterization.svg";
    {
      std::ofstream param_svg(final_output_name);
      input.param_svg(param_svg);
    }

    context.debug_viz = false;
  }

  bool in_width = context.widths;
  context.widths = false;
  auto fit_map_final = fea_fitting.fit(&input);

  if (fit_width) {
    for (auto &c_fit : fit_map_final) {
      input.clusters[c_fit.first].fit = c_fit.second;
    }
    context.widths = true;
    fit_map_final = fea_fitting.fit(&input);
  }

  for (auto &c_fit : fit_map_final)
    transform_fit_curve(c_fit.second, input.orig_center, input_thickness);

  std::ofstream fit_svg(output_filename);
  input.thickness = input_thickness;
  fea_fitting.fit_svg(fit_svg, input, fit_map_final, width, height);
  input.thickness = 1;

  if (!fit_width) {
    Capture capture;
    for (const auto &c : input.clusters) {
      capture.sketchedPolylines.emplace_back();
      capture.sketchedPolylines.back().stroke_ind = c.first;
      capture.sketchedPolylines.back().group_ind = c.first;
      for (const auto &p : fit_map_final.at(c.first).centerline)
        capture.sketchedPolylines.back().points.emplace_back(
          SketchUI::Point2D(p.x, p.y), 0);
    }
    capture.thickness = input_thickness;

    std::string scap_out =
      std::regex_replace(output_filename, std::regex("svg"), "scap");
    std::ofstream scap_ofs(scap_out);
    std::string out_buffer = capture.to_string();
    scap_ofs << "#" << width << "\t" << height << std::endl;
    scap_ofs.write(out_buffer.c_str(), out_buffer.size());
    scap_ofs.close();
  }

  context.widths = in_width;
}
