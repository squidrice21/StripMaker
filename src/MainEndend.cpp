#include "../stroke_strip_src/Serialization.h"
#include "Logger.h"
#include "Util.h"
#include "endpoint/Junction.h"
#include "endpoint/JunctionEndMerge.h"
#include "endpoint/JunctionMatching.h"
#include "endpoint/JunctionSerialization.h"
#include "feature/FeatureIO.h"
#include "solve/SolveEndMerge.h"
#include "solve/SolveEndSplit.h"
#include "solve/SolveLazy.h"
#include "solve/SolveUtil.h"
#include "uncutting/OutlierRemoval.h"

#include <CLI/CLI.hpp>
#include <set>
#include <string>

Input input;
Capture capture;
int width, height;
double input_thickness;

/**
 * The entrance function of post-processing (fitting continuation junctions and
 * removing outlier short strokes within larger strips).
 */
int main(int argc, char *argv[]) {
  CLI::App app{"StripMaker post-processing."};

  std::string scap_filename;
  std::string cache_filename;
  std::string split_filename;
  std::string fea_file;

  std::string vis_folder = "";
  std::string cut_file = "";
  std::string end_str = "";

  bool to_spiral = false;

  app.add_option("-i,--input", scap_filename, "The input scap.")->required();
  app.add_option("-v,--vis", vis_folder, "The visualization folder.");
  app.add_option("-c,--cache", cache_filename, "The cache parent folder.");
  app.add_option("-e,--end", end_str, "The parameterization cache ID string.");
  app.add_option("-u,--uncut", cut_file,
                 "TO generate stroke-stroke pairs for measurement.");
  app.add_flag("-s,--spiral", to_spiral, "Fit spiral.");

  CLI11_PARSE(app, argc, argv);

  context.to_spiral = to_spiral;
  std::set<size_t> spiral_sids;

  // 1. Read. This also sorts the input strokes based on indices.
  std::string base_filename = scap_filename;
  if (scap_filename.find_last_of("/\\") != std::string::npos) {
    base_filename = scap_filename.substr(scap_filename.find_last_of("/\\") + 1);
  }
  // Remove the "_out" from the name
  base_filename = base_filename.erase(base_filename.find_last_of("_"));
  if (base_filename.find(".") != std::string::npos) {
    base_filename = base_filename.erase(base_filename.find("."));
  }
  base_filename = base_filename.erase(base_filename.find_last_of("_"));

  std::map<size_t, Cluster> clusters;

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
        if (s.spiral)
          spiral_sids.emplace(s.stroke_ind);
      }
    }
    stroke_sampling_size =
      sum_points / stroke_sampling_total_ref * stroke_sampling_size_ref;
    stroke_sampling_size =
      std::min(input_thickness * 10,
               std::max(stroke_sampling_size_ref, stroke_sampling_size));

    if (stroke_sampling_size != stroke_sampling_size_ref) {
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

  // Prepare the clusters
  context.tighter_fit = true;
  if (!vis_folder.empty())
    make_cache_folder(vis_folder, false);
  {
    // Parameterize saved clusters
    std::string cache_folder = (cache_filename.empty())
                                 ? vis_folder + "/cache/"
                                 : cache_filename + "/cache/";
    auto file_exists = [](const std::string &file_name) {
      std::ifstream f(file_name.c_str());
      return f.good();
    };

    for (auto const &cc : input.clusters) {
      clusters[cc.first] = cc.second;
    }
    for (auto &cc : clusters) {
      std::string hash_str = get_hash_str(cc.second);
      std::string json_filename_all =
        cache_folder + "/" + hash_str + "_" + end_str + ".json";

      if ((!vis_folder.empty() || !cache_filename.empty()) &&
          !end_str.empty() && file_exists(json_filename_all)) {
        read_json(json_filename_all, cc.second);

        bool seen_spiral = false;
        for (auto &s : cc.second.strokes) {
          if (spiral_sids.count(s.stroke_ind)) {
            s.spiral = true;
            seen_spiral = true;
          }
        }
        for (auto &s : cc.second.original_input_strokes) {
          if (spiral_sids.count(s.stroke_ind)) {
            s.spiral = true;
          }
        }

        // Fit with thickness in the tighter setting
        if (vis_folder != cache_filename) {
          // Reparameterize spirals
          if (seen_spiral) {
            Input input_final;
            input_final.clusters[0] = cc.second;
            input_final.clusters[0].xsecs.clear();
            fea_param.parameterize(&input_final);
            cc.second = input_final.clusters[0];
          }

          if (cc.second.xsecs.size() < 2) {
            cc.second.xsecs.clear();
          }

          // Fit centerline
          {
            Input input_final;
            input_final.clusters[0] = cc.second;
            input_final.clusters[0].fit.centerline.clear();
            input_final.clusters[0].fit.widths.clear();
            input_final.clusters[0].fit.fit_sample_xsecs.clear();
            cc.second.fit = fea_fitting.fit(&input_final).begin()->second;
            cc.second.fit.cluster_idx = cc.first;
          }

          // Fit width
          {
            Input input_final;
            bool in_width = context.widths;
            context.widths = true;
            input_final.clusters[0] = cc.second;
            cc.second.fit = fea_fitting.fit(&input_final).begin()->second;
            cc.second.fit.cluster_idx = cc.first;
            context.widths = in_width;
          }

          transform_fit_curve(cc.second.fit, input.orig_center,
                              input_thickness);
        }
      } else {
        bool to_orient = true;
        bool test_time = true;
        fit_cluster(solve_context.norm_thickness, input.width, input.height,
                    cc.second, to_orient, test_time,
                    &solve_context.matching_pair);

        // Fit with thickness
        bool in_width = context.widths;
        Input input_final;
        context.widths = true;
        input_final.clusters[0] = cc.second;
        cc.second.fit = fea_fitting.fit(&input_final).begin()->second;
        cc.second.fit.cluster_idx = cc.first;

        transform_fit_curve(cc.second.fit, input.orig_center, input_thickness);

        context.widths = in_width;
      }
      assert(!cc.second.fit.centerline.empty());

      // Write to the solve cache
      solve_context.merged_cluster_cache[hash_str] = cc.second;
    }
  }
  solve_context.width = input.width;
  solve_context.height = input.height;
  solve_context.norm_thickness = 1;
  solve_context.input_thickness = input_thickness;

  ////////

  solve_context.stage = SolveContext::SolveStage::Endpoint;
  std::vector<Junction> junctions;
  junction_end_merge(input, clusters, junctions);

  std::set<size_t> s_has_junc;
  for (auto const &junc : junctions) {
    s_has_junc.emplace(junc.from.first);
    s_has_junc.emplace(junc.to.first);
  }

  std::map<int, FittedCurve> all_fits;
  for (auto const &cc : clusters) {
    all_fits[cc.first] = cc.second.fit;
  }

  // Continuation
  std::map<int, FittedCurve> cont_fits;
  std::vector<Junction> cont_junctions;
  std::vector<double> cont_angles;
  for (auto const &junc : junctions) {
    double angle;
    if (is_continuation(all_fits[junc.from.first], all_fits[junc.to.first],
                        junc, angle)) {
      cont_junctions.emplace_back(junc);
      cont_angles.emplace_back(angle);
    }
  }
  // Resolve branching
  resolve_continuation(cont_angles, cont_junctions);
  for (auto const &junc : cont_junctions) {
    cont_fits[junc.from.first] = all_fits[junc.from.first];
    cont_fits[junc.to.first] = all_fits[junc.to.first];
  }

  std::map<size_t, Cluster> merged_clusters = clusters;
  std::map<size_t, Cluster> merged_only_clusters;
  merge_ends(input, merged_clusters, merged_only_clusters, cont_junctions);

  {
    std::map<int, FittedCurve> final_fits;
    for (auto const &cc : merged_clusters) {
      if (cc.second.fit.centerline.empty())
        continue;
      assert(!cc.second.fit.widths.empty());
      final_fits[cc.first] = cc.second.fit;
    }
    // Individual fit output files
    std::ofstream ind_width_fit_svg(vis_folder +
                                    "/before_outlier_fit_width.svg");
    context.widths = true;
    fea_fitting.fit_svg(ind_width_fit_svg, input, final_fits, width, height);
    context.widths = false;
    double tmp_thickness = input.thickness;
    input.thickness = input_thickness;
    std::ofstream ind_fit_svg(vis_folder + "/before_outlier_fit.svg");
    fea_fitting.fit_svg(ind_fit_svg, input, final_fits, width, height);
    input.thickness = tmp_thickness;
  }

  // Remove outliers
  bool has_cut_left;
  std::map<int, FittedCurve> removed_fits;
  std::map<size_t, size_t> remvoed_sid;
  {
    std::vector<Junction> outlier_junctions;
    junction_end_merge(input, merged_clusters, outlier_junctions);

    std::map<size_t, Cluster> merged_clusters_single_removed;
    std::map<int, FittedCurve> merged_fits;
    for (auto const &cc : merged_clusters) {
      merged_fits[cc.first] = cc.second.fit;
    }
    bool skip_cut = false;
    has_cut_left = remove_outliers(input, merged_clusters, merged_fits,
                                   merged_clusters_single_removed,
                                   outlier_junctions, skip_cut);

    for (auto cc : merged_clusters) {
      if (!merged_clusters_single_removed.count(cc.first)) {
        removed_fits[cc.first] = cc.second.fit;
        for (size_t i = 0; i < removed_fits[cc.first].widths.size(); ++i) {
          removed_fits[cc.first].widths[i] *= 3;
        }
        for (auto const &s : cc.second.strokes) {
          remvoed_sid[s.stroke_ind] = cc.first;
        }
      }
    }
    merged_clusters = merged_clusters_single_removed;
    if (!has_cut_left) {
      std::ofstream removed_fit_svg(vis_folder + "/removed_fit.svg");
      context.widths = true;
      fea_fitting.fit_colored_svg(removed_fit_svg, input, removed_fits, width,
                                  height);
      context.widths = false;
      for (auto &s : capture.sketchedPolylines) {
        if (remvoed_sid.count(s.stroke_ind))
          s.group_ind = remvoed_sid[s.stroke_ind];
        else
          s.group_ind = -1;
      }
      std::ofstream scap_ofs(vis_folder + "/removed_strokes.scap");
      std::string out_buffer = capture.to_string();
      scap_ofs << "#" << width << "\t" << height << std::endl;
      scap_ofs.write(out_buffer.c_str(), out_buffer.size());
      scap_ofs.close();
    }
  }

  if (has_cut_left) {
    std::vector<Junction> outlier_junctions;
    junction_end_merge(input, merged_clusters, outlier_junctions);

    std::map<size_t, Cluster> merged_clusters_single_removed;
    std::map<int, FittedCurve> merged_fits;
    for (auto const &cc : merged_clusters) {
      merged_fits[cc.first] = cc.second.fit;
    }
    bool skip_cut = true;
    has_cut_left = remove_outliers(input, merged_clusters, merged_fits,
                                   merged_clusters_single_removed,
                                   outlier_junctions, skip_cut);

    for (auto cc : merged_clusters) {
      if (!merged_clusters_single_removed.count(cc.first)) {
        removed_fits[cc.first] = cc.second.fit;
        for (size_t i = 0; i < removed_fits[cc.first].widths.size(); ++i) {
          removed_fits[cc.first].widths[i] *= 3;
        }

        for (auto const &s : cc.second.strokes) {
          remvoed_sid[s.stroke_ind] = cc.first;
        }
      }
    }
    merged_clusters = merged_clusters_single_removed;
    std::ofstream removed_fit_svg(vis_folder + "/removed_fit.svg");
    context.widths = true;
    fea_fitting.fit_colored_svg(removed_fit_svg, input, removed_fits, width,
                                height);
    context.widths = false;

    for (auto &s : capture.sketchedPolylines) {
      if (remvoed_sid.count(s.stroke_ind))
        s.group_ind = remvoed_sid[s.stroke_ind];
      else
        s.group_ind = -1;
    }
    capture.thickness = input_thickness;
    std::ofstream scap_ofs(vis_folder + "/removed_strokes.scap");
    std::string out_buffer = capture.to_string();
    scap_ofs << "#" << width << "\t" << height << std::endl;
    scap_ofs.write(out_buffer.c_str(), out_buffer.size());
    scap_ofs.close();
  }

  // Output
  context.widths = true;
  std::map<int, FittedCurve> src_fits;
  for (auto const &cc : clusters) {
    if (cc.second.fit.centerline.empty() || !s_has_junc.count(cc.first))
      continue;
    assert(!cc.second.fit.widths.empty());
    src_fits[cc.first] = cc.second.fit;
  }

  input.thickness = input_thickness;
  std::ofstream refit_svg(vis_folder + "/junction_split_fit_width.svg");
  fea_fitting.fit_svg(refit_svg, input, all_fits, width, height);

  {
    context.widths = false;
    std::ofstream refit_svg(vis_folder + "/junction_split_fit.svg");
    fea_fitting.fit_svg(refit_svg, input, all_fits, width, height);
    context.widths = true;
  }

  std::ofstream fit_svg(vis_folder + "/junction_end_fit.svg");
  fea_fitting.fit_colored_svg(fit_svg, input, src_fits, width, height);

  std::ofstream cont_fit_svg(vis_folder + "/junction_cont_fit.svg");
  fea_fitting.fit_colored_svg(cont_fit_svg, input, cont_fits, width, height);

  std::map<int, FittedCurve> merged_fits;
  for (auto const &cc : merged_only_clusters) {
    if (cc.second.fit.centerline.empty())
      continue;
    assert(!cc.second.fit.widths.empty());
    merged_fits[cc.first] = cc.second.fit;
  }
  std::ofstream merged_fit_svg(vis_folder + "/junction_merged_fit.svg");
  fea_fitting.fit_colored_svg(merged_fit_svg, input, merged_fits, width,
                              height);

  // Final result
  std::map<int, FittedCurve> final_fits;
  for (auto const &cc : merged_clusters) {
    if (cc.second.fit.centerline.empty())
      continue;
    assert(!cc.second.fit.widths.empty());
    final_fits[cc.first] = cc.second.fit;
  }
  std::ofstream final_fit_svg(vis_folder + "/final_merged_fit.svg");
  fea_fitting.fit_svg(final_fit_svg, input, final_fits, width, height);

  // Final clustering
  std::map<size_t, size_t> s2c;
  for (const auto &cid2c : merged_clusters) {
    for (const auto &s : cid2c.second.strokes) {
      s2c[s.stroke_ind] = cid2c.first;
    }
  }

  for (auto &s : capture.sketchedPolylines) {
    assert(s2c.count(s.stroke_ind));
    s.group_ind = s2c[s.stroke_ind];
  }
  std::ofstream scap_ofs(vis_folder + "/" + base_filename + "_final_out.scap");
  std::string out_buffer = capture.to_string();
  scap_ofs << "#" << width << "\t" << height << std::endl;
  scap_ofs.write(out_buffer.c_str(), out_buffer.size());
  scap_ofs.close();

  // Individual fit output files
  std::ofstream ind_width_fit_svg(vis_folder + "/" + base_filename +
                                  "_fit_width.svg");
  fea_fitting.fit_svg(ind_width_fit_svg, input, final_fits, width, height);
  context.widths = false;
  std::ofstream ind_fit_svg(vis_folder + "/" + base_filename + "_fit.svg");
  fea_fitting.fit_svg(ind_fit_svg, input, final_fits, width, height);

  {
    Capture capture;
    for (const auto &c : merged_clusters) {
      capture.sketchedPolylines.emplace_back();
      capture.sketchedPolylines.back().stroke_ind = c.first;
      capture.sketchedPolylines.back().group_ind = c.first;
      for (const auto &p : final_fits.at(c.first).centerline)
        capture.sketchedPolylines.back().points.emplace_back(
          SketchUI::Point2D(p.x, p.y), 0);
    }
    capture.thickness = input_thickness;

    {
      std::string scap_out = vis_folder + "/" + base_filename + "_fit.scap";
      std::ofstream scap_ofs(scap_out);
      std::string out_buffer = capture.to_string();
      scap_ofs << "#" << width << "\t" << height << std::endl;
      scap_ofs.write(out_buffer.c_str(), out_buffer.size());
      scap_ofs.close();
    }
  }

  return 0;
}
