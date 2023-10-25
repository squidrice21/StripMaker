#define _USE_MATH_DEFINES
#include <cmath>

#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

#include <algorithm>
#include <deque>
#include <filesystem>
#include <iostream>

#include "../stroke_strip_src/Cluster.h"
#include "../stroke_strip_src/FittingEigenSparse.h"
#include "../stroke_strip_src/Parameterization.h"
#include "../stroke_strip_src/Polyline2D.h"
#include "../stroke_strip_src/StrokeCutting.h"
#include "../stroke_strip_src/StrokeOrientation.h"

#include "Util.h"

#include <CLI/CLI.hpp>

void to_scap(std::string const &filename, Capture const &capture, int width,
             int height) {
  std::ofstream scap_ofs(filename);
  std::string out_buffer;
  out_buffer += capture.to_string();

  scap_ofs << "#" << width << "\t" << height << std::endl;
  scap_ofs.write(out_buffer.c_str(), out_buffer.size());
  scap_ofs.close();
}

Input from_capture(Capture capture) {
  Input input;
  auto &clusters = input.clusters;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  input.thickness = capture.thickness;

  double rate = 2.75;

  for (auto &polyline : capture.sketchedPolylines) {
    polyline.reparameterize(
      std::min(rate * capture.thickness, polyline.totalLen() / 3));

    for (auto &point : polyline.points) {
      min_x = std::min(min_x, point.first.x);
      max_x = std::max(max_x, point.first.x);
      min_y = std::min(min_y, point.first.y);
      max_y = std::max(max_y, point.first.y);
    }
  }
  glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
  input.orig_center = center;

  input.width = (max_x - min_x) / capture.thickness;
  input.height = (max_y - min_y) / capture.thickness;

  for (auto &polyline : capture.sketchedPolylines) {
    if (polyline.points.empty())
      continue;

    clusters[polyline.group_ind].strokes.emplace_back();
    auto &stroke = clusters[polyline.group_ind].strokes.back();
    stroke.points.reserve(polyline.points.size());
    stroke.u.reserve(polyline.points.size());
    clusters[polyline.group_ind].strokes.back().stroke_ind =
      polyline.stroke_ind;
    clusters[polyline.group_ind].strokes.back().cluster_ind =
      polyline.group_ind;
    for (size_t i = 0; i < polyline.points.size(); ++i) {
      auto &point = polyline.points[i];

      stroke.points.push_back(
        (glm::dvec2(point.first.x, point.first.y) - center) /
        capture.thickness);
      stroke.u.push_back(0);
    }
  }

  return input;
}

Capture fit_capture(const Input &input, const auto fits) {
  Capture capture;
  auto &polylines = capture.sketchedPolylines;

  capture.thickness = input.thickness;

  double rate = 2.75;

  // Iterate over input's clusters and strokes
  size_t num_clusters = input.clusters.size();
  int stroke_count = 0;
  for (const auto &fit : fits) {
    SketchUI::Polyline2D polyline;
    auto centerline = fit.second.centerline;
    // Give the stroke an index and group ID
    polyline.stroke_ind = 1;
    polyline.group_ind = fit.second.cluster_idx;
    for (const auto &point : centerline) {
      // Add stroke width back and recenter
      glm::dvec2 original_point = (point * input.thickness) + input.orig_center;
      // Add the point to the polyline
      std::pair<SketchUI::Point2D, std::int64_t> original_point_object(
        SketchUI::Point2D(original_point.x, original_point.y), 0);
      polyline.points.emplace_back(original_point_object);
    }
    // Add the polyline to the capture
    capture.sketchedPolylines.push_back(polyline);
  }

  return capture;
}

/**
 * The entrance function for running StrokeStrip (an improved version that uses
 * higher sampling rate and runs faster) given the annotated input scap file.
 */
int main(int argc, char *argv[]) {
  CLI::App app{"Call StrokeStrip on clustered input."};

  bool enable_cut;

  std::string scap_filepath;
  std::string scap_filename;
  std::string output_dirname;

  app.add_option("-i,--input", scap_filepath, "The ground truth scap input.")
    ->required();
  app.add_option("-o,--outputdir", output_dirname, "The output directory.");

  CLI11_PARSE(app, argc, argv);

  std::filesystem::path filepath = scap_filepath;
  scap_filename = filepath.filename().string();
  scap_filename.erase(scap_filename.length() - 5, 5); // remove .scap

  if (output_dirname.empty()) {
    output_dirname = filepath.parent_path().string();
  } else {
    std::filesystem::create_directory(output_dirname);
  }

  Context context; // cut is set to false by default

  Input input;

  std::map<size_t, size_t> inv_reindexing;

  // 1. Preprocess
  {
    std::ifstream scap_ifs(scap_filepath);
    std::stringstream buffer;
    buffer << scap_ifs.rdbuf();
    scap_ifs.close();

    Capture capture;
    capture.from_string(buffer.str());

    std::string canvas_size_line =
      buffer.str().substr(0, buffer.str().find_first_of("\n"));
    int width;
    int height;
    sscanf(canvas_size_line.c_str(), "#%d\t%d", &width, &height);

    // Assign missing ids
    int max_ind = -1;
    for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
      max_ind = std::max(max_ind, capture.getSketchedPolyline(i).stroke_ind);
    }

    if (context.cut) {
      inv_reindexing =
        preprocess_cluster(1, width, height, &capture, context.cut);
    }

    input = from_capture(capture);
  }

  // 2. Orientation
  {
    StrokeOrientation orientation(context);
    orientation.orient_strokes(&input);
    if (context.debug_viz) {
      std::string final_output_name = scap_filename + "_param_debug.svg";
      final_output_name = output_dirname + "/" + final_output_name;
      std::ofstream orientation_svg(final_output_name);
      orientation.orientation_debug(orientation_svg, input);
    }
    std::map<int, std::vector<int>> orientations =
      orientation.orient_strokes(&input);
    orientation.flip_strokes(&input, orientations);
  }

  // 3. Parameterization
  {
    Parameterization param(context);
    param.parameterize(&input);
    if (context.debug_viz) {
      std::string final_output_name = scap_filename + "_param_debug.svg";
      final_output_name = output_dirname + "/" + final_output_name;
      std::ofstream debug_svg(final_output_name);
      param.debug_svg(debug_svg, input);
    }
    {
      std::string final_output_name = scap_filename + "_isolines.svg";
      final_output_name = output_dirname + "/" + final_output_name;
      std::ofstream isolines_svg(final_output_name);
      param.isolines_svg(isolines_svg, -1, input);
    }
  }

  // 4. Fitting
  {
    FittingEigenSparse fitting(context);
    auto fits = fitting.fit(&input);
    for (auto &c_fit : fits) {
      input.clusters[c_fit.first].fit = c_fit.second;
    }
    {
      std::string final_output_name = scap_filename + "_fit.svg";
      final_output_name = output_dirname + "/" + final_output_name;
      std::ofstream fit_svg(final_output_name);
      fitting.fit_svg(fit_svg, input, fits);
    }
    {
      std::string final_output_name = scap_filename + "_fit.scap";
      final_output_name = output_dirname + "/" + final_output_name;
      Capture result = fit_capture(input, fits);
      to_scap(final_output_name, result, input.width * input.thickness,
              input.height * input.thickness);
    }
  }

  {
    context.widths = true;
    FittingEigenSparse fitting(context);
    auto fits = fitting.fit(&input);
    for (auto &c_fit : fits) {
      input.clusters[c_fit.first].fit = c_fit.second;
    }
    {
      std::string final_output_name = scap_filename + "_fit_width.svg";
      final_output_name = output_dirname + "/" + final_output_name;
      std::ofstream fit_svg(final_output_name);
      fitting.fit_svg_center(fit_svg, input, fits, input.orig_center);
    }
    context.widths = false;
  }

  {
    std::string final_output_name = scap_filename + "_params.svg";
    final_output_name = output_dirname + "/" + final_output_name;
    std::ofstream param_svg(final_output_name);
    input.param_svg(param_svg, context.rainbow);
  }

  {
    std::string final_output_name = scap_filename + "_orientation.svg";
    final_output_name = output_dirname + "/" + final_output_name;
    std::ofstream orientation_svg(final_output_name);
    input.orientation_svg(orientation_svg);
  }

  return 0;
}
