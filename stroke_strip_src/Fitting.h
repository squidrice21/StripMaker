#pragma once

#include <glm/glm.hpp>

#include <map>

#include "Cluster.h"
#include "Context.h"

struct Fitting {
  Fitting(Context &context);
  std::map<int, FittedCurve> fit(Input *input);
  FittedCurve fit_cluster(Cluster cluster);

  void fit_svg(std::ostream &os, const Input &input,
               const std::map<int, FittedCurve> &fits, double width = -1,
               double height = -1);

private:
  Context &context;
  // Default value:
  // const double K_WEIGHT = 1e3;
  const double K_WEIGHT = 1e1;
  const double POS_WEIGHT = 5e-5;
  const double MIN_DIST = 1e-4;

  struct Sample {
    glm::dvec2 point;
    glm::dvec2 tangent;
    double k;
    bool no_k;
    bool gap;
    double width;
  };
  std::vector<Sample> samples_from_xsec(const Cluster &cluster,
                                        size_t xsec_idx);

  std::vector<glm::dvec2> fit_tangents(const std::vector<Sample> &samples,
                                       bool periodic);
  std::vector<glm::dvec2> fit_positions(const std::vector<Sample> &samples,
                                        const std::vector<glm::dvec2> &tangents,
                                        bool periodic);
  std::vector<glm::dvec2> fit_positions4one_stroke(Cluster &cluster);
  std::vector<double> fit_widths(const std::vector<Sample> &samples,
                                 bool periodic);
};
