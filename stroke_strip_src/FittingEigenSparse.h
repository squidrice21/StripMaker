#pragma once

#include "Cluster.h"
#include "Fitting.h"
#include <vector>

struct FittingEigenSparse {
  FittingEigenSparse(Context &context);
  std::map<int, FittedCurve> fit(Input *input);
  FittedCurve fit_cluster(Cluster &cluster);

  void fit_svg(std::ostream &os, const Input &input,
               const std::map<int, FittedCurve> &fits, double width = -1,
               double height = -1);
  void fit_svg_center(std::ostream &os, const Input &input,
                      const std::map<int, FittedCurve> &fits,
                      glm::dvec2 const &center, double width = -1,
                      double height = -1);
  void fit_colored_svg(std::ostream &os, const Input &input,
                       const std::map<int, FittedCurve> &fits,
                       double width = -1, double height = -1);
  void fit_svg_cluster(std::ostream &os, const Input &input,
                       const FittedCurve &kv,
                       const std::map<int, FittedCurve> &fits,
                       double width = -1, double height = -1);
  void fit_plain_text(std::ostream &os, const Input &input,
                      const std::map<int, FittedCurve> &fits);

private:
  Context &context;

  const double K_WEIGHT = 1e1;
  const double TIGHT_K_WEIGHT = 2;
  const double POS_WEIGHT = 5e-5;
  const double TIGHT_ENDPOINT_WEIGHT = 8.;
  const double TIGHT_POS_WEIGHT = 4e-4;
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
                                        const std::vector<Cluster::XSec> &xsecs,
                                        size_t xsec_idx);

  std::vector<glm::dvec2> fit_tangents(const std::vector<Sample> &samples,
                                       bool periodic);
  std::vector<glm::dvec2> fit_positions(const std::vector<Sample> &samples,
                                        const std::vector<glm::dvec2> &tangents,
                                        bool periodic);
  std::vector<double> fit_widths(const std::vector<Sample> &samples,
                                 bool periodic);
};
