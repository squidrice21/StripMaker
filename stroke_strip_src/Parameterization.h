#pragma once
#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

#include "Cluster.h"
#include "Context.h"

class Parameterization {
public:
  // HEAD
  Parameterization(Context &context);
  std::map<std::string, double> parameterize(Input *input);
  void isolines_svg(std::ostream &os, int new_stroke_ind, const Input &input);
  void orthogonal_isolines_svg(std::ostream &os, int new_stroke_ind,
                               const Input &input, size_t itr);
  void isolines_cluster_svg(std::ostream &os, int new_stroke_ind,
                            const Input &input, int cluster_num);
  void save_svg(std::ostream &os, const Input &input);
  void save_non_isoline_svg(std::ostream &os, const Input &input);

  // ddataset
  void isolines_only_cluster_svg(std::ostream &os, int new_stroke_ind,
                                 const Input &input, int cluster_num,
                                 std::map<int, int> stroke_to_cluster);
  void isolines_only_divided_cluster_svg(std::ostream &os, int new_stroke_ind,
                                         const Input &input, int cluster_num,
                                         std::vector<int> base_stroke_ind,
                                         std::vector<int> target_stroke_ind);

  void debug_svg(std::ostream &os, const Input &input);
  void debug_svg(std::ostream &os, int new_stroke_ind, const Input &input);

  // Exposing these functions for debugging parameterization failure
  std::vector<Cluster::XSec> orthogonal_xsecs(const Cluster &cluster,
                                              bool periodic_check = false);
  void ensure_connected(Cluster *cluster, Cluster::XSec *cut = nullptr);

  std::vector<Cluster::XSec>
  xsecs_from_params(const Cluster &cluster, bool ensure_sample_reference,
                    double sample_rate = 1.0, bool nonlinear = false);
  std::vector<Cluster::XSec>
  xsecs_from_single_stroke_params(const Cluster &cluster,
                                  bool ensure_sample_reference,
                                  double sample_rate = 1.0,
                                  bool nonlinear = false);

private:
  bool viz;

  Context &context;
  std::mutex viz_lock;

  struct DebugLine {
    glm::dvec2 from;
    glm::dvec2 to;
    std::string color;
  };
  void add_debug_line(DebugLine line);
  std::vector<DebugLine> debug_lines;

  std::map<std::string, double> parameterize_cluster(Cluster *cluster);
  std::map<std::string, double>
  params_from_xsecs(Cluster *cluster, bool relaxed, bool initial = false,
                    Cluster::XSec *cut = nullptr, double lambda = 0);
  std::map<std::string, double>
  params_from_xsecs_regularized(Cluster *cluster, bool relaxed,
                                bool initial = false,
                                Cluster::XSec *cut = nullptr);
  std::map<std::string, double>
  params_from_xsecs_retry(Cluster *cluster, double &final_lambda, bool relaxed,
                          bool initial = false, Cluster::XSec *cut = nullptr);
  void check_periodic(Cluster *cluster);

  Cluster::XSec xsec_at_u(const Cluster &cluster, double u);
  Cluster::XSec xsec_at_single_stroke_u(const Cluster &cluster, double u);
  Cluster::XSec orthogonal_xsec_at(const Cluster &cluster, size_t stroke,
                                   double i, bool is_fast = false);
};
