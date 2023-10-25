#include <cmath>

#include "FittingEigenSparse.h"
#include "Parameterization.h"
#include "SvgUtils.h"
#include "Utils.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <ostream>

// Note this kernel size for the tangent and curvature computation is tied with
// input stroke sampling step size
double d_i = 0.25;

FittingEigenSparse::FittingEigenSparse(Context &context)
  : context(context) {}

void FittingEigenSparse::fit_svg(std::ostream &os, const Input &input,
                                 const std::map<int, FittedCurve> &fits,
                                 double width, double height) {
  double padding = input.thickness;

  // Center
  glm::dvec2 center(0, 0);
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &c_fit : fits) {
    for (const auto &point : c_fit.second.centerline) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }
  }
  center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
  // double w = (max_x - min_x) * input.thickness + 2 * padding;
  // double h = (max_y - min_y) * input.thickness + 2 * padding;
  double w = (max_x - min_x) + 2 * padding;
  double h = (max_y - min_y) + 2 * padding;
  double x = -w / 2;
  double y = -h / 2;

  if (width < 0 || height < 0)
    SVG::begin(os, x, y, w, h);
  else
    SVG::begin(os, 0, 0, width, height);
  for (auto &kv : fits) {
    std::vector<glm::dvec2> path = kv.second.centerline;
    std::vector<double> w = kv.second.widths;
    if (path.empty())
      continue;
    if (width < 0 || height < 0) {
      for (auto &pt : path) {
        pt -= center;
        pt *= input.thickness;
      }
      for (auto &val : w) {
        val *= input.thickness;
      }
    }
    if (context.widths && path.size() > 1) {
      SVG::variable_polyline(os, path, w);
    } else {
      SVG::polyline(os, path, input.thickness);
    }
  }
  SVG::end(os);
}

void FittingEigenSparse::fit_svg_center(std::ostream &os, const Input &input,
                                        const std::map<int, FittedCurve> &fits,
                                        glm::dvec2 const &center, double width,
                                        double height) {
  double padding = input.thickness;

  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &c_fit : fits) {
    for (const auto &point : c_fit.second.centerline) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }
  }
  double w = (max_x - min_x) * input.thickness + 2 * padding;
  double h = (max_y - min_y) * input.thickness + 2 * padding;

  SVG::begin(os, 0, 0, w, h);
  for (auto &kv : fits) {
    std::vector<glm::dvec2> path = kv.second.centerline;
    std::vector<double> w = kv.second.widths;
    if (path.empty())
      continue;
    if (width < 0 || height < 0) {
      for (auto &pt : path) {
        pt = (pt * input.thickness) + center;
      }
      for (auto &val : w) {
        val *= input.thickness;
      }
    }
    if (context.widths && path.size() > 1) {
      SVG::variable_polyline(os, path, w);
    } else {
      SVG::polyline(os, path, input.thickness);
    }
  }
  SVG::end(os);
}

void FittingEigenSparse::fit_plain_text(
  std::ostream &os, const Input &input,
  const std::map<int, FittedCurve> &fits) {
  double padding = input.thickness;

  // Center
  glm::dvec2 center(0, 0);
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &c_fit : fits) {
    for (const auto &point : c_fit.second.centerline) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }
  }
  center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
  // double w = (max_x - min_x) * input.thickness + 2 * padding;
  // double h = (max_y - min_y) * input.thickness + 2 * padding;
  double w = (max_x - min_x) + 2 * padding;
  double h = (max_y - min_y) + 2 * padding;
  double x = -w / 2;
  double y = -h / 2;

  os << "0,0," << w << "," << h << std::endl;
  for (auto &kv : fits) {
    std::vector<glm::dvec2> path = kv.second.centerline;
    std::vector<double> w = kv.second.widths;
    if (path.empty())
      continue;
    for (size_t i = 0; i < path.size(); ++i) {
      os << path[i].x << "," << path[i].y << ","
         << ((context.widths) ? w[i] : 1.0) << ";";
    }
    os << std::endl;
  }
}

void FittingEigenSparse::fit_colored_svg(std::ostream &os, const Input &input,
                                         const std::map<int, FittedCurve> &fits,
                                         double width, double height) {
  const std::vector<std::string> colors = {
    "#00a391", "#a52a2a", "#ffd700", "#45769e", "#f77ff2", "#ff0000",
    "#db7093", "#1e90ff", "#00ee00", "#ff8c00", "#9046cc", "#6a7b53",
    "#ff1493", "#7fffd4", "#9acd32", "#ba55d3", "#ff00ff", "#d2691e",
    "#458c45", "#954595", "#f0e68c", "#00bfff", "#ffa07a",
  };
  double padding = input.thickness;

  // Center
  glm::dvec2 center(0, 0);
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &c_fit : fits) {
    for (const auto &point : c_fit.second.centerline) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }
  }
  center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
  // double w = (max_x - min_x) * input.thickness + 2 * padding;
  // double h = (max_y - min_y) * input.thickness + 2 * padding;
  double w = (max_x - min_x) + 2 * padding;
  double h = (max_y - min_y) + 2 * padding;
  double x = -w / 2;
  double y = -h / 2;

  if (width < 0 || height < 0)
    SVG::begin(os, x, y, w, h);
  else
    SVG::begin(os, 0, 0, width, height);
  for (auto &kv : fits) {
    std::vector<glm::dvec2> path = kv.second.centerline;
    std::vector<double> w = kv.second.widths;
    if (path.empty())
      continue;
    if (width < 0 || height < 0) {
      for (auto &pt : path) {
        pt -= center;
        pt *= input.thickness;
      }
      for (auto &val : w) {
        val *= input.thickness;
      }
    }

    std::string color = colors[kv.first % colors.size()];
    if (context.widths && path.size() > 1) {
      SVG::variable_polyline(os, path, w, color);
    } else {
      SVG::polyline(os, path, input.thickness, color);
    }
  }
  SVG::end(os);
}

void FittingEigenSparse::fit_svg_cluster(std::ostream &os, const Input &input,
                                         const FittedCurve &kv,
                                         const std::map<int, FittedCurve> &fits,
                                         double width, double height) {
  double padding = input.thickness;

  // Center
  glm::dvec2 center(0, 0);
  double min_x = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();

  for (const auto &c_fit : fits) {
    for (const auto &point : c_fit.second.centerline) {
      min_x = std::min(min_x, point.x);
      max_x = std::max(max_x, point.x);
      min_y = std::min(min_y, point.y);
      max_y = std::max(max_y, point.y);
    }
  }
  center = glm::dvec2((max_x + min_x) / 2, (max_y + min_y) / 2);
  double w = (max_x - min_x) * input.thickness + 2 * padding;
  double h = (max_y - min_y) * input.thickness + 2 * padding;
  double x = -w / 2;
  double y = -h / 2;

  if (width < 0 || height < 0)
    SVG::begin(os, x, y, w, h);
  else
    SVG::begin(os, 0, 0, width, height);
  {
    std::vector<glm::dvec2> path = kv.centerline;
    std::vector<double> w = kv.widths;
    if (width < 0 || height < 0) {
      for (auto &pt : path) {
        pt -= center;
        pt *= input.thickness;
      }
      for (auto &val : w) {
        val *= input.thickness;
      }
    }
    if (context.widths) {
      SVG::variable_polyline(os, path, w);
    } else {
      SVG::polyline(os, path, input.thickness);
    }
  }
  SVG::end(os);
}

std::map<int, FittedCurve> FittingEigenSparse::fit(Input *input) {
  return map_clusters<FittedCurve>(*input, [&](Cluster &c) {
    auto begin = std::chrono::high_resolution_clock::now();

    // This means parameterization failed
    if (c.strokes.empty())
      return FittedCurve();
    glm::dvec2 center;
    center_cluster(c, center);
    for (auto &p : c.fit.centerline)
      p -= center;
    for (auto &s : c.fit.fit_sample_xsecs) {
      for (auto &p : s.points) {
        p.point -= center;
      }
    }

    FittedCurve fit = fit_cluster(c);

    decenter_cluster(c, center);
    for (auto &p : fit.centerline)
      p += center;
    for (auto &s : fit.fit_sample_xsecs) {
      for (auto &p : s.points) {
        p.point += center;
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    context.fitting_timer +=
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin)
        .count() /
      1000.0;

    return fit;
  });
}

FittedCurve FittingEigenSparse::fit_cluster(Cluster &cluster) {
  FittedCurve curve;
  curve.fit_sample_xsecs = cluster.fit.fit_sample_xsecs;

  if (cluster.xsecs.size() == 1) {
    curve.centerline.emplace_back(cluster.xsecs.front().avg_point());

    if (context.widths) {
      curve.widths.emplace_back(cluster.xsecs.front().xsec_width());
    }

    curve.cluster_idx = cluster.strokes[0].cluster_ind;
    return curve;
  }
  if (cluster.strokes.size() == 1 &&
      !(context.to_spiral && cluster.strokes.front().spiral)) {
    double s_len = 0;
    cluster.strokes.front().u[0] = 0;
    for (size_t i = 1; i < cluster.strokes.front().points.size(); ++i) {
      s_len += glm::distance(cluster.strokes.front().points[i - 1],
                             cluster.strokes.front().points[i]);
      cluster.strokes.front().u[i] = s_len;
    }
    Parameterization parameterization(context);
    // double rate = cluster.periodic ? 0.025 : 0.1;
    bool ensure_sample_reference = false;
    double rate = cluster.periodic ? 0.5 : 2;
    curve.fit_sample_xsecs = parameterization.xsecs_from_single_stroke_params(
      cluster, ensure_sample_reference, rate, !cluster.periodic);

    for (auto const &xsec : curve.fit_sample_xsecs) {
      curve.centerline.emplace_back(xsec.avg_point());
      if (context.widths) {
        curve.widths.emplace_back(1);
      }
    }

    curve.cluster_idx = cluster.strokes[0].cluster_ind;
    return curve;
  }

  if (curve.fit_sample_xsecs.empty())
    curve.fit_sample_xsecs = cluster.xsecs;

  bool fit_width_with_init =
    (context.widths && !curve.fit_sample_xsecs.empty());

  // 1. Get xsec samples more densely sampled near endpoints
  if (!fit_width_with_init) {
    Parameterization parameterization(context);
    // double rate = cluster.periodic ? 0.025 : 0.1;
    bool ensure_sample_reference = false;
    double rate = cluster.periodic ? 0.5 : 2;

    curve.fit_sample_xsecs = parameterization.xsecs_from_params(
      cluster, ensure_sample_reference, rate, !cluster.periodic);
    for (auto &xsec : curve.fit_sample_xsecs) {
      sort_xsec(xsec);
    }
  }

  // 2. Make one fitting sample per xsec
  std::vector<FittingEigenSparse::Sample> samples;
  for (size_t i = 0; i < curve.fit_sample_xsecs.size(); ++i) {
    auto res = samples_from_xsec(cluster, curve.fit_sample_xsecs, i);
    samples.insert(samples.end(), res.begin(), res.end());
  }
  // std::cout << "step2 finish" << std::endl;

  bool seen_spiral = false;
  for (auto const &s : cluster.strokes) {
    if (s.spiral) {
      seen_spiral = true;
      break;
    }
  }

  if (!fit_width_with_init) {
    // 3. Sovle for tangents
    auto tangents = fit_tangents(samples, cluster.periodic);

    // std::cout << "step3 finish" << std::endl;
    // std::cout << cluster.strokes.size() << std::endl;

    // 4. Solve for positions
    curve.centerline = fit_positions(samples, tangents, cluster.periodic);

    // Enforce closure if the inputs contain spiral strokes and the endpoints
    // are not too far away
    if (seen_spiral && !cluster.periodic &&
        glm::distance(curve.centerline.front(), curve.centerline.back()) <
          0.1 * total_length(curve.centerline)) {
      cluster.periodic = true;
      tangents = fit_tangents(samples, cluster.periodic);
      curve.centerline = fit_positions(samples, tangents, cluster.periodic);
    }
  } else {
    curve.centerline = cluster.fit.centerline;
  }

  // std::cout << "step3.5 finish" << std::endl;
  if (context.widths) {
    curve.widths = fit_widths(samples, cluster.periodic);

    // Debug: Directly output init widths
    // curve.widths.clear();
    // for (auto const &s : samples)
    //   curve.widths.emplace_back(s.width);
    // if (cluster.periodic) {
    //   curve.widths.push_back(curve.widths.front());
    // }
  }
  // std::cout << "step4 finish" << std::endl;
  curve.cluster_idx = cluster.strokes[0].cluster_ind;
  return curve;
}

std::vector<double>
FittingEigenSparse::fit_widths(const std::vector<Sample> &samples,
                               bool periodic) {
  GRBModel model(context.grb);
  const unsigned int N = samples.size();
  // std::cout << "model finish 1" << std::endl;
  //  Create width variables
  std::vector<GRBVar> vars;
  vars.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    GRBVar w = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
    vars.push_back(w);
  }
  // std::cout << "model finish 2" << std::endl;
  std::vector<double> dists;
  for (size_t i = 1; i < N; ++i) {
    // Use 1e-4 to avoid Nan created by division of 0
    dists.push_back(
      std::max(1e-4, glm::distance(samples[i].point, samples[i - 1].point)));
  }

  // Distances used in Laplacian matrix
  std::vector<double> laplacian_dists;
  {
    double total = 0.0;
    laplacian_dists.push_back(dists.front());
    total += dists.front();

    for (auto d : dists) {
      laplacian_dists.push_back(d);
      total += d;
    }

    laplacian_dists.push_back(dists.back());
    total += dists.back();

    for (auto &d : laplacian_dists) {
      d /= total;
    }
  }

  // Weight each distance matching term by the distance to the next point
  std::vector<double> dist_weights;
  {
    double total = 0.0;
    for (auto d : dists) {
      dist_weights.push_back(d);
      total += d;
    }

    dist_weights.push_back(dists.back());
    total += dists.back();

    for (auto &d : dist_weights) {
      d /= total;
    }
  }

  // Laplacian matrix
  std::vector<std::vector<double>> L;
  for (size_t i = 0; i < N; ++i) {
    std::vector<double> col(N, 0.0);
    L.push_back(col);
  }
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      L[i][j] = 0.0;
    }
  }
  for (size_t i = 0; i < N; ++i) {
    if (context.taper_widths) {
      // Zero Dirichlet conditions
      L[i][i] = 1.0 / laplacian_dists[i] + 1.0 / laplacian_dists[i + 1];
      if (i > 0) {
        L[i][i - 1] = -1.0 / laplacian_dists[i];
      }
      if (i + 1 < N) {
        L[i][i + 1] = -1.0 / laplacian_dists[i + 1];
      }
    } else {
      // Zero Neumann conditions
      if (i > 0) {
        L[i][i - 1] = -1.0 / laplacian_dists[i - 1];
      }
      if (i > 0 && i + 1 < N) {
        // Middle
        L[i][i] = 1.0 / laplacian_dists[i - 1] + 1.0 / laplacian_dists[i];
      } else if (i == 0) {
        // First
        L[i][i] = 1.0 / laplacian_dists[i];
      } else {
        // Last
        L[i][i] = 1.0 / laplacian_dists[i - 1];
      }
      if (i + 1 < N) {
        L[i][i + 1] = -1.0 / laplacian_dists[i];
      }
    }
  }

  // Laplacian objective: w^T * L * w;
  GRBQuadExpr wT_L_w;
  //// std::cout << "model finish 3" << std::endl;
  for (size_t i = 0; i < N; ++i) {
    GRBLinExpr L_w;
    for (size_t j = 0; j < N; ++j) {
      L_w += L[i][j] * vars[j];
    }
    wT_L_w += vars[i] * L_w;
  }
  // std::cout << "model finish 4" << std::endl;
  // Width matching objectives
  std::vector<GRBLinExpr> width_matches;
  width_matches.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    width_matches.push_back(dist_weights[i] * (vars[i] - samples[i].width));
  }
  // std::cout << "model finish 5" << std::endl;
  try {
    model.setObjective(wT_L_w + 200.0 * l2_norm_sq(&model, width_matches),
                       GRB_MINIMIZE);
  } catch (GRBException e) {
    std::cout << "Error code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;

    throw e;
  }

  // Add constraints
  int taper_len = std::min<int>(N / 6, 500);
  for (size_t i = 0; i < N; ++i) {
    bool has_minimum = !(context.taper_widths && !periodic &&
                         (i < taper_len || i > N - taper_len));
    if (has_minimum) {
      model.addConstr(vars[i] >= std::max(1.0, 0.45 * samples[i].width));
      // model.addConstr(vars[i] >= std::max(1.0, 0.8 * samples[i].width));
    } else {
      model.addConstr(vars[i] >= 1.0);
    }
  }
  for (size_t i : {size_t(0), size_t(N - 1)}) {
    model.addConstr(vars[i] <= std::max(1.0, 1.5 * samples[i].width));
    // model.addConstr(vars[i] <= std::max(1.0, 2.0 * samples[i].width));
  }

  context.optimize_model(&model);

  std::vector<double> opt_widths;
  opt_widths.reserve(vars.size());
  for (auto &var : vars) {
    opt_widths.push_back(var.get(GRB_DoubleAttr_X));
  }

  if (periodic) {
    opt_widths.push_back(opt_widths.front());
  }

  return opt_widths;
}

std::vector<glm::dvec2>
FittingEigenSparse::fit_positions(const std::vector<Sample> &samples,
                                  const std::vector<glm::dvec2> &tangents,
                                  bool periodic) {
  size_t dim = samples.front().tangent.length();
  size_t rows = (periodic) ? 2 * samples.size() : 2 * samples.size() - 1;
  Eigen::SparseMatrix<double> A(rows, samples.size());

  // A
  // Forward difference matrix D
  double endpoint_weight =
    (context.tighter_fit) ? std::sqrt(TIGHT_ENDPOINT_WEIGHT) : std::sqrt(2.);
  const double sqrt_gamma =
    (context.tighter_fit) ? std::sqrt(TIGHT_POS_WEIGHT) : std::sqrt(POS_WEIGHT);
  std::vector<double> W_diag;
  W_diag.reserve(samples.size());
  for (size_t j = 0; j < samples.size(); j++) {
    if (j + 1 < samples.size())
      W_diag.emplace_back(
        std::max(MIN_DIST, glm::dot(samples[j + 1].point - samples[j].point,
                                    tangents[j])));
    else if (periodic)
      W_diag.emplace_back(std::max(
        MIN_DIST, glm::dot(samples[0].point - samples[j].point, tangents[j])));
  }

  // A
  // Forward difference matrix D
  for (size_t i = 0; i < tangents.size(); i++) {
    if (i + 1 != samples.size()) {
      A.coeffRef(i, i) = -1 / W_diag[i];
      A.coeffRef(i, i + 1) = 1 / W_diag[i];
    } else if (periodic) {
      A.coeffRef(i, i) = -1 / W_diag[i];
      A.coeffRef(i, 0) = 1 / W_diag[i];
    }
  }

  // Position term
  for (size_t i = 0; i < samples.size(); i++) {
    A.coeffRef(tangents.size() + i, i) = sqrt_gamma;

    if (!periodic && (i == 0 || i == samples.size() - 1))
      A.coeffRef(tangents.size() + i, i) *= endpoint_weight;
  }
  A.makeCompressed();

  std::vector<Eigen::VectorXd> b;
  b.resize(tangents.front().length(), Eigen::VectorXd(rows));

  // b
  // tangent
  for (size_t i = 0; i < tangents.size(); i++) {
    for (size_t j = 0; j < tangents.front().length(); j++) {
      b[j][i] = tangents[i][j];
    }
  }

  // Initial sample positions
  for (size_t i = 0; i < samples.size(); i++) {
    for (size_t j = 0; j < tangents.front().length(); j++) {
      b[j][tangents.size() + i] = sqrt_gamma * samples[i].point[j];

      if (!periodic && (i == 0 || i == samples.size() - 1))
        b[j][tangents.size() + i] *= endpoint_weight;
    }
  }

  // Solve for the normal equation
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                        Eigen::COLAMDOrdering<int>>
    solver;
  solver.compute(A.transpose() * A);
  if (solver.info() != Eigen::Success) {
    abort();
  }

  std::vector<Eigen::VectorXd> x;
  x.resize(samples.front().tangent.length());
  for (size_t j = 0; j < samples.front().tangent.length(); j++) {
    x[j] = solver.solve(A.transpose() * b[j]);
    if (solver.info() != Eigen::Success) {
      abort();
    }
  }

  std::vector<glm::dvec2> results;
  results.resize(samples.size());

  for (size_t i = 0; i < samples.size(); i++) {
    for (size_t j = 0; j < tangents.front().length(); j++) {
      results[i][j] = x[j][i];
    }
  }

  if (periodic) {
    results.push_back(results.front());
  }

  return results;
}

std::vector<glm::dvec2>
FittingEigenSparse::fit_tangents(const std::vector<Sample> &samples,
                                 bool periodic) {
  size_t dim = samples.front().tangent.length();
  size_t rows = (periodic) ? 2 * samples.size() : 2 * samples.size() - 1;
  Eigen::SparseMatrix<double> A(rows, samples.size());

  // A
  // Forward difference matrix D
  double k_weight =
    (context.tighter_fit) ? std::sqrt(TIGHT_K_WEIGHT) : std::sqrt(K_WEIGHT);
  for (size_t i = 0; i < samples.size(); i++) {
    double scale = samples[i].no_k ? 0.1 : 1.;
    scale *= k_weight;
    if (i + 1 != samples.size()) {
      A.coeffRef(i, i) = -1 * scale;
      A.coeffRef(i, i + 1) = 1 * scale;
    } else if (periodic) {
      A.coeffRef(i, i) = -1 * scale;
      A.coeffRef(i, 0) = 1 * scale;
    }
  }
  // Positional term (matching the input tangents)
  double endpoint_weight = std::sqrt(2.);
  for (size_t i = 0; i < samples.size(); i++) {
    size_t j = (periodic) ? samples.size() + i : samples.size() - 1 + i;
    A.coeffRef(j, i) +=
      (!periodic && (i == 0 || i == samples.size() - 1)) ? endpoint_weight : 1;
  }
  A.makeCompressed();

  std::vector<Eigen::VectorXd> b;
  b.resize(samples.front().tangent.length(), Eigen::VectorXd(rows));

  // b
  // Curvature
  for (size_t i = 0; i < samples.size(); i++) {
    for (size_t k = 0; k < samples.front().tangent.length(); k++) {
      size_t pos_i = (periodic) ? samples.size() + i : samples.size() - 1 + i;
      b[k][pos_i] = samples[i].tangent[k];
    }

    if (periodic || i < samples.size() - 1) {
      double scale = samples[i].no_k ? 0.1 : 1.;
      scale *= k_weight;
      int j = (i + 1) % samples.size();
      glm::dvec2 avg_tangent =
        glm::normalize(samples[i].tangent + samples[j].tangent);
      glm::dvec2 ortho = normal(avg_tangent);

      for (size_t k = 0; k < samples.front().tangent.length(); k++) {
        b[k][i] = scale * samples[i].k * ortho[k];
      }
    }
  }

  // Solve for the normal equation
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Lower,
                        Eigen::COLAMDOrdering<int>>
    solver;
  auto AtA = A.transpose() * A;
  solver.compute(AtA);
  if (solver.info() != Eigen::Success) {
    abort();
  }

  std::vector<Eigen::VectorXd> x;
  x.resize(samples.front().tangent.length());
  for (size_t j = 0; j < samples.front().tangent.length(); j++) {
    x[j] = solver.solve(A.transpose() * b[j]);
    if (solver.info() != Eigen::Success) {
      abort();
    }
  }

  std::vector<glm::dvec2> tangents;
  tangents.resize(samples.size());

  for (size_t i = 0; i < samples.size(); i++) {
    for (size_t j = 0; j < samples.front().tangent.length(); j++) {
      tangents[i][j] = x[j][i];
    }
    tangents[i] = glm::normalize(tangents[i]);
  }
  if (!periodic)
    tangents.pop_back();
  return tangents;
}

std::vector<FittingEigenSparse::Sample>
FittingEigenSparse::samples_from_xsec(const Cluster &cluster,
                                      const std::vector<Cluster::XSec> &xsecs,
                                      size_t xsec_idx) {
  std::vector<FittingEigenSparse::Sample> result;
  auto &xsec = xsecs[xsec_idx];

  if (xsec.connector &&
      glm::distance(xsec.points[0].point, xsec.points[1].point) > 1e-4) {
    glm::dvec2 pt1 = xsec.points[0].point;
    glm::dvec2 pt2 = xsec.points[1].point;
    glm::dvec2 tan1 = xsec.points[0].tangent;
    glm::dvec2 tan2 = xsec.points[1].tangent;

    glm::dvec2 tan = xsec.avg_tangent();
    glm::dvec2 ortho = normal(tan);

    // Find out which point leads into the next
    size_t stroke1 = xsec.points[0].stroke_idx;
    if (xsec_idx > 0) {
      auto &prev = xsecs[xsec_idx - 1];
      if (!std::any_of(prev.points.begin(), prev.points.end(),
                       [=](const Cluster::XSecPoint &p) {
                         return p.stroke_idx == stroke1;
                       })) {
        std::swap(pt1, pt2);
        std::swap(tan1, tan2);
      }
    } else {
      auto &next = xsecs[xsec_idx + 1];
      if (std::any_of(next.points.begin(), next.points.end(),
                      [=](const Cluster::XSecPoint &p) {
                        return p.stroke_idx == stroke1;
                      })) {
        std::swap(pt1, pt2);
        std::swap(tan1, tan2);
      }
    }

    double dist = glm::distance(pt1, pt2);
    size_t num_samples = std::max(1., std::ceil(dist / 0.05));
    double sign = glm::dot(tan2 - tan1, ortho) > 0 ? 1 : -1;
    double k = glm::length(tan2 - tan1) / double(num_samples) * sign;
    for (size_t i = 0; i < num_samples; ++i) {
      double t = double(i) / double(num_samples - 1);
      result.push_back({(1. - t) * pt1 + t * pt2, tan, k, false, true, 0.0});
    }
  } else {
    result.push_back(
      {xsec.avg_point(), xsec.avg_tangent(), 0., false, false, 0.0});
    auto &sample = result.back();

    double change_mags = 0.;
    int num_changes = 0;
    for (auto &p : xsec.points) {
      // Don't add tangent change at endpoints
      if (p.i < d_i ||
          p.i >= cluster.strokes[p.stroke_idx].points.size() - 1 - d_i)
        continue;

      glm::dvec2 ortho = normal(p.tangent);

      glm::dvec2 prev_pt(
        point(cluster.strokes[p.stroke_idx].points, p.i - d_i));
      glm::dvec2 pt(p.point);
      glm::dvec2 next_pt(
        point(cluster.strokes[p.stroke_idx].points, p.i + d_i));

      glm::dvec2 prev_tan(
        tangent(cluster.strokes[p.stroke_idx].points, p.i - d_i));
      glm::dvec2 tan(p.tangent);
      glm::dvec2 next_tan(
        tangent(cluster.strokes[p.stroke_idx].points, p.i + d_i));

      double k1 = glm::dot(ortho, tan - prev_tan) /
                  std::max(1e-4, glm::distance(pt, prev_pt));
      double k2 = glm::dot(ortho, next_tan - tan) /
                  std::max(1e-4, glm::distance(pt, next_pt));
      change_mags += (k1 + k2) / 2.;
      ++num_changes;
    }

    if (num_changes > 0) {
      change_mags /= double(num_changes);
    } else {
      sample.no_k = true;
    }
    double dist;
    if (xsec_idx < xsecs.size() - 1) {
      dist = glm::dot(normal(sample.tangent),
                      xsecs[xsec_idx + 1].avg_point() - sample.point);
    } else {
      dist = glm::dot(normal(sample.tangent),
                      sample.point - xsecs[xsec_idx - 1].avg_point());
    }
    sample.k = std::max(MIN_DIST, dist) * change_mags;

    if (xsec.points.size() > 0) {
      sample.width =
        std::abs(glm::dot(normal(sample.tangent), xsec.points.front().point -
                                                    xsec.points.back().point)) +
        1.0;
    }
  }

  return result;
}
