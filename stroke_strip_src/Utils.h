#pragma once
#include <cmath>

#include "Cluster.h"
#include <deque>
#include <future>
#include <glm/glm.hpp>
#include <gurobi_c++.h>
#include <mutex>
#include <vector>

extern size_t spiral_offset;

struct Intersection {
  double i;
  glm::dvec2 pt;
};
std::vector<Intersection> intersections(const std::vector<glm::dvec2> &polyline,
                                        glm::dvec2 from, glm::dvec2 to);

template <typename T>
T map(T in, T in_from, T in_to, T out_from, T out_to) {
  T mixed = (in - in_from) / (in_to - in_from);
  T res = mixed * (out_to - out_from) + out_from;
  if (res < std::min(out_from, out_to))
    res = std::min(out_from, out_to);
  if (res > std::max(out_from, out_to))
    res = std::max(out_from, out_to);
  return res;
}

double gaussian(double x, double sigma, double mu);

double poly_in_out(double t, double k);

GRBLinExpr l1_norm(GRBModel *model, const std::vector<GRBLinExpr> &x);
GRBQuadExpr l2_norm_sq(GRBModel *model, const std::vector<GRBLinExpr> &x);

struct LinExpr {
  std::map<size_t, double> var_idx_coeff;
  double const_coeff = 0;

  double getValue(const std::vector<GRBVar> &param_vars) const;

  std::map<size_t, double> quad_coeff;
  std::map<std::pair<size_t, size_t>, double> cross_coeff;
  std::map<size_t, double> lin_coeff;

  void expand() {
    std::vector<size_t> indices;
    indices.reserve(var_idx_coeff.size());
    for (const auto &v_coef : var_idx_coeff) {
      quad_coeff[v_coef.first] = v_coef.second * v_coef.second;
      lin_coeff[v_coef.first] = 2 * v_coef.second * const_coeff;
      indices.emplace_back(v_coef.first);
    }

    for (size_t i = 0; i + 1 < indices.size(); ++i) {
      for (size_t j = i + 1; j < indices.size(); ++j) {
        cross_coeff[std::make_pair(indices[i], indices[j])] =
          2 * var_idx_coeff[indices[i]] * var_idx_coeff[indices[j]];
      }
    }
  }
};
GRBLinExpr l1_norm(GRBModel *model, const std::vector<LinExpr> &x,
                   const std::vector<GRBVar> &param_vars,
                   bool add_length = false);
GRBLinExpr l1_norm_constr(GRBModel *model, const std::vector<LinExpr> &x,
                          const std::vector<GRBVar> &param_vars,
                          bool add_length = false);
GRBQuadExpr l2_norm_sq(GRBModel *model, const std::vector<LinExpr> &x,
                       const std::vector<GRBVar> &param_vars);

template <typename RetType>
std::map<int, RetType>
map_clusters(Input &input, std::function<RetType(Cluster &)> process) {
  const int NUM_TASKS = 1;
  if (NUM_TASKS == 1) {
    std::map<int, RetType> vals;
    for (auto &kv : input.clusters) {
      vals[kv.first] = process(kv.second);
    }
    return vals;
  } else {
    std::deque<int> keys;
    std::vector<std::future<void>> futures;
    std::map<int, RetType> vals;
    std::mutex queue_lock;
    std::mutex vals_lock;
    for (auto &kv : input.clusters) {
      keys.push_back(kv.first);
      vals[kv.first] = {};
    }

    for (size_t task = 0; task < NUM_TASKS; ++task) {
      futures.push_back(std::async(std::launch::async, [&]() -> void {
        while (true) {
          int id = -1;
          {
            std::lock_guard<std::mutex> lock(queue_lock);
            if (keys.empty())
              break;
            id = keys.front();
            keys.pop_front();
          }
          auto val = process(input.clusters[id]);
          {
            std::lock_guard<std::mutex> lock(vals_lock);
            vals[id] = val;
          }
        }
      }));
    }

    for (auto &future : futures) {
      future.get();
    }

    return vals;
  }
}

void center_cluster(Cluster &cluster, glm::dvec2 &center);
void decenter_cluster(Cluster &cluster, const glm::dvec2 &center);
void center_input(const Input &input, double &w, double &h, double &x,
                  double &y, glm::dvec2 &center);

bool is_roughly_circular(const Cluster::Stroke &s);
void compute_spiral_angles(Cluster::Stroke &spiral_s);
void cut_spiral_stroke(const Cluster::Stroke &spiral_s,
                       std::vector<Cluster::Stroke> &cut_ss,
                        const double cut_curv = 2 * M_PI / 4);
                       //  const double cut_curv = M_PI);
                      //  const double cut_curv = 2 * M_PI);
void assemble_spiral_stroke(const std::vector<Cluster::Stroke> &cut_ss,
                            std::vector<Cluster::Stroke> &spiral_ss,
                            std::vector<Cluster::XSec> &xsecs);

void connect_stroke(
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    &matching_pair,
  std::vector<Cluster::Stroke> &strokes,
  std::map<size_t, size_t> *connect_map_ptr = nullptr);

void sort_xsec(Cluster::XSec &xsec);
