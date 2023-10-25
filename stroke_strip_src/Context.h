#pragma once

#include <gurobi_c++.h>
#include <oneapi/tbb/global_control.h>

extern int model_number;

struct Context {
  Context(bool grb_log = false);
  void optimize_model(GRBModel *model) const;

  GRBEnv grb;
  // bool debug_viz = false;
  bool debug_viz = true;
  bool ignore_endpoint = false;
  // bool cut = false;
  bool cut = false;
  bool rainbow = false;
  bool widths = false;
  bool taper_widths = false;
  // bool taper_widths = true;

  bool tighter_fit = false;
  oneapi::tbb::global_control MAXTHREADS;

  // Fitting
  bool to_spiral = false;

  // Timer
  double orientation_timer = 0;
  double parameterization_timer = 0;
  double fitting_timer = 0;
};
