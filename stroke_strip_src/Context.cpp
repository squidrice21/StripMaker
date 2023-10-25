#include <cmath>

#include "Context.h"

// size_t max_allowed_parallelism = 4;
// size_t max_allowed_parallelism = 8;
size_t max_allowed_parallelism = 1;
int model_number = 0;

namespace {
void optimize_non_convex_model(GRBModel *model) {
  try {
    model->set(GRB_IntParam_NonConvex, 2);
    model->optimize();
  } catch (const GRBException &e2) {
    std::cout << "Error code = " << e2.getErrorCode() << std::endl;
    std::cout << e2.getMessage() << std::endl;
    throw;
  }
}
} // namespace

Context::Context(bool grb_log)
  : grb(true)
  , MAXTHREADS(tbb::global_control::max_allowed_parallelism,
               max_allowed_parallelism) {
  if (!grb_log) {
    grb.set(GRB_IntParam_LogToConsole, 0);
  }
  grb.start();
}

void solve_setup_sampling_2_75(GRBModel *model) {
  model->set(GRB_DoubleParam_IterationLimit, 10000);
  model->set(GRB_DoubleParam_NodeLimit, 10000);
  model->set(GRB_IntParam_SolutionLimit, 5);

  // model->set(GRB_DoubleParam_TimeLimit, 180);
  // model->set(GRB_DoubleParam_TimeLimit, 60);
}

void solve_setup_sampling_1_2_non_convex(GRBModel *model) {
  model->set(GRB_DoubleParam_IterationLimit, 1000);
  model->set(GRB_DoubleParam_NodeLimit, 1000);
  model->set(GRB_IntParam_SolutionLimit, 5);
}

void solve_setup(GRBModel *model) {
  model->set(GRB_DoubleParam_IterationLimit, 1000);
  model->set(GRB_DoubleParam_NodeLimit, 1000);
  model->set(GRB_IntParam_SolutionLimit, 5);
}

void Context::optimize_model(GRBModel *model) const {
  // Timelimit and thread limit
  // model->set(GRB_DoubleParam_TimeLimit, 60);
  solve_setup(model);

  // model->set(GRB_IntParam_Threads, 4);
  model->set(GRB_IntParam_Threads, 1);

  try {
    // model.set(GRB_DoubleParam_FeasibilityTol, 1e-2);
    model->optimize();
    // {
    //   int optimstatus = model->get(GRB_IntAttr_Status);
    //   int iterations = model->get(GRB_DoubleAttr_IterCount);

    //   if (optimstatus != GRB_OPTIMAL) {
    //     std::cout << "@Time-iteration: " << iterations << ", " <<
    //     optimstatus
    //               << "@" << std::endl;
    //   }
    // }
  } catch (GRBException e) {
    try {
      std::cout << "Error code = " << e.getErrorCode() << std::endl;
      // model.set(GRB_IntParam_DualReductions, 0);
      model->set(GRB_IntParam_BarHomogeneous, 1);
      model->set(GRB_DoubleParam_PSDTol, 1e-4);
      // model->set(GRB_IntParam_NonConvex, 2);
      model->optimize();

      // {
      //   int optimstatus = model->get(GRB_IntAttr_Status);
      //   int iterations = model->get(GRB_DoubleAttr_IterCount);

      //   if (optimstatus != GRB_OPTIMAL) {
      //     std::cout << "@Time-iteration 2: " << iterations << "@" <<
      //     std::endl;
      //   }
      // }
    } catch (const GRBException &e2) {
      std::cout << "Error code = " << e2.getErrorCode() << std::endl;
      std::cout << e2.getMessage() << std::endl;
      // We don't throw here since this would return empty result which would
      // be caught and handled properly elsewhere
    }
  }
  if (model->get(GRB_IntAttr_Status) == GRB_NUMERIC) {
    // model.set(GRB_IntParam_DualReductions, 0);
    model->set(GRB_IntParam_BarHomogeneous, 1);
    model->optimize();
  }
  if (model->get(GRB_IntAttr_Status) == GRB_NUMERIC) {
    // model->set(GRB_IntParam_BarHomogeneous, -1);
    model->set(GRB_IntParam_ScaleFlag, 2);
    model->set(GRB_DoubleParam_ObjScale, -0.5);
    model->optimize();
  }

  if (model->get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
    int optimstatus = model->get(GRB_IntAttr_Status);

    if (optimstatus == GRB_INF_OR_UNBD) {
      std::cout << "Model is infeasible or unbounded" << std::endl;
    } else if (optimstatus == GRB_INFEASIBLE) {
      std::cout << "Model is infeasible" << std::endl;
    } else if (optimstatus == GRB_UNBOUNDED) {
      std::cout << "Model is unbounded" << std::endl;
    } else if (optimstatus == GRB_TIME_LIMIT) {
      std::cout << "Solve reached time limit" << std::endl;
    }
    /*else {
        std::cout << "Optimization was stopped with status = " << optimstatus
                  << std::endl;
      }*/
  }
}
