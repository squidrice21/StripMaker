#pragma once

#include "types.h"

namespace sketching {
enum DecompositionType : std::uint32_t { //
  Connected = 0,
  Probability = 1,
};

enum HardConstraintType : std::uint32_t { //
  JunctionDist = 0,
  StrokeWidth = 1,
};

static const size_t snap_candidate_count = 3;
static const Float junction_cutoff = 0.5;
static const std::vector<Float> region_junc_threshold{0.7, 0.8};
static const Float min_final_prob = 0.05;

extern bool to_bridge;
extern bool to_include_prev_connections;

// Variables
extern bool include_ts;
extern bool include_corners;
extern JunctionType::Type recosider_type;

extern Float corner_projection_ratio;

// Stroke pairing related
static const Float end_alignment_angle_offset = 45.0 / 180 * M_PI;
static const Float stroke_alignment_angle = 25.0 / 180 * M_PI;
static const Float far_end_ratio = 2.0;
static const Float far_stroke_width_ratio = 2.0;
static const Float dist_change_ratio = 2.0;

// State related
static const size_t unlimited_last_steps = 2;
// static const size_t max_num_sub_region_vars = 10;
static const size_t max_num_sub_region_vars =
  std::numeric_limits<size_t>::max(); // For exact decomposition only
static const size_t max_num_sub_region_sols = 100;
extern size_t max_num_states_region;
extern Float min_decomp_prob;
extern DecompositionType decomposition_type;
extern Float largest_gap_ratio;
extern Float largest_positive_region;
extern Float serial_largest_non_region_gap;

// Objective
extern Float region_term_ratio;
extern Float outside_term_ratio;
extern Float max_pos_weight;
extern HardConstraintType hard_constraint_type;
extern Float hard_region_junc_ratio;
static const Float accurate_radius_check_ratio = 3;
extern Float region_linear_slope;

// Baseline parameters
extern Float baseline_accept_threshold;
static const Float baseline_accept_threshold_final = 0.5;

// Feature computation
extern FeatureType prediction_feature_type;
extern bool high_valence_prob_update;
} // namespace sketching
