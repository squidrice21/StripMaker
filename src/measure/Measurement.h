#pragma once

#include "../../stroke_strip_src/Cluster.h"
#include "../../stroke_strip_src/Fitting.h"
#include "glm/detail/type_vec.hpp"

#include <map>
#include <vector>

/**
 * @brief Calculates the total signed curvature of a stroke.
 *
 * @param s The stroke to calculate the total signed curvature for.
 * @return The total signed curvature of the stroke.
 */
double stroke_total_signed_curvature(const Cluster::Stroke &s);

/**
 * @brief Calculates the maximum gap between two consecutive points if the input
 * stroke were a spiral stroke.
 *
 * @param s The spiral stroke to calculate the maximum gap for.
 * @return The maximum gap between two consecutive points on the spiral stroke.
 */
double max_spiral_gap(const Cluster::Stroke &s);

/**
 * @brief Calculates the minimum gap between two consecutive points if the input
 * stroke were a spiral stroke.
 *
 * @param s The spiral stroke to calculate the minimum gap for.
 * @return The minimum gap between two consecutive points on the spiral stroke.
 */
double min_spiral_gap(const Cluster::Stroke &s);

/**
 * @brief Calculates the length of the closed fitting curve of the stroke.
 *
 * @param s The spiral stroke to calculate the closed length for.
 * @return The closed length of the spiral stroke.
 */
double spiral_closed_length(const Cluster::Stroke &s);

/**
 * @brief Calculates the mass centers of a spiral stroke.
 *
 * @param s The spiral stroke to calculate the mass centers for.
 * @return A vector of mass centers of the spiral stroke.
 */
std::vector<glm::dvec2> spiral_mass_centers(const Cluster::Stroke &s);

/**
 * @brief Cuts the input stroke with trials of several different turning angle
 * thresholds.
 *
 * @param s The spiral stroke to set the cut angle for.
 * @param spiral_cluster_ptr A pointer to the spiral cluster.
 */
void set_spiral_cut_angle(Cluster::Stroke &s,
                          Cluster *spiral_cluster_ptr = nullptr);

/**
 * @brief Determines if a stroke is a spiral stroke (thus needs to be pre-cut
 * for parameterization and fitting).
 *
 * @param s The stroke to detect if it is a spiral stroke.
 */
void detect_spiral_stroke(Cluster::Stroke &s);
