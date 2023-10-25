#pragma once

#include "../stroke_strip_src/Cluster.h"
#include "Junction.h"

/**
 * Checks if two strips are separated in the forward direction.
 *
 * @param c1 The first strip to check.
 * @param c2 The second strip to check.
 * @param merged_cluster The merged strip of c1 and c2.
 * @param fits The fitted curves of the strips.
 * @param junctions The vector of junctions.
 * @param decision_junctions The vector of decision junctions.
 * @return True if the strips are separated in the forward direction, false
 * otherwise.
 */
bool is_separate_forward_check(const Cluster &c1, const Cluster &c2,
                               const Cluster &merged_cluster,
                               const std::vector<FittedCurve> &fits,
                               const std::vector<Junction> &junctions,
                               std::vector<Junction> &decision_junctions);
