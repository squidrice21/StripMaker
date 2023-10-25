#pragma once

#include <cmath>

#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>

#include "../stroke_strip_src/Fitting.h"
#include "../stroke_strip_src/SketchInfo.h"
#include <glm/glm.hpp>

/**
 * @brief A class that provides methods to convert a Capture object (the input
 * type used by SketchAggregator) to an Input object (the input type used by
 * StrokeStrip).
 *
 */
class Capture2Input {
public:
  /**
   * @brief Convert a Capture object to an Input object.
   *
   * @param capture The Capture object to convert.
   * @param recenter A boolean indicating whether to recenter the sketch.
   * @return Input The resulting Input object.
   */
  Input from_capture(Capture capture, bool recenter = true);

  /**
   * @brief Convert a Capture object to an Input object for set operations.
   *
   * @param capture The Capture object to convert.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param ct The center point of the sketch.
   * @param stroke_num The number of strokes in the sketch.
   * @param num_cluster The number of clusters in the sketch.
   * @return Input The resulting Input object.
   */
  Input from_capture4set(Capture capture, std::map<int, bool> &rev_list,
                         glm::dvec2 &ct, int stroke_num, int num_cluster);

  /**
   * @brief Convert a Capture object to an Input object for clustering.
   *
   * @param capture The Capture object to convert.
   * @param base_num The base number.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param target_num The target number.
   * @param flip A boolean indicating whether to flip the sketch.
   * @return Input The resulting Input object.
   */
  Input from_capture4cluster(Capture capture, int base_num,
                             std::map<int, bool> &rev_list, int target_num,
                             bool flip);

  /**
   * @brief Convert a Capture object to an Input object for orientation.
   *
   * @param capture The Capture object to convert.
   * @param fits A map of stroke indices to FittedCurve objects.
   * @param ct The center point of the sketch.
   * @param stroke_num The number of strokes in the sketch.
   * @return Input The resulting Input object.
   */
  Input from_capture4orientation(Capture capture,
                                 std::map<int, FittedCurve> fits, glm::dvec2 ct,
                                 int stroke_num);

  /**
   * @brief Convert a Capture object to an Input object for input.
   *
   * @param capture The Capture object to convert.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param crossing_cluster A map of cluster indices to boolean values
   * indicating whether the cluster crosses the centerline.
   * @param rev_clusters A map of cluster indices to boolean values indicating
   * whether the cluster is reversed.
   * @param stroke_num The number of strokes in the sketch.
   * @param num_cluster The number of clusters in the sketch.
   * @return Input The resulting Input object.
   */
  Input from_capture4input(Capture capture, std::map<int, bool> &rev_list,
                           std::map<int, bool> crossing_cluster,
                           std::map<int, bool> rev_clusters, int stroke_num,
                           int num_cluster);

  /**
   * @brief Add a stroke to an Input object.
   *
   * @param capture The Capture object to convert.
   * @param input The Input object to add the stroke to.
   * @param stroke_num The stroke number.
   * @param group_num The group number.
   */
  void add_stroke(Capture capture, Input &input, int stroke_num, int group_num);

  /**
   * @brief Convert a Capture object to an Input object for two clusters.
   *
   * @param capture The Capture object to convert.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param base_num The base number.
   * @param base_stroke The base stroke.
   * @param upper_ind The upper index.
   * @return Input The resulting Input object.
   */
  Input from_capture4two_cluster(Capture capture, std::map<int, bool> &rev_list,
                                 int base_num, int base_stroke, int upper_ind);

  /**
   * @brief Convert a Capture object to an Input object for set merging.
   *
   * @param capture The Capture object to convert.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param base_num The base number.
   * @param base_stroke The base stroke.
   * @param upper_ind The upper index.
   * @return Input The resulting Input object.
   */
  Input from_capture4set_merge(Capture capture, std::map<int, bool> &rev_list,
                               int base_num, int base_stroke, int upper_ind);

  /**
   * @brief Convert a Capture object to an Input object for negative set
   * merging.
   *
   * @param capture The Capture object to convert.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param base_num The base number.
   * @param target_num The target number.
   * @return Input The resulting Input object.
   */
  Input from_capture4set_merge_negative(Capture capture,
                                        std::map<int, bool> &rev_list,
                                        int base_num, int target_num);

  /**
   * @brief Convert a Capture object to an Input object for negative two
   * clusters.
   *
   * @param capture The Capture object to convert.
   * @param rev_list A map of stroke indices to boolean values indicating
   * whether the stroke is reversed.
   * @param base_num The base number.
   * @param target_num The target number.
   * @return Input The resulting Input object.
   */
  Input from_capture4two_cluster_negative(Capture capture,
                                          std::map<int, bool> &rev_list,
                                          int base_num, int target_num);
};
