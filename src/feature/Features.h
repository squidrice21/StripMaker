/**
 * @file Features.h
 * @brief This file contains the declarations of various feature extraction
 * functions and classes.
 */
#pragma once

#include "../sample/SampleGeneration.h"
#include "../stroke_strip_src/Cluster.h"

#include <memory>
#include <string>
#include <vector>

/**
 * @brief Calculates a hash value (used for cache files and debug visualization)
 * for the given stroke indices.
 */
std::string example_hash(const std::vector<size_t> &stroke_indices);
std::string example_hash(const std::vector<int> &stroke_indices);
std::string example_hash(const std::vector<size_t> &stroke_indices1,
                         const std::vector<size_t> &stroke_indices2);
std::string get_hash_str(const Cluster &cluster);
std::string get_hash_str_comb(const Cluster &cluster1, const Cluster &cluster2);

struct FeatureVector;

/**
 * @brief Abstract base class for feature extraction.
 *
 * This class defines an interface for feature extraction. It provides two
 * virtual functions:
 * - description(): returns a string describing the feature.
 * - csv(): returns the feature name used in CSV.
 *
 * The class also has a protected member variable num_binned_ which specifies
 * the number of bins used for binning the feature values.
 */
struct Feature {
  /**
   * @brief Returns a string describing the feature.
   *
   * @return std::string
   */
  virtual std::string description() const = 0;

  /**
   * @brief Returns a string representation of the feature in CSV format.
   *
   * @return std::string
   */
  virtual std::string csv() const { return "undefined"; }

  /**
   * @brief Default destructor.
   *
   */
  virtual ~Feature() = default;

  /**
   * @brief Number of bins used for binning the feature values.
   *
   */
  size_t num_binned_;

protected:
  /**
   * @brief Constructor.
   *
   * @param num_binned Number of bins used for binning the feature values.
   */
  Feature(size_t num_binned = 1)
    : num_binned_(num_binned) {}
};

// Local feature
struct StrokeClusterFeature : Feature {
  virtual ~StrokeClusterFeature() = default;

  virtual std::vector<double> operator()(const Cluster &merged_cluster,
                                         const FittedCurve &fit1,
                                         const FittedCurve &fit2) const = 0;

protected:
  StrokeClusterFeature(size_t num_binned = 1)
    : Feature(num_binned) {}
};

// Global feature
struct SecondaryClusterFeature : Feature {
  virtual ~SecondaryClusterFeature() = default;

  virtual std::vector<std::string> csvs() const {
    return std::vector<std::string>{"undefined"};
  }
  virtual std::vector<double>
  operator()(const FeatureVector &base_feature,
             const std::vector<FeatureVector> &features) const = 0;

protected:
  SecondaryClusterFeature(size_t num_binned = 1)
    : Feature(num_binned) {}
};

struct TypicalClusterFeature : Feature {
  virtual ~TypicalClusterFeature() = default;

  virtual std::vector<double> operator()(const Cluster &cluster) const = 0;

protected:
  TypicalClusterFeature(size_t num_binned = 1)
    : Feature(num_binned) {}
};

/**
 * @brief These functions define the features that are computed. These computed
 * features are later picked by the prediction code based on predefined or given
 * feature definitions.
 */
std::vector<std::unique_ptr<StrokeClusterFeature>>
get_stroke_cluster_features();
std::vector<std::unique_ptr<StrokeClusterFeature>>
get_stroke_cluster_probability_features();
std::vector<std::unique_ptr<SecondaryClusterFeature>>
get_secondary_cluster_features();
std::vector<std::unique_ptr<TypicalClusterFeature>>
get_typical_cluster_features();

/**
 * @brief Struct providing feature computation functions and containing the
 * computed feature values.
 */
struct FeatureVector {
  enum Type : std::uint8_t {
    StrokeCluster = 0,
    ClusterCluster = 1,
  };
  FeatureVector()
    : width_(-1)
    , height_(-1) {}
  FeatureVector(int width, int height)
    : width_(width)
    , height_(height) {}
  void compute(Cluster &, Cluster &, double, Cluster &, bool test_time,
               const std::vector<
                 std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
                 *matching_pair_ptr,
               bool reorient = false, std::string vis_output = "");
  void
  compute_angle(Cluster &, Cluster &, double, Cluster &,
                const std::vector<
                  std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
                  *matching_pair_ptr,
                bool reorient = false, std::string vis_output = "");
  void compute_distance(
    Cluster &cluster1, Cluster &cluster2, double thickness,
    Cluster &merged_cluster,
    const std::vector<std::pair<std::pair<size_t, bool>,
                                std::pair<size_t, bool>>> *matching_pair_ptr,
    bool reorient, std::string vis_output = "");
  void compute_probability(
    Cluster &, Cluster &, double, Cluster &,
    const std::vector<std::pair<std::pair<size_t, bool>,
                                std::pair<size_t, bool>>> *matching_pair_ptr,
    bool reorient = false, std::string vis_output = "");
  void compute_secondary(
    Cluster &, Cluster &, double, Cluster &,
    const std::map<size_t, FeatureVector> &cluster_features,
    const std::vector<std::pair<std::pair<size_t, bool>,
                                std::pair<size_t, bool>>> *matching_pair_ptr,
    bool reorient = false, std::string vis_output = "");
  void
  compute_typical(Cluster &, double,
                  const std::vector<
                    std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
                    *matching_pair_ptr,
                  bool reorient = false, std::string vis_output = "");

  struct ClusterInfo {
    size_t cluster_idx;
    std::vector<size_t> stroke_indices;
  };

  int width_, height_;
  std::vector<double> features_;
  std::vector<std::string> descriptions_;
  std::vector<ClusterInfo> cluster_info_vec_;
  std::map<std::string, double> obj_term_values_;
  std::map<std::string, std::string> vis_files_;
};

/**
 * Fills the inf values with the calculated features from the provided
 * feature calculators.
 *
 * @param feature_calculators A vector of unique pointers to
 * StrokeClusterFeature objects used to calculate features.
 * @param descriptions A vector of strings to store the descriptions of each
 * calculated feature.
 * @param features A vector of doubles to store the calculated features.
 */
void fill_inf(
  const std::vector<std::unique_ptr<StrokeClusterFeature>> &feature_calculators,
  std::vector<std::string> &descriptions, std::vector<double> &features);
/**
 * Fills the inf values with the information contained in
 * the SecondaryClusterFeature vector.
 *
 * @param SecondaryClusterFeature A vector of unique pointers to
 * SecondaryClusterFeature objects.
 * @param descriptions A vector of strings to be filled with the descriptions of
 * the features.
 * @param features A vector of doubles to be filled with the values of the
 * features.
 */
void fill_inf(const std::vector<std::unique_ptr<SecondaryClusterFeature>>
                &SecondaryClusterFeature,
              std::vector<std::string> &descriptions,
              std::vector<double> &features);

/**
 * Computes the feature vectors for a given filename and set of samples.
 *
 * @param filename The name of the file to compute features for.
 * @param samples The set of samples to use for computing features.
 * @param matching_pair_ptr A pointer to a vector of matching pairs, where each
 * pair is represented as a pair of pairs. The first pair represents the indices
 * of the matching samples, and the second pair represents whether the samples
 * are from the same sketch or not. If this pointer is null, then all pairs of
 * samples will be considered as matching pairs.
 * @param test_time Whether or not to output timing information for the feature
 * computation.
 * @param output_vis_folder The folder to output visualization files to. If this
 * is an empty string, then no visualization files will be output.
 * @param cache_folder The folder to use for caching intermediate results. If
 * this is an empty string, then no caching will be used.
 *
 * @return A vector of feature vectors, where each feature vector corresponds to
 * a sample in the input set.
 */
std::vector<FeatureVector> compute_features(
  const std::string &filename, const std::vector<Sample> &samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool test_time, const std::string &output_vis_folder = "",
  const std::string &cache_folder = "");

/**
 * Computes the feature vectors for the given input and samples.
 *
 * @param input The input to compute features for.
 * @param samples The samples to use for computing features.
 * @param matching_pair_ptr A pointer to the matching pairs to use for computing
 * features.
 * @param test_time Whether or not this is a test run.
 * @param output_vis_folder The folder to output visualizations to.
 * @param cache_folder The folder to cache results in.
 * @return A vector of feature vectors.
 */
std::vector<FeatureVector> compute_features(
  const Input &input, const std::vector<Sample> &samples,
  const std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>
    *matching_pair_ptr,
  bool test_time, const std::string &output_vis_folder = "",
  const std::string &cache_folder = "");
