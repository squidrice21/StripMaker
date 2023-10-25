#pragma once

#include "../sample/SampleGeneration.h"
#include "../stroke_strip_src/Cluster.h"
#include "Features.h"

#include <string>
#include <vector>

struct FilteredSample;

/**
 * @brief Reads sample data from a file and each sample is represented as a pair
 * of strips.
 *
 * @param input The input data.
 * @param example_filename The name of the file containing the sample data.
 * @param clusters The vector of pairs of clusters to store the sample data.
 */
void read_samples(const Input &input, const std::string &example_filename,
                  std::vector<std::pair<Cluster, Cluster>> &clusters);

/**
 * @brief Reads sample data from a file and each sample is represented as a
 * strip.
 *
 * @param input The input data.
 * @param example_filename The name of the file containing the sample data.
 * @param clusters The vector of clusters to store the sample data.
 */
void read_samples(const Input &input, const std::string &example_filename,
                  std::vector<Cluster> &clusters);

/**
 * @brief Reads feature definitions from a file and stores them in a vector of
 * strings.
 *
 * @param fea_filename The name of the file containing the feature definitions.
 * @param fea_csv_names The vector of strings to store the feature definitions.
 */
void read_feature_definitions(const std::string &fea_filename,
                              std::vector<std::string> &fea_csv_names);

/**
 * @brief Writes filtered sample data to a file.
 *
 * @param filtered_example_filename The name of the file to write the filtered
 * sample data to.
 * @param filtered_samples The vector of filtered samples to write to the file.
 */
void write_filtered_samples(
  const std::string &filtered_example_filename,
  const std::vector<FilteredSample> &filtered_samples);

/**
 * @brief Writes sample data to a file.
 *
 * @param example_filename The name of the file to write the sample data to.
 * @param pos_samples The vector of sample data to write to the file.
 */
void write_samples(const std::string &example_filename,
                   const std::vector<Sample> &pos_samples);

/**
 * @brief Writes feature vectors to a CSV file.
 *
 * @param feature_filename The name of the CSV file to write the feature vectors
 * to.
 * @param fea The vector of feature vectors to write to the CSV file.
 */
void write_feature_csv(const std::string &feature_filename,
                       const std::vector<FeatureVector> &fea);

/**
 * @brief Writes probability feature vectors to a CSV file.
 *
 * @param feature_filename The name of the CSV file to write the probability
 * feature vectors to.
 * @param fea The vector of probability feature vectors to write to the CSV
 * file.
 */
void write_probability_feature_csv(const std::string &feature_filename,
                                   const std::vector<FeatureVector> &fea);

/**
 * @brief Writes secondary feature vectors to a CSV file.
 *
 * @param feature_filename The name of the CSV file to write the secondary
 * feature vectors to.
 * @param fea The vector of secondary feature vectors to write to the CSV file.
 */
void write_secondary_feature_csv(const std::string &feature_filename,
                                 const std::vector<FeatureVector> &fea);

/**
 * @brief Fits curves to a clustered input scap and save the fitting curves to a
 * file.
 *
 * @param scap_filename The name of the file containing the clustered data.
 * @param output_filename The name of the file to write the fitting data to.
 * @param fit_width Whether to fit the width of the data.
 * @param vis_folder The name of the folder to write visualization data to.
 * @param cut_filename The name of the file containing cut data.
 * @param fit_single Whether to fit the data as a single cluster.
 * @param to_spiral Whether to fit the data as spirals.
 * @param disable_cut_orientation Whether to disable cut orientation.
 */
void save_fitting(const std::string &scap_filename,
                  const std::string &output_filename,
                  const bool fit_width = false,
                  const std::string &vis_folder = "",
                  const std::string &cut_filename = "",
                  const bool fit_single = false, const bool to_spiral = false,
                  const bool disable_cut_orientation = false);
