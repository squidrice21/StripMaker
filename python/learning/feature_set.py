#!/usr/bin/env python

"""
Feature definition for the local and global classifier.
"""

import yaml

# Local classifier features
feature_set = [
    'average_angle',
    'median_angle',
    'average_distance',
    'mean_distance',

    'overlapping_ratio',

    'max_avg_width_to_fitting_ratio_cluster',
    'max_median_width_to_fitting_ratio_cluster',
    'min_avg_width_to_fitting_ratio_cluster',
    'min_median_width_to_fitting_ratio_cluster',
    'avg_width_to_fitting_ratio_overall',
    'median_width_to_fitting_ratio_overall',

    'average_cwise_distance_to_max_xsection',
    'median_cwise_distance_to_max_xsection',
    'average_cwise_distance_to_min_xsection',
    'median_cwise_distance_to_min_xsection',

    'average_cwise_distance_to_overall_xsection',
    'median_cwise_distance_to_overall_xsection',

    'max_p90th_cluster_gap_to_cwise_distance',
    'min_p90th_cluster_gap_to_cwise_distance',
    'max_p10th_cluster_gap_to_cwise_distance',
    'min_p10th_cluster_gap_to_cwise_distance',

    'max_overlapping_xsection_change_ratio',
    'min_overlapping_xsection_change_ratio',
    'min_avg_overlapping_xsection_change_ratio',
    'max_avg_overlapping_xsection_change_ratio',

    'velocity',
    'alignment',
]

# Global classifier features
secondary_feature_set = [
    'average_angle',
    'median_angle',
    'average_distance',
    'mean_distance',

    'overlapping_ratio',

    'max_avg_width_to_fitting_ratio_cluster',
    'max_median_width_to_fitting_ratio_cluster',
    'min_avg_width_to_fitting_ratio_cluster',
    'min_median_width_to_fitting_ratio_cluster',
    'avg_width_to_fitting_ratio_overall',
    'median_width_to_fitting_ratio_overall',

    ##
    'average_cwise_distance_to_max_xsection',
    'median_cwise_distance_to_max_xsection',
    'average_cwise_distance_to_min_xsection',
    'median_cwise_distance_to_min_xsection',

    'average_cwise_distance_to_overall_xsection',
    'median_cwise_distance_to_overall_xsection',

    'max_p90th_cluster_gap_to_cwise_distance',
    'min_p90th_cluster_gap_to_cwise_distance',
    'max_p10th_cluster_gap_to_cwise_distance',
    'min_p10th_cluster_gap_to_cwise_distance',

    ##
    'max_overlapping_xsection_change_ratio',
    'min_overlapping_xsection_change_ratio',
    'min_avg_overlapping_xsection_change_ratio',
    'max_avg_overlapping_xsection_change_ratio',

    'velocity',
    'alignment',

    'average_distance2global_average_average_distance',
    'median_distance2global_average_median_distance',
    'average_distance2global_median_average_distance',
    'median_distance2global_median_median_distance',
    'average_distance2global_p90th_average_distance',
    'median_distance2global_p90th_median_distance',
]

typical_feature_set = []


def save_feature_set(yml_file, in_feature_set):
    with open(yml_file, 'w') as f:
        yaml.dump(in_feature_set, f, default_flow_style=False)


def load_feature_set(yml_file):
    global feature_set

    with open(yml_file, 'r') as f:
        feature_set = yaml.safe_load(f)


def load_secondary_feature_set(yml_file):
    global secondary_feature_set

    with open(yml_file, 'r') as f:
        secondary_feature_set = yaml.safe_load(f)
