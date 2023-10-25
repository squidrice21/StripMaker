#!/usr/bin/env python

"""
Utilities for feature renaming, timing and pandas data frame operations.
"""

import os
import re
import sys
import time
import yaml

import pandas as pd

import _sketching as s

feature_rename_mapping = {
    'mean_distance': 'median_distance',
    'avg_width_to_fitting_ratio_cluster1': 'avg_narrowness_cluster1',
    'median_width_to_fitting_ratio_cluster1': 'median_narrowness_cluster1',
    'avg_width_to_fitting_ratio_cluster2': 'avg_narrowness_cluster2',
    'median_width_to_fitting_ratio_cluster2': 'median_narrowness_cluster2',
    'avg_width_to_fitting_ratio_overall': 'avg_narrowness_combined',
    'median_width_to_fitting_ratio_overall': 'median_narrowness_combined',
    'max_avg_width_to_fitting_ratio_cluster': 'maxC12_avg_narrowness',
    'max_median_width_to_fitting_ratio_cluster': 'maxC12_median_narrowness',
    'min_avg_width_to_fitting_ratio_cluster': 'minC12_avg_narrowness',
    'min_median_width_to_fitting_ratio_cluster': 'minC12_median_narrowness',
    'average_cwise_distance_to_overall_xsection': 'avg_distance2width_combined',
    'median_cwise_distance_to_overall_xsection': 'median_distance2width_combined',
    'average_cwise_distance_to_max_xsection': 'avg_distance2LocalMaxWidth',
    'median_cwise_distance_to_max_xsection': 'median_distance2LocalMaxWidth',
    'average_cwise_distance_to_min_xsection': 'avg_distance2LocalMinWidth',
    'median_cwise_distance_to_min_xsection': 'median_distance2LocalMinWidth',
    'max_cluster_gap1_to_cwise_distance': 'max_LocalMaxGap2width_cluster1',
    'max_cluster_gap2_to_cwise_distance': 'max_LocalMaxGap2width_cluster2',
    'min_cluster_gap1_to_cwise_distance': 'min_LocalMaxGap2width_cluster1',
    'min_cluster_gap2_to_cwise_distance': 'min_LocalMaxGap2width_cluster2',
    'first_overlapping_xsection_change_ratio': 'first_local_nonoverlapping2overlapping',
    'last_overlapping_xsection_change_ratio': 'last_local_nonoverlapping2overlapping',
    'overlapping_ratio': 'overlapping_length2combined_length',
    'min_overlapping_xsection_change_ratio': 'minC12_local_nonoverlapping2overlapping',
    'max_overlapping_xsection_change_ratio': 'maxC12_local_nonoverlapping2overlapping',
    'min_avg_overlapping_xsection_change_ratio': 'maxC12_AvgNonoverlapping2AvgOverlapping',
    'max_avg_overlapping_xsection_change_ratio': 'minC12_AvgNonoverlapping2AvgOverlapping',
    'max_p90th_cluster_gap_to_cwise_distance': 'maxC12_90th_LocalMedianGap2width',
    'min_p90th_cluster_gap_to_cwise_distance': 'minC12_90th_LocalMedianGap2width',
    'max_p10th_cluster_gap_to_cwise_distance': 'maxC12_10th_LocalMedianGap2width',
    'min_p10th_cluster_gap_to_cwise_distance': 'minC12_10th_LocalMedianGap2width',
    'max_min_cluster_gap_to_cwise_distance':  'maxC12_min_LocalMedianGap2width',
    'max_max_cluster_gap_to_cwise_distance':  'maxC12_max_LocalMedianGap2width',
    'min_min_cluster_gap_to_cwise_distance':  'minC12_min_LocalMedianGap2width',
    'min_max_cluster_gap_to_cwise_distance':  'minC12_max_LocalMedianGap2width',
}


def safe_drop(data, col_name):
    try:
        data.drop(col_name, axis=1, inplace=True)
    except KeyError:
        print(f'Missing column: {col_name}')


def add_table_to_df(table, df=[]):
    header_table = {col[0].replace('| ', ''): [col[1]] for col in table}
    if len(df) == 0:
        df = pd.DataFrame(columns=[k for k, v in header_table.items()])
    new_df = pd.DataFrame.from_dict(header_table)
    df = pd.concat([df, new_df], ignore_index=True)
    return df


class Timing:
    def __init__(self, msg=''):
        self.msg = msg

    def __enter__(self):
        if self.msg:
            print(self.msg, end='... ', file=sys.stderr, flush=True)
        self.start = time.perf_counter_ns()

    def __exit__(self, _type, _value, _traceback):
        elapsed = 1e-6 * (time.perf_counter_ns() - self.start)
        if elapsed < 10000:
            print(f'{int(elapsed)} ms', file=sys.stderr)
        else:
            print(f'{int(elapsed / 1000)} s', file=sys.stderr)


def filter_df(df, param_failure_yml):
    if df.shape[0] == 0:
        return df

    # Get info from df
    input_str = str(df['input'].values[0])
    if os.sep in input_str:
        id_list = [f for f in df['input'].apply(
            lambda x: x.split(os.sep)[-2]).to_list() if f != 'input']
        if len(id_list) > 0:
            df['id'] = df['input'].apply(lambda x: str(x).split(os.sep)[-2])
        else:
            df['id'] = df['input'].apply(lambda x: str(x).split(os.sep)[-3])
    else:
        df['id'] = df['input']
    # In case that the cluster 1 and cluster 2 are swapped
    df['id1'] = df.apply(lambda x: x['id'] + '/' + ('pos_' if x['ans'] else 'neg_') + s.example_hash([int(ind)
                         for ind in str(x['stroke_indices1']).split(' ')],
        [int(ind) for ind in str(x['stroke_indices2']).split(' ')]), axis=1)
    df['id2'] = df.apply(lambda x: x['id'] + '/' + ('pos_' if x['ans'] else 'neg_') + s.example_hash([int(ind)
                         for ind in str(x['stroke_indices2']).split(' ')],
        [int(ind) for ind in str(x['stroke_indices1']).split(' ')]), axis=1)

    # df.to_csv('filter_df.csv', index=False)

    # Load failure cases
    with open(param_failure_yml) as f:
        failure_list = yaml.safe_load(f)['failure_examples']
    failure_list = [re.sub(f'_[0-9]+_[0-9]+_', '_', hash_str)
                    for hash_str in failure_list]
    # print(failure_list)
    # df = df[df['id1'].isin(failure_list) | df['id2'].isin(failure_list)]
    df.drop(df[df['id1'].isin(failure_list) |
            df['id2'].isin(failure_list)].index, inplace=True)

    # df.to_csv('filter_df_after.csv', index=False)
    # exit()

    safe_drop(df, 'id')
    safe_drop(df, 'id1')
    safe_drop(df, 'id2')

    return df


def rename_features(df):
    rename_map = {}
    for col in df.columns:
        match_text = [(col, k) for k, _ in feature_rename_mapping.items()
                      if k == re.sub(r'_b[0-9]*$', '', col)]
        for (col, k) in match_text:
            rename_map.update({col: col.replace(k, feature_rename_mapping[k])})

    if len(rename_map) > 0:
        df = df.rename(columns=rename_map)

    return df, len(rename_map) > 0


def rename_feature_set(feature_set):
    rename_map = {}
    for col in feature_set:
        match_text = [(col, k) for k, _ in feature_rename_mapping.items()
                      if k == re.sub(r'_b[0-9]*$', '', col) and 'global_' not in col]
        for (col, k) in match_text:
            rename_map.update({col: col.replace(k, feature_rename_mapping[k])})

    renamed_feature_set = []
    for f in feature_set:
        if f in rename_map:
            renamed_feature_set.append(rename_map[f])
        else:
            renamed_feature_set.append(f)

    return renamed_feature_set
