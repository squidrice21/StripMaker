#!/usr/bin/env python

"""
Train RandomForestClassifier on generated data.
"""

import argparse
import glob
import os
import pickle
import re
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils.class_weight import compute_class_weight

from utils.util import (safe_drop, filter_df,
                        rename_features, rename_feature_set)
try:
    from learning.measure import (print_stats, print_forest_stats,
                                  print_feature_importance, cutoff)
except ModuleNotFoundError:
    from measure import (print_stats, print_forest_stats,
                         print_feature_importance, cutoff)

FAILURES_YML = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/param_failures.yml')


def preprocess(path, to_filter=False, in_feature_set=[]):
    if len(in_feature_set) == 0:
        from learning.feature_set import feature_set, typical_feature_set, secondary_feature_set

    if 'csv' in path:
        in_data = pd.read_csv(path)
    else:
        in_data = pd.read_excel(path)

    is_typical = 'stroke_indices2' not in in_data.columns
    is_global = len([col for col in in_data.columns if 'global_' in col]) > 0
    if is_typical and len(in_feature_set) == 0:
        feature_set = typical_feature_set
    if is_global and len(in_feature_set) == 0:
        feature_set = secondary_feature_set
    if len(in_feature_set) > 0:
        feature_set = in_feature_set

    if 'ans' not in in_data.columns and 'label' in in_data.columns:
        label = in_data['label'].to_list()
        in_data.insert(in_data.columns.get_loc("label") + 1, "ans", label)
        safe_drop(in_data, 'label')
    if 'prod' in in_data.columns:
        in_data['pred'] = in_data['prod']
        safe_drop(in_data, 'prod')

    if not is_typical:
        # Rename the features. If any of them is changed, overwrite the csv file.
        in_data, to_save = rename_features(in_data)
        if to_save:
            in_data.to_csv(path, index=False)
        renamed_feature_set = rename_feature_set(feature_set)
        # print(renamed_feature_set)

        # Filter out cases with known parameterization failures
        if to_filter:
            in_data = filter_df(in_data, FAILURES_YML)
    else:
        renamed_feature_set = feature_set

    in_data.replace([np.inf, -np.inf], np.nan, inplace=True)
    if in_data[in_data.isnull().any(axis=1)].shape[0] > 0:
        debug_nan_csv = './nan.csv'
        nan_df = in_data[in_data.isnull().any(axis=1)].copy()
        if os.path.exists(debug_nan_csv):
            df = pd.read_csv(path)
            nan_df = pd.concat([df, nan_df], ignore_index=True)
        nan_df.to_csv(debug_nan_csv, index=False)
        del nan_df

    data = pd.DataFrame()
    # Pick features
    if 'pred' in in_data.columns:
        data['pred'] = in_data['pred']
    data['ans'] = in_data['ans']

    # Order based on renamed_feature_set
    in_fea = []
    for col in in_data.columns:
        for fi, fea in enumerate(renamed_feature_set):
            if fea == re.sub(r'_b[0-9]*$', '_b', col):
                in_fea.append((fi, col))
                break
    in_fea = sorted(in_fea)
    in_fea = [f[1] for f in in_fea]
    for fea in in_fea:
        data[fea] = in_data[fea]
    # print(data.columns)

    # Drop infinity (implying no overlapping)
    data.dropna(how='any', inplace=True)
    # data.replace(np.inf, 0.0, inplace=True)

    full_data = in_data.loc[data.index, :]
    return data, full_data


def read_data(path_list: list, folder=None, to_filter=False, feature_set=[]):
    df, merged_df = preprocess(path_list[0], to_filter, feature_set)
    for i in range(len(path_list)-1):
        add_df, add_merged_df = preprocess(
            path_list[i + 1], to_filter, feature_set)
        df = pd.concat([df, add_df], ignore_index=True)
        merged_df = pd.concat([merged_df, add_merged_df], ignore_index=True)
    if folder:
        merged_df.to_csv(folder + '/merged_df.csv', index=False)
        df.to_csv(folder + '/correct_df.csv', index=False)

    return df, merged_df


def process4learning(df, to_transform=True):
    Y = df['ans'].values
    safe_drop(df, 'ans')
    safe_drop(df, 'pred')
    safe_drop(df, 'path')
    X = df.values

    X = np.array(X, dtype=float)
    Y = np.array(Y, dtype=float)
    sc = None

    descriptions = df.columns

    return X, Y, sc, descriptions


def balanced_class_weights(y):
    (w1, w2) = compute_class_weight(
        'balanced', classes=np.array([0, 1]), y=y[y < 2])
    return {0: w1, 1: w2}


def learn(X, Y, random_state=42):
    print(f'random_state={random_state}')
    rfc = RandomForestClassifier(
        max_depth=20, min_samples_leaf=1, n_estimators=150, n_jobs=-1, random_state=random_state)
    rfc.fit(X, Y)

    return rfc


def predict_custom_cutoff(rfc, X):
    prob = rfc.predict_proba(X)
    YY = (rfc.predict_proba(X)[:, 1] >= cutoff).astype(
        bool)
    return YY


def learn_from_data(path_list: list,
                    model_save_file=None, sc_save_file=None, feature_set=[],
                    random_state=42):
    df, merged_df = read_data(path_list, folder=None,
                              feature_set=feature_set, to_filter=True)
    to_transform = True
    X, Y, sc, descriptions = process4learning(df, to_transform=to_transform)

    start = time.time()
    rfc = learn(X, Y, random_state=random_state)

    end = time.time()
    print(f'Learning time: {end - start} s')

    # Measurements
    YY = predict_custom_cutoff(rfc, X)
    prob = rfc.predict_proba(X)
    if prob.shape[1] > 1:
        prob = prob[:, 1]
    print_stats('Overfitting', Y, YY)
    print('')
    print_forest_stats(rfc)
    print('')
    print_feature_importance(rfc, descriptions, X, Y)

    if model_save_file:
        pickle.dump(rfc, open(model_save_file, 'wb'))
    if sc_save_file:
        pickle.dump(sc, open(sc_save_file, 'wb'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='Input file(s)')
    parser.add_argument('-o', '--output', help='Output model path')
    parser.add_argument('--debug', action='store_true',
                        help='Place a breakpoint before execution')
    args = parser.parse_args()

    if args.debug:
        breakpoint()

    model_save_file = './model.sav'
    sc_save_file = './sc.sav'
    if args.output:
        model_save_file = args.output
        sc_save_file = os.path.dirname(
            args.output) + '/sc_' + os.path.basename(args.output)

    input_regex = 'examples/**/*.csv'
    if args.input:
        input_regex = args.input
    path_list = glob.glob(input_regex)
    learn_from_data(path_list, model_save_file=model_save_file,
                    sc_save_file=sc_save_file)
