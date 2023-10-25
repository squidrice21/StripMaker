#!/usr/bin/env python

"""
Run data generation, training and cross validation for the local classifier.
"""

import argparse
import os
import shutil
import subprocess
import time
import yaml

import pandas as pd

from learning.cross_validate import cv_from_data
from learning.learning import learn_from_data
from learning.feature_set import load_feature_set


def generate_data(path, sample_exe, feature_exe):
    data_list = []

    # Generate data
    os.makedirs(path.replace('.scap', ''), exist_ok=True)

    # Generate samples
    output_file = path.replace('.scap', '') + '/samples'
    cache_folder = os.path.abspath(os.path.dirname(output_file) + '/cache')
    os.makedirs(cache_folder, exist_ok=True)
    cmd = [sample_exe, '-i', path, '-o', output_file,
           '-c', cache_folder]
    # Use cut point as parameterization hint
    cmd += ['--uncut']

    if not os.path.exists(output_file + '_filtered.csv'):
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)

    # Compute features based on the samples
    sample_input_files = [
        output_file + '_pos_stroke.csv',
        output_file + '_neg_stroke.csv',
        output_file + '_pos_cluster.csv',
        output_file + '_neg_cluster.csv',
    ]
    feature_output_files = []
    for sample_f in sample_input_files:
        output_file = sample_f.replace('/samples', '/features')

        if not os.path.exists(sample_f):
            continue
        feature_output_files.append(output_file)
        cache_folder = os.path.abspath(os.path.dirname(output_file) + '/cache')
        os.makedirs(cache_folder, exist_ok=True)
        cmd = [feature_exe, '-i', path, '-e', sample_f, '-o', output_file,
               '-v', os.path.abspath(os.path.dirname(output_file)),
               '-c', cache_folder]
        cmd += ['--uncut']
        if not os.path.exists(output_file):
            subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True)
        # Remove cache files
        if os.path.exists(cache_folder):
            shutil.rmtree(cache_folder)

    # Merge the feature files
    df = pd.DataFrame()
    for p in feature_output_files:
        add_df = pd.read_csv(p)
        df = pd.concat([df, add_df], ignore_index=True)
    df.drop_duplicates(
        subset=['stroke_indices1', 'stroke_indices2'], keep='first', inplace=True)

    output_file = path.replace('.scap', '') + '/features_full.csv'
    df.to_csv(output_file, index=False)

    # Add results to the list
    data_list.append(output_file)

    return data_list


def train_data(data_list, model_folder,
               model_save=None, feature_set=[],
               random_state=42):
    start = time.time()
    # Learn
    model_file = model_folder + '/model.sav'
    if model_save:
        model_file = model_save
    sc_file = os.path.dirname(model_file) + '/sc_' + \
        os.path.basename(model_file)
    print(f'data_list: {data_list}')
    print(f'model_save_file={model_file}, sc_save_file={sc_file}')
    learn_from_data(data_list, model_save_file=model_file,
                    sc_save_file=sc_file, feature_set=feature_set,
                    random_state=random_state)
    end = time.time()
    # print(f'{model_folder} time: {end - start} s')

    return model_file


def cv_data(data_list, scap_folder,
            model_folder, model_save=None, feature_set=[],
            random_state=42):
    # CV
    model_file = model_folder + '/model.sav'
    if model_save:
        model_file = model_save
    else:
        os.makedirs(model_folder, exist_ok=True)
    sc_file = os.path.dirname(model_file) + '/sc_' + \
        os.path.basename(model_file)
    print(f'data_list: {data_list}')
    print(f'cv_model_save_file={model_file}, cv_sc_save_file={sc_file}')
    cv_from_data(data_list, scap_folder, model_save_file=model_file,
                 sc_save_file=sc_file, feature_set=feature_set,
                 random_state=random_state)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='Input file')
    parser.add_argument('--cv', action='store_true',
                        help='Run cross-validation')
    parser.add_argument('-m', '--model', help='Input model file')
    parser.add_argument('-e', '--exe', help='data generation exe')
    parser.add_argument('-s', '--snapshot', help='snapshot path')
    parser.add_argument('-c', '--scap', help='scap file path')
    parser.add_argument(
        '-f', '--fea', help='YAML file of the feature definitions')
    parser.add_argument('--seed', type=int, default=42,
                        help='Seed number for random forest')
    parser.add_argument('--debug', action='store_true',
                        help='place a breakpoint before execution')
    args = parser.parse_args()

    if args.debug:
        # Give time to attach a debugger.
        breakpoint()

    if args.snapshot:
        os.makedirs(args.snapshot, exist_ok=True)
    else:
        os.makedirs("snapshot/", exist_ok=True)
        args.snapshot = 'snapshot'

    if args.fea:
        load_feature_set(args.fea)

    from learning.feature_set import feature_set

    # Generate data and learn
    if '.scap' in args.input and args.exe:
        sample_exe, feature_exe = args.exe.split(',')
        data_list = generate_data(args.input, sample_exe, feature_exe)
    elif '.yml' in args.input:
        with open(args.input) as f:
            data_list = yaml.safe_load(f)

        if not args.cv:
            model_file = train_data(
                data_list, args.snapshot,
                model_save=args.model, feature_set=feature_set,
                random_state=args.seed)
            print(f'Model: {model_file}')
        else:
            assert args.scap
            cv_data(data_list, args.scap, args.snapshot,
                    model_save=args.model, feature_set=feature_set,
                    random_state=args.seed)
    else:
        assert False, 'Unknown input file type'


if __name__ == '__main__':
    main()
