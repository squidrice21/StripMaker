#!/usr/bin/env python

"""
Round robin cross validation and statistics.
"""

import glob
import os
import shutil
import subprocess

from learning.learning import read_data, process4learning, learn, predict_custom_cutoff
from learning.measure import (
    print_stats, print_forest_stats, print_feature_importance, sum_stats)
from utils.util import Timing, add_table_to_df
from launching.run_prediction import temporary_directory

import numpy as np
from tabulate import tabulate

VIS_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/visualize_training.py')
PRED_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/predict.py')

TEMPLATE = r'''
\pdfsuppresswarningpagegroup=1
\documentclass[10pt,letterpaper]{article}
\usepackage[paperheight=10in, paperwidth=5in, margin=10mm]{geometry}
\usepackage{graphicx}
\usepackage{listings}
\usepackage{longtable}
\setlength{\parindent}{0em}
\setlength{\parskip}{0em}
\begin{document}
'''


def svg_to_pdf(svg_file):
    pdf_file = svg_file.replace('.svg', '.pdf')
    cmd = ['rsvg-convert', '-f', 'pdf', '-o',
           pdf_file, svg_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    return pdf_file


def all_svg_to_pdf(path):
    path_list = glob.glob(path + "/*.svg")
    # print(path_list)
    for file_path in path_list:
        svg_to_pdf(file_path)
    return True


def cv_from_data(path_list: list, scap_folder,
                 model_save_file=None, sc_save_file=None, feature_set=[],
                 random_state=42):
    df, merged_df = read_data(
        path_list, feature_set=feature_set, to_filter=True)

    # Convert the input field for the hold-out filtering
    id_list = [f for f in merged_df['input'].apply(
        lambda x: x.split(os.sep)[-2]).to_list() if f != 'input']
    if len(id_list) > 0:
        merged_df['input'] = merged_df['input'].apply(
            lambda x: x.split(os.sep)[-2])
    else:
        merged_df['input'] = merged_df['input'].apply(
            lambda x: x.split(os.sep)[-3])

    cv_output_folder = os.path.dirname(model_save_file)
    os.makedirs(cv_output_folder, exist_ok=True)

    table_latex = []
    training_sum_num = [len(df[df['ans'] == 1]), len(df[df['ans'] == 0])]

    for hold_file in path_list:
        input_file = hold_file.split(os.sep)[-2]
        print('-------------------------------------------------------------------------------')
        with Timing(f'Fitting with holdout {input_file}'):
            training_df = df[merged_df['input'] != input_file].copy()
            to_transform = True
            X, Y, sc, descriptions = process4learning(
                training_df, to_transform=to_transform)
            rfc = learn(X, Y, random_state=random_state)

            # Measurements on the validation set
            val_df = df[merged_df['input'] == input_file].copy()
            output_df = merged_df[merged_df['input'] == input_file].copy()
            X_val, Y_val, _, _ = process4learning(
                val_df, to_transform=False)
            if X_val.shape[0] == 0:
                continue
            if sc:
                X_val = sc.transform(X_val)

            YY_val = predict_custom_cutoff(rfc, X_val)
            prob = rfc.predict_proba(X_val)
            if prob.shape[1] > 1:
                prob = prob[:, 1]

            training_pos = np.count_nonzero(Y == 1)
            training_neg = np.count_nonzero(Y == 0)
            table_list = print_stats(f'{input_file}', Y_val, YY_val,
                                     num_train=[training_pos, training_neg])
            table_latex.append((input_file, table_list))
            print('')
            print_forest_stats(rfc)
            print('')
            print_feature_importance(rfc, descriptions, X, Y)

            output_df.insert(output_df.columns.get_loc(
                "ans") + 1, "pred", prob)
            failure_output_df = output_df[output_df['ans'] != YY_val]

            del training_df, output_df, val_df

    # Sum up the statistics
    sum_table_list = sum_stats(
        [table_list for input_file, table_list in table_latex])
    sum_table_list = [
        ['Training count', training_sum_num[0] + training_sum_num[1]],
        ['| Positive training',
         training_sum_num[0]],
        ['| Negative training',
         training_sum_num[1]]
    ] + sum_table_list
    df = add_table_to_df([['Input', 'Total']] + sum_table_list, df=[])

    with temporary_directory() as tmpdir:
        print(tmpdir)

        # Generate the latex summary
        tex_str = TEMPLATE

        tex_str += str(tabulate(sum_table_list, tablefmt='latex') + '\n\n')
        tex_str += r'\clearpage' + '\n\n'

        for input_file, table_list in table_latex:
            df = add_table_to_df([['Input', input_file]] + table_list, df)
            tex_str += input_file.replace('_', r'\_') + '\n\n'
            tex_str += str(tabulate(table_list, tablefmt='latex') + '\n\n')
            tex_str += r'\vspace{3em}' + '\n\n'
        tex_str += r'\end{document}'
        tex_name = tmpdir + '/summary.tex'
        with open(tex_name, 'w') as f:
            f.write(tex_str)
        subprocess.check_call(['pdflatex', '--interaction=batchmode', tex_name],
                              cwd=tmpdir)
        # Since the temporary folder will be removed,
        # we move the pdf to the current directory
        shutil.move(os.path.join(tmpdir, tex_name.replace('.tex', '.pdf')),
                    cv_output_folder + '/summary.pdf')
        df.to_csv(cv_output_folder + '/summary.csv', index=False)


def save(fig, fname: str):
    if os.path.dirname(fname):
        os.makedirs(os.path.dirname(fname), exist_ok=True)
    fig.savefig(fname, transparent=False, bbox_inches=0, pad_inches=0)
