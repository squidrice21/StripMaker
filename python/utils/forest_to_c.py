#!/usr/bin/env python

"""
Convert the sklearn endpoint-endpoint and endpoint-stroke classifiers into C++
code.
"""

import argparse
import pickle

from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.tree import _tree


def tree_to_c(tree, function_name: str) -> str:
    tree_ = tree.tree_

    string_builder = []
    string_builder.append(
        f"static double {function_name}(const double* feature_vec) {{")

    def recurse(node, depth):
        indent = "  " * depth
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            threshold = tree_.threshold[node]
            feature = tree_.feature[node]
            string_builder.append(
                f"{indent}if (feature_vec[{feature}] <= {threshold}) {{")
            recurse(tree_.children_left[node], depth + 1)
            string_builder.append(f"{indent}}} else {{")
            recurse(tree_.children_right[node], depth + 1)
            string_builder.append(f"{indent}}}")
        else:
            val = tree_.value[node, 0, 1] / tree_.value[node, 0].sum()
            string_builder.append(f"{indent}return {val};")

    recurse(0, 1)

    string_builder.append('}\n')

    return '\n'.join(string_builder)


def forest_to_c(clf: RandomForestClassifier, function_name: str) -> str:
    string_builder = []
    for i, est in enumerate(clf.estimators_):
        string_builder.append(tree_to_c(est, f'{function_name}_tree_{i}'))
        # string_builder.append(export_text(est))
    param_type = 'std::vector<double>'
    string_builder.append(
        f'double {function_name}(const {param_type}& feature_vec) {{')
    string_builder.append(
        f'  assert(feature_vec.size() == {clf.n_features_in_});')
    string_builder.append(f'  const double* feat = feature_vec.data();')
    string_builder.append('  double pred = 0.0;')
    for i in range(len(clf.estimators_)):
        string_builder.append(f'  pred += {function_name}_tree_{i}(feat);')
    string_builder.append(f'  return pred / {len(clf.estimators_)}.0;\n}}\n')
    return '\n'.join(string_builder)


def scaler_to_c(scaler: StandardScaler, function_name: str) -> str:
    string_builder = []
    param_type = 'std::vector<double>'
    string_builder.append(
        f'void {function_name}({param_type}& feature_vec) {{')

    if scaler is not None and (scaler.mean_ is not None or scaler.scale_ is not None):
        string_builder.append(
            f'  assert(feature_vec.size() == {scaler.n_features_in_});')
        if scaler.mean_ is not None:
            string_builder.append(f'  std::vector<double> mean{{')
            for v in scaler.mean_:
                string_builder.append(f'    {v},')
            string_builder.append(f'  }};')
        if scaler.scale_ is not None:
            string_builder.append(f'  std::vector<double> var{{')
            for v in scaler.scale_:
                string_builder.append(f'    {v},')
            string_builder.append(f'  }};')
        string_builder.append(
            f'  for (size_t i = 0; i < feature_vec.size(); ++i) {{')
        if scaler.mean_ is not None:
            string_builder.append(f'    feature_vec[i] -= mean[i];')
        if scaler.scale_ is not None:
            string_builder.append(f'    feature_vec[i] /= var[i];')
        string_builder.append(f'  }}\n')

    string_builder.append(f'}}\n')
    return '\n'.join(string_builder)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('classifier', help='path to pickled classifier')
    parser.add_argument('out', help='path to c++ output')
    parser.add_argument('-n', '--norm', help='normalizer')
    parser.add_argument('-s', '--sec', action='store_true',
                        help='Whether input is a secondary cluster classifier')
    args = parser.parse_args()

    with open(args.classifier, 'rb') as f:
        clf = pickle.load(f)

    string_builder = []
    string_builder.append('// Automatically generated using forest_to_c.py.')
    string_builder.append('// Do not edit.\n')
    if args.sec:
        string_builder.append('#include "forest_sec.h"\n')
    else:
        string_builder.append('#include "forest.h"\n')
    string_builder.append('#include <cassert>\n')
    if args.sec:
        string_builder.append(forest_to_c(clf, 'clf_sec'))
    else:
        string_builder.append(forest_to_c(clf, 'clf'))

    with open(args.out, 'w') as f:
        f.write('\n'.join(string_builder))

    if args.norm:
        with open(args.norm, 'rb') as f:
            scaler = pickle.load(f)

        string_builder = []
        string_builder.append(
            '// Automatically generated using forest_to_c.py.')
        string_builder.append('// Do not edit.\n')
        if args.sec:
            string_builder.append('#include "forest_sec.h"\n')
        else:
            string_builder.append('#include "forest.h"\n')

        string_builder.append('#include <cassert>\n')
        if args.sec:
            string_builder.append(scaler_to_c(scaler, 'normalize_sec'))
        else:
            string_builder.append(scaler_to_c(scaler, 'normalize'))

        with open(args.out.replace('.cpp', '_norm.cpp'), 'w') as f:
            f.write('\n'.join(string_builder))
    else:
        string_builder = []
        string_builder.append(
            '// Automatically generated using forest_to_c.py.')
        string_builder.append('// Do not edit.\n')
        string_builder.append('#include "forest.h"\n')
        string_builder.append('#include <cassert>\n')
        string_builder.append(scaler_to_c(None, 'normalize'))
        with open(args.out.replace('.cpp', '_norm.cpp'), 'w') as f:
            f.write('\n'.join(string_builder))

    print('Converted')


if __name__ == '__main__':
    main()
