#!/usr/bin/env python

"""
Run local inference and visualize resulting files.
Optionally, convert given model files into C++ code.
"""

import argparse
import codecs
import os
import subprocess
import sys
import time

CONVERT_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../utils/forest_to_c.py')
VIS_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/scap_to_svg.py')
MODEL_CPP = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../../src/classifier/forest.cpp')
SEC_MODEL_CPP = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../../src/classifier/forest_sec.cpp')

# Change this to your own MSBuild path
MSBUILD_PATH = 'C:/Program Files (x86)/Microsoft Visual Studio/2017/BuildTools/MSBuild/15.0/Bin/MSBuild.exe'
SLN_PATH = os.path.abspath(os.path.split(os.path.abspath(__file__))[
                           0]+'/../../../sketch_clustering-build/strokestrip.sln')


def prepare_forest(model, type='solve'):
    type_cmd = []
    model_cpp = MODEL_CPP
    if type == 'solve':
        model_cpp = MODEL_CPP
    elif type == 'secondary':
        model_cpp = SEC_MODEL_CPP
        type_cmd = ['-s']
    else:
        assert False

    model_file = model
    sc_file = os.path.dirname(model_file) + '/sc_' + \
        os.path.basename(model_file)

    assert os.path.exists(model_file) and os.path.exists(sc_file)
    cmd = [sys.executable, CONVERT_PY,
           model_file, model_cpp, '-n', sc_file]
    cmd += type_cmd
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    confirm_line = [l for l in result.stdout.split(
        '\n') if 'Converted' in l]
    if len(confirm_line) == 0:
        print(f'Error: Cannot convert model {model_file}')
        print(result.stdout)
        print(result.sterr)
        exit(-1)
    print(f'Converted: {model_file} => {model_cpp}')


def prepare_classifiers(exe, model, exe_type='solve', sec_model=None):
    exe_name = exe

    # Convert and recompile
    to_recompile = False
    if not exe:
        assert model
    if model:
        assert exe

        to_recompile = True

        model_file = model
        sc_file = os.path.dirname(model_file) + '/sc_' + \
            os.path.basename(model_file)

        assert os.path.exists(model_file) and os.path.exists(sc_file)
        cmd = [sys.executable, CONVERT_PY,
               model_file, MODEL_CPP, '-n', sc_file]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        confirm_line = [l for l in result.stdout.split(
            '\n') if 'Converted' in l]
        if len(confirm_line) == 0:
            print(f'Error: Cannot convert model {model_file}')
            print(result.stdout)
            print(result.sterr)
            exit(-1)
        print(f'Converted: {model_file} => {MODEL_CPP}')

    if sec_model:
        assert exe

        to_recompile = True

        model_file = sec_model
        sc_file = os.path.dirname(model_file) + '/sc_' + \
            os.path.basename(model_file)

        assert os.path.exists(model_file) and os.path.exists(sc_file)
        cmd = [sys.executable, CONVERT_PY, '-s',
               model_file, SEC_MODEL_CPP, '-n', sc_file]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True)
        confirm_line = [l for l in result.stdout.split(
            '\n') if 'Converted' in l]
        if len(confirm_line) == 0:
            print(f'Error: Cannot convert model {model_file}')
            print(result.stdout)
            print(result.sterr)
            exit(-1)
        print(f'Converted: {model_file} => {SEC_MODEL_CPP}')

    if to_recompile:
        # In this case, the input exe should be the build directory
        if sys.platform == "darwin":
            if os.path.isfile(exe) or not os.path.exists(exe):
                exe = os.path.dirname(exe)
            if model:
                model_name = os.path.basename(model).replace('.sav', '')
            elif sec_model:
                model_name = os.path.basename(sec_model).replace('.sav', '')

            if exe_type == 'solve':
                if os.path.exists(exe + '/clustering-solve'):
                    os.remove(exe + '/clustering-solve')
                exe_name = exe + f'/clustering-solve-{model_name}'
                if os.path.exists(exe_name):
                    os.remove(exe_name)

                result = subprocess.run(['make'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=False, cwd=exe)
                result.stdout = codecs.decode(result.stdout, 'utf-8', 'ignore')
                result.stderr = codecs.decode(result.stderr, 'utf-8', 'ignore')
                if not os.path.exists(exe + '/clustering-solve'):
                    print(f'Error: Cannot compile model {model_file}')
                    print(result.stdout)
                    print(result.stderr)
                    exit(-1)

                os.rename(exe + '/clustering-solve', exe_name)
                print(f'Compiled model: {model_file} into {exe_name}')
            elif exe_type == 'secondary':
                exe_target = ['clustering-secondary-sample',
                              'clustering-secondary-feature']
                exe_name = []
                for target in exe_target:
                    if os.path.exists(exe + f'/{target}'):
                        os.remove(exe + f'/{target}')
                    exe_name.append(exe + f'/{target}-{model_name}')
                    if os.path.exists(exe_name[-1]):
                        os.remove(exe_name[-1])

                result = subprocess.run(['make'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=False, cwd=exe)
                result.stdout = codecs.decode(result.stdout, 'utf-8', 'ignore')
                result.stderr = codecs.decode(result.stderr, 'utf-8', 'ignore')
                if not os.path.exists(exe + f'/{target}'):
                    print(f'Error: Cannot compile model {model_file}')
                    print(result.stdout)
                    print(result.stderr)
                    exit(-1)

                for i in range(len(exe_target)):
                    target = exe_target[i]
                    os.rename(exe + f'/{target}', exe_name[i])
                    print(f'Compiled model: {model_file} into {exe_name[i]}')
                exe_name = ','.join(exe_name)
            elif exe_type == 'junction':
                if os.path.exists(exe + '/clustering-junction'):
                    os.remove(exe + '/clustering-junction')
                exe_name = exe + f'/clustering-junction-{model_name}'
                if os.path.exists(exe_name):
                    os.remove(exe_name)

                result = subprocess.run(['make'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=False, cwd=exe)
                result.stdout = codecs.decode(result.stdout, 'utf-8', 'ignore')
                result.stderr = codecs.decode(result.stderr, 'utf-8', 'ignore')
                if not os.path.exists(exe + '/clustering-junction'):
                    print(f'Error: Cannot compile model {model_file}')
                    print(result.stdout)
                    print(result.stderr)
                    exit(-1)

                os.rename(exe + '/clustering-junction', exe_name)
                print(f'Compiled model: {model_file} into {exe_name}')
        elif sys.platform == 'win32':
            if os.path.isfile(exe) or not os.path.exists(exe):
                exe = os.path.dirname(exe)
            if os.path.exists(exe + '/clustering-solve.exe'):
                os.remove(exe + '/clustering-solve.exe')

            model_name = os.path.basename(model).replace('.sav', '')
            exe_name = exe + f'/clustering-solve-{model_name}.exe'
            if os.path.exists(exe_name):
                os.remove(exe_name)

            build_cmd = [MSBUILD_PATH, SLN_PATH, '/t:clustering-solve',
                         '/p:Configuration=Release', '/p:Platform=x64']
            try:
                subprocess.check_call(build_cmd)
            except subprocess.CalledProcessError:
                print(f'Error: Cannot compile model {model_file}')
                print(result.stdout)
                print(result.stderr)
                exit(-1)
            os.rename(exe + '/clustering-solve.exe', exe_name)
            print(f'Compiled model: {model_file} into {exe_name}')

    return exe_name


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='Input file')
    parser.add_argument('-m', '--model', help='Input model file')
    parser.add_argument('-e', '--exe', help='data generation exe')
    parser.add_argument('-o', '--output', help='Output path')
    parser.add_argument(
        '-f', '--fea', help='YAML file of the feature definitions')
    parser.add_argument('-v', '--vis', action='store_true',
                        help='Save intermediate visualizations')
    parser.add_argument('--sec', type=str,
                        help='Input secondary model file')
    parser.add_argument('-u', '--uncut', type=str,
                        help='To merge cut points under the hood')
    parser.add_argument(
        '-c', '--cached', action='store_true', help='the input scap files store first step results')
    parser.add_argument('-s', '--solver', type=str,
                        help='solver command line options')
    parser.add_argument('--debug', action='store_true',
                        help='Place a breakpoint before execution')
    args = parser.parse_args()

    if args.debug:
        breakpoint()

    predict_file = args.input

    args.exe = prepare_classifiers(
        args.exe, args.model, sec_model=args.sec)

    if not args.output:
        args.output = args.input.replace('.scap', '_pred.scap')

    if os.path.isfile(args.output):
        out_folder = os.path.dirname(args.output)
    else:
        out_folder = args.output
    os.makedirs(out_folder, exist_ok=True)

    # Run C++ pipeline
    print(
        f'Solving for {os.path.dirname(predict_file)} at {args.output}...')
    cmd = [args.exe, '-i', predict_file, '-o', args.output]
    if args.fea:
        for fea in args.fea.split(','):
            if '_sec_' in fea:
                cmd += ['--sec', fea]
            else:
                cmd += ['-f', fea]
    if args.vis:
        cmd += ['--interm', args.output]

    # For under-hood merging of cut positions
    if args.uncut:
        cmd += ['--uncut', args.uncut]

    if args.cached:
        cmd += ['--cached']
    if args.solver:
        solver_options = args.solver.split(',')
        for opt in solver_options:
            cmd += ['--' + opt]

    print(' '.join(cmd))

    # Commented out to just run the GT fit
    start = time.time()
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)
    end = time.time()

    print(result.stdout)

    predict_file = os.path.basename(predict_file)
    out_predict_file = predict_file.replace('.scap', '_out.scap')
    cmd = [sys.executable, VIS_PY, '-c', '-f', args.output + '/' + out_predict_file,
           (args.output + '/' + out_predict_file).replace('.scap', '.svg')]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True)

    print(f'{predict_file} time: {end - start} s')
