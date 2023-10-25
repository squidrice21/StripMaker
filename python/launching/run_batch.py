#!/usr/bin/env python

"""
Run data generation (optional), training and cross validation in batch.
"""

import argparse
import asyncio
import codecs
import locale
import os
from pathlib import Path
import subprocess
import sys
import threading
import yaml

from learning.feature_set import (
    feature_set, save_feature_set, secondary_feature_set)
from utils.util import rename_feature_set

lock = threading.Lock()
processes = []

INSTANCE_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/run_pipeline.py')
INSTANCE_SEC_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/run_pipeline_secondary.py')
HIST_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/plot_feature_distribution.py')

CONVERT_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../utils/forest_to_c.py')
MODEL_CPP = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../../src/classifier/forest.cpp')
SEC_MODEL_CPP = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../../src/classifier/forest_sec.cpp')


async def read_stream(stream):
    lines = []
    while not stream.at_eof():
        data = await stream.readline()
        line = data.decode(locale.getpreferredencoding(False))
        if len(line):
            lines.append(line)
    return lines


async def run_subprocess_in_bg(case, cmd, timeout):
    # print("Running {} `{}`".format(case, " ".join([str(c) for c in cmd])))
    print("Running {}".format(case))
    # https://stackoverflow.com/questions/45769985/asyncio-create-subprocess-exec-console-window-opening-for-each-call
    # This flag avoids popup window after crashing
    if sys.platform == "win32":
        DETACHED_PROCESS = 0x00000008
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            creationflags=DETACHED_PROCESS,
        )
    else:
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )

    with lock:
        processes.append(process)

    task_code = asyncio.ensure_future(process.wait())
    task_out = asyncio.ensure_future(read_stream(process.stdout))
    task_err = asyncio.ensure_future(read_stream(process.stderr))

    done, pending = await asyncio.wait([task_code, task_out, task_err], timeout=timeout)
    if pending:
        # timeout
        if process.returncode is None:
            # kill the subprocess, then `await future` will return soon
            try:
                print(f'* Killing {case} after {timeout} s...')
                process.kill()
            except ProcessLookupError:
                pass

    code = await task_code
    out = await task_out
    err = await task_err

    if code == 0:
        print(f'* Success: {case}.')
    else:
        print(f'* Failure: {case} {code}...')

    return code, out, err, case


async def run_test(case, cmd, timeout, semaphore):
    # prevent more than certain number things to run at the same time
    async with semaphore:
        return await run_subprocess_in_bg(case, cmd, timeout)


async def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'train', help='YAML file of inputs to train or a model file')
    parser.add_argument(
        '-f', '--fea', help='YAML file of the feature definitions')
    parser.add_argument('-e', '--exe', help='data generation exe')
    parser.add_argument('-s', '--snapshot', help='snapshot path')
    parser.add_argument('--secondary', action='store_true',
                        help='indicate the secondary model procedure')
    parser.add_argument('--save', nargs='?', const=Path(''),
                        type=Path, help='save model as a cpp file')
    parser.add_argument('--debug', action='store_true',
                        help='place a breakpoint before execution')
    parser.add_argument('--timeout', type=int, default=1800,
                        help='timeout threshold in seconds')
    parser.add_argument('--max-procs', type=int, default=4,
                        help='max number of tests to run at the same time')
    args = parser.parse_args()

    if args.debug:
        # Give time to attach a debugger.
        breakpoint()

    if args.save:
        if str(args.save) == '.':
            args.save = Path(MODEL_CPP) if not args.secondary else Path(
                SEC_MODEL_CPP)

    if args.secondary:
        args.secondary = args.snapshot

    if args.snapshot:
        os.makedirs(args.snapshot, exist_ok=True)
    else:
        os.makedirs("snapshot/", exist_ok=True)
        args.snapshot = 'snapshot'

    # Train or specify the model file if already trained
    model_file = args.train
    scap_folder = None
    if 'yml' in args.train:
        print('================================')
        print('Run data generation.')
        assert args.exe

        semaphore = asyncio.Semaphore(args.max_procs)

        coroutines = []
        feature_list_all = []
        with open(args.train) as f:
            input_ = yaml.safe_load(f)
        root = input_.get('root')
        # Create tasks
        for path in input_['drawings']:
            if isinstance(path, dict):
                path = path['path']

            path = root + "/" + path
            if scap_folder == None:
                scap_folder = os.path.dirname(path)
            if not args.secondary:
                cmd = [sys.executable, INSTANCE_PY, path,
                       '-e', args.exe, '-s', args.snapshot]
            else:
                # Match the model and input
                if '.sav' not in args.secondary:
                    model_file = args.secondary + '/cv/model_' + \
                        os.path.basename(path).replace('.scap', '.sav')
                else:
                    model_file = args.secondary
                cmd = [sys.executable, INSTANCE_SEC_PY, path,
                       '-e', args.exe, '-s', args.snapshot]

            coroutines.append(
                run_test(path, cmd, args.timeout, semaphore=semaphore))
            if not args.secondary:
                feature_list_all.append(os.path.abspath(
                    path).replace('.scap', '/features_full.csv'))
            else:
                feature_list_all.append(os.path.abspath(
                    path).replace('.scap', '/features_sec_full.csv'))

        # Run data generation
        results = await asyncio.gather(*coroutines)
        for result in results:
            time_line = [l for l in result[1] if ' time: ' in l]
            if len(time_line) > 0:
                print(''.join(time_line))

        feature_list = []
        for fea_file in feature_list_all:
            if os.path.exists(fea_file):
                feature_list.append(fea_file)

        if not args.secondary:
            feature_list_yml = args.snapshot + '/feature_list.yml'
        else:
            feature_list_yml = args.snapshot + '/feature_sec_list.yml'
        with open(feature_list_yml, 'w') as outfile:
            yaml.dump(feature_list, outfile, default_flow_style=False)

        if not args.fea:
            if not args.secondary:
                # Save the defined feature set for this run
                renamed_feature_set = rename_feature_set(feature_set)
                save_feature_set(args.snapshot + '/feature_set.yml',
                                 renamed_feature_set)
            else:
                renamed_feature_set = rename_feature_set(secondary_feature_set)
                save_feature_set(args.snapshot + '/feature_sec_set.yml',
                                 renamed_feature_set)

        print('================================')
        print('Run training.')
        if args.secondary:
            cmd = [sys.executable, INSTANCE_SEC_PY,
                   feature_list_yml, '-s', args.snapshot]
        else:
            cmd = [sys.executable, INSTANCE_PY,
                   feature_list_yml, '-s', args.snapshot]

        # Load feature definitions
        if args.fea:
            cmd += ['-f', args.fea]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=False)
        result.stdout = codecs.decode(result.stdout, 'utf-8', 'ignore')
        print(result.stdout)
        model_line = [l for l in result.stdout.split('\n') if 'Model: ' in l]
        time_line = [l for l in result.stdout.split('\n') if ' time: ' in l]
        if len(model_line) == 0:
            print(f'Training failed. Returned code: {result.returncode}...')
            raise RuntimeError('Training failed.')
        model_file = model_line[0].replace('Model: ', '')
        if len(time_line) > 0:
            print(time_line[-1])

        # Save the model as a C++ file
        if args.save:
            type_cmd = []
            model_cpp = str(args.save)
            if args.secondary:
                type_cmd = ['-s']

            sc_file = os.path.dirname(model_file) + '/sc_' + \
                os.path.basename(model_file)
            model_file
            cmd = [sys.executable, CONVERT_PY,
                   model_file, model_cpp, '-n', sc_file]
            cmd += type_cmd

            confirmed = True
            if os.path.exists(model_cpp):
                while True:
                    confirmation = input(
                        f"The C++ file '{model_cpp}' already exists. Do you want to overwrite? (yes/no): ").strip().lower()
                    if confirmation in ["yes", "no"]:
                        break
                    else:
                        print("Invalid response. Please enter 'yes' or 'no'.")

                if confirmation == "no":
                    print("Operation aborted.")
                    confirmed = False

            if confirmed:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)

        print('================================')
        print('Run CV.')
        assert scap_folder
        if args.secondary:
            cmd = [sys.executable, INSTANCE_SEC_PY,
                   feature_list_yml, '--cv', '-c', scap_folder]
        else:
            cmd = [sys.executable, INSTANCE_PY,
                   feature_list_yml, '--cv', '-c', scap_folder]
        cmd += ['-s', args.snapshot + '/cv/', ]

        # Load feature definitions
        if args.fea:
            cmd += ['-f', args.fea]

        subprocess.run(cmd)


if __name__ == '__main__':
    if sys.platform == "win32":
        loop = asyncio.ProactorEventLoop()
        asyncio.set_event_loop(loop)
    try:
        asyncio.get_event_loop().run_until_complete(main())
    except KeyboardInterrupt:
        with lock:
            for proc in processes:
                try:
                    proc.kill()
                except ProcessLookupError:
                    pass
