#!/usr/bin/env python
import yaml
import sys
import subprocess
import os
import glob
import asyncio
import argparse

from launching.run_prediction import temporary_directory
from utils.async_utils import async_entrance, run_with_limit

sys.path.insert(0, os.path.abspath(
    os.path.split(os.path.abspath(__file__))[0] + '/../'))


VIS_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/scap_to_svg.py')
UNCUT_VIS_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/visualize_uncutting.py')


async def run_preprocessing(data_yml, preprocessing_exe, timeout, max_procs, output_directory):
    semaphore = asyncio.Semaphore(max_procs)
    coroutines = []

    with open(data_yml) as f:
        input_ = yaml.safe_load(f)
    root = input_.get('root')

    # Create tasks
    preprossed_files = []
    for path in input_['drawings']:
        if isinstance(path, dict):
            path = path['path']

        simple_path = path
        path = root + "/" + path
        output_path = path.replace('.scap', '_prep.scap')
        if output_directory:
            output_path = output_directory + "/" + simple_path
            if os.path.abspath(path) == os.path.abspath(output_path):
                output_path = path.replace('.scap', '_prep.scap')

        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        if preprocessing_exe:
            cmd = [preprocessing_exe, '-u', path, '-o', output_path]
            # print(' '.join(cmd))
            # exit()
            coroutines.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))
        preprossed_files.append(output_path)

    # Run
    if preprocessing_exe:
        results = await asyncio.gather(*coroutines)

    coroutines_vis = []
    for cluster_f in preprossed_files:
        cmd = [sys.executable, VIS_PY, '-c', cluster_f,
               cluster_f.replace('.scap', '.svg')]
        # print(' '.join(cmd))
        # exit()
        coroutines_vis.append(
            run_with_limit(path, cmd, timeout, semaphore=semaphore))

    results = await asyncio.gather(*coroutines_vis)

    return preprossed_files


async def run_order(data_yml, order_exe, ref_folders, timeout, max_procs, output_directory):
    semaphore = asyncio.Semaphore(max_procs)
    coroutines = []

    with open(data_yml) as f:
        input_ = yaml.safe_load(f)
    root = input_.get('root')

    # Create tasks
    preprossed_files = []
    for path in input_['drawings']:
        if isinstance(path, dict):
            path = path['path']

        simple_path = path
        path = root + "/" + path
        output_path = path.replace('.scap', '_prep.scap')
        if output_directory:
            output_path = output_directory + "/" + simple_path
            if os.path.abspath(path) == os.path.abspath(output_path):
                output_path = path.replace('.scap', '_prep.scap')

        if order_exe:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            cmd = [order_exe, '-i', path, '-o', output_path]
            name_match = [os.path.basename(path).replace('.scap', ''), os.path.basename(
                path).replace('.scap', '').replace('_SA_cluster', '')]
            ref_files = []
            for folder in ref_folders:
                ref_files = glob.glob(folder + "/*.scap")
                ref_files = [f for f in ref_files if os.path.basename(
                    f).replace('.scap', '') in name_match]
                if len(ref_files) > 0:
                    break
            if len(ref_files) == 0:
                continue

            cmd += ['-d', ref_files[0]]
            # print(' '.join(cmd))
            # exit()
            coroutines.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))
            preprossed_files.append(output_path)

    # Run
    if order_exe:
        results = await asyncio.gather(*coroutines)

    coroutines_vis = []
    for cluster_f in preprossed_files:
        if os.path.exists(cluster_f):
            cmd = [sys.executable, VIS_PY, '-c', cluster_f,
                   cluster_f.replace('.scap', '.svg')]
            # print(' '.join(cmd))
            # exit()
            coroutines_vis.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))

    results = await asyncio.gather(*coroutines_vis)

    return preprossed_files


async def run_uncutting(data_yml, uncutting_exe, timeout, max_procs, output_directory):
    semaphore = asyncio.Semaphore(max_procs)
    coroutines = []

    with open(data_yml) as f:
        input_ = yaml.safe_load(f)
    root = input_.get('root')

    # Create tasks
    preprossed_files = []
    for path in input_['drawings']:
        use_gt = False
        use_angle = False
        raw_limited = False
        original_input = ''
        if isinstance(path, dict):
            if 'use_gt' in path:
                use_gt = path['use_gt']
            if 'use_angle' in path:
                use_angle = path['use_angle']
            if 'raw_limited' in path:
                raw_limited = path['raw_limited']
            if 'original_input' in path:
                original_input = root + "/" + path['original_input']
            # if 'original_input' not in path:
            #     continue
            path = path['path']
        # use_gt = False
        # use_angle = False

        simple_path = path
        path = root + "/" + path
        output_path = path.replace('.scap', '_prep.scap')
        if output_directory:
            output_path = output_directory + "/" + simple_path
            if os.path.abspath(path) == os.path.abspath(output_path):
                output_path = path.replace('.scap', '_prep.scap')

        if uncutting_exe:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)

            cmd = [uncutting_exe, '-i', path, '-o', output_path,
                   '--vis', output_path.replace('.scap', '_uncut.csv'),
                   '--csv', output_path.replace('.scap', '_cut.csv')]
            if use_gt:
                cmd += ['--gt']
            if use_angle:
                cmd += ['--angle']
            if len(original_input) > 0 and os.path.exists(original_input):
                cmd += ['--raw', original_input]
                if raw_limited:
                    cmd += ['--limit']
            # print(' '.join(cmd))
            # exit()
            coroutines.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))
            preprossed_files.append(output_path)

    # Run
    if uncutting_exe:
        results = await asyncio.gather(*coroutines)

    coroutines_vis = []
    for cluster_f in preprossed_files:
        if os.path.exists(cluster_f):
            cmd = [sys.executable, VIS_PY, '-c', cluster_f,
                   cluster_f.replace('.scap', '.svg')]
            # print(' '.join(cmd))
            # exit()
            coroutines_vis.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))

    results = await asyncio.gather(*coroutines_vis)

    return preprossed_files


async def run_sa(preprossed_files, sa_exe, timeout, max_procs, output_directory):
    semaphore = asyncio.Semaphore(max_procs)
    coroutines = []

    # Create tasks
    fit_scap = []
    cluster_scap = []
    for path in preprossed_files:
        if isinstance(path, dict):
            path = path['path']

        # -t: Run post processing; -m: Output intermediate files
        # -c: Conservative setting; -p: Skip preprocessing
        # abs_f = os.path.abspath(path)
        # des_f = abs_f.replace(os.path.abspath(output_directory), output_directory + '/sa_out/')
        des_f = output_directory + '/' + os.path.basename(path)
        os.makedirs(os.path.dirname(des_f), exist_ok=True)

        # cmd = [sa_exe, '-t', '-m', '-c', '-p', path, '-o', des_f]
        # cmd = [sa_exe, '-t', '-m', '-p', path, '-o', des_f]
        cmd = [sa_exe, '-p', path, '-o', des_f]
        # print(' '.join(cmd))
        if not os.path.exists(des_f.replace('.scap', '_cluster.scap')):
            coroutines.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))
        else:
            print(f'Skip: {os.path.basename(des_f)}')

        fit_scap.append(des_f)
        cluster_scap.append(des_f.replace('.scap', '_cluster.scap'))

    # Run
    sa_results = await asyncio.gather(*coroutines)

    # Visualize as SVG files
    coroutines_vis = []
    for i in range(len(fit_scap)):
        fit_f = fit_scap[i]
        cluster_f = cluster_scap[i]
        if os.path.exists(fit_f) and not os.path.exists(fit_f.replace('.scap', '.svg')):
            cmd = [sys.executable, VIS_PY, fit_f,
                   fit_f.replace('.scap', '.svg')]
            coroutines_vis.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))
        if os.path.exists(cluster_f) and not os.path.exists(cluster_f.replace('.scap', '.svg')):
            cmd = [sys.executable, VIS_PY, '-c', cluster_f,
                   cluster_f.replace('.scap', '.svg')]
            coroutines_vis.append(
                run_with_limit(path, cmd, timeout, semaphore=semaphore))

    results = await asyncio.gather(*coroutines_vis)

    # Read out
    print('======================================')
    for code, out, err, case in sa_results:
        overall_out = err + out
        entry_line = [l for l in overall_out if 'time: ' in l]

        if len(entry_line) > 0:
            print(entry_line[0].replace('\n', ''))
    print('======================================')


async def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', help='Input file')
    parser.add_argument('-e', '--exe', help='data preprocessing exe')
    parser.add_argument('-o', '--output', help='output directory')
    parser.add_argument('-p', '--pre', action='store_true',
                        help='StrokeAggregator preprocessing exe')
    parser.add_argument('-d', '--order', help='data ordering exe')
    parser.add_argument(
        '-r', '--ref', help='reference folder(s) for data ordering')
    parser.add_argument('-i', '--intermediate',
                        help='intermediate folder for data ordering')
    parser.add_argument('-u', '--uncut', action='store_true',
                        help='run uncutting on preprocessed inputs')
    parser.add_argument('-v', '--vis',
                        help='visualization output path')
    parser.add_argument('-c', '--comparison',
                        action='store_true', help='StrokeAggregator exe')
    parser.add_argument('--timeout', type=int, default=2400,
                        help='timeout threshold in seconds')
    parser.add_argument('--max-procs', type=int, default=5,
                        help='max number of tests to run at the same time')
    parser.add_argument('--debug', action='store_true',
                        help='place a breakpoint before execution')
    args = parser.parse_args()

    if args.debug:
        # Give time to attach a debugger.
        breakpoint()

    # Run preprocessing
    preprossed_files = None
    if args.pre:
        preprossed_files = await run_preprocessing(args.input, args.exe, args.timeout,
                                                   args.max_procs, args.output)

    with temporary_directory() as tmpdir:
        if args.order:
            assert args.ref
            ref_folders = args.ref.split(',')
            order_folder = args.intermediate
            if not order_folder:
                order_folder = tmpdir
            preprossed_files = await run_order(args.input, args.order, ref_folders, args.timeout,
                                               args.max_procs, order_folder)
            if args.uncut:
                with open(args.input) as f:
                    input_ = yaml.safe_load(f)
                input_['root'] = os.path.abspath(order_folder)
                with open(order_folder + '/' + os.path.basename(args.input), 'w') as outfile:
                    yaml.dump(input_, outfile,
                              default_flow_style=False, sort_keys=False)
                args.input = order_folder + '/' + os.path.basename(args.input)
            if len(preprossed_files) > 0:
                order_folder = os.path.dirname(preprossed_files[0])

        # Run uncutting on preprocessed data (to remove visually incorrect cuts already exist
        # in our input to the preprocessing)
        if args.uncut:
            preprossed_files = await run_uncutting(args.input, args.exe, args.timeout,
                                                   args.max_procs, args.output)
            uncut_folder = ''
            if len(preprossed_files) > 0:
                uncut_folder = os.path.dirname(preprossed_files[0])
            if args.vis and len(uncut_folder) > 0:
                cmd = [sys.executable, UNCUT_VIS_PY, '-i',
                       uncut_folder, '-s', order_folder, '-o', args.vis + '_uncut', '--csv', '_uncut.csv']
                # print(' '.join(cmd))
                subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)

                cmd = [sys.executable, UNCUT_VIS_PY, '-i',
                       uncut_folder, '-s', uncut_folder, '-o', args.vis + '_cut', '--csv', '_cut.csv']
                # print(' '.join(cmd))
                subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)

    # Optional: Run comparison method (StrokeAggregator)
    if args.comparison:
        if not preprossed_files:
            with open(args.input) as f:
                input_ = yaml.safe_load(f)
            root = input_.get('root')

            # Create tasks
            preprossed_files = []
            for path in input_['drawings']:
                if isinstance(path, dict):
                    path = path['path']
                preprossed_files.append(root + '/' + path)

        await run_sa(preprossed_files, args.exe,
                     args.timeout, args.max_procs, args.output)


if __name__ == '__main__':
    async_entrance(main)
