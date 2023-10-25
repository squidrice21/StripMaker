#!/usr/bin/env python

"""
Run consolidation and visualize in batch or on a single input.
"""

import argparse
import codecs
from contextlib import contextmanager
import glob
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import yaml

import jinja2

sys.path.insert(0, os.path.abspath(
    os.path.split(os.path.abspath(__file__))[0] + '/../'))


INSTANCE_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../learning/predict_cpp.py')
VIS_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/scap_to_svg.py')

TEMPLATE = r'''
\pdfsuppresswarningpagegroup=1
\documentclass[12pt,letterpaper]{article}
\usepackage[paperheight=22in, paperwidth=20in, margin=10mm]{geometry}
\usepackage{graphicx}
\usepackage{listings}
\usepackage{longtable}
\usepackage{tikz}
\setlength{\parindent}{0em}
\setlength{\parskip}{0em}
\begin{document}

\title{Results of ``StripMaker: Perception-driven Learned Vector Sketch Consolidation"}
\date{}
\maketitle

% This is a final result visualization page
\begin{centering}

{% for row in rows %}
\begin{longtable}{ccc}
     {{row['input']}} & {{row['refined_fit']}} & {{row['post_fit']}}  \\
     {{row['input_name']}} & Refined Consolidation & Post-processed Consolidation \\
     & {{row['refined_clustering']}} & {{row['post_clustering']}}  \\
     & Refined Clustering & Post-processed Clustering \\
     & {{row['refined_width_fit']}} & {{row['post_width_fit']}}  \\
     & Refined Consolidation with Width & Post-processed Consolidation with Width \\
\end{longtable}
\clearpage
{% endfor %}

\end{centering}

\end{document}
'''


@contextmanager
def temporary_directory():
    tmpdir = tempfile.mkdtemp()
    try:
        yield tmpdir
    except:
        print('\n*\n* An error occurred. The intermediate files are left in ' +
              f'\n* "{tmpdir}".\n*\n', file=sys.stderr)
        raise
    # Unlike tempfile.TemporaryDirectory, we do not want to remove the temporary
    # directory when an exception occurs because we want to be able to manually
    # inspect the directory afterwards to see what went wrong.
    shutil.rmtree(tmpdir)


def tex_graphics(path: str, size_ratio=1, check_existence=True) -> str:
    if (path == None) or (check_existence and not os.path.exists(path)):
        return '\\small(missing)'
    height = 6 * size_ratio
    width = 6 * size_ratio
    abs_path = os.path.abspath(path)
    return fr'\includegraphics[height={height:.1f}in, width={width:.1f}in, keepaspectratio]{{{abs_path}}}'


def svg_to_pdf(svg_file):
    pdf_file = svg_file.replace('.svg', '.pdf')
    cmd = ['rsvg-convert', '-f', 'pdf', '-o',
           pdf_file, svg_file]
    subprocess.run(cmd)
    pdf_file = os.path.abspath(pdf_file)
    return pdf_file


def scap_to_pdf(tmpdir, f, to_color=False):
    if not os.path.exists(f):
        return None
    svg_filename = tmpdir + '/' + \
        os.path.basename(f.replace('.scap', '.svg'))
    cmd = [sys.executable, VIS_PY, f, svg_filename]
    if to_color:
        cmd += ['-c']
    subprocess.run(cmd)

    pdf_filename = svg_to_pdf(svg_filename)

    return pdf_filename


def prepare_cmd(path, output_dir, use_junc, exe_str):
    exes = exe_str.split(',')

    cmd = [sys.executable, INSTANCE_PY, path,
           '--exe', exes[0], '-o', output_dir]
    cmd += ['--vis']

    if use_junc:
        cmd += ['-s=split,merge,junc']
    else:
        cmd += ['-s=split,merge']

    return cmd


def clean_intermediate_files(output_dir, row):
    if Path(output_dir, 'cache').exists():
        shutil.rmtree(Path(output_dir, 'cache'))

    keep_files = [Path(file).name for key,
                  file in row.items() if 'input' not in key]
    for filename in os.listdir(output_dir):
        if os.path.isfile(os.path.join(output_dir, filename)):
            if filename not in keep_files:
                # print(f'Removing {filename}...')
                os.remove(os.path.join(output_dir, filename))


def check_pdf_extension(value):
    _, extension = os.path.splitext(value)
    if extension.lower() == '.pdf':
        return value
    else:
        raise argparse.ArgumentTypeError(f"'{value}' is not a valid PDF file.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        'input', help='A single scap input file or a YAML file of inputs')
    parser.add_argument('-e', '--exe', required=True,
                        help='solver exe')
    parser.add_argument('-s', '--snapshot',
                        default='./snapshot', help='snapshot path')
    parser.add_argument('--nojunc', action='store_true',
                        help='disable junction split')
    parser.add_argument('--spiral', action='store_true',
                        help='detect and fit spirals')
    parser.add_argument(
        '-o', '--output', type=check_pdf_extension, help='output PDF path')
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

    assert len(args.exe.split(
        ',')) > 1, "Need to specify both consolidation and post-processing exes"

    # Parse inputs
    paths = []
    if args.input.endswith('.scap'):
        paths.append(args.input)
    else:
        with open(args.input) as f:
            input_ = yaml.safe_load(f)
        root = input_.get('root')

        for path in input_['drawings']:
            if isinstance(path, dict):
                path = path['path']

            path = root + "/" + path
            paths.append(path)

    # Run temporal consolidation and refinement
    rows = []
    for path in paths:
        cut_filename = path.split(os.sep)[-2] + '_cut.csv'
        file_id = path.split(os.sep)[-1].replace('.scap', '')

        python_path = INSTANCE_PY
        output_dir = args.snapshot + '/' + \
            os.path.basename(path).replace('.scap', '')

        os.makedirs(output_dir, exist_ok=True)

        if not Path(output_dir + '/' + Path(path).stem + '_out.svg').exists():
            cmd = prepare_cmd(path, output_dir, not args.nojunc, args.exe)

            # Find cut point csv
            cut_files = glob.glob(str(Path(path).parent / "*.csv"))
            cut_filename = [f for f in cut_files if cut_filename in f]
            if len(cut_filename) > 0:
                cut_filename = cut_filename[0]
            else:
                cut_filename = path.split(
                    os.sep)[-1].replace('.scap', '') + '_cut.csv'
                cut_filename = [f for f in cut_files if cut_filename in f]
                if len(cut_filename) > 0:
                    cut_filename = cut_filename[0]
            if Path(cut_filename).exists():
                cmd += ['--uncut', cut_filename]

            # Run inference
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=False)
            result.stdout = codecs.decode(result.stdout, 'utf-8', 'ignore')

            # Run post-processing: Enforce junctions; Merge outlier strokes
            # Find refined scap output
            split_filename = Path(output_dir, Path(
                path).stem + '_merge_out.scap')
            split_filename = str(split_filename)

            end_str = 'end'

            post_exe = args.exe.split(',')[1]
            cmd = [post_exe, '-i', split_filename, '--end', end_str,
                   '--vis', output_dir]
            if Path(cut_filename).exists():
                cmd += ['--uncut', cut_filename]
            if args.spiral:
                cmd += ['--spiral']

            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=False)

            # Rename intermediate files
            os.rename(Path(output_dir, 'merge_prediction_fit.svg'),
                      Path(output_dir, Path(path).stem + '_rfit.svg'))
            os.rename(Path(output_dir, 'merge_prediction_fit_width.svg'),
                      Path(output_dir, Path(path).stem + '_rfit_width.svg'))

            print(f'Finished processing {path}...')
        else:
            print(f'Found existing output for {path}. Skipping...')

        row = {
            'input': path,
            'refined_fit': output_dir + '/' + Path(path).stem + '_rfit.svg',
            'post_fit': output_dir + '/' + Path(path).stem + '_fit.svg',
            'input_name': file_id.replace('_', r'\_'),
            'refined_clustering': output_dir + '/' + Path(path).stem + '_merge_out.svg',
            'post_clustering': output_dir + '/' + Path(path).stem + '_final_out.scap',
            'refined_width_fit': output_dir + '/' + Path(path).stem + '_rfit_width.svg',
            'post_width_fit': output_dir + '/' + Path(path).stem + '_fit_width.svg'
        }

        rows.append(row)

        # Need to clean intermediate results and cache or different inputs would use the same cache
        clean_intermediate_files(output_dir, row)

    # Visualize
    if args.output:
        with temporary_directory() as tmpdir:
            for row in rows:
                for key, file in row.items():
                    if file.endswith('.scap'):
                        row[key] = tex_graphics(scap_to_pdf(
                            tmpdir, file, to_color='clustering' in key), 0.8)
                    elif file.endswith('.svg'):
                        row[key] = tex_graphics(svg_to_pdf(file), 0.8)

            # Render the tex given the template and values
            template = jinja2.Template(TEMPLATE)
            tex_name = os.path.basename(args.output).replace('.pdf', '.tex')
            with open(os.path.join(tmpdir, tex_name), 'w') as f:
                # Note we bind the variables here
                f.write(template.render(rows=rows))

            # Compile via latexmk
            print(tmpdir)
            subprocess.check_call(['latexmk', '-quiet', '-pdf', tex_name],
                                  cwd=tmpdir)
            shutil.move(os.path.join(
                tmpdir, tex_name.replace('.tex', '.pdf')), args.output)
