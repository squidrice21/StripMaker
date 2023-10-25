#!/usr/bin/env python

"""
Helper script to convert a scap file into a png via inkscape and imagemagick.
"""

import sys
import subprocess
import os
import argparse
import xml.etree.ElementTree as ET

from PIL import Image

sys.path.insert(0, os.path.abspath(
    os.path.split(os.path.abspath(__file__))[0] + '/../'))

VIS_PY = os.path.abspath(os.path.split(
    os.path.abspath(__file__))[0]+'/../visualization/scap_to_svg.py')


# Modified from https://github.com/Nauhcnay/A-Benchmark-for-Rough-Sketch-Cleanup
def get_sketch_size(sketch, root=None):
    '''
    Given:
        sketch, pixel sketch(.png) or vector sketch(.svg) name
        path, path to search folder
        root, xml root element
    Return:
        the size of sketch in give path with same name
        if path is given, it could be either png file or svg file
        if root is given, it could only be a xml root element of a svg file
    Caution:
        The unit of svg artboart should only be "pixel"
    '''
    if sketch is None:
        if root is None:
            print("xml mode detected, but no root given!")
            raise ValueError
        else:
            extension = '.svg'
    else:
        _, extension = os.path.splitext(sketch)

    if '.svg' == extension:
        if root is None:
            tree = ET.parse(sketch)
            root = tree.getroot()
        if "viewBox" in root.attrib:
            if ',' in root.attrib['viewBox']:
                width = float(root.attrib['viewBox'].split(',')[2])
                height = float(root.attrib['viewBox'].split(',')[3])
            else:
                width = float(root.attrib['viewBox'].split(' ')[2])
                height = float(root.attrib['viewBox'].split(' ')[3])
        elif "width" in root.attrib and "height" in root.attrib:
            width = float(root.attrib['width'].strip('px'))
            height = float(root.attrib['height'].strip('px'))
        else:
            print("Error:\tparsing svg failed")
        return width, height
    elif '.png' == extension:
        img = Image.open(sketch)
        width, height = img.size
        return (float(width), float(height))
    else:
        print(extension)
        print(sketch)
        print("Error:\tunspport input")
        raise ValueError


def scap_to_png(scap_path, inkscape_path, imagemagick_path, png_path, padding=None):
    svg_path = os.path.dirname(
        png_path) + '/' + os.path.basename(scap_path).replace('.scap', '.svg')
    if not os.path.exists(svg_path):
        cmd = [sys.executable, VIS_PY, scap_path, svg_path]
        subprocess.run(cmd)

    inkscape_args = [inkscape_path, svg_path, "--export-filename=" + png_path,
                     "--export-background=white", '--export-background-opacity=255']
    # Inkscape 0.92:
    # inkscape_args = [inkscape_path, "--without-gui", "--file=" + join(os.getcwd(),svg), "--export-png=" + join(os.getcwd(), save_dir, name+'.png'),  "--export-background=white"]
    dpi = int(96.0)
    subprocess.run(inkscape_args + ["--export-dpi=%d" % dpi])

    # Pad generated png
    if padding != None:
        padding_trbl = padding
        if not isinstance(padding, list):
            padding_trbl = [padding, padding, padding, padding]

        width, height = get_sketch_size(svg_path)
        width2 = int(width + padding_trbl[1] + padding_trbl[3])
        height2 = int(height + padding_trbl[0] + padding_trbl[2])

        output_png = png_path.replace('.png', '2.png')
        output_png = png_path
        cmd = [imagemagick_path, f'{os.path.abspath(png_path)}', '-background', 'white', '-extent',
               f'{width2}x{height2}-{padding_trbl[3]}-{padding_trbl[0]}', f'{os.path.abspath(output_png)}']

        subprocess.run(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('scap', help='scap input')
    parser.add_argument('png', help='png input')
    parser.add_argument('--ink', required=True, help='inkscape path')
    parser.add_argument('--img', required=True, help='imagemagick path')
    parser.add_argument('--padding', help='padding')
    parser.add_argument('--debug', action='store_true',
                        help='place a breakpoint before execution')
    args = parser.parse_args()

    if args.debug:
        # Give time to attach a debugger.
        breakpoint()

    padding = None
    if args.padding:
        padding = [int(s) for s in args.padding.split(',')]

    scap_to_png(args.scap, args.ink, args.img, args.png, padding)
