#!/usr/bin/env python

"""
Helper script to convert a scap file into a SVG.
"""

import svgwrite
import sys
import os
import re
import math
import statistics

try:
    import scapio, transfer_scap
except ImportError:
    sys.path.insert(0, os.path.abspath(
        os.path.split(os.path.abspath(__file__))[0]+'/../'))
    import visualization.scapio as scapio
    import visualization.transfer_scap as transfer_scap

read_scap = scapio.read_scap
rotate_capture = transfer_scap.rotate_capture
reflect_capture = transfer_scap.reflect_capture

if svgwrite.version < (1, 0, 1):
    print("This script requires svgwrite 1.0.1 or newer for internal stylesheets.")
    sys.exit()

infile = '-'
outfile = '-'
offset_file = '-'
to_fit_size = False
to_viz_cut = False
to_viz_short = False
is_closure = False
is_mls = False
vis_shadow = False
viz_short_ref = '-'
board_size = (-1, -1)
CSS_STYLES = """
	.background { fill: white; }
	.line { stroke: firebrick; stroke-width: .1mm; }
	.blacksquare { fill: indigo; }
	.whitesquare { fill: white; }
"""
to_color = False

# A palette with optimally distinct colors
# http://tools.medialab.sciences-po.fr/iwanthue/
palette = ['#c583d9', '#9de993', '#ad5ba7', '#58f2c3', '#c94e76',
           '#57b870', '#d271bb', '#3d8e43', '#f375a0', '#00d4db', '#d55969',
           '#64e4ff', '#e68e54', '#54aeff', '#eed973', '#008bd4', '#fdb968',
           '#6f73bc', '#bee481', '#8f6aaf', '#8bb658', '#bbabff', '#658424',
           '#f1b9ff', '#989c37', '#bb568c', '#cfe395', '#c45372', '#008c6f',
           '#ff94af', '#488656', '#ffc0e1', '#6c8139', '#a9d0ff', '#a98c2c',
           '#0294bb', '#d7b555', '#4a7ea6', '#ffa372', '#96ddff', '#b06735',
           '#bae0f9', '#b76246', '#b0e4e9', '#ff9b9c', '#378673', '#ffa98e',
           '#6a79a2', '#ffcca9', '#53827e', '#ffcdcf', '#648168', '#e7d3fa',
           '#93744c', '#cde1bc', '#9b6a8a', '#efd8b1', '#928098', '#7c7b5c', '#a6686c']
palette = ["#00a391", "#a52a2a", "#ffd700", "#45769e", "#f77ff2", "#ff0000",
           "#db7093", "#1e90ff", "#00ee00", "#ff8c00", "#9046cc", "#6a7b53",
           "#ff1493", "#7fffd4", "#9acd32", "#ba55d3", "#ff00ff", "#d2691e",
           "#458c45", "#954595", "#f0e68c", "#00bfff", "#ffa07a", ]
highlight_palette = ["#00ffff", "#00ff00"]

to_scap = False
highlight_offset = 10000


def eprint(args):
    sys.stderr.write(args + '\n')


def draw_capture(capture, dwg, colormap={}):
    global to_color, palette, highlight_palette, to_viz_cut, \
        to_viz_short, viz_short_ref, vis_shadow

    is_short = float('inf')
    if to_viz_short:
        is_short = viz_short_ref

    out_s_width = 0.5

    for stroke in capture:
        d_str = ''
        s_width = out_s_width
        for point in stroke:
            control_str = ' L '
            if d_str == '':
                control_str = 'M '
            d_str += control_str + '%f %f' % (point[0], point[1])
        color_str = 'black'
        s_width = stroke.thickness

        if to_color:
            color_str = palette[stroke.group_ind % len(palette)]
            if stroke.group_ind >= highlight_offset:
                stroke.group_ind -= highlight_offset
                color_str = highlight_palette[stroke.group_ind % len(
                    highlight_palette)]

            if is_closure and stroke.group_ind == 0:
                color_str = '#000000'

        if stroke.group_ind == -1 or (vis_shadow and stroke[0][2] == 0):
            color_str = '#e0e0e0'
            if not to_color:
                color_str = '#aeaeae'

            #s_width = out_s_width

        if to_viz_short and stroke.length() < is_short:
            color_str = 'crimson'
            s_width = 2

        if stroke.stroke_ind in colormap:
            color_str = colormap[stroke.stroke_ind]
        elif len(colormap) > 0:
            color_str = '#e0e0e0'

        dwg.add(dwg.path(d=d_str, fill='none', stroke=color_str, stroke_width=s_width,
                id=('%d' % stroke.stroke_ind)))
        # dwg.save()

        # if to_viz_cut and len(stroke) > 0:
        #     dwg.add(dwg.circle(center=(
        #         stroke[0][0], stroke[0][1]), r='5px', fill='red', stroke='none', stroke_width=0))
        #     dwg.add(dwg.circle(center=(
        #         stroke[-1][0], stroke[-1][1]), r='5px', fill='red', stroke='none', stroke_width=0))


def draw_cut(capture, dwg):
    endpoints = [s[0] for s in capture if len(
        s) > 0] + [s[-1] for s in capture if len(s) > 0]
    endpoints_count = {}
    str2p = {}
    for p in endpoints:
        p_str = "({:.5f},{:.5f})".format(p[0], p[1])
        if p_str not in endpoints_count:
            endpoints_count[p_str] = 0
        endpoints_count[p_str] += 1
        str2p[p_str] = p

    for s, c in endpoints_count.items():
        if c > 1:
            p = str2p[s]
            dwg.add(dwg.circle(center=(
                p[0], p[1]), r='5px', fill='red', stroke='none', stroke_width=0))


def draw_mls(capture, dwg):
    out_radius = 2
    out_line_w = 1

    dots = []
    lines = []

    for stroke in capture:
        if stroke.stroke_ind == -2:
            dots.append(stroke)
        elif stroke.stroke_ind == -4:
            lines.append(stroke)

    # we want dots on top of lines
    for stroke in lines:
        start = (stroke[0][0], stroke[0][1])
        end = (stroke[1][0], stroke[1][1])
        s_width = out_line_w
        color_str = 'gray'
        dwg.add(dwg.line(start=start, end=end, fill='none',
                stroke=color_str, stroke_width=s_width))

    for stroke in dots:
        center = (stroke[0][0], stroke[0][1])
        s_width = out_radius
        color_str = 'green'

        dwg.add(dwg.circle(center=center, r=s_width, fill=color_str,
                stroke=color_str, stroke_width=s_width))


def parse_args():
    global infile, outfile, board_size, to_fit_size, to_viz_cut, \
        to_viz_short, viz_short_ref, to_color, is_closure, is_mls, \
        offset_file, to_scap, vis_shadow

    for arg in sys.argv:
        if '-s=' in arg:
            to_viz_short = True
            viz_short_ref = arg.replace('-s=', '')
        elif '.scap' in arg:
            infile = arg
        elif '.svg' in arg:
            outfile = arg
        elif '.offset' in arg:
            offset_file = arg
        elif '-f' == arg:
            to_fit_size = True
        elif '-c' == arg:
            to_color = True
        elif '-v' == arg:
            to_viz_cut = True
        elif '-l' == arg:
            is_closure = True
        elif '-m' == arg:
            is_mls = True
        elif '-s' == arg:
            to_scap = True
        elif '-h' == arg:
            vis_shadow = True
        else:
            m = re.match(r'[0-9]+x[0-9]+', arg)

            if m != None:
                size_str = m.group(0)
                board_size = (int(size_str.split('x')[0]), int(
                    size_str.split('x')[1]))

    #print('offset_file: {} - {}'.format(offset_file, sys.argv))

    if outfile == '-':
        sys.exit(
            "Usage:\n\tscap-to-svg.py [-f|-c] [infile.scap] [outfile.svg]\nTo omit infile.scap and use stdio, replace file name with '-'")


def color_svg(capture, width, height, outfile, colormap={}):
    # Fit sketch
    if len(colormap) == 0:
        (width, height, capture) = scapio.fit_capture(capture)
    thickness = 0
    for s in capture:
        thickness = max(thickness, s.thickness)

    width += thickness
    height += thickness

    transformed_capture = []
    for stroke in capture:
        transformed_capture.append(
            scapio.Stroke(thickness=stroke.thickness))
        transformed_capture[-1].stroke_ind = stroke.stroke_ind
        transformed_capture[-1].group_ind = stroke.group_ind

        for point in stroke:
            transformed_capture[-1].append((point[0] + thickness / 2,
                                            (point[1] + thickness / 2), 0))

    capture = transformed_capture

    dwg = svgwrite.Drawing(outfile, size=(width, height))
    dwg.viewbox(0, 0, width, height)

    # always use css for styling
    dwg.defs.add(dwg.style(CSS_STYLES))

    # set background
    #dwg.add(dwg.rect(size=('100%','100%'), class_='background'))

    draw_capture(capture, dwg, colormap=colormap)
    dwg.save()


if __name__ == '__main__':
    parse_args()

    if infile == '-':
        # eprint("Reading scap from stdin.")
        infile = sys.stdin
    else:
        # eprint("Reading scap from '" + infile + "'")
        infile = open(infile, 'rt')

    if to_viz_short:
        eprint("Reading ref scap from '" + viz_short_ref + "'")
        viz_short_ref = open(viz_short_ref, 'rt')
        ref_capture, width, height = read_scap(viz_short_ref.read())

        s_len = []
        for stroke in ref_capture:
            s_len.append(stroke.length())

        sorted(s_len)
        ninty_ind = math.floor(len(s_len) * 0.9)
        is_short = 0.1 * statistics.median(s_len)  # s_len[ninty_ind]
        print('Short threshold: %f' % is_short)
        viz_short_ref.close()
        viz_short_ref = is_short

    assert(outfile != '-')
    # eprint("Writing svg to '" + outfile + "'")

    capture, width, height = read_scap(infile.read())

    if offset_file != '-':
        to_fit_size = False
        out_scale = 1
        translate_x = 0
        translate_y = 0

        with open(offset_file, 'rt') as offset_f:
            offset_str = offset_f.read().replace('\n', '').split('\t')
            # print(offset_str)
            translate_x = float(offset_str[0])
            translate_y = float(offset_str[1])
            if len(offset_str) >= 3:
                out_scale = float(offset_str[2])

        capture = transfer_scap.scale_capture(capture, 1 / out_scale)
        width *= 1 / out_scale
        height *= 1 / out_scale
        width += translate_x
        height += translate_y
        capture = transfer_scap.translate_capture(
            capture, (translate_x, translate_y))
        for s in capture:
            s.thickness /= out_scale

        if len(offset_str) >= 5:
            width = float(offset_str[3])
            height = float(offset_str[4])

    if board_size[0] >= 0 and board_size[1] >= 0:
        width = board_size[0]
        height = board_size[1]
    elif to_fit_size:
        (width, height, capture) = scapio.fit_capture(capture)
        thickness = 0
        for s in capture:
            thickness = max(thickness, s.thickness)

        width += thickness
        height += thickness

        transformed_capture = []
        for stroke in capture:
            transformed_capture.append(
                scapio.Stroke(thickness=stroke.thickness))
            transformed_capture[-1].stroke_ind = stroke.stroke_ind
            transformed_capture[-1].group_ind = stroke.group_ind

            for point in stroke:
                transformed_capture[-1].append((point[0] + thickness / 2,
                                                (point[1] + thickness / 2), 0))

        capture = transformed_capture

    if is_closure:
        capture = reflect_capture(capture, 1, -1)
        capture = rotate_capture(capture, 90)
        (width, height, capture) = scapio.fit_capture(capture)

    dwg = svgwrite.Drawing(outfile, size=(width, height))
    dwg.viewbox(0, 0, width, height)
    # checkerboard has a size of 10cm x 10cm;
    # defining a viewbox with the size of 80x80 means, that a length of 1
    # is 10cm/80 == 0.125cm (which is for now the famous USER UNIT)
    # but I don't have to care about it, I just draw 8x8 squares, each 10x10 USER-UNITS

    # always use css for styling
    dwg.defs.add(dwg.style(CSS_STYLES))

    # set background
    #dwg.add(dwg.rect(size=('100%','100%'), class_='background'))

    if not is_mls:
        draw_capture(capture, dwg)
    else:
        draw_mls(capture, dwg)
    if to_viz_cut:
        draw_cut(capture, dwg)
    dwg.save()

    if to_scap:
        outfile_scap = outfile.replace('.svg', '_s.scap')
        with open(outfile_scap, 'wt') as outf:
            outf.write(scapio.scap_to_string(
                capture, width=width, height=height))

    if infile != sys.stdin:
        infile.close()
