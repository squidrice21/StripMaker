#!/usr/bin/env python

"""
Transfer a scap (translation and scaling).
"""

import sys
import math
import os

try:
    import scapio
except ImportError:
    sys.path.insert(0, os.path.abspath(
        os.path.split(os.path.abspath(__file__))[0]+'/../'))
    import visualization.scapio as scapio

translation = (0, 0)
scale = 1
scap_file = ''
transferred_file = '-'


def parse_args():
    global translation, scale, scap_file, transferred_file

    error = False
    for i, arg in enumerate(sys.argv):
        if '-t' == arg:
            if i + 2 >= len(sys.argv):
                error = True
                break
            try:
                translation = (float(sys.argv[i + 1]), float(sys.argv[i + 2]))
            except ValueError:
                error = True
                break
        elif '-s' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            try:
                scale = float(sys.argv[i + 1])
            except ValueError:
                error = True
                break
        elif '-' == arg or '.scap' in arg:
            if scap_file == '':
                scap_file = arg
            else:
                transferred_file = arg

    if scap_file == '':
        scap_file = '-'

    if error:
        sys.exit("Usage:\n\tcanvas-calibrate.py [-n canvas.png [_widthx_height] [margin] | \
            -c to_calibrate.png [calibration.cal] [_widthx_height] [margin]]")


def scale_capture(capture, scale):
    scaled = []
    for sketch in capture:
        ss = scapio.Stroke(thickness=sketch.thickness)
        ss.stroke_ind = sketch.stroke_ind
        ss.group_ind = sketch.group_ind

        for p in sketch:
            ss.append((p[0] * scale, p[1] * scale, p[2]))
        scaled.append(ss)

    return scaled


def translate_capture(capture, translation):
    translated = []
    for sketch in capture:
        ts = scapio.Stroke(thickness=sketch.thickness)
        ts.stroke_ind = sketch.stroke_ind
        ts.group_ind = sketch.group_ind

        for p in sketch:
            ts.append((p[0] + translation[0], p[1] + translation[1], p[2]))
        translated.append(ts)

    return translated


def reflect_capture(capture, x_ref, y_ref):
    reflected = []
    for sketch in capture:
        ts = scapio.Stroke(thickness=sketch.thickness)
        ts.stroke_ind = sketch.stroke_ind
        ts.group_ind = sketch.group_ind

        for p in sketch:
            ts.append((p[0] * x_ref, p[1] * y_ref, p[2]))
        reflected.append(ts)

    return reflected


def rotate_capture(capture, angle_deg):
    rotated = []
    angle_rad = angle_deg / 180. * math.pi
    r_mat = [math.cos(angle_rad), -math.sin(angle_rad),
             math.sin(angle_rad), math.cos(angle_rad)]

    for sketch in capture:
        ts = scapio.Stroke(thickness=sketch.thickness)
        ts.stroke_ind = sketch.stroke_ind
        ts.group_ind = sketch.group_ind

        for p in sketch:
            ts.append((p[0] * r_mat[0] + p[1] * r_mat[1], p[0]
                      * r_mat[2] + p[1] * r_mat[3], p[2]))
        rotated.append(ts)

    return rotated


if __name__ == '__main__':
    parse_args()

    if scap_file == '-':
        scap_file = sys.stdin
    else:
        scap_file = open(scap_file, 'r')
    if transferred_file == '-':
        transferred_file = sys.stdout
    else:
        transferred_file = open(transferred_file, 'w')

    print('Transferring... %f, %f - %f.' %
          (scale, translation[0], translation[1]))

    (capture, width, height) = scapio.read_scap(scap_file.read())
    capture = translate_capture(scale_capture(capture, scale), translation)

    transferred_file.write(scapio.scap_to_string(capture, width, height))

    if scap_file != sys.stdin:
        scap_file.close()
    if transferred_file != sys.stdout:
        transferred_file.close()
