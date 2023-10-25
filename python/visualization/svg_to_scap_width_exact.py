#!/usr/bin/env python

"""
Convert a SVG file into a scap (a set of polyline strokes). 
It resamples non polyline paths and keeps the polylines the same.
"""

import sys
from svgpathtools import Path, Line, Arc, svg2paths, wsvg
import xml.etree.ElementTree as ET

from visualization.transfer_scap import scale_capture
import visualization.scapio as scapio

sample_rate_mm = 5
path_filename = []
translation = (0, 0)
out_filename = '-'
out_scale = 1
y_flip = 1

to_matlab = False

thickness_greater = 0
thickness_less = 1000

max_dimension = -1


def parse_args():
    global path_filename, sample_rate_mm, \
        translation, out_filename, out_scale, y_flip, \
        thickness_greater, thickness_less, \
        to_matlab, max_dimension

    error = False
    for i, arg in enumerate(sys.argv):
        if '.svg' in arg:
            path_filename.append(arg)
        elif '-o' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            out_filename = sys.argv[i + 1]
        elif '-r' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            try:
                sample_rate_mm = float(sys.argv[i + 1])
            except ValueError:
                error = True
                break
        elif '-s' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            try:
                out_scale = float(sys.argv[i + 1])
            except ValueError:
                error = True
                break
        elif '-t' == arg:
            if i + 2 >= len(sys.argv):
                error = True
                break
            try:
                translation = (float(sys.argv[i + 1]), float(sys.argv[i + 2]))
            except ValueError:
                error = True
                break
        elif '-f' == arg:
            y_flip = -1
        elif '-g' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            try:
                thickness_greater = float(sys.argv[i + 1])
            except ValueError:
                error = True
                break
        elif '-m' == arg:
            to_matlab = True
        elif '-l' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            try:
                thickness_less = float(sys.argv[i + 1])
            except ValueError:
                error = True
                break
        elif '-d' == arg:
            if i + 1 >= len(sys.argv):
                error = True
                break
            try:
                max_dimension = float(sys.argv[i + 1])
            except ValueError:
                error = True
                break

    if len(path_filename) == 0 or error:
        sys.exit(
            "Usage:\n\tsvg-to-scap.py [-r sample_rate_mm] [-s out_scale] [-t translation_x translation_y] [-g greater_thickness] [-l less_thickness] [-f] path_filename.svg -o out.svg|out.scap")


def get_svg_size(svg_file):
    tree = ET.parse(svg_file)
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


def get_svg_thickness(svg_file):
    _, attributes = svg2paths(svg_file)
    w = -1
    for attr_dict in attributes:
        if 'stroke-width' in attr_dict:
            w = max(w, float(attr_dict['stroke-width']))
    return w


def resample_path(path, sample_rate):
    resampled_path = []

    path_len = [0]
    for p in path:
        path_len.append(path_len[-1] + p.length())

    sample_len = sample_rate * path_len[-1]
    # print(sample_len)

    acc_sample_len = 0
    c = 0
    local_poly = None
    while acc_sample_len < path_len[-1]:
        large_ind = [ind for (l, ind) in zip(
            path_len, range(0, len(path_len))) if l > acc_sample_len]

        #print(('\t%d: ' % c) + str(large_ind))

        if len(large_ind) < 1:
            break

        curr_ind = large_ind[0] - 1
        next_ind = large_ind[0]

        local_t = (acc_sample_len - path_len[curr_ind]) / \
            (path_len[next_ind] - path_len[curr_ind])

        if type(path[curr_ind]) is Arc:
            local_poly = path[curr_ind]
            resampled_path.append(
                (local_poly.point(local_t).real, local_poly.point(local_t).imag))
        else:
            local_poly = path[curr_ind].poly()
            resampled_path.append(
                (local_poly(local_t).real, local_poly(local_t).imag))

        acc_sample_len = acc_sample_len + sample_len
        c = c + 1

    if len(resampled_path) == 0:
        if type(path[-1]) is not Arc:
            resampled_path.append(
                (path[0].poly()(0).real, path[0].poly()(0).imag))
        else:
            resampled_path.append(
                (path[0].point(0).real, path[0].point(0).imag))

    if type(path[-1]) is not Arc:
        local_poly = path[-1].poly()
        resampled_path.append((local_poly(1).real, local_poly(1).imag))
    else:
        local_poly = path[-1]
        resampled_path.append(
            (local_poly.point(1).real, local_poly.point(1).imag))

    return resampled_path


def convert_path(path):
    converted_path = []
    for i in range(len(path)):
        if i == 0:
            converted_path.append((path[i].start.real, path[i].start.imag))
        converted_path.append((path[i].end.real, path[i].end.imag))

    return converted_path


def parse_path(paths, attributes, sample_rate_mm):
    resampled_paths = []
    thickness_arr = []
    color_arr = []

    name_map = []

    for i in range(0, len(paths)):
        # It is a stroke not a clip path
        # if 'style' in attributes[i]:
        if True:
            if len(paths[i]) == 0 or paths[i].length() == 0:
                continue

            is_opaque = True

            # filter opacity if exists

            if 'style' in attributes[i]:
                attr_str = attributes[i]['style'].split(';')
                attr_dict = {}
                for str in attr_str:
                    temp_str = str.replace(' ', '').split(':')
                    attr_dict[temp_str[0]] = temp_str[1]
            else:
                attr_dict = attributes[i]

            # if 'stroke-opacity' in attr_dict:
            #    if float(attr_dict['stroke-opacity'] ) < 0.4:
            #        is_opaque = False
            #        if float(attr_dict['stroke-width'] ) > 0.35:
            #            is_opaque = True

            # if not is_opaque or len(paths[i]) == 0 or paths[i].length() == 0:
            # if len(paths[i]) == 0 or paths[i].length() < 5:
            #    continue

            if type(paths[i][0]) is not Line:
                resampled_paths.append(resample_path(
                    paths[i], sample_rate_mm / paths[i].length()))
            else:
                resampled_paths.append(convert_path(paths[i]))
            if 'stroke-width' in attr_dict:
                thickness_arr.append(
                    float(attr_dict['stroke-width'].replace('px', '')))
            elif 'stroke-width' in attributes[i]:
                thickness_arr.append(
                    float(attributes[i]['stroke-width'].replace('px', '')))
            else:
                thickness_arr.append(1)

            if 'class' in attr_dict:
                color_arr.append(attr_dict['class'])
            elif 'class' in attributes[i]:
                color_arr.append(attributes[i]['class'])
            else:
                if 'stroke' in attr_dict:
                    color_arr.append(attr_dict['stroke'])
                elif 'stroke' in attributes[i]:
                    color_arr.append(attributes[i]['stroke'])
                else:
                    color_arr.append('cls-empty')

            # print(attributes[i])
            # if 'number' in attributes[i]:
            #    name_map.append(int(attributes[i]['number']))
            if 'id' in attributes[i]:
                try:
                    name_map.append(int(attributes[i]['id'].replace('path', '').replace(
                        'line', '').replace('polygon', '').replace('poly', '').replace('ellipse', '')))
                except ValueError:
                    print('Invalid id str: ' + attr[j]['id'])

            # Too short set an invalid thickness
            if len(paths[i]) == 0:
                # if len(paths[i]) == 0 or paths[i].length() < 5:
                thickness_arr[-1] = -1

    return resampled_paths, thickness_arr, color_arr, name_map


def build_bbox(resampled_paths):
    bbox = [(float('inf'), float('inf')), (-float('inf'), -float('inf'))]
    for i in range(len(resampled_paths)):
        path = resampled_paths[i]

        #print('Thickness: %f' % thickness_arr[i])

        if thickness_arr[i] < 0:
            continue

        if not(thickness_arr[i] >= thickness_greater and
               thickness_arr[i] < thickness_less):
            continue

        for point in path:
            bbox[0] = (min(bbox[0][0], point[0]), min(bbox[0][1], point[1]))
            bbox[1] = (max(bbox[1][0], point[0]), max(
                bbox[1][1], y_flip * point[1]))

    if bbox[0][0] == float('inf'):
        bbox = [(0, 0), (0, 0)]

    return bbox


if __name__ == '__main__':
    parse_args()

    paths = []
    attributes = []
    sid_offset = 1000
    for i, p_name in enumerate(path_filename):
        print('Converting: {}'.format(p_name))
        width, height = get_svg_size(p_name)
        p, attr = svg2paths(p_name)

        for j in range(len(attr)):
            if 'id' in attr[j]:
                try:
                    attr[j]['id'] = str(int(attr[j]['id'].replace('path', '').replace('line', '').replace(
                        'polygon', '').replace('poly', '').replace('ellipse', '')) + sid_offset * i)
                except ValueError:
                    print('Invalid id str: ' + attr[j]['id'])
        paths += p
        attributes += attr

    # filter 0-length path
    long_paths = []
    long_attributes = []
    for i, p in enumerate(paths):
        try:
            if p.length() > 0.5:
                long_paths.append(p)
                long_attributes.append(attributes[i])
        except ZeroDivisionError:
            print('Error path: ' + str(p))

    paths = long_paths
    attributes = long_attributes

    assert len(paths) == len(attributes), 'Paths and attributes don\'t match.'

    resampled_paths, thickness_arr, color_arr, name_map = parse_path(
        paths, attributes, sample_rate_mm)

    # Re-adjust sampling distance
    if max_dimension > 0:
        # Build bbox
        resampled_bbox = build_bbox(resampled_paths)
        print(resampled_bbox)
        sampling_scale = max_dimension / \
            max(resampled_bbox[1][0] - resampled_bbox[0][0],
                resampled_bbox[1][1] - resampled_bbox[0][1])
        print("Readjust sampling rate: {} -> {}".format(sample_rate_mm,
              sample_rate_mm / sampling_scale))
        sample_rate_mm /= sampling_scale
        resampled_paths, thickness_arr, color_arr, name_map = parse_path(
            paths, attributes, sample_rate_mm)

    # print(name_map)

    if '.svg' in out_filename:
        out_paths = []
        out_attributes = []
        style_str = 'fill:none;stroke-width:0.2;stroke-opacity:1;stroke:#000000'

        for path in resampled_paths:
            svg_path = Path()

            for i in range(0, len(path) - 1):
                svg_path.append(Line(complex(path[i][0], path[i][1]),
                                     complex(path[i + 1][0], path[i + 1][1])))

            out_paths.append(svg_path)
            out_attributes.append({'style': style_str})

        wsvg(out_paths, attributes=out_attributes, filename=out_filename)
    else:
        name_file = '-'
        if out_filename == '-':
            out_filename = sys.stdout
        else:
            name_file = open(out_filename.replace('scap', 'name'), 'w')
            offset_file = out_filename.replace('scap', 'offset')
            out_filename = open(out_filename, 'w')

        # Clusters
        clusters = list(set(color_arr))
        clusters = {c: i for i, c in enumerate(clusters)}

        capture = []
        s_ind = 0
        g_ind = 0

        name_str = ''
        for i in range(len(resampled_paths)):
            path = resampled_paths[i]

            #print('Thickness: %f' % thickness_arr[i])

            if thickness_arr[i] < 0:
                continue

            if not(thickness_arr[i] >= thickness_greater and
                   thickness_arr[i] < thickness_less):
                continue

            stroke = scapio.Stroke(thickness=thickness_arr[i])
            stroke.thickness = thickness_arr[i]

            #stroke.stroke_ind = s_ind
            #s_ind = s_ind + 1

            # Use fixed indexing to match the input indexing
            stroke.stroke_ind = i + 1

            # stroke.group_ind = g_ind
            stroke.group_ind = clusters[color_arr[i]]

            for p in path:
                stroke.append((p[0] * out_scale, y_flip * p[1] * out_scale, 0))

            # Use the index in the original svg
            # if i < len(name_map):
            #     stroke.stroke_ind = name_map[i]
            capture.append(stroke)
            if i < len(name_map):
                name_str += '%s\t%s\n' % (name_map[i], str(stroke.stroke_ind))

            #print('%s\t%s\n' % (name_map[i], str(stroke.stroke_ind)))

        bbox = [(float('inf'), float('inf')), (-float('inf'), -float('inf'))]

        for stroke in capture:
            for point in stroke:
                bbox[0] = (min(bbox[0][0], point[0]),
                           min(bbox[0][1], point[1]))
                bbox[1] = (max(bbox[1][0], point[0]),
                           max(bbox[1][1], point[1]))

        # Normalize size
        if max_dimension > 0:
            (_, _, capture) = scapio.fit_capture(capture)
            out_scale = max_dimension / \
                max(bbox[1][0] - bbox[0][0], bbox[1][1] - bbox[0][1])
            capture = scale_capture(capture, out_scale)
            for s in capture:
                s.thickness *= out_scale

        with open(offset_file, 'w') as offset_f:
            offset_f.write('{}\t{}\t{}\t{}\t{}'.format(
                bbox[0][0], bbox[0][1], out_scale, width, height))

        if to_matlab:
            out_filename.write(scapio.matlab_to_string(capture))
        else:
            out_filename.write(scapio.scap_to_string(capture))
            name_file.write(name_str)

        if out_filename != sys.stdout:
            out_filename.close()
            name_file.close()
