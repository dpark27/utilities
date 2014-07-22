#!/usr/bin/python

import sys
import os
import optparse
import numpy as np


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input_dir', dest='input_dir', action='store', type='string', default='')
    parser.add_option('-f', '--file_type', dest='file_type', action='store', type='string', default='.lanc')
    parser.add_option('-t', '--trios', dest='trios', action='store', type='string', default='1')
    parser.add_option('-l', '--idv_list', dest='idv_list', action='store', type='string', default='idvs.txt')
    parser.add_option('-o', '--output_file', dest='output_file', action='store', type='string', default='.')

    (options, args) = parser.parse_args()

    trios = int(options.trios)

    idx_idv_map = load_idv_order(options.idv_list)
    ganc = compute_ganc(options.input_dir, idx_idv_map, options.file_type, trios)
    output_ganc(ganc, idx_idv_map, options.output_file)


def output_ganc(ganc, idx_idv_map, output_file):
    o_file = open(output_file, 'wb')
    for i in ganc:
        idv_data = idx_idv_map[i]
        idv_ganc = ganc[i]
        total = sum(idv_ganc)
        anc0 = idv_ganc[0] / float(total)
        anc1 = idv_ganc[1] / float(total)
        anc2 = idv_ganc[2] / float(total)
        row_to_write = idv_data.split() + [anc0, anc1, anc2]
        row_to_write = ' '.join(map(str, row_to_write)) + '\n'
        o_file.write(row_to_write)
    o_file.close()


def compute_ganc(input_dir, idx_idv_map, file_type=".lanc", trios=1):
    ganc = {}
    for lanc_file in os.listdir(input_dir):
        if lanc_file.endswith(file_type):
            lanc_file_path = input_dir + '/' + lanc_file
            lanc_array = load_lanc(lanc_file_path, trios)
            for i in xrange(0, len(idx_idv_map.keys())):
                row_idx = i * 2
                lanc1 = lanc_array[row_idx,]
                lanc2 = lanc_array[row_idx+1, ]
                anc0 = np.where(lanc1 == 0)[0].shape[0] + np.where(lanc2 == 0)[0].shape[0]
                anc1 = np.where(lanc1 == 1)[0].shape[0] + np.where(lanc2 == 1)[0].shape[0]
                anc2 = np.where(lanc1 == 2)[0].shape[0] + np.where(lanc2 == 2)[0].shape[0]
                if i not in ganc:
                    ganc[i] = [0.0, 0.0, 0.0]
                ganc[i][0] += anc0
                ganc[i][1] += anc1
                ganc[i][2] += anc2
    return ganc


def load_lanc(lanc_file, trios=1):
    hap_array = []
    i_file = open(lanc_file, 'rb')
    for line in i_file:
        data = list(line.strip().replace('?', '3'))
        hap_array.append(data)
    i_file.close()
    hap_array = np.array(hap_array, dtype=np.int8)
    lanc_array = np.empty([1, 1])
    if trios:
        n_rows = (hap_array.shape[0]/4 * 6)
        n_cols = hap_array.shape[1]
        row_idx = 0
        lanc_array = np.empty([n_rows, n_cols])
        for m1, m2, f1, f2, c1, c2 in zip(hap_array[::4], hap_array[1::4], hap_array[2::4], hap_array[3::4], hap_array[::4], hap_array[2::4]):
            lanc_array[row_idx, ] = m1
            row_idx += 1
            lanc_array[row_idx, ] = m2
            row_idx += 1
            lanc_array[row_idx, ] = f1
            row_idx += 1
            lanc_array[row_idx, ] = f2
            row_idx += 1
            lanc_array[row_idx, ] = c1
            row_idx += 1
            lanc_array[row_idx, ] = c2
            row_idx += 1
    else:
        lanc_array = hap_array
    return lanc_array


def load_idv_order(idv_list):
    idx_idv_map = {}
    idv_idx = 0
    i_file = open(idv_list, 'rb')
    for line in i_file:
        data = line.strip()
        idx_idv_map[idv_idx] = data
        idv_idx += 1
    i_file.close()
    return idx_idv_map


if __name__ == '__main__':
    main()
