#!/usr/bin/python

import sys
import optparse
import numpy as np
import pandas as pd


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input_dir', dest='input_dir', action='store', type='string', default='')
    parser.add_option('-ftype', '--file_type', dest='file_type', action='store', type='string', default='.lanc')
    parser.add_option('-t', '--trios', dest='trios', action='store', type='string', default='1')

    (options, args) = parser.parse_args()


def compute_ganc(input_dir, ftype=".lanc", trios=1, idx_idv_map):
    ganc = {}
    for lanc_file in os.listdir(input_dir):
        if lanc_file.endswith(ftype):
            lanc_file_path = input_dir + '/' + lanc_file
            lanc_array = load_lanc(lanc_file_path, trios)
            for i in xrange(0, len(idx_idv_map.keys())):
                idv_lanc = 


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
    lanc_array = pd.DataFrame(lanc_array)
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
