#!/usr/bin/python

import sys
import os
import optparse
import numpy as np


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input_file', dest = 'input_file', action = 'store', type = 'string', default = None)
    parser.add_option('-o', '--output_file', dest = 'output_file', action = 'store', type = 'string', default  = '')

    (options, args) = parser.parse_args()

    lanc_array = load_lanc_file(options.input_file)
    avgs = compute_lanc_avg(lanc_array)
    output_lanc_avg(avgs, options.output_file)


def output_lanc_avg(avgs, output_file):
    o_file = open(output_file, 'wb')
    for i in sorted(avgs.keys()):
        avg = avgs[i]
        row_to_write = [i] + avg
        row_to_write = ' '.join(row_to_write) + '\n'
        o_file.write(row_to_write)
    o_file.close()


def compute_lanc_avg(lanc_array):
    avgs = {}
    snp_idx = 0
    for col in lanc_array.T:
        anc0 = np.equal(col, 0).sum()
        anc1 = np.equal(col, 1).sum()
        anc2 = np.equal(col, 2).sum()
        total = float(anc0 + anc1 + anc2)
        avgs[snp_idx] = [anc0 / total, anc1 / total, anc2 / total]
        snp_idx += 1

    return avgs

def load_lanc_file(input_file):
    lanc_array = []
    i_file = open(input_file, 'rb')
    for line in i_file:
        data = list(line.strip())
        lanc_array.append(data)
    i_file.close()
    lanc_array = np.array(lanc_array, dtype=np.int8)

    return lanc_array


if __name__ == '__main__':
    main()
