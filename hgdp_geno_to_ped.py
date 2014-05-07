#!/usr/bin/python

import os
import sys
import optparse
import numpy as np


def main():
    parser = optparse.OptionParser()
    parser.add_option('-l', '--hgdp_geno_file', dest='hgdp_geno_file', action='store', type='string', default=None)
    parser.add_option('-m', '--hgdp_map_file', dest='hgdp_map_file', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    parse_geno_file(options.hgdp_geno_file)
    parse_map_file(options.hgdp_map_file)


def parse_map_file(hgdp_map_file):
    o_file = open(hgdp_map_file.replace('.map', '.plink.map'), 'wb')

    i_file = open(hgdp_map_file, 'rb')
    for line in i_file:
        data = line.strip().split()

        rsid = data[0]
        chrom = data[1]
        pos = data[2]

        plink_row = map(str, [chrom, rsid, 0, pos])
        plink_row = '\t'.join(plink_row) + '\n'

        o_file.write(plink_row)

    o_file.close()


def parse_geno_file(hgdp_geno_file):
    geno_matrix = []

    i_file = open('chr1.geno', 'rb')
    iids = i_file.readline().strip().split()
    geno_matrix.append(iids)

    for line in i_file:
        gt_data = line.strip().split()[1:]
        geno_matrix.append(gt_data)

    i_file.close()

    o_file = open(hgdp_geno_file.replace('.geno', '.ped'), 'wb')
    gt_array = np.array(geno_matrix, dtype=np.str)
    vfunc = np.vectorize(create_plink_format_gt)
    for d in gt_array.T:
        idv_data = map(str, [d[0], d[0], 0, 0, 0, 0])
        gt_data = ' '.join(vfunc(d[1:]))

        plink_row = idv_data + ' ' + gt_data + '\n'

        o_file.write(plink_row)

    o_file.close()


def create_plink_format_gt(gt):
    gt = gt[0] + ' ' + gt[1]

    return gt

if __name__ == '__main__':
    main()
