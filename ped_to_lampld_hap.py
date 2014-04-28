#!/usr/bin/python

import sys
import subprocess
import optparse
from itertools import tee, izip


def main():
    parser = optparse.OptionParser()
    parser.add_option('-m', '--ped_file', dest='ped_file', action='store', type='string', default=None)
    parser.add_option('-p', '--map_file', dest='map_file', action='store', type='string', default=None)
    parser.add_option('-o', '--output_file', dest='output_file', action='store', type='string', default='1.hap')

    (options, args) = parser.parse_args()

    snp_map = build_snp_map(options.map_file)
    build_haplotypes(options.ped_file, snp_map, options.output_file)


def build_haplotypes(ped_file, snp_map, output_file):
    o_file = open(output_file, 'wb')

    i_file = open(ped_file, 'rb')
    for line in i_file:
        data = line.strip().split()

        hap1 = ''
        hap2 = ''

        gts = data[6:]
        gt_idx = 0
        for gt in grouped(gts, 2): 
            gt1 = gt[0]
            gt2 = gt[1]

            snp_data = snp_map[gt_idx]
            ref = snp_data[3]

            if gt1 == '0':
                hap1 += '?'
            elif gt1 == ref:
                hap1 += '1'
            else:
                hap1 += '0'

            if gt2 == '0':
                hap2 += '?'
            elif gt2 == ref:
                hap2 += '1'
            else:
                hap2 += '0'

            gt_idx += 1

        o_file.write(hap1 + '\n')
        o_file.write(hap2 + '\n')

    i_file.close()
    o_file.close()


def build_snp_map(map_file):
    snp_map = {}

    line_idx = 0
    i_file = open(map_file, 'rb')
    for line in i_file:
        data = filter(None, line.strip().split())

        rsid = data[0]
        chrom = data[1]
        pos = data[3]
        ref = data[4]
        alt = data[5]

        snp_map[line_idx] = [rsid, chrom, pos, ref, alt]

        line_idx += 1

    i_file.close()

    return snp_map


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


if __name__ == '__main__':
    main()

