#!/usr/bin/python -u 

import os
import sys
import optparse
import re
from itertools import tee, izip


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--hapmap_phased_file', dest='hapmap_phased_file', action='store', type='string', default='hapmap.phased')
    parser.add_option('-o', '--output_ped_file', dest='output_ped_file', action='store', type='string', default='hapmap.ped')

    (options, args) = parser.parse_args()

    # get chromosome from hapmap filename
    g = re.search('chr\d+', options.hapmap_phased_file)
    chrom = g.group(0).replace('chr', '')

    idv_data, gt_data, map_data = parse_phased_file(options.hapmap_phased_file)
    output_ped_file(options.output_ped_file, chrom, idv_data, gt_data, map_data)


def output_ped_file(output_ped_file, chrom, idv_data, gt_data, map_data):
    col_indexes = sorted(idv_data.keys())
    row_indexes = sorted(map_data.keys())

    map_file = output_ped_file.replace('.ped', 'map')
    o_file = open(map_file, 'wb')
    for row_idx in row_indexes:
        d = map_data[row_idx]
        rsid = d[0]
        pos = d[1]

        row_to_write = [chrom, rsid, 0, pos]
        row_to_write = '\t'.join(map(str, row_to_write)) + '\n'

        o_file.write(row_to_write)

    o_file.close()

    o_file = open(output_ped_file, 'wb')
    for col_idx in col_indexes:
        idv_d = idv_data[col_idx]
        idv_gts = gt_data[col_idx]

        row_to_write = idv_d + idv_gts
        row_to_write = ' '.join(map(str, row_to_write)) + '\n'

        o_file.write(row_to_write)

    o_file.close()


def parse_phased_file(hapmap_phased_file):
    idv_data = {}
    gt_data = {}
    map_data = {}

    row_idx = 0
    i_file = open(hapmap_phased_file, 'rb')
    for line in i_file:
        data = line.strip().split()
        col_idx = 0
        if row_idx == 0:
            for id1, id2 in grouped(data[2:], 2):
                iid = id1.replace('_A', '')
                idv_data[col_idx] = [iid, iid, 0, 0, 0, 0]

                col_idx += 1
        else:
            rsid = data[0]
            pos = data[1]
            map_data[row_idx] = [rsid, pos]

            for allele1, allele2 in grouped(data[2:], 2):
                if col_idx in gt_data:
                    gt_data[col_idx] += allele1 + ' ' + allele2
                else:
                    gt_data[col_idx] = ''
                    gt_data[col_idx] += allele1 + ' ' + allele2

                col_idx += 1


        row_idx += 1

    i_file.close()

    return idv_data, gt_data, map_data


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)


if __name__ == '__main__':
    main()
