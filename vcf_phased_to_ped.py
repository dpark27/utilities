#!/usr/bin/python -u 

import os
import sys
import optparse
import re
from itertools import tee, izip


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--vcf_phased_file', dest='vcf_phased_file', action='store', type='string', default='hapmap.phased')
    parser.add_option('-o', '--output_ped_file', dest='output_ped_file', action='store', type='string', default='hapmap.ped')

    (options, args) = parser.parse_args()

    idv_data, gt_data, map_data = parse_vcf(options.vcf_phased_file)
    output_ped_file(options.output_ped_file, idv_data, gt_data, map_data)


def output_ped_file(output_ped_file, idv_data, gt_data, map_data):
    col_indexes = sorted(idv_data.keys())
    row_indexes = sorted(map_data.keys())

    output_map_file = output_ped_file.replace('.ped', '.map')
    o_file = open(output_map_file, 'wb')
    for row_idx in row_indexes:
        row_to_write = map_data[row_idx]
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


def parse_vcf(vcf_phased_file):
    idv_data = {}
    gt_data = {}
    map_data = {}

    row_idx = 0
    i_file = open(vcf_phased_file, 'rb')
    for line in i_file:
        if line.startswith('##'):
            continue
        else:
            data = line.strip().split()
            col_idx = 0
            if row_idx == 0:
                iids = data[9:]
                for iid in iids:
                    idv_data[col_idx] = [iid, iid, 0, 0, 0, 0]

                    col_idx += 1
            else:
                chrom = data[0]
                pos = data[1]
                rsid = data[2]
                map_data[row_idx] = [chrom, rsid, 0, pos]

                ref = data[3]
                alt = data[4]
                gts = data[9:]
                for gt_d in gts:
                    gt_d = gt_d.split(':')[0]
                    gts = gt_d.split('|')
                    allele1 = gts[0].replace('0', ref).replace('1', alt)
                    allele2 = gts[1].replace('0', ref).replace('1', alt)

                    if col_idx in gt_data:
                        gt_data[col_idx].append(allele1)
                        gt_data[col_idx].append(allele2)
                    else:
                        gt_data[col_idx] = []
                        gt_data[col_idx].append(allele1)
                        gt_data[col_idx].append(allele2)

                    col_idx += 1

        row_idx += 1

    i_file.close()

    return idv_data, gt_data, map_data


if __name__ == '__main__':
    main()
