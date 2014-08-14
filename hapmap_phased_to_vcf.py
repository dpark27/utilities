#!/usr/bin/python -u 

import os
import sys
import optparse
import re
import numpy as np
from itertools import izip


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--hapmap_phased_file', dest='hapmap_phased_file', action='store', type='string', default='hapmap.phased')
    parser.add_option('-b', '--bim_file', dest='bim_file', action='store', type='string', default='test.bim')
    parser.add_option('-o', '--output_vcf_file', dest='output_vcf_file', action='store', type='string', default='hapmap.ped')

    (options, args) = parser.parse_args()

    # get chromosome from hapmap filename
    g = re.search('chr\d+', options.hapmap_phased_file)
    chrom = g.group(0).replace('chr', '')

    rsid_positions, rsid_genotypes, idvs = parse_phased_file(options.hapmap_phased_file)
    # rsid_ref_alt_map = create_ref_alt_map(options.bim_file)
    output_vcf_file(idvs, rsid_positions, rsid_genotypes, chrom, options.output_vcf_file)


def output_vcf_file(idvs, rsid_positions, rsid_genotypes, chrom, output_vcf_file):
    o_file = open(output_vcf_file, 'wb')
    o_file.write('##fileformat=VCFv4.0' + '\n')
    o_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">' + '\n')
    header_row = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + map(lambda l: l.replace('_A', ''), idvs[::2])
    print len(header_row)
    raw_input()
    header_row = '\t'.join(header_row) + '\n'
    o_file.write(header_row)
    sorted_positons = sorted(rsid_positions.keys())
    for pos in sorted_positons:
        genotype_row = []
        rsid = rsid_positions[pos]
        genotype_row += [chrom, pos, rsid]
        genotypes = rsid_genotypes[rsid]
        uniques = np.unique(genotypes)
        ref = uniques[0]
        if uniques.shape[0] > 1:
            alt = uniques[1]
        else:
            alt = '.'
        gt_map = {}
        gt_map[ref] = 0
        gt_map[alt] = 1
        genotype_row += [ref, alt, '.', 'PASS', '.', 'GT']
        for idv_gt in grouped(genotypes, 2):
            genotype_row += ['|'.join(map(str, [gt_map[idv_gt[0]], gt_map[idv_gt[1]]]))]
        print len(genotype_row)
        raw_input()
        genotype_row = '\t'.join(genotype_row) + '\n'
        o_file.write(genotype_row)
    o_file.close()


def parse_phased_file(hapmap_phased_file):
    i_file = open(hapmap_phased_file, 'rb')
    header = next(i_file).strip().split()
    idvs = header[2:]
    rsid_positions = {}
    rsid_genotypes = {}
    for line in i_file:
        data = line.strip().split()
        rsid = data[0]
        pos = data[1]
        rsid_positions[pos] = rsid
        rsid_genotypes[rsid] = data[2:]
    i_file.close()
    return rsid_positions, rsid_genotypes, idvs


# def create_ref_alt_map(bim_file):
#     rsid_ref_alt_map = {}
#     i_file = open(bim_file, 'rb')
#     for line in bim_file:
#         data = line.strip().split()
#         rsid = data[1]
#         ref = data[4]
#         alt = data[5]
#         rsid_ref_alt_map[rsid] = [ref, alt]
#     i_file.close()
#     return rsid_ref_alt_map


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)


if __name__ == '__main__':
    main()
