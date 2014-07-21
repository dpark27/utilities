#!/usr/bin/python

import optparse
import numpy as np 


def main():
    parser = optparse.OptionParser()
    parser.add_option('-r', '--ref_map_file', dest='ref_map_file', action='store', type='string', default=None)
    parser.add_option('-m', '--ped_map_file', dest='ped_map_file', action='store', type='string', default=None)
    parser.add_option('-p', '--ped_file', dest='ped_file', action='store', type='string', default=None)
    parser.add_option('-o', '--output_prefix', dest='output_prefix', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    snp_map = load_map_file(options.map_file)
    ref_map, ref_snp_order = load_ref_map_file(options.ref_map_file)
    ped_data = load_ped_file(options.ped_file, len(snp_map.keys()))
    output_haps(ref_snp_order, snp_map, ref_map, ped_data, options.output_prefix)


def output_haps(ref_snp_order, snp_map, ref_map, ped_data, output_prefix):
    o_file = open(output_prefix + '.haps', 'wb')
    for row_idx in xrange(0, ped_data.shape[0]):
        hap1 = ped_data[row_idx, ::2]
        hap2 = ped_data[row_idx, 1::2]
        lampld_hap1 = []
        lampld_hap2 = []
        for rsid in ref_snp_order:
            if rsid in snp_map:
                snp_idx = snp_map[rsid]
                allele1 = hap1[snp_idx]
                allele2 = hap2[snp_idx]
                ref, alt = ref_map[rsid]
                if allele1 == ref:
                    lampld_hap1.append('0')
                elif allele1 == alt:
                    lampld_hap1.append('1')
                else:
                    lampld_hap1.append('?')
                if allele2 == ref:
                    lampld_hap2.append('0')
                elif allele2 == alt:
                    lampld_hap2.append('1')
                else:
                    lampld_hap2.append('?')
        row_to_write = ''.join(map(str, lampld_hap1)) + '\n'
        o_file.write(row_to_write)
        row_to_write = ''.join(map(str, lampld_hap2)) + '\n'
        o_file.write(row_to_write)
    o_file.close()


def load_ped_file(ped_file, n_col):
    n_row = 0
    i_file = open(ped_file, 'rb')
    for line in i_file:
        n_row += 1
    i_file.close()
    # load the genotype data
    ped_data = np.empty([n_row, n_col], dtype='S1')
    row_idx = 0
    i_file = open(ped_file, 'rb')
    for line in i_file:
        g_data = np.array(line.strip().split()[6:], dtype='S1')
        ped_data[row_idx,] = g_data
        row_idx += 1
    i_file.close()
    return ped_data


def load_map_file(map_file):
    snp_idx = 0
    snp_map = {}
    i_file = open(map_file, 'rb')
    for line in i_file:
        data = line.strip().split()
        rsid = data[1]
        snp_map[rsid] = snp_idx
        snp_idx += 1
    i_file.close()
    return snp_map


def load_ref_map_file(ref_map_file):
    ref_map = {}
    ref_snp_order = []
    i_file = open(ref_map_file, 'rb')
    for line in i_file:
        data = line.strip().split()
        rsid = data[0]
        ref = data[4]
        alt = data[5]
        ref_map[rsid] = [ref, alt]
        ref_snp_order.append(rsid)
    i_file.close()
    return ref_map, ref_snp_order


if __name__ == '__main__':
    main()

