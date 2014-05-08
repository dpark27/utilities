#!/usr/bin/python

import os
import sys
import optparse


def main():
    parser = optparse.OptionParser()
    parser.add_option('-m', '--reference_bim_file', dest='reference_bim_file', action='store', type='string', default=None)
    parser.add_option('-p', '--bim_file_to_fix', dest='bim_file_to_fix', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    rsid_allele_map = build_rsid_allele_map(options.reference_bim_file)
    fix_bim_file(options.bim_file_to_fix, rsid_allele_map)


def fix_bim_file(bim_file_to_fix, rsid_allele_map):
    o_file = open(bim_file_to_fix + '.fixed', 'wb')

    i_file = open(bim_file_to_fix, 'rb')
    for line in i_file:
        data = line.strip().split()
        rsid = data[1]

        if data[4] != '0' and data[5] != '0':
            data[4] = rsid_allele_map[rsid][data[4]]
            data[5] = rsid_allele_map[rsid][data[5]]
        else:
            if data[5] != '0':
                if data[5] == 'A' or data[5] == 'T':
                    data[4] = rsid_allele_map[rsid]['C']
                    data[5] = rsid_allele_map[rsid][data[5]]
                else:
                    data[4] = rsid_allele_map[rsid]['A']
                    data[5] = rsid_allele_map[rsid][data[5]]               
            else:
                if data[4] == 'A' or data[4] == 'T':
                    data[5] = rsid_allele_map[rsid]['C']
                    data[4] = rsid_allele_map[rsid][data[4]]
                else:
                    data[5] = rsid_allele_map[rsid]['A']
                    data[4] = rsid_allele_map[rsid][data[4]] 

        row_to_write = '\t'.join(map(str, data)) + '\n'
        o_file.write(row_to_write)

    i_file.close()
    o_file.close()

    os.remove(bim_file_to_fix)
    os.rename(bim_file_to_fix + '.fixed', bim_file_to_fix)


def build_rsid_allele_map(reference_bim_file):
    rsid_allele_map = {}

    i_file = open(reference_bim_file, 'rb')
    for line in i_file:
        data = line.strip().split()

        chrom = data[0]
        rsid = data[1]
        ref = data[4]
        alt = data[5]

        rsid_allele_map[rsid] = {}
        if ref == 'A' or ref == 'T':    
            rsid_allele_map[rsid]['A'] = ref
            rsid_allele_map[rsid]['T'] = ref
            rsid_allele_map[rsid]['C'] = alt
            rsid_allele_map[rsid]['G'] = alt
        else:
            rsid_allele_map[rsid]['C'] = ref
            rsid_allele_map[rsid]['G'] = ref
            rsid_allele_map[rsid]['A'] = alt
            rsid_allele_map[rsid]['T'] = alt

    i_file.close()

    return rsid_allele_map


if __name__ == '__main__':
    main()
