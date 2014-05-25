#!/usr/bin/python

import os
import sys
import optparse
import numpy as np


def main():
    parser = optparse.OptionParser()
    parser.add_option('-l', '--lanc_file', dest='lanc_file', action='store', type='string', default=None)
    parser.add_option('-p', '--phased_genos_file', dest='phased_genos_file', action='store', type='string', default=None)
    parser.add_option('-m', '--map_file', dest='map_file', action='store', type='string', default=None)
    parser.add_option('-o', '--output_file', dest='output_file', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    lanc, genos = load_lanc_genos(options.lanc_file, options.phased_genos_file)
    allele_freqs = compute_allele_frequencies(lanc, genos)
    output_freqs(allele_freqs, options.map_file, options.output_file)


def output_freqs(allele_freqs, map_file, output_file):
    rsids = {}
    i_file = open(map_file, 'rb')
    idx = 0
    for line in i_file:
        data = line.strip().split()

        rsid = data[0]

        rsids[idx] = rsid

        idx += 1

    i_file.close()

    o_file = open(output_file, 'wb')
    for i in xrange(0, len(rsids.keys())):
        allele_freq = allele_freqs[i]
        rsid = rsids[i]

        row_to_write = [rsid] + allele_freq
        row_to_write = ' '.join(map(str, row_to_write)) + '\n'

        o_file.write(row_to_write)

    o_file.close()


def load_lanc_genos(lanc_file, phased_genos_file):
    lanc = []
    i_file = open(lanc_file, 'rb')
    for line in i_file:
        lanc.append(list(line.strip()))

    i_file.close()

    lanc = np.array(lanc, dtype=np.uint8)

    genos = []
    i_file = open(phased_genos_file, 'rb')
    for line in i_file:
        genos.append(list(line.strip().replace('?', '3')))

    i_file.close()

    genos = np.array(genos, dtype=np.uint8)

    return lanc, genos


def compute_allele_frequencies(lanc, genos):
    allele_freqs = {}

    # compute at each snp
    n_snps = lanc.shape[1]
    for i in xrange(0, n_snps):
        # determine which is minor allele
        minor_allele = 0
        if np.where(genos[:, i]== 0)[0].shape[0] > np.where(genos[:, i]== 1)[0].shape[0]:
            minor_allele = 1

        eur_snps = np.where(lanc[:, i]==0)
        eur_ref_count = np.where(genos[eur_snps[0], i]==minor_allele)[0].shape[0]
        eur_snps = eur_snps[0].shape[0]

        eur_freq = 0
        if eur_snps != 0:
            eur_freq = float(eur_ref_count) / eur_snps

        nam_snps = np.where(lanc[:, i]==1)
        nam_ref_count = np.where(genos[nam_snps[0], i]==minor_allele)[0].shape[0]
        nam_snps = nam_snps[0].shape[0]
        nam_freq = 0
        if nam_snps != 0:
            nam_freq = float(nam_ref_count) / nam_snps

        afr_snps = np.where(lanc[:, i]==2)
        afr_ref_count = np.where(genos[afr_snps[0], i]==minor_allele)[0].shape[0]
        afr_snps = afr_snps[0].shape[0]
        afr_freq = 0
        if afr_snps != 0:
            afr_freq = float(afr_ref_count) / afr_snps

        allele_freqs[i] = [eur_freq, nam_freq, afr_freq]

    return allele_freqs


if __name__ == '__main__':
    main()
