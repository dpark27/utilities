#!/usr/bin/python -u

import os
import sys
import optparse
import linecache
import numpy as np
import math
import scipy.stats
from itertools import izip


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', '--lampld_interval_file', dest='lampld_interval_file', action='store', type='string', default='.')
    parser.add_option('-p', '--ped_file', dest='ped_file', action='store', type='string', default='.')
    parser.add_option('-d', '--spouse_list', dest='spouse_list', action='store', type='string', default='.')
    parser.add_option('-m', '--freq_file', dest='freq_file', action='store', type='string', default='.')
    parser.add_option('-o', '--output_file', dest='output_file', action='store', type='string', default='.')

    (options, args) = parser.parse_args()

    spouses = load_spouses(options.spouse_list)
    idv_line_map = get_idv_line_map(options.ped_file)
    idv_anc = load_idv_anc(options.lampld_interval_file)
    snp_order = load_snp_order(options.ped_file)
    anc_spec_means = load_anc_specific_means(options.freq_file)
    final_stats_by_snp_index = perform_analysis_on_spouses(spouses, idv_line_map, anc_spec_means, idv_anc, options.ped_file, snp_order)
    output_final_stats(final_stats_by_snp_index, options.output_file)


def output_final_stats(final_stats_by_snp_index, output_file):
    o_file = open(output_file, 'wb')
    for index in final_stats_by_snp_index:
        row_to_write = [index] + final_stats_by_snp_index[index]
        row_to_write = '\t'.join(map(str, row_to_write)) + '\n'

        o_file.write(row_to_write)

    o_file.close()


def perform_analysis_on_spouses(spouses, idv_line_map, anc_spec_means, idv_anc, ped_file, snp_order):
    aggregated_stats_by_index = {}
    for pat_id in spouses:
        mat_id = spouses[pat_id]

        if (pat_id not in idv_line_map) or (mat_id not in idv_line_map):
            print "here"
            continue

        pat_genotypes = linecache.getline(ped_file, idv_line_map[pat_id]).strip().split(' ')[6:]
        print linecache.getline(ped_file, idv_line_map[pat_id]).strip().split(' ')[1:6]
        pat_anc_map = get_idv_anc_map(idv_anc, pat_id)

        mat_genotypes = linecache.getline(ped_file, idv_line_map[mat_id]).strip().split(' ')[6:]
        mat_anc_map = get_idv_anc_map(idv_anc, mat_id)
        print linecache.getline(ped_file, idv_line_map[mat_id]).strip().split(' ')[1:6]

        index = 0
        for i, j in pairwise(zip(pat_genotypes, mat_genotypes)):
            # "i=(pg1, mg1), j=(pg2, mg2)"
            if i[0] == '0' or j[0] == '0' or i[1] == '0' or j[1] == '0':
                index += 1
                continue

            pat_freq = 0
            if i[0] == '1':
                pat_freq += 1
            if j[0] == '1':
                pat_freq += 1

            rsid = snp_order[index]

            pat_anc = map(int, pat_anc_map[index])            
            pat_anc_mean1 = anc_spec_means[rsid][pat_anc[0]]
            pat_anc_mean2 = anc_spec_means[rsid][pat_anc[1]]
            pat_anc_mean = (pat_anc_mean1 + pat_anc_mean2) / 2.0
            if pat_anc_mean < 0.05 or pat_anc_mean > 0.95:
                if index not in aggregated_stats_by_index:
                    aggregated_stats_by_index[index] = {}
                    aggregated_stats_by_index[index]["mat_stats"] = []
                    aggregated_stats_by_index[index]["pat_stats"] = []
                index += 1
                continue

            pat_stat = (pat_freq - (2.0 * pat_anc_mean)) / np.sqrt((2.0 * pat_anc_mean) * (1.0 - pat_anc_mean))

            mat_freq = 0
            if i[1] == '1':
                mat_freq += 1
            if j[1] == '1':
                mat_freq += 1

            mat_anc = map(int, mat_anc_map[index])            
            mat_anc_mean1 = anc_spec_means[rsid][mat_anc[0]]
            mat_anc_mean2 = anc_spec_means[rsid][mat_anc[1]]
            mat_anc_mean = (mat_anc_mean1 + mat_anc_mean2) / 2.0
            if mat_anc_mean < 0.05 or mat_anc_mean > 0.95:
                if index not in aggregated_stats_by_index:
                    aggregated_stats_by_index[index] = {}
                    aggregated_stats_by_index[index]["mat_stats"] = []
                    aggregated_stats_by_index[index]["pat_stats"] = []
                index += 1
                continue

            mat_stat = (mat_freq - (2.0 * mat_anc_mean)) / np.sqrt((2.0 * mat_anc_mean) * (1.0 - mat_anc_mean))

            if index not in aggregated_stats_by_index:
                aggregated_stats_by_index[index] = {}
                aggregated_stats_by_index[index]["mat_stats"] = []
                aggregated_stats_by_index[index]["pat_stats"] = []

            aggregated_stats_by_index[index]["mat_stats"].append(mat_stat)
            aggregated_stats_by_index[index]["pat_stats"].append(pat_stat)

            index += 1

    final_stats_by_snp_index = {}
    for index in aggregated_stats_by_index:
        mat_stats = aggregated_stats_by_index[index]["mat_stats"]
        pat_stats = aggregated_stats_by_index[index]["pat_stats"]
        rsid = snp_order[index]
        corr = scipy.stats.pearsonr(pat_stats, mat_stats)[0]
        final_stat = (corr**2) * (len(mat_stats))
        if not np.isnan(final_stat):
            final_stats_by_snp_index[index] = [rsid, final_stat, corr, len(mat_stats)]


    return final_stats_by_snp_index


def load_snp_order(ped_file):
    snp_order = {}
    idx = 0
    i_file = open(ped_file.replace('.ped', '.map'), 'rb')
    for line in i_file:
        data = line.strip().split('\t')
        snp_order[idx] = data[1] # this is the rsid
        idx += 1

    i_file.close()

    return snp_order


def get_idv_anc_map(idv_anc, idv):
    idv_anc_map = {}

    idv_anc_data = idv_anc[idv]
    for l in idv_anc_data:
        anc = l.split(':')[0]

        start = int(l.split(':')[1].split('-')[0])
        stop = int(l.split(':')[1].split('-')[1])
        for i in xrange(start, stop+1):
            idv_anc_map[i] = anc

    return idv_anc_map


def load_idv_anc(lampld_interval_file):
    idv_anc = {}

    i_file = open(lampld_interval_file, 'rb')
    for line in i_file:
        data = line.strip().split('\t')

        idv = data[0]
        ancestry = data[1:]

        idv_anc[idv] = ancestry

    i_file.close()

    return idv_anc


def load_anc_specific_means(freq_file):
    anc_spec_means = {}

    i_file = open(freq_file, 'rb')

    index = 0
    for line in i_file:
        data = line.strip().split()
       
        rsid = data[0]
        anc0_mean = float(data[1])
        anc1_mean = float(data[2])
        anc2_mean = float(data[3])

        anc_spec_means[rsid] = [anc0_mean, anc1_mean, anc2_mean]

        index += 1

    i_file.close()

    return anc_spec_means


def load_spouses(spouse_list):
    spouses = {}

    i_file = open(spouse_list, 'rb')
    for line in i_file:
        data = line.strip().split(' ')
        pat_id = data[0]
        mat_id = data[1]

        spouses[pat_id] = mat_id

    i_file.close()

    return spouses


def get_idv_line_map(ped_file):
    idv_line_map = {}

    line_index = 1
    i_file = open(ped_file, 'rb')
    for line in i_file:
        data = line.strip().split(' ')

        idv_id = data[1]
        idv_line_map[idv_id] = line_index

        line_index += 1

    i_file.close()

    return idv_line_map


def pairwise(iterable):
    "s -> (s0,s1), (s2,s3), (s4, s5), ..."
    a = iter(iterable)

    return izip(a, a)


if __name__ == '__main__':
    main()



