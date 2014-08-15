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
    parser.add_option('-g', '--geno_file', dest='geno_file', action='store', type='string', default='.')
    parser.add_option('-c', '--chunk', dest='chunk', action='store', type='string', default='.')
    parser.add_option('-s', '--snp_file', dest='snp_file', action='store', type='string', default='.')
    parser.add_option('-l', '--idv_file', dest='idv_file', action='store', type='string', default='.')
    parser.add_option('-o', '--ped_file_output', dest='ped_file_output', action='store', type='string', default='.')

    (options, args) = parser.parse_args()

    chunk = int(options.chunk)
    chunk_size = 100
    idv_ped_info = get_idv_data(options.idv_file)
    snp_data = load_snp_data(options.snp_file)
    genotypes = get_idv_genotypes(options.geno_file, snp_data, len(idv_ped_info.keys()), chunk, chunk_size)
    output_ped_file(genotypes, idv_ped_info, options.ped_file_output, snp_data, chunk, chunk_size)


def output_ped_file(genotypes, idv_ped_info, ped_file_output, snp_data, chunk, chunk_size):
    chunk_start = ((chunk-1) * chunk_size)
    chunk_end = ((chunk-1) * chunk_size) + chunk_size
    chunk_idxs = range(chunk_start, chunk_end)
    o_file = open(ped_file_output, 'wb')
    for idx in xrange(0, genotypes.shape[1]/2):
        print chunk_idxs[idx], idx
        idv_data = idv_ped_info[chunk_idxs[idx]]
        gt_idxs = [(2 * idx) , ((2 * idx)+1)]
        idv_gts = genotypes[:, gt_idxs]
        row_to_write = idv_data + idv_gts.flatten().tolist()
        row_to_write = ' '.join(map(str, row_to_write)) + '\n'
        o_file.write(row_to_write)
    o_file.close()
    map_file_output = ped_file_output.replace('.ped', '.map')
    o_file = open(map_file_output, 'wb')
    sorted_snp_idx = sorted(snp_data.keys())
    for idx in sorted_snp_idx:
        snp_info = snp_data[idx]
        chrom = snp_info[0]
        rsid = snp_info[1]
        pos = snp_info[2]
        row_to_write = map(str, [chrom, rsid, 0, pos])
        row_to_write = '\t'.join(row_to_write) + '\n'
        o_file.write(row_to_write)
    o_file.close()


def get_idv_genotypes(geno_file, snp_data, n_idvs, chunk, chunk_size):
    genotypes = np.empty([len(snp_data.keys()), 2*chunk_size], dtype='S1')
    chunk_start = ((chunk-1) * chunk_size)
    chunk_end = ((chunk-1) * chunk_size) + chunk_size
    idx = 0
    i_file = open(geno_file, 'rb')
    for line in i_file:
        data = np.array(list(line.strip()), dtype='S1')[chunk_start:chunk_end]
        snp_info = snp_data[idx]
        ref = snp_info[3]
        alt = snp_info[4]
        plink_gts = []
        for i in xrange(0, data.shape[0]):
            gt = data[i]
            plink_gt = []
            if gt == '9':
                plink_gt = ['0', '0']
            elif gt == '0':
                plink_gt = [alt, alt]
            elif gt == '1':
                plink_gt = [alt, ref]
            else:
                plink_gt = [ref, ref]
            plink_gts += plink_gt
        genotypes[idx, ] = plink_gts
        idx += 1
    i_file.close()
    return genotypes


def load_snp_data(snp_file):
    snp_data = {}
    idx = 0
    i_file = open(snp_file, 'rb')
    for line in i_file:
        data = filter(None, line.strip().split())
        rsid = data[0]
        chrom = data[1]
        pos = data[3]
        ref = data[4]
        alt = data[5]
        snp_data[idx] = [chrom, rsid, pos, ref, alt]
        idx += 1
    i_file.close()
    return snp_data


def get_idv_data(idv_file):
    idv_ped_info = {}
    idx = 0
    i_file = open(idv_file, 'rb')
    for line in i_file:
        data = line.strip().split()
        fid = ''
        iid = ''
        if len(data[0].split(':')) > 1:
            fid = data[0].split(':')[0]
            iid = data[0].split(':')[1]
        else:
            fid = data[0].split(':')[0]
            iid = fid
        sex = data[1]
        if sex == 'M':
            sex = 1
        else:
            sex = 2
        idv_data = [fid, iid, 0, 0, sex, 0]
        idv_ped_info[idx] = idv_data
        idx += 1
    i_file.close()
    return idv_ped_info


if __name__ == '__main__':
    main()
