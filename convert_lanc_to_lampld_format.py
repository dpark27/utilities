#!/usr/bin/python 

import os
import sys
import optparse
import numpy as np


def main():
    parser = optparse.OptionParser()
    parser.add_option('-l', '--lanc_file', dest='lanc_file', action='store', type='string', default='.')
    parser.add_option('-o', '--out_file', dest='out_file', action='store', type='string', default='.')

    (options, args) = parser.parse_args()

    convert_to_lampld_out(options.lanc_file, options.out_file)


def convert_to_lampld_out(lanc_file, out_file):
    lanc = []
    i_file = open(lanc_file, 'rb')
    for line in i_file:
        lanc.append(list(line.strip()))

    i_file.close()
    lanc = np.array(lanc, dtype=np.uint8)

    #do this in a loop
    o_file = open(out_file, 'wb')
    for i in xrange(0, lanc.shape[0], 2):
        maternal_breakpoints = np.where(lanc[i, :-1] != lanc[i, 1:])[0]
        paternal_breakpoints = np.where(lanc[i+1, :-1] != lanc[i+1, 1:])[0]

        all_breakpoints = sorted(unique(list(maternal_breakpoints) + list(paternal_breakpoints)))
        idv_breaks = []
        for break_end in all_breakpoints:
            anc1 = lanc[i, break_end]
            anc2 = lanc[i+1, break_end]

            anc = sorted([anc1, anc2])
            anc = ''.join(map(str, anc))

            interval = [anc, ':', break_end]
            interval = ''.join(map(str, interval))

            idv_breaks.append(interval)

        anc1 = lanc[i, lanc.shape[1] - 1]
        anc2 = lanc[i+1, lanc.shape[1] - 1]

        anc = sorted([anc1, anc2])
        anc = ''.join(map(str, anc))

        interval = [anc, ':', lanc.shape[1]]
        interval = ''.join(map(str, interval))

        idv_breaks.append(interval)

        o_file.write(' '.join(idv_breaks) + '\n')


def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]


if __name__ == '__main__':
    main()
