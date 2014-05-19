#!/usr/bin/python

import sys
import optparse
import subprocess


def main():
    parser = optparse.OptionParser()
    parser.add_option('-m', '--map_file', dest='map_file', action='store', type='string', default=None)
    parser.add_option('-p', '--ped_file', dest='ped_file', action='store', type='string', default=None)
    parser.add_option('-f', '--first_six_file', dest='fsix_file', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    split_map_files = []

    current_chrom = 1
    split_map_output_file = options.ped_file.replace('.ped', str(current_chrom) + '.map')
    split_map_files.append(split_map_output_file)
    o_file = open(split_map_output_file, 'wb')
    i_file = open(options.map_file, 'rb')
    for line in i_file:
        data = line.strip().split('\t')
        chrom = int(data[0])
        if chrom != current_chrom:
            o_file.close()
            current_chrom = chrom
            split_map_output_file = options.ped_file.replace('.ped', str(current_chrom) + '.map')
            split_map_files.append(split_map_output_file)
            o_file = open(split_map_output_file, 'wb')

        o_file.write(line)

    o_file.close()
    i_file.close()

    offset = 6
    for split_map_file in split_map_files:
        i_file = open(split_map_file, 'rb')
        rows = len(i_file.readlines())
        i_file.close()

        first_cut = offset + 1
        last_cut = offset + (2 * rows)

        cut_command = """cut -d ' ' -f """ + str(first_cut) + """-""" + str(last_cut) + """ """ + options.ped_file + """ > """ + split_map_file.replace('.map', '.pped')
        result = subprocess.call(cut_command, shell=True)
        if result != 0:
            print "Problem with cut "
            sys.exit()

        offset = last_cut

    for split_map_file in split_map_files:
        paste_command = """paste -d ' ' """ + options.first_six_file + """ """ + split_map_file.replace('.map', '.pped') + """ > """ + split_map_file.replace('.map', '.ped')
        result = subprocess.call(paste_command, shell=True)
        if result != 0:
            print "Problem with paste"
            sys.exit()

    subprocess.call("rm " + options.ped_file.replace('.ped', '') + "*.pped", shell=True)


if __name__ == '__main__':
    main()
