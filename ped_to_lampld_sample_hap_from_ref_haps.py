#!/usr/bin/python

import optparse
from itertools import izip


def main():
    parser = optparse.OptionParser()
    parser.add_option('-p', '--haps_map_file', dest='haps_map_file', action='store', type='string', default=None)
    parser.add_option('-m', '--ped_map_file', dest='ped_map_file', action='store', type='string', default=None)
    parser.add_option('-d', '--ped_file', dest='ped_file', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    ped_rsid_list = load_ped_map_file(options.ped_map_file)
    haps_rsid_list, rsid_dictionary = load_haps_map_file(options.haps_map_file, ped_rsid_list)
    convert_ped_to_lampld_format(options.ped_file, ped_rsid_list, haps_rsid_list, rsid_dictionary)
    # check_rsid_order(ped_rsid_list, haps_rsid_list)


def check_rsid_order(ped_rsid_list, haps_rsid_list):
    print "checking"
    intersection = []

    for rsid in haps_rsid_list:
        if rsid in ped_rsid_list:
            intersection.append(rsid)

    if set(intersection) != set(ped_rsid_list):
        print "Error with intersection"

    for i in range(0, len(intersection)):
        hap_rsid = intersection[i]
        ped_rsid = ped_rsid_list[i]
        if hap_rsid != ped_rsid:
            print "Problem with ordering"


def convert_ped_to_lampld_format(ped_file, ped_rsid_list, haps_rsid_list, rsid_dictionary):
    print "converting genotype data"
    i_file = open(ped_file, 'rb')
    o_file = open('sample.haps', 'wb')
    for line in i_file:
        data = line.strip().split(' ')
        idv_data = data[0:6]
        print "converting individual " + idv_data[1]

        idv_dict = {}
        for rsid in haps_rsid_list:
            idv_dict[rsid] = '?'

        rsid_index = 0
        genotype_data = data[6:]
        for f, s in grouped(genotype_data, 2):
            rsid = ped_rsid_list[rsid_index]

            genotype = rsid_dictionary[rsid]
            ref = genotype[0]
            # alt = genotype[1]

            reference_allele_count = 0
            if f == ref:
                reference_allele_count += 1
            if s == ref:
                reference_allele_count += 1

            if f == '0' or s == '0':
                reference_allele_count = '?'

            idv_dict[rsid] = str(reference_allele_count)

            rsid_index += 1

        row_to_write = []
        for rsid in haps_rsid_list:
            row_to_write.append(idv_dict[rsid])

        o_file.write(''.join(row_to_write) + '\n')

    i_file.close()
    o_file.close()


def load_haps_map_file(haps_map_file, ped_rsid_list):
    print "loading haps map file"
    rsid_dictionary = {}
    rsid_list = []

    i_file = open(haps_map_file, 'rb')
    for line in i_file:
        data = line.strip().split('\t')

        rsid = data[0]

        # TODO: change this
        ref = data[4]
        alt = data[5]
        # print ref, alt

        rsid_list.append(rsid)

        rsid_dictionary[rsid] = [ref, alt]

    i_file.close()

    return rsid_list, rsid_dictionary


def load_ped_map_file(ped_map_file):
    print "loading ped map file"
    rsid_list = []

    i_file = open(ped_map_file, 'rb')
    for line in i_file:
        data = line.strip().split('\t')

        rsid = data[1]
        rsid_list.append(rsid)

    i_file.close()

    return rsid_list


def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)


if __name__ == '__main__':
    main()

