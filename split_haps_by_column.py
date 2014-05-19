#!/usr/bin/python

import optparse

def main():
    parser = optparse.OptionParser()
    parser.add_option('-c', '--column_file', dest='column_file', action='store', type='string', default=None)
    parser.add_option('-i', '--input_file', dest='input_file', action='store', type='string', default=None)

    (options, args) = parser.parse_args()

    column_right_boundaries = []
    i_file = open(options.column_file, 'rb')
    for line in i_file:
        boundary = int(line.strip())
        column_right_boundaries.append(boundary)

    i_file.close()

    # open up file handles
    writers = []
    for i in range(0, len(column_right_boundaries)):
        writers.append(open(str(i+1) + '.haps', 'wb'))

    i_file = open(options.input_file, 'rb')
    for line in i_file:
        data = list(line.strip())

        left_boundary = 0
        tmp_right_boundary = 0
        for right_boundary, writer in zip(column_right_boundaries, writers):
            tmp_right_boundary = left_boundary + right_boundary

            if tmp_right_boundary != len(data):
                row_to_write = data[left_boundary:tmp_right_boundary]
                writer.write(''.join(row_to_write) + '\n')
            else:
                row_to_write = data[left_boundary:]
                writer.write(''.join(row_to_write) + '\n')

            print left_boundary, tmp_right_boundary

            left_boundary = tmp_right_boundary

    i_file.close()

    for writer in writers:
        writer.close()


if __name__ == '__main__':
    main()
