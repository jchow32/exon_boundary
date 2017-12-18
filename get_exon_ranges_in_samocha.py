# get_exon_ranges_in_samocha.py
# to print out the regions that were found to be significantly different from each other according to Samocha,

import sys
import getopt
import csv


def main(argv):
    input_file = ''
    boundary_file = ''
    output_file = ''

    transcript_dictionary = {}
    blank_ends = []
    writing_lines = ""

    try:
        opts, args = getopt.getopt(argv, "hi:b:o:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python get_exon_ranges_in_samocha.py -i <input_file> -b <boundary_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python get_exon_ranges_in_samocha.py -i <input_file> -b <boundary_file> -o <output_file>"
            sys.exit()
        elif opt in ("-i", "--inputfile"):
            input_file = arg
        elif opt in ("-b", "--boundaryfile"):
            boundary_file = arg
        elif opt in ("-o", "--outputfile"):
            output_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()
    if boundary_file == '':
        print "No exon boundary file specified."
        sys.exit()
    if output_file == '':
        print "No output file specified."
        sys.exit()

    with open(boundary_file, 'r') as b_file:
        exon_line = csv.reader(b_file, delimiter=" ")

        for line in exon_line:
            start_end = "%s-%s" % (line[1], line[2])
            if line[0] in transcript_dictionary:
                transcript_dictionary[line[0]].append(start_end)
            else:
                transcript_dictionary[line[0]] = [start_end]

    with open(input_file, 'r') as i_file:
        sig_line = csv.reader(i_file, delimiter="\t")

        for line in sig_line:
            exon_number = 0
            start = line[2]
            end = line[3]
            for start_end_pair in transcript_dictionary[line[0]]:  # looping through array of head-tail, head-tail,...
                exon_number += 1
                head = start_end_pair.split('-')[0]
                tail = start_end_pair.split('-')[1]
                if head <= start <= tail:  # the start belongs in this exon
                    blank_ends.append(exon_number)
                if head <= end <= tail:  # the end belongs in this exon
                    blank_ends.append(exon_number)
            for number in range(blank_ends[0], blank_ends[1]+1):
                writing_lines += '%s_%s\t%s\t%s\t%s\n' % (line[0], number, line[4], line[5], line[6])
            blank_ends = []

    with open(output_file, 'w') as write_file:
        write_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
