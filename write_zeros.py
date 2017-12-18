# write_zeros.py
# for a transcript, if there are no mutations in an exon, write in zero.

import sys
import getopt
import csv


def main(argv):
    input_file = ''
    boundary_file = ''
    output_file = ''

    transcript_dictionary = {}
    record_transcript_exon = ""

    try:
        opts, args = getopt.getopt(argv, "hi:b:o:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python write_zeros.py -i <input_file> -b <boundary_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python write_zeros.py -i <input_file> -b <boundary_file> -o <output_file>"
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
        print "No boundary file specified."
        sys.exit()
    if output_file == '':
        print "No output file specified."
        sys.exit()

    with open(boundary_file, 'r') as b_file:
        exon_line = csv.reader(b_file, delimiter=" ")

        for line in exon_line:
            if line[0] not in transcript_dictionary:  # haven't seen this transcript before
                transcript_dictionary[line[0]] = 1
            elif line[0] in transcript_dictionary:  # have seen this transcript before
                transcript_dictionary[line[0]] += 1
                # storing the total number of exons per transcript in dictionary

    previous_transcript = ""
    transcript_name = ""
    exon_range = []

    with open(input_file, 'r') as sum_file:
        sums = csv.reader(sum_file, delimiter="\t")

        for line in sums:
            transcript_name = line[0]
            if previous_transcript == transcript_name:  # still on the same transcript
                exon_range.remove(int(line[1]))
            else:  # onto a new transcript, reset the exon range
                # immediately write in the missing exons
                if exon_range:  # is not empty
                    for item in exon_range:
                        record_transcript_exon += "%s\t%s\t%s\n" % (previous_transcript, item, 0)
                exon_range = range(1, int(transcript_dictionary[transcript_name])+1)
                exon_range.remove(int(line[1]))  # remove this first match in the exon range
                previous_transcript = transcript_name

            record_transcript_exon += "%s\t%s\t%s\n" % (transcript_name, line[1], line[2])

        # don't forget to write the last exons of the last transcript
        for item in exon_range:
            record_transcript_exon += "%s\t%s\t%s\n" % (transcript_name, item, 0)

    with open(output_file, 'w') as write_file:
        write_file.write(record_transcript_exon)


if __name__ == "__main__":
    main(sys.argv[1:])
