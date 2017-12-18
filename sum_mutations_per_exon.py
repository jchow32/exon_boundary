# sum_mutations_per_exon.py
# determine the total number of missense OR synonymous mutations in an exon given transcript, position, exon #
# ENST00000000233.5 127228572 1
# ENST00000000233.5 127228586 1
# ENST00000000233.5 127228606 1

import sys
import getopt
import csv


def main(argv):
    input_file = ''
    output_file = ''
    previous_exon = 0
    previous_transcript = ""
    exon_sum = 0
    record_transcript_exon = ""
    transcript_name = ""

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python sum_mutations_per_exon.py -i <input_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python sum_mutations_per_exon.py -i <input_file> -o <output_file>"
            sys.exit()
        elif opt in ("-i", "--inputfile"):
            input_file = arg
        elif opt in ("-o", "--outputfile"):
            output_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()
    if output_file == '':
        print "No output file specified."
        sys.exit()

    with open(input_file, 'r') as t_file:
        transcript = csv.reader(t_file, delimiter=" ")

        for line in transcript:
            transcript_name = line[0]
            if line[2] == previous_exon and transcript_name == previous_transcript:  # same exon, same transcript
                exon_sum += 1
            elif line[2] != previous_exon and transcript_name == previous_transcript:  # same transcript, but different exon
                # immediately record the previous exon
                record_transcript_exon += "%s\t%s\t%s\n" % (transcript_name, previous_exon, exon_sum)
                exon_sum = 1  # reset the sum count to 1
                previous_exon = line[2]  # previous will now be current
                previous_transcript = line[0]
            else:  # different exon, different transcript
                # immediately record the last exon of previous transcript
                if previous_exon != 0:  # don't write before very first exon done counting
                    record_transcript_exon += "%s\t%s\t%s\n" % (previous_transcript, previous_exon, exon_sum)
                # what if there are no mutations in a given exon

                exon_sum = 1
                previous_exon = line[2]
                previous_transcript = line[0]

        # don't forget to write the last exon of the last transcript
        record_transcript_exon += "%s\t%s\t%s\n" % (transcript_name, previous_exon, exon_sum)

    with open(output_file, 'w') as write_file:
        write_file.write(record_transcript_exon)


if __name__ == "__main__":
    main(sys.argv[1:])
