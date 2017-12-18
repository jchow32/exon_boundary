# exon_boundary_from_1_for_S.py
# From exon boundaries from S, convert chromosomal positions to relative positions (starting from the number 1)
# Get 'region' formatted file. Retrieve minimum start position from transcript, subtract this value from
# positions (not exon specific) in S-identified regions of significantly different missense depletion.

import sys
import getopt
import csv


def main(argv):
    input_file = ''
    output_file = ''
    chi_file = ''

    try:
        opts, args = getopt.getopt(argv, "hi:c:o:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python exon_boundary_from_1_for_S.py -i <input_file> -c <chi-squared file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python exon_boundary_from_1_for_S.py -i <input_file> -c <chi-squared file> -o <output_file>"
            sys.exit()
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-c", "--chi"):
            chi_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()
    if chi_file == '':
        print "No input chi-squared file specified."
        sys.exit()
    if output_file == '':
        print "No output file specified."
        sys.exit()

    writing_lines = ""
    previous_transcript = ""
    transcript = ""
    transcript_dict = {}
    transcript_array = []
    start_array = []

    with open(input_file, 'r') as i_file:
        exon_line = csv.reader(i_file, delimiter=" ")

        for line in exon_line:
            transcript = line[0]
            start = int(line[1])
            if transcript not in transcript_dict:  # have not seen this transcript before

                if len(transcript_dict) == 0:  # the very first transcript of the file
                    start_array.append(start)

                    transcript_dict[transcript] = 0  # initiate with something to make len != 0

                else:  # NOT the very first transcript of the file
                    # immediately determine the minimum start of the previous transcript
                    min_start = min(start_array) - 1
                    transcript_dict[previous_transcript] = min_start
                    # reset for the new transcript
                    transcript_dict[transcript] = 0  # initiate with something
                    start_array = [start]
                previous_transcript = transcript

            else:  # have seen this transcript before
                # just keep adding the starts
                start_array.append(start)
                previous_transcript = transcript

        # don't forget about the very last transcript
        min_start = min(start_array) - 1
        transcript_dict[transcript] = min_start

    with open(chi_file, 'r') as c_file:
        chi_line = csv.reader(c_file, delimiter="\t")

        for line in chi_line:
            transcript = line[0]
            gene = line[1]
            min_start = transcript_dict[transcript]  # get the stored minimum from the dictionary
            start = int(line[2]) - min_start
            end = int(line[3]) - min_start

            writing_lines += "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, transcript, start, end, line[4], line[5], line[6])

    with open(output_file, 'w') as out_file:
        out_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
