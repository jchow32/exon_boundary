# get_avg_score_per_transcript.py
# calculate the average score over all exons for a transcript

import sys
import getopt
import csv


def main(argv):
    input_file = ''
    output_file = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python get_avg_score_per_transcript.py -i <input_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python get_avg_score_per_transcript.py -i <input_file> -o <output_file>"
            sys.exit()
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()
    if output_file == '':
        print "No output file specified."
        sys.exit()

    writing_lines = ""
    transcript_dict = {}
    non_nan_dict = {}

    with open(input_file, 'r') as i_file:
        g_line = csv.reader(i_file, delimiter="\t")

        for line in g_line:
            transcript = line[0]
            score = line[6]
            if transcript not in transcript_dict:
                transcript_dict[transcript] = 0
                non_nan_dict[transcript] = 0
            if score != "NaN":
                transcript_dict[transcript] += float(score)  # only add non NaN scores
                non_nan_dict[transcript] += 1.0  # keep track of how many non nan scores

        for key in transcript_dict:
            if non_nan_dict[key] != 0:
                avg_score = transcript_dict[key] / non_nan_dict[key]
            else:
                avg_score = "NaN"
            writing_lines += "%s\t%s\n" % (key, avg_score)

    with open(output_file, 'w') as out_file:
        out_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
