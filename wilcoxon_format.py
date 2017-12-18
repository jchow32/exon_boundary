# wilcoxon_format.py
# get desired format:
# for every transcript (kruskal p-value < 0.05), get all possible pairwise combinations of scores, print each to line

import sys
import getopt
import csv
import itertools


def main(argv):
    input_file = ''
    output_file = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python wilcoxon_format.py -i <input_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python wilcoxon_format.py -i <input_file> -o <output_file>"
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
    string_pair0 = ""
    string_pair1 = ""
    score_array = []

    with open(input_file, 'r') as i_file:
        k_line = csv.reader(i_file, delimiter="\t")

        for line in k_line:
            transcript = line[0]
            gene = line[1]
            scores = line[2].split(",")
            exon_range = line[4]

            for e_range in exon_range.split(","):  # every individual exon range
                if "-" in e_range:  # it's a range
                    start, end = e_range.split("-")
                    score_array.append(scores[int(start)-1:int(end)])
                else:  # it's a single exon
                    score_array.append(scores[int(e_range)-1])

            for item in score_array:
                if isinstance(item, list):  # check if it's a list of NaN
                    if item.count("NaN") == len(item):
                        score_array.remove(item)
                else:  # check if it's a single NaN
                    if item == "NaN":
                        score_array.remove(item)

            for pair in itertools.combinations(score_array, 2):  # all possible pairwise comparisons for scores
                # convert each part of pair from list to comma delimited string

                if isinstance(pair[0], list):  # there are multiple scores
                    for a in pair[0]:
                        string_pair0 += "%s," % a
                    string_pair0 = string_pair0[:-1]  # remove the trailing comma
                else:
                    string_pair0 += str(pair[0])

                if isinstance(pair[1], list):
                    for b in pair[1]:
                        string_pair1 += "%s," % b
                    string_pair1 = string_pair1[:-1]
                else:
                    string_pair1 += str(pair[1])

                # Correction for when one of the pairs is entirely NaN
                # Do not write lines where NaN count == number of commas + 1
                writing_lines += "%s\t%s\t%s\t%s\n" % (transcript, gene, string_pair0, string_pair1)
                string_pair0 = ""
                string_pair1 = ""

            score_array = []  # reset for every transcript

    with open(output_file, 'w') as out_file:
        out_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
