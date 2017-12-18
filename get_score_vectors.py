# get_score_vectors_per_transcript.py
# input file of transcripts, exon number, gene name, if identified by S or not.
# return transcript, region (exon-exon), comma separated scores in region, gene

import sys
import getopt
import csv


def main(argv):
    input_file = ''
    output_file = ''
    output_no_nan_file = ''

    try:
        opts, args = getopt.getopt(argv, "hi:o:n:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python new_get_score_vectors.py -i <input_file> -o <output_file> -n <no_NaN>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "new_get_score_vectors.py -i <input_file> -o <output_file> -n <no_NaN>"
            sys.exit()
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-o", "--output"):
            output_file = arg
        elif opt in ("-n", "--noNaN"):
            output_no_nan_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()
    if output_file == '':
        print "No output file specified."
        sys.exit()
    if output_no_nan_file == '':
        print "No output without NaN values file specified."
        sys.exit()

    region_exons = []
    region_scores = []
    writing_lines = ""
    no_nan_lines = ""
    previous_state = ""
    previous_transcript = ""
    previous_gene = ""
    score_vector = ""
    score_vector_no_nan = ""

    with open(input_file, 'r') as i_file:
        trans_line = csv.reader(i_file, delimiter=" ")

        for line in trans_line:

            if len(line) == 4:  # this is not a S identified sig region
                state = 0
            else:
                state = line[4]  # this was identified as a sig region by S

            transcript = line[0]

            if state == previous_state and transcript == previous_transcript:  # same region and transcript
                region_exons.append(int(line[1]))
                region_scores.append(line[2])

            elif state != previous_state and transcript == previous_transcript:  # different region in same transcript
                # record the previous region immediately
                min_exon = min(region_exons)
                max_exon = max(region_exons)

                for score in region_scores:
                    score_vector += "%s," % score
                    if score != "NaN":
                        score_vector_no_nan += "%s," % score
                score_vector = score_vector[:-1]  # remove the trailing comma at the end of vector
                score_vector_no_nan = score_vector_no_nan[:-1]

                if min_exon == max_exon:
                    writing_lines += "%s\t%s\t%s\t%s\n" % (previous_transcript, min_exon, score_vector, previous_gene)
                    if score_vector_no_nan != "":
                        no_nan_lines += "%s\t%s\t%s\t%s\n" % (previous_transcript, min_exon, score_vector_no_nan,
                                                              previous_gene)
                else:
                    writing_lines += "%s\t%s-%s\t%s\t%s\n" % (previous_transcript, min_exon, max_exon, score_vector,
                                                              previous_gene)
                    if score_vector_no_nan != "":  # sometimes it's blank because all scores were NaN
                        no_nan_lines += "%s\t%s-%s\t%s\t%s\n" % (previous_transcript, min_exon, max_exon,
                                                                 score_vector_no_nan, previous_gene)

                # reset for new transcript, new gene, state, new region exon, scores
                previous_transcript = transcript
                previous_state = state
                previous_gene = line[3]
                region_exons = [int(line[1])]
                region_scores = [line[2]]
                score_vector = ""
                score_vector_no_nan = ""

            elif transcript != previous_transcript:  # different transcript
                # handle just starting out the file
                if previous_transcript == "":  # just started the file
                    region_exons = [int(line[1])]
                    region_scores = [line[2]]
                    previous_transcript = transcript
                    previous_state = state
                else:  # just a different transcript
                    # write the last region of the previous transcript immediately
                    min_exon = min(region_exons)
                    max_exon = max(region_exons)

                    for score in region_scores:
                        score_vector += "%s," % score
                        if score != "NaN":
                            score_vector_no_nan += "%s," % score
                    score_vector = score_vector[:-1]  # remove the trailing comma at the end of vector
                    score_vector_no_nan = score_vector_no_nan[:-1]

                    if min_exon == max_exon:
                        writing_lines += "%s\t%s\t%s\t%s\n" % (previous_transcript, min_exon, score_vector,
                                                               previous_gene)
                        if score_vector_no_nan != "":
                            no_nan_lines += "%s\t%s\t%s\t%s\n" % (previous_transcript, min_exon,
                                                                  score_vector_no_nan, previous_gene)
                    else:
                        writing_lines += "%s\t%s-%s\t%s\t%s\n" % (previous_transcript, min_exon, max_exon, score_vector,
                                                                  previous_gene)
                        if score_vector_no_nan != "":
                            no_nan_lines += "%s\t%s-%s\t%s\t%s\n" % (previous_transcript, min_exon, max_exon,
                                                                     score_vector_no_nan, previous_gene)

                # reset for the new transcript, new gene, state, new region exon, scores
                previous_transcript = transcript
                previous_gene = line[3]
                previous_state = state
                region_exons = [int(line[1])]
                region_scores = [line[2]]
                score_vector = ""
                score_vector_no_nan = ""

        # don't forget about the very last exon of the very last transcript
        min_exon = min(region_exons)
        max_exon = max(region_exons)

        for score in region_scores:
            score_vector += "%s," % score
            if score != "NaN":
                score_vector_no_nan += "%s," % score
        score_vector = score_vector[:-1]
        score_vector_no_nan = score_vector_no_nan[:-1]

        if min_exon == max_exon:
            writing_lines += "%s\t%s\t%s\t%s\n" % (previous_transcript, min_exon, score_vector, previous_gene)
            if score_vector_no_nan != "":
                no_nan_lines += "%s\t%s\t%s\t%s\n" % (previous_transcript, min_exon, score_vector_no_nan, previous_gene)
        else:
            writing_lines += "%s\t%s-%s\t%s\t%s\n" % (previous_transcript, min_exon, max_exon, score_vector,
                                                      previous_gene)
            if score_vector_no_nan != "":
                no_nan_lines += "%s\t%s-%s\t%s\t%s\n" % (previous_transcript, min_exon, max_exon, score_vector_no_nan,
                                                         previous_gene)

    with open(output_file, 'w') as write_file:
        write_file.write(writing_lines)

    with open(output_no_nan_file, 'w') as nan_file:
        nan_file.write(no_nan_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
