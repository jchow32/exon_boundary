# kruskal_wallis_format.py
# get desired format:
# transcript, gene, all scores in transcript (,), # of scores contributed per region(,), exon ranges(,)

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
        print "python kruskal_wallis_format.py -i <input_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python kruskal_wallis_format.py -i <input_file> -o <output_file>"
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

    region_exons = []
    region_scores = []
    region_num_scores = []

    writing_lines = ""
    previous_transcript = ""
    previous_gene = ""
    score_vector = ""
    exon_vector = ""
    final_vector = ""
    num_scores_per_region_vector = ""
    all_nan = 0

    with open(input_file, 'r') as i_file:
        trans_line = csv.reader(i_file, delimiter="\t")

        for line in trans_line:

            transcript = line[0]

            if transcript == previous_transcript:  # same transcript
                region_exons.append(line[1])
                region_scores.append(line[2])

                # reset for new transcript, new gene, state, new region exon, scores
                previous_transcript = transcript
                previous_gene = line[3]

            elif transcript != previous_transcript:  # different transcript
                # handle just starting out the file
                if previous_transcript == "":  # just started the file
                    region_exons = [line[1]]
                    region_scores = [line[2]]
                    previous_transcript = transcript
                else:  # just a different transcript
                    # write the last region of the previous transcript immediately
                    for score in region_scores:
                        score_vector += "%s," % score
                    score_vector = score_vector[:-1]

                    for exon_range in region_exons:
                        exon_vector += "%s," % exon_range
                    exon_vector = exon_vector[:-1]

                    exons_per_region_array = exon_vector.split(",")
                    scores_per_region_array = score_vector.split(",")

                    for e_range in exons_per_region_array:
                        if "-" in e_range:  # a range between 2 exons
                            start, end = e_range.split("-")
                            num_scores = int(end) - int(start) + 1
                            #  one_score_range = scores_per_region_array[int(start)-1:int(end)]
                            #  if one_score_range.count("NaN") == len(one_score_range):  # then they were all NaNs
                            #    all_nan += 1

                        else:  # is just a single number exon
                            num_scores = 1
                            # if scores_per_region_array[int(e_range)-1].count("NaN") == 1:
                            #     all_nan += 1

                        region_num_scores.append(num_scores)
                        num_scores_per_region_vector += "%s," % num_scores
                    num_scores_per_region_vector = num_scores_per_region_vector[:-1]  # remove the trailing comma

                    for i, j in zip(region_exons, region_num_scores):
                        final_vector += ("%s," % i) * j
                    final_vector = final_vector[:-1]

                    # Correction for when last exon was split into 2 separate regions
                    if len(region_exons) == 1:  # there is only 1 region, the last exon was split
                        score_vector += ",%s" % score_vector.rsplit(',', 1)[1]  # add the final score again
                        num_scores_per_region_vector += ',1'
                        exon_vector += ",%s" % exon_vector.rsplit('-', 1)[1]  # add the last element only
                        final_vector += ",%s" % final_vector.rsplit('-', 1)[1]  # add the last element only
                        # append "1" to exons_per_region_array
                        exons_per_region_array.append("1")

                    # Correction for when n-1 of n regions are entirely NaN and should not be written
                    for e_range in exon_vector.split(","):
                        if "-" in e_range:
                            start, end = e_range.split("-")
                            one_score_range = scores_per_region_array[int(start)-1:int(end)]
                            if one_score_range.count("NaN") == len(one_score_range):
                                all_nan += 1
                        else:
                            if scores_per_region_array[int(e_range)-1].count("NaN") == 1:
                                all_nan += 1

                    if all_nan < len(exons_per_region_array)-1:  # then write the transcript
                        writing_lines += "%s\t%s\t%s\t%s\t%s\t%s\n" % (previous_transcript, previous_gene, score_vector,
                                                                       num_scores_per_region_vector, exon_vector,
                                                                       final_vector)

                # reset for the new transcript, new gene, state, new region exon, scores
                previous_transcript = transcript
                previous_gene = line[3]
                region_exons = [line[1]]
                region_scores = [line[2]]
                region_num_scores = []
                all_nan = 0
                score_vector = ""
                exon_vector = ""
                final_vector = ""
                num_scores_per_region_vector = ""

        # don't forget about the very last transcript

        for score in region_scores:
            score_vector += "%s," % score
        score_vector = score_vector[:-1]

        for exon_range in region_exons:
            exon_vector += "%s," % exon_range
        exon_vector = exon_vector[:-1]

        exons_per_region_array = exon_vector.split(",")
        for e_range in exons_per_region_array:
            if "-" in e_range:  # a range between 2 exons
                start, end = e_range.split("-")
                num_scores = int(end) - int(start) + 1
            else:  # is just a single number exon
                num_scores = 1
            region_num_scores.append(num_scores)
            num_scores_per_region_vector += "%s," % num_scores
        num_scores_per_region_vector = num_scores_per_region_vector[:-1]

        for i, j in zip(region_exons, region_num_scores):
            final_vector += ("%s," % i) * j
        final_vector = final_vector[:-1]

        # Correction for when last exon was split into 2 separate regions
        if len(region_exons) == 1:  # there is only 1 region, the last exon was split
            score_vector += ",%s" % score_vector.rsplit(',', 1)[1]  # add the final score again
            num_scores_per_region_vector += ',1'
            exon_vector += ",%s" % exon_vector.rsplit('-', 1)[1]  # add the last element only
            final_vector += ",%s" % final_vector.rsplit('-', 1)[1]  # add the last element only

        writing_lines += "%s\t%s\t%s\t%s\t%s\t%s\n" % (previous_transcript, previous_gene, score_vector,
                                                       num_scores_per_region_vector, exon_vector, final_vector)

    with open(output_file, 'w') as out_file:
        out_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
