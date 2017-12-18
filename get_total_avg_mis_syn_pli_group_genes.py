# get_total_avg_mis_syn_pli_per_group_genes.py
# given averages per gene, get the total average #miss, #syn, #miss/syn, pLI for ALL genes (return 4 values total)

import sys
import getopt
import csv


def main(argv):
    input_file = ''

    try:
        opts, args = getopt.getopt(argv, "hi:")
        # requires input
    except getopt.GetoptError:
        print "See command line format:"
        print "python get_total_avg_mis_syn_pli_per_group_genes.py -i <input_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python get_total_avg_mis_syn_pli_per_group_genes.py -i <input_file> "
            sys.exit()
        elif opt in ("-i", "--input"):
            input_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()

    sum_miss = 0
    sum_syn = 0
    sum_miss_syn = 0
    sum_pli = 0
    sum_score = 0
    total_gene_count = 0
    non_nan_miss_syn_count = 0
    non_nan_pli_count = 0
    non_nan_score_count = 0

    with open(input_file, 'r') as i_file:
        g_line = csv.reader(i_file, delimiter="\t")

        for line in g_line:
            miss = float(line[1])
            syn = float(line[2])
            miss_syn = line[3]
            score = line[4]
            pli = line[5]

            sum_miss += miss
            sum_syn += syn
            if miss_syn != "NaN":
                sum_miss_syn += float(miss_syn)
                non_nan_miss_syn_count += 1.0

            if pli != "NaN":
                sum_pli += float(pli)
                non_nan_pli_count += 1.0

            if score != "NaN":
                sum_score += float(score)
                non_nan_score_count += 1.0

            total_gene_count += 1.0

        # get the average total #miss, #syn, #miss/syn, pLI for all genes lumped into 1

        total_miss_avg = sum_miss / total_gene_count
        total_syn_avg = sum_syn / total_gene_count
        total_miss_syn_avg = sum_miss_syn / non_nan_miss_syn_count
        total_pli_avg = sum_pli / non_nan_pli_count
        total_score_avg = sum_score /non_nan_score_count

        print("Average # missense for all genes: %s" % total_miss_avg)
        print("Average # synonymous for all genes: %s" % total_syn_avg)
        print("Average # missense/synonymous for all genes: %s" % total_miss_syn_avg)
        print("Average score for all genes: %s" % total_score_avg)
        print("Average pLI for all genes: %s" % total_pli_avg)


if __name__ == "__main__":
    main(sys.argv[1:])
