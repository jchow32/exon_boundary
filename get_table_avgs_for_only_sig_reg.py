# get_table_avgs_for_only_sig_reg.py
# given file with columns
# transcript0, exon1, # missense2, # synonymous3, total missense4, total synonymous5, score6, p-value7, q-value8,
# gene9

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
        print "python get_table_avgs_for_only_sig_reg.py -i <input_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python get_table_avgs_for_only_sig_reg.py -i <input_file> "
            sys.exit()
        elif opt in ("-i", "--input"):
            input_file = arg

    if input_file == '':
        print "No input file specified."
        sys.exit()

    gene_dict_miss = {}
    gene_dict_syn = {}
    gene_dict_miss_syn = {}
    gene_dict_score = {}
    gene_dict_pli = {}

    gene_dict_count = {}
    gene_dict_non_nan_miss_syn_count = {}
    gene_dict_non_nan_score_count = {}
    gene_dict_non_nan_pli_count = {}

    sum_miss_avg = 0
    sum_syn_avg = 0
    sum_miss_syn_avg = 0
    sum_pli_avg = 0
    sum_score_avg = 0

    non_nan_miss_syn_gene = 0
    non_nan_pli_gene = 0
    non_nan_score_gene = 0
    genes_total = 0

    with open(input_file, 'r') as i_file:
        g_line = csv.reader(i_file, delimiter="\t")

        for line in g_line:
            gene = line[0]
            miss = float(line[3])
            syn = float(line[4])
            if syn == 0:
                miss_syn = "NaN"
            else:
                miss_syn = miss / syn

            score = line[7]
            if len(line) == 11:  # means that pLI exists
                pli = line[10]
            else:
                pli = "NaN"

            # initialize gene dictionaries
            if gene not in gene_dict_miss:
                gene_dict_miss[gene] = 0
                gene_dict_syn[gene] = 0
                gene_dict_miss_syn[gene] = 0
                gene_dict_score[gene] = 0
                gene_dict_pli[gene] = 0

                gene_dict_count[gene] = 0
                gene_dict_non_nan_miss_syn_count[gene] = 0
                gene_dict_non_nan_score_count[gene] = 0
                gene_dict_non_nan_pli_count[gene] = 0

            gene_dict_miss[gene] += miss
            gene_dict_syn[gene] += syn

            if miss_syn != "NaN":
                gene_dict_miss_syn[gene] += float(miss_syn)
                gene_dict_non_nan_miss_syn_count[gene] += 1.0

            if score != "NaN":
                gene_dict_score[gene] += float(score)
                gene_dict_non_nan_score_count[gene] += 1.0

            if pli != "NaN":
                gene_dict_pli[gene] += float(pli)
                gene_dict_non_nan_pli_count[gene] += 1.0

            gene_dict_count[gene] += 1.0

    # get the average total #miss, #syn, #miss/syn, pLI for all genes lumped into 1

    for gene_key in gene_dict_miss:
        sum_miss_avg += gene_dict_miss[gene_key] / gene_dict_count[gene_key]
        sum_syn_avg += gene_dict_syn[gene_key] / gene_dict_count[gene_key]
        if gene_dict_non_nan_miss_syn_count[gene_key] != 0:
            sum_miss_syn_avg += gene_dict_miss_syn[gene_key] / gene_dict_non_nan_miss_syn_count[gene_key]
            non_nan_miss_syn_gene += 1.0
        if gene_dict_non_nan_pli_count[gene_key] != 0:
            sum_pli_avg += gene_dict_pli[gene_key] / gene_dict_non_nan_pli_count[gene_key]
            non_nan_pli_gene += 1.0
        if gene_dict_non_nan_score_count[gene_key]:
            sum_score_avg += gene_dict_score[gene_key] / gene_dict_non_nan_score_count[gene_key]
            non_nan_score_gene += 1.0

        genes_total += 1.0

    print("Average # missense for all sig regions: %s" % (sum_miss_avg / genes_total))
    print("Average # synonymous for all sig regions: %s" % (sum_syn_avg / genes_total))
    print("Average # missense/synonymous for all sig regions: %s" % (sum_miss_syn_avg / non_nan_miss_syn_gene))
    print("Average score for all sig regions: %s" % (sum_score_avg / non_nan_score_gene))
    print("Average pLI for all sig regions: %s" % (sum_pli_avg / non_nan_pli_gene))


if __name__ == "__main__":
    main(sys.argv[1:])
