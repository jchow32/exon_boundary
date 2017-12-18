# get_average_mis_syn_per_gene.py
# calculate average #missense, #synonymous, #mis/syn per gene over several transcripts.

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
        print "python get_average_mis_syn_per_gene.py -i <input_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python get_average_mis_syn_per_gene.py -i <input_file> -o <output_file>"
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
    gene_dict_miss = {}
    gene_dict_syn = {}
    gene_dict_miss_syn = {}
    gene_dict_count = {}
    gene_dict_miss_syn_count = {}
    gene_dict_score = {}
    gene_dict_score_count = {}
    gene_dict_pli = {}

    sum_avg_miss = 0
    sum_avg_syn = 0
    sum_avg_miss_syn = 0
    sum_avg_score = 0
    sum_avg_pli = 0

    miss_syn_not_nan_count = 0
    score_not_nan_count = 0
    pli_not_nan_count = 0

    with open(input_file, 'r') as i_file:
        g_line = csv.reader(i_file, delimiter="\t")

        for line in g_line:
            gene = line[0]
            miss = float(line[2])
            syn = float(line[3])
            if line[4] != "NaN":
                miss_syn = float(line[4])
            else:
                miss_syn = 0

            if line[5] != "NaN":
                score = float(line[5])
            else:
                score = 0

            if line[6] != "":  # pLI value exists
                pli = line[6]
            else:
                pli = "NaN"

            # initialize dictionaries only once
            if gene not in gene_dict_miss:
                gene_dict_miss[gene] = 0
                gene_dict_syn[gene] = 0
                gene_dict_miss_syn[gene] = 0
                gene_dict_count[gene] = 0
                gene_dict_miss_syn_count[gene] = 0
                gene_dict_score_count[gene]= 0
                gene_dict_score[gene] = 0

            gene_dict_miss[gene] += miss
            gene_dict_syn[gene] += syn
            gene_dict_miss_syn[gene] += miss_syn
            gene_dict_count[gene] += 1.0
            gene_dict_score[gene] += score
            if line[4] != "NaN":  # only increment this m/s score if not NaN. Don't include NaNs in avg
                gene_dict_miss_syn_count[gene] += 1.0
            if line[5] != "NaN":
                gene_dict_score_count[gene] += 1.0
            gene_dict_pli[gene] = pli

        for gene_key in gene_dict_miss:
            avg_miss = gene_dict_miss[gene_key] / gene_dict_count[gene_key]
            avg_syn = gene_dict_syn[gene_key] / gene_dict_count[gene_key]
            if gene_dict_miss_syn_count[gene_key] != 0:
                avg_miss_syn = gene_dict_miss_syn[gene_key] / gene_dict_miss_syn_count[gene_key]
            else:
                avg_miss_syn = "NaN"

            if gene_dict_score_count[gene_key] != 0:
                avg_score = gene_dict_score[gene_key] / gene_dict_score_count[gene_key]
            else:
                avg_score = "NaN"
            pli = gene_dict_pli[gene_key]

            writing_lines += "%s\t%s\t%s\t%s\t%s\t%s\n" % (gene_key, avg_miss, avg_syn, avg_miss_syn, avg_score, pli)

            # get the total sum of the average (#miss, #syn, #miss/syn, or pLI) per gene for all genes
            sum_avg_miss += gene_dict_miss[gene_key]
            sum_avg_syn += gene_dict_syn[gene_key]
            if gene_dict_miss_syn[gene_key] != "NaN":
                sum_avg_miss_syn += gene_dict_miss_syn[gene_key]  # handle NaNs
                miss_syn_not_nan_count += 1.0
            if gene_dict_score[gene_key] != "NaN":
                sum_avg_score += gene_dict_score[gene_key]
                score_not_nan_count += 1.0
            if gene_dict_pli[gene_key] != "NaN":
                sum_avg_pli += float(gene_dict_pli[gene_key])  # handle NaNs
                pli_not_nan_count += 1.0

        # get the average total #miss, #syn, #miss/syn, pLI for all genes lumped into 1
        total_miss_avg = sum_avg_miss / len(gene_dict_miss)
        total_syn_avg = sum_avg_syn / len(gene_dict_syn)
        total_miss_syn_avg = sum_avg_miss_syn / miss_syn_not_nan_count
        total_score_avg = sum_avg_score / score_not_nan_count
        total_pli_avg = sum_avg_pli/ pli_not_nan_count
        print("Average # missense for all genes: %s" % total_miss_avg)
        print("Average # synonymous for all genes: %s" % total_syn_avg)
        print("Average # missense/synonymous for all genes: %s" % total_miss_syn_avg)
        print("Average score for all genes: %s" % total_score_avg)
        print("Average pLI for all genes: %s" % total_pli_avg)

        with open(output_file, 'w') as out_file:
            out_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
