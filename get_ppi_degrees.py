# get_ppi_degrees.py
# Given a list of PPIs, determine if genes 1 and 2 are both identified by gene list A, gene list B,
# and/or are separately identified
# Gene lists are several lines, each line containing 1 gene name

import sys
import getopt
import csv


def main(argv):
    ppi_file = ''
    a_genes_file = ''
    b_genes_file = ''

    a_output_file = ''
    b_output_file = ''
    ab_output_file = ''

    try:
        opts, args = getopt.getopt(argv, "hp:a:b:x:y:z:")
        # requires input and output file
    except getopt.GetoptError:
        print "See command line format:"
        print "python get_ppi_degrees.py -p <ppi_file> -a <a_genes_file> -b <b_genes_file> -x <ppi_A> -y <ppi_B> -z <ppi_AB>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print "python get_ppi_degrees.py -p <ppi_file> -a <a_genes_file> -b <b_genes_file> -x <ppi_A> -y <ppi_B> -z <ppi_AB>"
            print "Outputs x, y, and z corresponding to if the PPI was found in gene list A, B, or partially in A or B."
            sys.exit()
        elif opt in ("-p", "--ppi"):
            ppi_file = arg
        elif opt in ("-a", "--a_gene"):
            a_genes_file = arg
        elif opt in ("-b", "--b_gene"):
            b_genes_file = arg
        elif opt in "-x":
            a_output_file = arg
        elif opt in "-y":
            b_output_file = arg
        elif opt in "-z":
            ab_output_file = arg

    if ppi_file == '':
        print "No ppi file specified."
        sys.exit()
    if a_genes_file == '':
        print "No a genes list file specified."
        sys.exit()
    if b_genes_file == '':
        print "No b genes list file specified."
        sys.exit()
    if a_output_file == '':
        print "No output file for A genes specified."
        sys.exit()
    if b_output_file == '':
        print "No output file for B genes specified."
        sys.exit()
    if ab_output_file == '':
        print "No output for AB genes specified."

    a_writing_lines = ""
    b_writing_lines = ""
    ab_writing_lines = ""
    a_genes_array = []
    b_genes_array = []

    # read all the genes from list A and put into array
    with open(a_genes_file, 'r') as a_file:

        for line in a_file:
            a_genes_array.append(line.rstrip())

    # read all the genes from list B and put into array
    with open(b_genes_file, 'r') as b_file:

        for line in b_file:
            b_genes_array.append(line.rstrip())

    # read the ppi file
    with open(ppi_file, 'r') as p_file:
        ppi_line = csv.reader(p_file, delimiter="\t")

        for line in ppi_line:
            gene1 = line[0]
            gene2 = line[1]

            if gene1 in a_genes_array and gene2 in a_genes_array:
                a_writing_lines += "%s\t%s\n" % (gene1, gene2)
            if gene1 in b_genes_array and gene2 in b_genes_array:
                b_writing_lines += "%s\t%s\n" % (gene1, gene2)
            if gene1 in a_genes_array and gene2 in b_genes_array:
                ab_writing_lines += "%s\t%s\n" % (gene1, gene2)
            if gene2 in a_genes_array and gene1 in b_genes_array:
                ab_writing_lines += "%s\t%s\n" % (gene2, gene1)

    with open(a_output_file, 'w') as a_out_file:
        a_out_file.write(a_writing_lines)

    with open(b_output_file, 'w') as b_out_file:
        b_out_file.write(b_writing_lines)

    with open(ab_output_file, 'w') as ab_out_file:
        ab_out_file.write(ab_writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
