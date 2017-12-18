# exon_boundary_from_1.py
# From exon boundaries, convert chromosomal positions to relative positions (starting from the number 1)

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
        print "python exon_boundary_from_1.py -i <input_file> -o <output_file>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print "python exon_boundary_from_1.py -i <input_file> -o <output_file>"
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
    transcript_array = []
    start_array = []
    end_array = []

    with open(input_file, 'r') as i_file:
        e_line = csv.reader(i_file, delimiter=" ")

        for line in e_line:
            transcript = line[0]
            start = int(line[1])
            end = int(line[2])

            if transcript not in transcript_array:  # moving onto new transcript
                transcript_array.append(transcript)

                if len(transcript_array) == 1:  # the very first transcript
                    start_array.append(start)
                    end_array.append(end)
                else:  # not the very first transcript
                    # immediately write the previous transcript
                    min_start = min(start_array) - 1
                    start_array[:] = [x - min_start for x in start_array]
                    end_array[:] = [y - min_start for y in end_array]
                    index = 1
                    for i, j in zip(start_array, end_array):
                        writing_lines += "%s\t%s\t%s\t%s\n" % (transcript_array[-2], i, j, index)
                        index += 1
                    start_array = [start]
                    end_array = [end]

            else:  # same transcript as before
                start_array.append(start)
                end_array.append(end)

        # don't forget about the very last transcript

        min_start = min(start_array) - 1
        start_array[:] = [x - min_start for x in start_array]  # subtract the min from every array element
        end_array[:] = [x - min_start for x in end_array]  # subtract the min from every array element
        index = 1
        for i, j in zip(start_array, end_array):
            writing_lines += "%s\t%s\t%s\t%s\n" % (transcript_array[-1], i, j, index)
            index += 1

    with open(output_file, 'w') as out_file:
        out_file.write(writing_lines)


if __name__ == "__main__":
    main(sys.argv[1:])
