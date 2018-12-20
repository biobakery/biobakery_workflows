#!/usr/bin/env python

""" Remove files that have read-count less than min read count."""

import argparse
import sys

def parse_arguments(args):
    """ Parse the arguments from the user """

    parser = argparse.ArgumentParser(
        description="Discarding samples with low number of reads\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input-fasta",
        help="input fasta file with all samples concatanated filtered and truncated\n[REQUIRED]",
        required=True)
    parser.add_argument(
        "--output-fasta",
        help="output fasta file with low read number samples removed\n[REQUIRED]",
        required=False)
    parser.add_argument(
        "--output-discarded",
        help="output fasta file with low read number samples\n[REQUIRED]",
        required=False)
    parser.add_argument(
        "--min-read-count", default=5000,
        help="minimal number of reads in sample threshold",
        required=False)

    return parser.parse_args()

def main():

        args = parse_arguments(sys.argv)
        file = args.input_fasta

        # count reads in samples
        allsamples = {}
        for line in open(file).read().split('>'):
            samplename = line.split(".")[0]

            if samplename in allsamples:
                allsamples[samplename] += 1
            else:
                allsamples[samplename] = 1

        # remove samples with read counts less than min_read_count
        file_handle_nolowreads = open(args.output_fasta, "w")
        file_handle_lowreads = open(args.output_discarded, "w")
        for line in open(file).read().split('>'):
            samplename = line.split(".")[0]
            if allsamples[samplename] > int(args.min_read_count):
                file_handle_nolowreads.write(">" + line)
            else:
                if len(line) > 0:
                    file_handle_lowreads.write(">" + line)

        file_handle_nolowreads.close()
        file_handle_lowreads.close()


if __name__ == "__main__":
    main()

