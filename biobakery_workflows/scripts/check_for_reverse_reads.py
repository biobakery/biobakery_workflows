#!/usr/bin/env python

# This script will check for many unclassified alignments which can indicate the reads are of the reverse strand.

import sys
import os
import argparse

UNCLASSIFIED="k__Eukaryota; p__ ; c__ ; o__ ; f__ ; g__ ; s__"
FAIL_RATIO=0.50

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Check for unaligned reads\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input", 
        help="the taxonomy file\n[REQUIRED]", 
        metavar="<input_folder>", 
        required=True)

    return parser.parse_args()

def main():
    # parse arguments from the user
    args = parse_arguments(sys)

    # count the total ASVs and those unknown
    total_asv=-1
    total_unknown=0
    for line in open(args.input):
        total_asv+=1
        if UNCLASSIFIED in line:
            total_unknown+=1

    ratio = total_unknown / (total_asv*1.0)
    if ratio > FAIL_RATIO:
        sys.exit("Total number of unknown ASVs, those classifed as '{0}', exceeds allowed RATIO {1} at {2}, it is likely your reads are from the reverse strand. Please try running again adding the option '--tryRC='TRUE'' to run DADA2 on the inverse complement of the reads.".format(UNCLASSIFIED, FAIL_RATIO, ratio))
    else:
        print("Total number of unknown ASVs, those classifed as '{0}', is less then the allowed RATIO {1} at {2}".format(UNCLASSIFIED, FAIL_RATIO, ratio))


if __name__ == "__main__":
    main()

