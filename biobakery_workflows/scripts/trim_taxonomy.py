#!/usr/bin/env python

import sys
import os
import argparse

try:
    from biobakery_workflows import utilities
except ImportError:
    sys.exit("Please install biobakery workflows.")

# This script will trim a taxonomy file to only include the most specific clade
# followed by any unknown taxonomy if included in the full name.

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Trim taxonomy\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the otu table\n[REQUIRED]",
        metavar="<otu.tsv>",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="file to write the trimmed otu table\n[REQUIRED]",
        metavar="<output.tsv>",
        required=True)

    return parser.parse_args()


def main():

    args=parse_arguments(sys)

    try:
        file_handle_write=open(args.output,"wt")
    except EnvironmentError:
        sys.exit("Error: Unable to open output file: " + args.output)

    try:
        file_handle_read=open(args.input,"rt")
    except EnvironmentError:
        sys.exit("Error: Unable to read input file: " + args.input)

    # write the header to the new file
    file_handle_write.write(file_handle_read.readline())
    
    # trim the taxonomy
    for line in file_handle_read:
        data=line.rstrip().split("\t")
        trimmed_taxon=utilities.taxonomy_trim([data[-1]])
        file_handle_write.write("\t".join(data[:-1]+trimmed_taxon)+"\n")

    file_handle_read.close()
    file_handle_write.close()

if __name__ == "__main__":
    main()


