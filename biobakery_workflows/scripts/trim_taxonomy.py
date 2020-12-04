#!/usr/bin/env python

import sys
import os
import argparse
import gzip

try:
    from biobakery_workflows import utilities
except ImportError:
    sys.exit("Please install biobakery workflows.")

# This script will trim a taxonomy file to only include the most specific clade
# followed by any unknown taxonomy if included in the full name.

BIOM_COMMENT = "# Constructed from biom file"

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
    parser.add_argument(
        "-t", "--taxonomy-column",
        help="the column number, zero based index, containing taxonomy information [auto-detected]",
        default=None,
        type=int)
    parser.add_argument(
        "-e", "--end-taxonomy-column",
        help="the column number, zero based index, to write the taxonomy information [default is original index]",
        default=None,
        type=int)

    return parser.parse_args()


def main():

    args=parse_arguments(sys)

    try:
        file_handle_write=open(args.output,"wt")
    except EnvironmentError:
        sys.exit("Error: Unable to open output file: " + args.output)

    try:
        if args.input.endswith(".gz"):
            file_handle_read=gzip.open(args.input,"rt")
        else:
            file_handle_read=open(args.input,"rt")
    except EnvironmentError:
        sys.exit("Error: Unable to read input file: " + args.input)

    # write the header to the new file
    header = file_handle_read.readline().rstrip().split("\t")
    
    # ignore comment if present
    if header[0].startswith(BIOM_COMMENT):
        header = file_handle_read.readline().rstrip().split("\t")

    # trim the taxonomy and sum species
    taxonomy_data = {}
    for line in file_handle_read:
        # ignore lines that are comments
        if line.startswith("#"):
            file_handle_write.write(line)
        else:
            data=line.rstrip().split("\t")
            if args.taxonomy_column is None:
                # try to figure out which column has the taxonomy data
                try:
                    args.taxonomy_column=[index for index, value in enumerate(data) if "k__" in value][0]
                except IndexError:
                    sys.exit("Error unable to find the taxonomy column. Please provide it with the option --taxonomy-column <0>.")
            if args.end_taxonomy_column is None:
                args.end_taxonomy_column=args.taxonomy_column

            new_taxonomy=utilities.taxonomy_trim([data[args.taxonomy_column]])[0]
            data.pop(args.taxonomy_column)
            if new_taxonomy in taxonomy_data:
                data = [data[0]]+[str(float(a)+float(b)) for a,b in zip(taxonomy_data[new_taxonomy][1:], data[1:])] 
            taxonomy_data[new_taxonomy]=data               
            
    # write the header
    old_taxon = header.pop(args.taxonomy_column)
    header[args.end_taxonomy_column]=old_taxon

    file_handle_write.write("\t".join(header)+"\n")

    # write the new data  
    for taxon,data in taxonomy_data.items():
        data[args.end_taxonomy_column]=taxon
        file_handle_write.write("\t".join(data)+"\n")

    file_handle_read.close()
    file_handle_write.close()

if __name__ == "__main__":
    main()


