#!/n/sw/eb/apps/centos7/Anaconda3/2019.10/bin/python

import sys
import os
import argparse
import gzip

# This script will take any type of tab-delimited table and reformat it as a feature table
# to be used as input for Maaslin2 and other downstream stats processing.

STRATIFIED_DELIMITER = "|"

BIOM_COMMENT = "# Constructed from biom file"

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Create feature table\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the count table\n[REQUIRED]",
        metavar="<input.tsv>",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="file to write the feature table\n[REQUIRED]",
        metavar="<output>",
        required=True)
    parser.add_argument(
        "--sample-tag-columns",
        help="remove this string from the sample names in columns")
    parser.add_argument(
        "--remove-stratified",
        help="remove stratified rows",
        action="store_true")
    parser.add_argument(
        "--reduce-stratified-species-only",
        help="reduce stratified rows to just species",
        action="store_true")
    parser.add_argument(
        "--max-length-feature-name",
        help="the max length for each feature name",
        default=65)

    return parser.parse_args()


def main():

    args=parse_arguments(sys)

    # determine the open function based on file type
    if args.input.endswith(".gz"):
        open_read = gzip.open
    else:
        open_read = open

    # read in the file and process depending on the arguments provided
    with open_read(args.input, mode='rt') as file_handle_read:
        with open(args.output,"w") as file_handle_write:
            # remove sample tags from column headers if present
            header = file_handle_read.readline()
            # remove biom comment
            if header.startswith(BIOM_COMMENT):
                header = file_handle_read.readline()
            if args.sample_tag_columns:
                header = "\t".join([i.split(args.sample_tag_columns)[0] for i in header.rstrip().split("\t")])
            file_handle_write.write(header+"\n")

            for line in file_handle_read:
                # ignore and do not write out commented lines
                if not line.startswith("#"):
                    filter=False
                    if args.remove_stratified and STRATIFIED_DELIMITER in line:
                        filter=True

                    # only print out species, reducing taxonomy information to just genus and species
                    if args.reduce_stratified_species_only:
                        filter=True
                        if "s__" in line and not "t__" in line:
                            filter=False
                            info = line.split("\t")
                            if STRATIFIED_DELIMITER in info[0]:
                                taxon = info[0].split(STRATIFIED_DELIMITER)
                                line = "\t".join([STRATIFIED_DELIMITER.join([taxon[-2],taxon[-1]])]+info[1:])

                    # if max length for feature name set then reduce if needed
                    info = line.split("\t")
                    if args.max_length_feature_name and len(info[0]) > args.max_length_feature_name:
                        info[0]=info[0][0:args.max_length_feature_name+1]+"..."
                        line="\t".join(info)

                    if not filter:
                        file_handle_write.write(line)
        
if __name__ == "__main__":
    main()


