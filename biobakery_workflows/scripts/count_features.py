#!/usr/bin/env python

import sys
import os
import argparse

# This script will count all features for each sample. It can ignore stratification by species
# and also "un"-features if these options are set by the user.

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Count features for each species\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the feature table\n[REQUIRED]",
        metavar="<input.tsv>",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="file to write counts table\n[REQUIRED]",
        metavar="<output>",
        required=True)
    parser.add_argument(
        "--reduce-sample-name",
        help="remove the extra strings from the sample names",
        action='store_true')
    parser.add_argument(
        "--ignore-un-features",
        help="do not count UNMAPPED, UNGROUPED, or UNINTEGRATED features",
        action='store_true')
    parser.add_argument(
        "--ignore-stratification",
        help="only count the main feature not the stratified instances",
        action='store_true')
    parser.add_argument(
        "--include",
        help="only count features with this string included")
    parser.add_argument(
        "--filter",
        help="do not count features with this string included")

    return parser.parse_args()


def main():

    args=parse_arguments(sys)

    data=[]
    samples=[]
    with open(args.input) as file_handle:
        # remove RPKs from sample name if included
        if args.reduce_sample_name:
            samples=[sample.replace("_Abundance-RPKs","").replace("_genefamilies_Abundance","").replace("_Abundance","").replace("_taxonomic_profile","") for sample in file_handle.readline().rstrip().split("\t")[1:]]
        else:
            samples=file_handle.readline().rstrip().split("\t")[1:]
             
        for line in file_handle:
            if "|" in line and args.ignore_stratification:
                store=False
            else:
                store=True

            if "UNMAPPED" in line or "UNGROUPED" in line or "UNINTEGRATED" in line and args.ignore_un_features:
                store=False

            if args.include and not args.include in line:
                store=False

            if args.filter and args.filter in line:
                store=False

            if store:
                data.append(line.rstrip().split("\t")[1:])

    try:
        file_handle=open(args.output,"w")
    except EnvironmentError:
        sys.exit("Error: Unable to open output file: " + args.output)

    # write out the header
    file_handle.write("\t".join(["# samples","total features"])+"\n")

    # count the total features for each sample
    for i, sample in enumerate(samples):
        features=0
        for row in data:
            if float(row[i]) > 0:
                features+=1
        file_handle.write(sample+"\t"+str(features)+"\n")

    file_handle.close()

if __name__ == "__main__":
    main()


