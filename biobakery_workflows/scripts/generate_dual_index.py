#!/usr/bin/env python

""" This script will generate dual index file from barcode files.
    It is used for demultiplexing dual indexed reads.
"""

import argparse
import sys
import os, fnmatch
import itertools


def parse_arguments(args):
    """ Parse the arguments from the user """

    parser = argparse.ArgumentParser(
        description="Generate dual index file for demultiplexing\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input",
        help="input folder with dual indexed reads and barcodes\n[REQUIRED]",
        required=True)

    return parser.parse_args()

def main():
    # parse arguments from the user
    args = parse_arguments(sys.argv)

    barcode_files = fnmatch.filter(os.listdir(args.input), '*barcode_*')
    barcode_files = [os.path.join(args.input, file) for file in barcode_files]
    allbarcodes = set()
    for barcode_file in barcode_files:
        try:
            file_handle = open(barcode_file)
            lines = file_handle.readlines()
            file_handle.close()
        except EnvironmentError:
            sys.exit("ERROR: Unable to read barcode: " + barcode_file)
        allbarcodes.update(lines[1::4])

    dual_indexes_all = list(itertools.combinations(allbarcodes, 2))
    dual_indexes = set(dual_indexes_all)

    dual_index_file = os.path.join(args.input,"barcodes_gen.fil")
    fh = open(dual_index_file,"w")
    i=0
    for ind in dual_indexes:
        i+= 1
        fh.write(str(i) + " " + ind[0].replace("\n","")+ "-" + ind[1].replace("\n","") + " Nextra\n")
    fh.close()

    print("Dual index file " + dual_index_file + " has been generated")

if __name__ == "__main__":
    main()