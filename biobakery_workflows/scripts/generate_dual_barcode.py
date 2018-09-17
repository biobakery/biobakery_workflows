#!/usr/bin/env python

""" This script will generate dual barcode file from index(barcode) files.
    It is used for demultiplexing dual indexed reads.
"""

import argparse
import sys
import os, fnmatch
from biobakery_workflows.tasks import general


def parse_arguments(args):
    """ Parse the arguments from the user """

    parser = argparse.ArgumentParser(
        description="Generate dual barcode file for demultiplexing\n",
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
    # name dual barcode file
    dual_barcode_file = os.path.join(args.input,"dual_barcode_file.txt")
    # generate dual barcode file
    general.generate_dual_barcode(barcode_files, dual_barcode_file)

if __name__ == "__main__":
    main()