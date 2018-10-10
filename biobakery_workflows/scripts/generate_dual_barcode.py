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
    parser.add_argument(
        "--barcode-pair-identifier", default="_1",
        help="pair identifier for barcodes(index) files\n[REQUIRED]",
        required=False)

    return parser.parse_args()

def main():
    # parse arguments from the user
    args = parse_arguments(sys.argv)
    pair_id1 = args.barcode_pair_identifier
    pair_id2 = pair_id1.replace("1","2")
    barcode_files1 = fnmatch.filter(os.listdir(args.input), "*barcode"+ pair_id1 + ".fastq*")
    barcode_files2 = fnmatch.filter(os.listdir(args.input), "*barcode" + pair_id2 + ".fastq*")
    barcode_files = barcode_files1 + barcode_files2
    barcode_files = [os.path.join(args.input, file) for file in barcode_files]
    # name dual barcode file
    dual_barcode_file = os.path.join(args.input,"dual_barcode_file.txt")
    # generate dual barcode file
    general.generate_dual_barcode(barcode_files, dual_barcode_file)

if __name__ == "__main__":
    main()