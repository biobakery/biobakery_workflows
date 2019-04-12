#!/usr/bin/env python

""" This script will remove primers from fastq files using cutadapt
    Primers and reverse compliments are stored in files
"""

import argparse
import sys
import os, fnmatch


def parse_arguments(args):
    """ Parse the arguments from the user """

    parser = argparse.ArgumentParser(
        description="Remove primers\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input_folder",
        help="input folder with N filtered files\n[REQUIRED]",
        required=True)
    parser.add_argument(
        "--fwd_primer_file",
        help="Path to file with fwd primers\n[REQUIRED]",
        required=True)
    parser.add_argument(
        "--rev_primer_file",
        help="Path to file with rev primers\n[REQUIRED]",
        required=True)
    parser.add_argument(
        "--pair_id", default="R1_0001",
        help="Pair identifier\n[REQUIRED]",
        required=False)

    return parser.parse_args()

def main():
    # parse arguments from the user
    args = parse_arguments(sys.argv)

    with open(args.fwd_primer_file) as f:
        FWD = f.read().splitlines()
    with open(args.rev_primer_file) as f:
        REV = f.read().splitlines()
    pair_id2 = args.pair_id.replace("1", "2")
    fwd_reads = fnmatch.filter(os.listdir(args.input_folder), "*" + args.pair_id + "*.fastq*")
    rev_reads = fnmatch.filter(os.listdir(args.input_folder), "*" + pair_id2 + "*.fastq*")
    fwd_reads_inp = sorted([os.path.join(args.input_folder, f) for f in fwd_reads])
    rev_reads_inp = sorted([os.path.join(args.input_folder, f) for f in rev_reads])
    cutadapt_folder = os.path.join(args.input_folder, "cutadapt")
    if not os.path.exists(cutadapt_folder):
        os.mkdir(cutadapt_folder)
    fwd_reads_out = sorted([os.path.join(cutadapt_folder, f) for f in fwd_reads])
    rev_reads_out = sorted([os.path.join(cutadapt_folder, f) for f in rev_reads])

    for i in range(0,len(rev_reads_inp)):
        command_string="cutadapt -g "+FWD[0]+" -a "+REV[1]+" -G "+REV[0]+" -A "+FWD[1]+" -n 2 \-o "+fwd_reads_out[i]+\
                       " -p "+rev_reads_out[i]+" "+fwd_reads_inp[i]+" "+ rev_reads_inp[i]
        os.system(command_string)

if __name__ == "__main__":
    main()


