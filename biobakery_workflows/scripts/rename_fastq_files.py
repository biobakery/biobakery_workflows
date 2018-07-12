#!/usr/bin/env python

import sys
import os
import subprocess
import string
import argparse

# rename the fastq files based on sample name instead of barcode
# to run: $ rename_fastq_files.py --input input_folder --input-barcodes barcodes.txt --output output_folder

def parse_arguments(args):
    """
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Rename files from barcodes to sample names\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the input folder containing the fastq files\n[REQUIRED]",
        metavar="<input_folder>",
        required=True)
    parser.add_argument(
        "-b", "--input-barcodes",
        help="the tab-delimited file of sample names and barcodes\n[REQUIRED]",
        metavar="<barcodes.tsv>",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="the output folder to write the new fastq files\n[REQUIRED]",
        metavar="<output_folder>",
        required=True)
    parser.add_argument(
        "--input-basename1",
        help="the basename of the input file with SEQ replacing the barcode for R1",
        default="1_ABCD.1.SEQ.unmapped.1.fastq.gz")
    parser.add_argument(
        "--input-basename2",
        help="the basename of the input file with SEQ replacing the barcode for R2",
        default="1_ABCD.1.SEQ.unmapped.2.fastq.gz")
    parser.add_argument(
        "--output-basename1",
        help="the basename of the output file with SAMPLE to be replaced with the sample name for R1",
        default="SAMPLE.R1.fastq.gz")
    parser.add_argument(
        "--output-basename2",
        help="the basename of the output file with SAMPLE to be replaced with the sample name for R2",
        default="SAMPLE.R2.fastq.gz")

    return parser.parse_args()


def main():

    args=parse_arguments(sys)

    # get the set of barcodes (and reverse complements) for each sample
    revcmp=string.maketrans("ATGC","TACG")
    barcodes={}
    for line in open(args.input_barcodes):
        data=line.rstrip().split("\t")
        newseq=data[1].translate(revcmp)[::-1]
        barcodes[newseq]=data[0]

    # create the output folder if it does not exist
    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    # rename files, copying to new output folder
    not_found=0
    for seq, sample in barcodes.items():
        command1=["cp",os.path.join(args.input,args.input_basename1.replace("SEQ",seq)),
            os.path.join(args.output,args.output_basename1.replace("SAMPLE",sample))]
        command2=["cp",os.path.join(args.input,args.input_basename2.replace("SEQ",seq)),
            os.path.join(args.output,args.output_basename2.replace("SAMPLE",sample))]
        try:
            print(command1)
            subprocess.check_output(command1)
            print(command2)
            subprocess.check_output(command2)
        except subprocess.CalledProcessError:
            print("Unable to find sample: {}".format(sample))
            not_found+=1

    print("Total not found: {}".format(not_found))

if __name__ == "__main__":
    main()

