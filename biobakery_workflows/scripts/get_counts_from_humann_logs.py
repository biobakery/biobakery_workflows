#!/usr/bin/env python

# This script will create a table of reads from humann log files. The table will have 
# total reads, unaligned after nucleotide search, and unaligned after translated search.
# The table will also include the total species number.

import sys
import os
import argparse

TOTAL_COUNT_TAG="reads; of these:"
NUCLEOTIDE_COUNT_TAG="Unaligned reads after nucleotide alignment:"
TRANSLATED_COUNT_TAG="Unaligned reads after translated alignment:"
SPECIES_COUNT_TAG="Total species selected from prescreen:"

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Reads the HUMAnN2 logs and prints a table of read counts\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input", 
        help="the folder of log files\n[REQUIRED]", 
        metavar="<input_folder>", 
        required=True)
    parser.add_argument(
        "-o", "--output", 
        help="file to write counts table\n[REQUIRED]", 
        metavar="<output>", 
        required=True)

    return parser.parse_args()

def main():
    # parse arguments from the user
    args = parse_arguments(sys)

    # search subfolders to find the log files
    log_files=[]
    for path, directories, files in os.walk(args.input):
        for file in filter(lambda x: x.endswith(".log"), files):
            log_files.append(os.path.join(path,file))

    try:
        file_handle=open(args.output,"w")
    except EnvironmentError:
        sys.exit("Error: Unable to open output file: " + args.output)

    # write out the header
    file_handle.write("\t".join(["# samples","total reads","total nucleotide aligned","total translated aligned","total species"])+"\n")

    for file in log_files:
        sample=os.path.basename(file).split(".log")[0]
        data=[sample,"NA","NA","NA","NA"]
        for line in open(file):
            if TOTAL_COUNT_TAG in line:
                data[1]=int(line.split()[7][2:])
            elif NUCLEOTIDE_COUNT_TAG in line:
                try:
                    data[2]=int(data[1]*((100-float(line.split()[-2]))/100.0))
                except (ValueError, TypeError):
                    print("Warning: Unable to compute nucleotide reads aligned from log line: " + line)
            elif TRANSLATED_COUNT_TAG in line:
                try:
                    data[3]=int(data[1]*((100-float(line.split()[-2]))/100.0))
                except (ValueError, TypeError):
                    print("Warning: Unable to compute translated reads aligned from log line: " + line)
            elif SPECIES_COUNT_TAG in line:
                data[4]=line.split()[-1]

        file_handle.write("\t".join([str(i) for i in data])+"\n")

    file_handle.close()
    print("Read table written to file: " + args.output)

if __name__ == "__main__":
    main()

