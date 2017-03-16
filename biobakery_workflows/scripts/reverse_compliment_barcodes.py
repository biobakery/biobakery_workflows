#!/usr/bin/env python

# This script will print out a bar codes file with reverse complement sequences.

import sys
import os
import string
import argparse

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Reverse complement barcodes (add another group of barcodes to the file)\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i","--input",
        help="the barcodes tsv file (format: <id>\t<barcode>)\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the new barcodes tsv file (format: <id>\t<barcode>\t<group>)\n",
        required=True)

    return parser.parse_args()

def main():
    # parse the arguments
    args=parse_arguments(sys.argv)
    
    input=os.path.abspath(args.input)
    output=os.path.abspath(args.output)
    
    if not os.path.isfile(input):
        sys.exit("The input file provided can not be found: " + 
            input+"  Please provide the full path or enter a new file.")
    # try to open the input file
    try:
        file_handle_read=open(input,"rt")
    except EnvironmentError:
        sys.exit("Unable to open input file.")

    # try to open the output file
    try:
        file_handle_write=open(output,"wt")
    except EnvironmentError:
        sys.exit("Unable to open output file.")

    rev_cmp=string.maketrans("ACGT","TGCA")

    newlines=[]
    # write the header
    header=file_handle_read.readline()
    file_handle_write.write(header)
    for line in file_handle_read:
        data=line.rstrip().split("\t")
        file_handle_write.write("\t".join([data[0],data[1],"original"])+"\n")
        # reverse complement the bar code
        newlines.append("\t".join([data[0],data[1].translate(rev_cmp)[::-1],"reverse_complement"]))
        
    # print the reverse complements
    for line in newlines:
        file_handle_write.write(line+"\n")

    file_handle_write.close()
    file_handle_read.close()

if __name__ == "__main__":
    main()

