#!/usr/bin/env python

""" This script will merge fastq files with the 
    given string in the folder provided """

import sys

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to at least python v2.7")
    
import os
    
def catch_open(file,write=None):
    """ Try to open the file, return file handle, or exit if unable to open """
    try:
        if write:
            file_handle=open(file,"w")
        else:
            file_handle=open(file)
    except EnvironmentError:
        sys.exit("ERROR: Unable to open file: " + file)
        
    return file_handle

def write_file(infile,outfile):
    """ Append the outfile with the line in infile """
    
    with catch_open(infile) as file_handle:
        for line in file_handle:
            outfile.write(line)
    
def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Merge fastq files.")
    parser.add_argument('input_folder',help="The input folder with the fastq files.")
    parser.add_argument('input_filenames',help="The search string for file names.")
    parser.add_argument('output_fastq',help="The merged fastq output file.")
    
    return parser.parse_args()

def main():
    # parse arguments
    args = parse_arguments(sys.argv)
    
    # get a list of all files in the input folder to merge
    input_merge_files=filter(lambda x: args.input_filenames in x, os.listdir(args.input_folder))

    # try to open the output file
    with catch_open(args.output_fastq, write=True) as file_handle_write:
        for file in input_merge_files:
            write_file(os.path.join(args.input_folder,file),file_handle_write)
 
if __name__ == "__main__":
    main()
