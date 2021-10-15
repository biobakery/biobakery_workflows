#!/usr/bin/env python

""" This script will subset a tab-delimited file to only include a subset of ids. """

import sys
import argparse
import subprocess
import os
    
def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Sub-set tab-delimited file")
    parser.add_argument('--input',help="The input file to subset",required=True)
    parser.add_argument('--input-ids',help="The file of sample ids to use",required=True)
    parser.add_argument('--output',help="The output subset file",required=True)
    parser.add_argument('--exclude-id',help="The ids to exclude (only for samples as columns)",default="")
    parser.add_argument('--exclude-column',help="The columns to exclude (only for input files with samples as rows)",default="")
    parser.add_argument('--exclude-row-name',help="The rows to exclude (only for input files with samples as columns)",default="")
    
    return parser.parse_args()

def main():
    # parse arguments
    args = parse_arguments(sys.argv)

    # create output dir if needed
    output_dir=os.path.dirname(os.path.abspath(args.output))
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # check if the file exists
    if not os.path.isfile(args.input):
        sys.exit("ERROR: Input file does not exist: "+args.input)

    # check to see if the file has samples as columns or rows
    try:
        subprocess.check_output("head -n 1 {0} | grep -f {1}".format(args.input,args.input_ids),shell=True)
        samples_as_columns=True
    except subprocess.CalledProcessError:
        samples_as_columns=False
    
    # add exclude of row name if provided
    exclude_row=""
    if args.exclude_row_name:
        exclude_row="grep -v '{0}' | "

    exclude_id_grep=""
    if args.exclude_id:
        exclude_id_grep=" | grep -v "+args.exclude_id+" "

    if samples_as_columns:
        cmmd="COLUMNS=$(head -n 1 "+args.input+" | " + exclude_row + "tr '\\t' '\\n' | cat -n | grep -f "+args.input_ids+exclude_id_grep+" | awk -F' ' '{ print $1}' | tr '\\n' ',' | sed 's/,$//' ) && cut -f '1,'$COLUMNS  "+args.input+" > "+args.output
    else:
        if args.exclude_column:
            cmmd="cut --complement -f {0} {1} | head -n 1 > {3} && cut --complement -f {0} {1} | grep -f {2} >> {3}".format(args.exclude_column,args.input,args.input_ids,args.output)
        else:
            cmmd="cat {0} | head -n 1 > {2} && cat {0} | grep -f {1} >> {2}".format(args.input,args.input_ids,args.output)

    print("Running command to subset file: "+args.input)
    print(cmmd)
    subprocess.check_call(cmmd,shell=True)
    print("Subset file created: "+args.output)
 
if __name__ == "__main__":
    main()
