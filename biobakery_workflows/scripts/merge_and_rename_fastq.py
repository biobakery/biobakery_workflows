#!/usr/bin/env python

""" This script will merge two fastq files and will rename the sequences to the
    basename of the first file and the read number.
    It will first sort the sequences by ids before renaming. """

import sys

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to at least python v2.7")
    
import os
import re
import gzip
    
FASTQ_LINE1_START="@"
FASTQ_LINE2_REGEXP="^[A|a|T|t|G|g|C|c|N|n]+$"
FASTQ_LINE3_START="+"
NEW_SEQUENCE_NAME_DELIMITER="."
    
def catch_open(file,write=None):
    """ Try to open the file, return file handle, or exit if unable to open """

    open_function=open
    if file.endswith(".gz"):
        open_function=gzip.open

    try:
        if write:
            file_handle=open_function(file,"w")
        else:
            file_handle=open_function(file)
    except EnvironmentError:
        sys.exit("ERROR: Unable to open file: " + file)
        
    return file_handle

def read_file_n_lines(file,n):
    """ Read a file n lines at a time """
    
    line_set=[]
    with catch_open(file) as file_handle:
        for line in file_handle:
            if len(line_set) == n:
                yield line_set
                line_set=[]
            line_set.append(line)
    
    # yield the last set
    if len(line_set) == n:
        yield line_set
    
def fastq_format_error_message(lines):
    """ Return all lines with the error message on formatting """
    
    sys.exit("ERROR: Issue in fastq file format\nSequence lines:\n"+"".join(lines))
    
def fastq_rename(new_sequence_name,file,file_handle_write,read_count=None):
    """
    Read the fastq files and return renamed sequences
    Check the fastq file is formatted as 4 lines of sequence id, sequence, '+' line,
    and then quality values 
    """
    
    # set read count to 1 if not already set
    if not read_count:
        read_count=1
        
    # read the file 4 lines at a time
    # store the sequences
    store_sequences={}
    for lines in read_file_n_lines(file, 4):
        # check formatting is correct
        if not lines[0][0] == FASTQ_LINE1_START:
            fastq_format_error_message(lines)
        if not lines[2][0] == FASTQ_LINE3_START:
            fastq_format_error_message(lines)
            
        # only store sequences with proper characters
        if re.search(FASTQ_LINE2_REGEXP,lines[1]):
            store_sequences[lines[0]]="".join(lines[1:])

    # print sequences sorted
    for sequence_id in sorted(store_sequences.keys()):
        # rename the sequence id
        new_sequence_id=FASTQ_LINE1_START+new_sequence_name+NEW_SEQUENCE_NAME_DELIMITER+str(read_count)+'\n'
        read_count+=1
        
        file_handle_write.write(new_sequence_id+store_sequences[sequence_id])
        
    return read_count

def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Merge and rename fastq files.")
    parser.add_argument('input_paired_fastq',help="Paired fastq input file.")
    parser.add_argument('input_unpaired_fastq',help="Unpaired fastq input file.")
    parser.add_argument('input_remove_string',help="Remove string from file basename.")
    parser.add_argument('output_fastq',help="The merged renamed fastq output file.")
    
    return parser.parse_args()

def check_file_nonempty(file):
    """ Check the file is non-empty """
    
    try:
        size = os.stat(file).st_size
    except EnvironmentError:
        sys.exit("ERROR: Unable to access input file: " + file)
        
    if size == 0:
        sys.exit("ERROR: Input file is empty: " + file)

def main():
    # parse arguments
    args = parse_arguments(sys.argv)
    
    # get the new sequence name as the basename of the first input file
    input_file_basename=os.path.basename(args.input_paired_fastq)
    # remove gzip extension if present
    if input_file_basename.endswith(".gz"):
        input_file_basename=input_file_basename[:-3]    

    # remove the extension
    new_sequence_name=input_file_basename.replace(os.path.splitext(input_file_basename)[-1],"")
 
    # remove non-alphanumeric characters from name
    new_sequence_name=re.sub(r"[^a-zA-Z\d]","_",new_sequence_name)
 
    # remove string from end of file name if present
    if args.input_remove_string:
        new_sequence_name=new_sequence_name.replace(args.input_remove_string,"")
 
    # try to open the output file
    file_handle_write=catch_open(args.output_fastq, write=True)
    
    # write the paired reads to the output file with sequences renamed
    read_count=fastq_rename(new_sequence_name,args.input_paired_fastq,file_handle_write)
    
    # now write the unpaired reads
    if args.input_unpaired_fastq:
        read_count=fastq_rename(new_sequence_name,args.input_unpaired_fastq,file_handle_write,read_count)
    
    # close the output file
    file_handle_write.close()

    
if __name__ == "__main__":
    main()
