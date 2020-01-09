#!/usr/bin/env python

""" This script will take as input paired files of sequences which are already 
demultiplexed and split them into files for each sample name using the 
corresponding index file and barcodes file.
"""

import argparse
import sys
import os
import string

def parse_arguments(args):
    """ Parse the arguments from the user """
    
    parser = argparse.ArgumentParser(
        description= "Demultiplex files using barcodes and index (sequences no longer include barcodes)\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input-index", 
        help="input index file \n[REQUIRED if barcodes are not in read ids]")
    parser.add_argument(
        "--input-read1", 
        help="input read1 file \n[REQUIRED]", 
        required=True)
    parser.add_argument(
        "--input-read2", 
        help="input read2 file \n[REQUIRED]", 
        required=True)
    parser.add_argument(
        "--input-barcodes", 
        help="input barcodes file \n[REQUIRED]", 
        required=True)
    parser.add_argument(
        "--reverse-complement", 
        help="use reverse complements of index sequences\n", 
        action="store_true")
    parser.add_argument(
        "-o", "--output", 
        help="directory to write output files\n[REQUIRED]", 
        metavar="<output>", 
        required=True)
    
    return parser.parse_args()

def read_fastq_file(file):
    """ Read a fastq file (not compressed) """
    
    line_set=[]
    with open(file) as file_handle:
        for line in file_handle:
            if len(line_set) == 4:
                yield line_set
                line_set=[]
            line_set.append(line)
            
def read_barcode_file(file):
    """ Read through the barcode file storing samples and sequences """
    
    barcodes={}
    with open(file) as file_handle:
        header=file_handle.readline()
        for line in file_handle:
            data=line.rstrip().split("\t")
            # replace periods with underscore in sample ids
            data[0]=data[0].replace(".","_")
            if data[1] in barcodes and data[0] != barcodes[data[1]]:
                sys.exit("ERROR: Duplicate barcode with different sample match: "+ data[1])
            barcodes[data[1]]=data[0]
            
    return barcodes
            
def main():
    # parse arguments from the user
    args = parse_arguments(sys.argv)
    
    # if the output folder does not exist then create it
    if not os.path.isdir(args.output):
        try:
            os.makedirs(args.output)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create output directory: " + args.output)
    
    # read through the barcode file
    barcodes=read_barcode_file(args.input_barcodes)
    
    # create rev cmp function
    rev_cmp=string.maketrans("ACGT","TGCA")
    
    # read through read1/2 files pulling a sequence set each time
    read1_lines = read_fastq_file(args.input_read1)
    read2_lines = read_fastq_file(args.input_read2)
    
    # read through the index and pair files matching up sequences to barcodes
    missing_barcodes=0
    total_barcodes=0
    read_counts_by_sample={}
    new_files=set()
    if args.input_index:
        all_index_lines = read_fastq_file(args.input_index)
    else:
        all_index_lines = read_fastq_file(args.input_read1)

    for index_lines in all_index_lines:
        total_barcodes+=1
       
        if args.input_index: 
            sequence=index_lines[1].rstrip()
        else:
            sequence=index_lines[0].rstrip().split(":")[-1]
        
        # reverse complement
        if args.reverse_complement:
            sequence=sequence.translate(rev_cmp)[::-1]
        
        try:
            sample_id=barcodes[sequence]
        except KeyError:
            missing_barcodes+=1
            sample_id="unknown"
            
        # check that read1/2 have the same sequence id
        read1=next(read1_lines)
        read1_id = read1[0].split(" ")[0]
        read2=next(read2_lines)
        read2_id = read2[0].split(" ")[0]

        # get sequence id depending on if index is provided
        if args.input_index:
            sequence_id=index_lines[0].split(" ")[0]
        else:
            sequence_id=read1_id
        
        if not ( (read1_id == read2_id) and (read2_id == sequence_id) ):
            sys.exit("Reads are not ordered: " + read1_id)
            
        # increase read counts by sample
        read_counts_by_sample[sample_id]=read_counts_by_sample.get(sample_id,0)+1
        
        # write pairs to new files
        new_read1_file=os.path.join(args.output,sample_id+"_R1.fastq")
        new_read2_file=os.path.join(args.output,sample_id+"_R2.fastq")
        
        # add to the list of new files
        new_files.add(new_read1_file)
        new_files.add(new_read2_file)
        
        # open each file in append mode
        # can't keep all files open as this might be too many open files
        # depending on the total number of samples
        with open(new_read1_file, "a") as fh:
            fh.write("".join(read1))
            
        with open(new_read2_file, "a") as fh:
            fh.write("".join(read2))
        
    # write the read counts to the file
    with open(os.path.join(args.output,"read_counts.txt"),"w") as fh:
        fh.write("# sample\ttotal_reads\n")
        for sample, count in read_counts_by_sample.items():
            fh.write("{}\t{}\n".format(sample,count))
        
    print("Total barcodes: {}".format(total_barcodes))
    print("Unknown barcodes: {} {}%".format(missing_barcodes,missing_barcodes/(total_barcodes*1.0)*100))
    print("Files written: "+"\n".join(sorted(list(new_files))))

if __name__ == "__main__":
    main()
