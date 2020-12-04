#!/usr/bin/env python

"""
To Run: 
$ ./join_taxonomic_profiles.py -i <input_dir> -o <data_table.tsv>

"""

import argparse
import sys
import os
import re

from humann.tools import util

TAXON_TABLE_DELIMITER="\t"
TAXON_TABLE_COMMENT="#"

def join_data_tables(data_tables,output,verbose=None):
    """
    Join the taxon tables to a single taxon table
    """
    
    sample_data_store={}
    file_basenames=[]
    taxa_store=set()
    for data_table in data_tables:
        
        if verbose:
            print("Reading file: " + data_table)
        
        # get the basename of the file
        file_basename='.'.join(os.path.basename(data_table).split('.')[:-1])
        file_basenames.append(file_basename)
        
        sample_data_store[file_basename]={}
        for line in open(data_table):
            if not line.startswith(TAXON_TABLE_COMMENT):
                data=line.rstrip().split(TAXON_TABLE_DELIMITER)
                taxon=data[0]
                if data[-1].replace(".","").replace("e-","").isdigit():
                    read_percent=data[-1]
                else:
                    read_percent=data[-2]

                taxa_store.add(taxon)
                sample_data_store[file_basename][taxon]=read_percent
            
    # write the joined table
    sample_header=["# taxonomy "]+file_basenames
    try:
        file_handle=open(output,"w")
        file_handle.write(TAXON_TABLE_DELIMITER.join(sample_header)+"\n")
    except EnvironmentError:
        sys.exit("Unable to write file: " + output)  
        
    for taxon in sorted(taxa_store, key=len):
        new_line=[taxon]
        for sample_name in file_basenames:
            new_line.append(sample_data_store[sample_name].get(taxon,"0"))
        file_handle.write(TAXON_TABLE_DELIMITER.join(new_line)+"\n")
    
    file_handle.close()

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    
    parser = argparse.ArgumentParser(
        description= "Join taxonomy tables\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-v","--verbose", 
        help="additional output is printed\n", 
        action="store_true",
        default=False)
    parser.add_argument(
        "-i","--input",
        help="the directory of tables\n",
        required=True)
    parser.add_argument(
        "-o","--output",
        help="the table to write\n",
        required=True)
    parser.add_argument(
        "--file_name",
        help="only join tables with this string included in the file name")

    return parser.parse_args()


def main():
    # Parse arguments from command line
    args=parse_arguments(sys.argv)
    
    # check for format of the taxon tables
    input_dir=os.path.abspath(args.input)
    
    # check the directory exists
    if not os.path.isdir(input_dir):
        sys.exit("The input directory provided can not be found." + 
            "  Please enter a new directory.")
    
    data_tables=[]
    file_list=os.listdir(input_dir)
   
    reduced_file_list=[]
    # filter out files which do not meet the name requirement if set
    if args.file_name:
        for file in file_list:
            if re.search(args.file_name,file):
                reduced_file_list.append(file)
    else: 
        for file in file_list:
            # ignore dot files, like ".DS_Store" on Apple OS X
            if file[0] == ".":
                if args.verbose:
                    print("Not including file in input folder: " + file)
            else:
                reduced_file_list.append(file)
    file_list=reduced_file_list
            
    args.output=os.path.abspath(args.output)
    output_dir=os.path.dirname(args.output)
    
    for file in file_list:
        data_tables.append(os.path.join(input_dir,file))
            
    # sort the data tables so they are in the same order on all platforms
    data_tables.sort()
        
    # split the data table
    if data_tables:
        if args.verbose:
            print("Joining taxonomy table")
            
        join_data_tables(data_tables,args.output,verbose=args.verbose)
                
        print("Taxonomy table created: " + args.output)
    else:
        print("Zero taxonomy tables were found to join.")

if __name__ == "__main__":
    main()
