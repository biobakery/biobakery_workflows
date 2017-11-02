#!/usr/bin/env python

"""
This script will create fasta files for each genus/species taxonomy from the alignments 
These fasta output files will then be used as input for the oligotyping workflow.

To run:
$ python create_fasta_per_taxonomy_from_alignments.py all_samples_taxonomy_closed_reference.tsv all_samples_otu_mapping_results.uc  all_samples_green_genes_mapping_results.uc all_samples_truncated.fasta output_folder

"""

import sys

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to at least python v2.7")

import os
import re
from operator import add

USEARCH_HIT="H"
USEARCH_HIT_INDEX=0
USEARCH_QUERY_INDEX=8
USEARCH_TARGET_INDEX=9
SAMPLE_READ_DELIMITER="."
OLIGOTYPING_READ_DELIMITER="_"
UNNAMED_TAXONOMY="Unclassified"
FASTA_SEQ_START=">"

TAXONOMY_IDS={"order":"o__","family":"f__","genus":"g__","species":"s__"}

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

def get_sample_id(sample):
    """ Get the sample id from a string with the read number """
    return SAMPLE_READ_DELIMITER.join(sample.split(SAMPLE_READ_DELIMITER)[0:-1])

def create_fasta_files(matches, all_taxa, trimmed_fasta, output):
    """ Create fasta files organized by taxonomy """

    # write out all of the sequences for a single taxonomy at a time   
    for taxa in all_taxa:
        # try to open the output file
        current_output=os.path.join(output,taxa+".fasta")
        print("Writing file: " + current_output)
        file_handle_write=catch_open(current_output,write=True)
        write_seq=False
        with catch_open(trimmed_fasta) as file_handle_read:
            for line in file_handle_read:
                if line[0] == FASTA_SEQ_START:
                    id=line.rstrip().replace(FASTA_SEQ_START,"")
                    id_taxonomy=matches.get(id,"")
                    if id_taxonomy == taxa:
                        file_handle_write.write(line.replace(SAMPLE_READ_DELIMITER, OLIGOTYPING_READ_DELIMITER))
                        write_seq=True
                    else:
                        write_seq=False
                else:
                    if write_seq:
                        file_handle_write.write(line)
        file_handle_write.close()

def process_alignments(input_otu, input_nonchimera_uc, input_greengenes_uc,rank):
    """ Read in the alignment information to get taxonomy for each read """

    # read in the taxonomy from the otu table
    taxonomy={}
    header=""
    print("Reading otu table")
    for data in read_tab_file(input_otu):
        if not header:
            header=data[1:]
        else:
            # process the taxonomy for each otu 
            taxonomy_data=data[-1].replace(" ","").split(";")
            # filter to the rank provided 
            max_taxonomy=list(filter(lambda x: x.startswith(TAXONOMY_IDS[rank]), taxonomy_data))[0]
            max_taxonomy_index=taxonomy_data.index(max_taxonomy)
            # filter out those missing taxonomy information (ie filter out those that are "g__"
            # as these are "unclassified")
            if not (TAXONOMY_IDS[rank] == max_taxonomy):
                taxonomy[data[0]]="_".join(taxonomy_data[:(max_taxonomy_index+1)]).replace("[","").replace("]","")

    # read the green genes alignment
    otu_to_ggid={}
    print("Reading greengenes uc file")
    for otu_id, target in read_uc_file(input_greengenes_uc):
        otu_to_ggid[otu_id]=target

    # read the cluster alignment
    read_to_otuid={}
    print("Reading otu uc file")
    for read_id, otu_id in read_uc_file(input_nonchimera_uc):
        read_to_otuid[read_id]=otu_id

    # create a data set of read id to taxonomy
    matches={}
    all_taxa=set()
    print("Matching read ids to taxonomy")
    for readid,otuid in read_to_otuid.items():
        # get the greengenes id and taxonomy
        ggid=otu_to_ggid.get(otuid,"")
        taxon=taxonomy.get(ggid,"")

        if taxon:
            matches[readid]=taxon
            all_taxa.add(taxon)

    return matches, all_taxa

def read_tab_file(file):
    """ Read a tab delimited file """
    with catch_open(file) as file_handle:
        for line in file_handle:
            yield line.rstrip().split("\t")
    
def read_uc_file(uc_file, all_alignments=None):
    """ Read the uc file and return the type, query, and target """
    for data in read_tab_file(uc_file):
        try:
            query = data[USEARCH_QUERY_INDEX]
        except IndexError:
            query = ""
                  
        try:
            target = data[USEARCH_TARGET_INDEX]
        except IndexError:
            target = ""
                
        try:
            type = data[USEARCH_HIT_INDEX]
        except IndexError:
            type = ""
                
        if type and query and target:
            if all_alignments:
                yield type, query, target
            elif type == USEARCH_HIT:
                yield query, target
    
def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Process usearch output to get taxonomy.")
    parser.add_argument('input_otu',help="The closed reference otu table.")
    parser.add_argument('input_nonchimera_uc',help="The uc file of alignment results from usearch to nonchimeras.")
    parser.add_argument('input_greengenes_uc',help="The uc file of alignment results from nonchimeras to greengenes.")
    parser.add_argument('input_trimmed_fasta',help="The fasta file of trimmed reads from the samples.")
    parser.add_argument('--rank',help="The rank to user for the output fasta files.",choices=TAXONOMY_IDS.keys(),default="species")
    parser.add_argument('output_folder',help="The folder to write the per genus/species output files.")
    
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
    
    # check the input files are non-empty
    check_file_nonempty(args.input_nonchimera_uc)

    # if the output folder does not exist, then create it
    if not os.path.isdir(args.output_folder):
        try:
            print("Creating output directory: " + args.output_folder)
            os.mkdir(args.output_folder)
        except EnvironmentError:
            sys.exit("ERROR: Unable to create output directory")

    matches, all_taxa=process_alignments(args.input_otu,args.input_nonchimera_uc,args.input_greengenes_uc,args.rank)
    create_fasta_files(matches,all_taxa,args.input_trimmed_fasta,args.output_folder)

if __name__ == "__main__":
    main()

    
    
