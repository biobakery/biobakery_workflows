#!/usr/bin/env python

import sys
import os
import argparse

try:
    from humann2 import config
    from humann2 import store
except ImportError:
    sys.exit("Please install humann2.")

# This script will take as input the bowtie2 alignment file from
# running humann2 along with a list of species and pathways. The selection list
# will be formatted as two tab delimited columns of species \t pathway. It can
# have headers or comments starting with "#". The script will output a 
# reduced fasta (or fastq if selected) file that only includes reads that 
# will map to those species and pathways when running the fastq file as input to humann2.

# Please note since some pathways have reactions that overlap with other pathways
# more than just the selected pathways may appear for a specific species. 

# If running humann2 with the output from the script does not identify any species
# in the metaphlan2 step, recreate the humann2 input file with this script adding
# the option "--add-unintegrated". This option will add reads that map to a species
# in the humann2 nucleotide alignment step but the gene families associated with these
# alignments are not included in any pathways (and are counted as "unintegrated" in 
# the humann2 pathway abundance output file).

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Create demo file\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--input-sam",
        help="the sam bowtie2 alignment temp file from running HUMAnN2\n[REQUIRED]",
        metavar="<input.sam>",
        required=True)
    parser.add_argument(
        "--input-selection",
        help="the species and pathways requested\n[REQUIRED]",
        metavar="<input.tsv>",
        required=True)
    parser.add_argument(
        "--add-unintegrated",
        help="add reads that do not map to pathways for each selected species\n[OPTIONAL]",
        action="store_true")
    parser.add_argument(
        "--output",
        help="file to write the sub-sampled file\n[REQUIRED]",
        metavar="<output.fastq>",
        required=True)
    parser.add_argument(
        "--output-format",
        help="the format to write the output file\n[OPTIONAL]",
        choices=["fastq","fasta"],
        default="fasta")
    parser.add_argument(
        "--gene-families",
        help="the gene families to select with\n[OPTIONAL]",
        choices=["UniRef50","UniRef90"],
        default="UniRef90")

    return parser.parse_args()

def write_sequence(file_handle, read_name, sequence, quality_scores, output_format):
    """ Write the sequence to the output file in the selected format """
    
    if output_format == "fasta":
        file_handle.write(">"+read_name+"\n")
        file_handle.write(sequence+"\n")
    else:
        file_handle.write("@"+read_name+"\n")
        file_handle.write(sequence+"\n")
        file_handle.write("+\n")
        file_handle.write(quality_scores+"\n")

def main():

    args=parse_arguments(sys)
    
    # read in the gene families to reactions database
    print("Reading gene families to reactions database")
    reactions_database=store.ReactionsDatabase(config.pathways_database_part1)
    
    # read in the reactions to pathways database
    print("Reading reactions to pathways database")
    pathways_database=store.PathwaysDatabase(config.pathways_database_part2, reactions_database)
    
    # read in the species and pathways selected
    print("Finding gene families for each species for the pathways selected")
    genefamilies={}
    for line in open(args.input_selection):
        # remove starting and ending spaces
        line=line.strip()
        if not line.startswith("#"):
            try:
                species, pathway = line.split("\t")
            except IndexError:
                print("Warning: Skipping selection line because of format: " + line)
                continue
            if not species in genefamilies:
                genefamilies[species]=set()
            # get the reactions and then gene families for the pathway
            for reaction in pathways_database.find_reactions(pathway):
                genefamilies[species].update(list(filter(lambda x: x.startswith(args.gene_families), reactions_database.find_genes(reaction))))
     
    # get the set of gene families in any pathway           
    all_genefamilies_in_pathways=set()
    for reaction in pathways_database.reaction_list():
        all_genefamilies_in_pathways.update(list(filter(lambda x: x.startswith(args.gene_families), reactions_database.find_genes(reaction))))

    # open the output file
    try:
        file_handle=open(args.output,"w")
    except EnvironmentError:
        sys.exit("ERROR: Unable to open output file: " + args.output)
                
    # find the reads for the species and pathways selected
    print("Reading sam file and writing output file")
    for line in open(args.input_sam):
        data=line.rstrip().split("\t")
        if len(data) > config.sam_read_quality:
            read_name=data[config.sam_read_name_index]
            reference_name=data[config.sam_reference_index]
            sequence=data[config.sam_read_index]
            quality_scores=data[config.sam_read_quality]

            # check for a requested species and gene family
            reference_data = reference_name.split(config.chocophlan_delimiter)
            try:
                uniref50=reference_data[-2]
                uniref90=reference_data[-3]
                species=reference_data[-4].split(".")[-1]
            except IndexError:
                # ignore reads that do not map to the reference
                continue
            
            selected_gene_family=uniref90
            if args.gene_families == "UniRef50":
                selected_gene_family=uniref50

            requested_genes_for_species=genefamilies.get(species,[])
            if selected_gene_family in requested_genes_for_species:
                # print this sequence to the output file
                write_sequence(file_handle, read_name, sequence, quality_scores, args.output_format)                  
            elif args.add_unintegrated and not(selected_gene_family in all_genefamilies_in_pathways):  
                write_sequence(file_handle, read_name, sequence, quality_scores, args.output_format)                  

    print("Output file written: " + args.output)
    
if __name__ == "__main__":
    main()
