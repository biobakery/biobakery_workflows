#!/usr/bin/env python

""" This script will create otu tables from usearch alignments """

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
UNNAMED_TAXONOMY="Unclassified"
FASTA_SEQ_START=">"

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

def create_fasta_files(targets, green_genes_fasta, nonchimera_fasta, output_open_ref_fasta, output_closed_ref_fasta):
    """ Created open and closed reference fasta files with the sequence ids matching those in the otu tables """

    # try to open the output files
    file_handle_write_open_ref=catch_open(output_open_ref_fasta,write=True)
    file_handle_write_closed_ref=catch_open(output_closed_ref_fasta,write=True)

    # read through the nonchimera fasta file, writing to the output files
    write_seq=False
    with catch_open(nonchimera_fasta) as file_handle_read:
        for line in file_handle_read:
            if line[0] == FASTA_SEQ_START:
                # check if this sequence should be written to either output file
                id=line.rstrip().replace(FASTA_SEQ_START,"")
                if id in targets:
                    # write this to the open reference file
                    file_handle_write_open_ref.write(line)
                    write_seq=True
                else:
                    write_seq=False
            else:
                if write_seq:
                    file_handle_write_open_ref.write(line)

    # read through the green genes fasta file, writing to the output files
    write_seq=False
    with catch_open(green_genes_fasta) as file_handle_read:
        for line in file_handle_read:
            if line[0] == FASTA_SEQ_START:
                # check if this sequence should be written to either output file
                id=line.rstrip().replace(FASTA_SEQ_START,"")
                if id in targets:
                    # write this to the open and closed reference files
                    file_handle_write_open_ref.write(line)
                    file_handle_write_closed_ref.write(line)
                    write_seq=True
                else:
                    write_seq=False
            else:
                if write_seq:
                    file_handle_write_open_ref.write(line)
                    file_handle_write_closed_ref.write(line)

def create_otu_table(taxonomy_file, samples, denovo_otu_table, green_genes_uc, out_tsv, filtered_out_tsv, reads_to_otus):
    """ Create open and closed reference otu tables """

    # read in the green genes taxonomy ids to names file
    green_genes_taxonomy_ids={}
    for otu_id, taxonomy_name in read_tab_file(taxonomy_file):
        green_genes_taxonomy_ids[otu_id] = taxonomy_name

    otu_table = {}
    targets = set()
    known_otus = set()
    for alignment_type, otu_id, target in read_uc_file(green_genes_uc, all_alignments=True):
        if otu_id in denovo_otu_table:
            hits = denovo_otu_table[otu_id]
            # check if this is an alignment hit
            if alignment_type == USEARCH_HIT:
                if target in green_genes_taxonomy_ids:
                    taxonomy_name = green_genes_taxonomy_ids[target]
                    known_otus.add(otu_id)
                else:
                    sys.exit("ERROR: Alignment to green genes id not included in taxonomy file:" + target)
            else:
                taxonomy_name = UNNAMED_TAXONOMY
                target = otu_id

            # store the otu mapping to target
            targets.add(target)

            # record hits by the otu id and also the taxonomy name
            full_id = (target, taxonomy_name)
            if not full_id in otu_table:
                otu_table[full_id] = hits
            else:
                # add the current hits to those already recorded
                otu_table[full_id] = list(map(add, otu_table[full_id], hits))

    # write the two output files
    with catch_open(out_tsv, write=True) as file_handle:
        with catch_open(filtered_out_tsv, write=True) as filtered_file_handle:
            header="\t".join(["# OTU"]+samples+["taxonomy"])+"\n"
            file_handle.write(header)
            filtered_file_handle.write(header)
            for (taxonomy_id, taxonomy_name), hits in otu_table.items():
                output_line="\t".join([taxonomy_id]+list(map(str, hits))+[taxonomy_name])+"\n"
                file_handle.write(output_line)
                # do not write unclassified output to filtered file
                if not taxonomy_name == UNNAMED_TAXONOMY:
                    filtered_file_handle.write(output_line) 

    # compute the counts of reads mapping to known (included in green genes) and unknown otus 
    sample_known_read_counts={sample:0 for sample in samples}
    sample_unclassified_read_counts={sample:0 for sample in samples}
    for read, otus in reads_to_otus.items():
        sample = get_sample_id(read)
        # check if any of the otus are known
        if known_otus.intersection(otus):
            sample_known_read_counts[sample]+=1
        else:
            sample_unclassified_read_counts[sample]+=1

    return targets, sample_known_read_counts, sample_unclassified_read_counts

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
    
def write_denovo_otu_table(nonchimera_uc_mapping, output_file):
    """ Read the chimera uc mapping file and write a denovo table """
    
    # read the uc mapping file
    table = {}
    samples = set()
    queries_to_otus = {}
    for query, otu in read_uc_file(nonchimera_uc_mapping):
        sample = get_sample_id(query)
        samples.add(sample)
        if not otu in table:
            table[otu]={}
        table[otu][sample] = table[otu].get(sample,0) + 1
        if not query in queries_to_otus:
            queries_to_otus[query]=set()
        queries_to_otus[query].add(otu)
    
    # write the denovo output file
    sample_list=list(samples)
    denovo_otu_table={}
    with catch_open(output_file, write=True) as file_handle:
        # write the header
        file_handle.write("\t".join(["# OTU"]+list(map(str,sample_list)))+"\n")
        for otu, hits_to_otu in table.items():
            hits_by_sample = [hits_to_otu.get(sample, 0) for sample in sample_list]
            file_handle.write("\t".join([otu]+list(map(str,hits_by_sample)))+"\n")
            denovo_otu_table[otu]=hits_by_sample
            
    return sample_list, denovo_otu_table, queries_to_otus
                    
def count_reads_per_sample(file):
    """ Count the reads for each sample from the original fasta file """

    samples={}
    for line in catch_open(file):
        if line.startswith(">"):
            try:
                sample, read = line.replace(">","").split(SAMPLE_READ_DELIMITER)
            except ValueError:
                print("Warning: Sequence id has unexpected format, not included in total read count: " + line)
                continue
            if not sample in samples:
                samples[sample]=set()
            samples[sample].add(read)

    sample_counts={sample:len(list(reads)) for sample, reads in samples.items()}

    return sample_counts

def write_read_count_table(output_file, original_counts, known_counts, unknown_counts):
    """ Write a table of read counts per sample """

    samples=original_counts.keys()
    with catch_open(output_file, write=True) as file_handle:
        file_handle.write("\t".join(["# sample","original read count","reads mapping to OTU with taxonomy","reads mapping to unclassifed OTU"])+"\n")
        for sample in samples:
            file_handle.write("\t".join([sample]+[str(i) for i in [original_counts[sample],known_counts.get(sample,0),unknown_counts.get(sample,0)]])+"\n")

def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Process usearch output to get taxonomy.")
    parser.add_argument('input_taxonomy',help="Taxonomy file from green genes.")
    parser.add_argument('input_greengenes_fasta',help="Fasta file from green genes.")
    parser.add_argument('input_greengenes_uc',help="The uc file of the alignment results from usearch to green genes.")
    parser.add_argument('input_nonchimera_uc',help="The uc file of alignment results from usearch to nonchimeras.")
    parser.add_argument('input_nonchimera_fasta',help="The uc nonchimeras fasta file.")
    parser.add_argument('input_original_fasta',help="The original fasta file containing all reads for all samples.")
    parser.add_argument('output_open_ref_tsv',help="The open ref otu table written.")
    parser.add_argument('output_open_ref_fasta',help="The open ref fasta written.")
    parser.add_argument('output_closed_ref_tsv',help="The closed reference otu table written.")
    parser.add_argument('output_closed_ref_fasta',help="The closed reference fasta written.")
    parser.add_argument('output_denovo_otu',help="The denovo otu table written.")
    parser.add_argument('output_read_count',help="The read count table written.")   
 
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
    check_file_nonempty(args.input_greengenes_uc)

    # get the original read counts
    original_read_counts = count_reads_per_sample(args.input_original_fasta)

    samples, denovo_otu_table, reads_to_otus = write_denovo_otu_table(args.input_nonchimera_uc,args.output_denovo_otu)
    targets, sample_known_read_counts, sample_unknown_read_counts=create_otu_table(args.input_taxonomy,
        samples,denovo_otu_table,args.input_greengenes_uc,args.output_open_ref_tsv,args.output_closed_ref_tsv, reads_to_otus)
    create_fasta_files(targets,args.input_greengenes_fasta,args.input_nonchimera_fasta,args.output_open_ref_fasta,
        args.output_closed_ref_fasta)

    write_read_count_table(args.output_read_count, original_read_counts, sample_known_read_counts, sample_unknown_read_counts)

if __name__ == "__main__":
    main()

    
    
