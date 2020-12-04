#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
import random

try:
    from humann import config
    from humann import store
except ImportError:
    sys.exit("Please install humann.")

# This script will take as input the bowtie2 alignment file from
# running humann along with a list of species and pathways. The selection list
# will be formatted as two tab delimited columns of species \t pathway. 
# To include all pathways for a species use the key "ALL". It can
# have headers or comments starting with "#". The script will output a 
# reduced fasta (or fastq if selected) file that only includes reads that 
# will map to those species and pathways when running the fastq file as input to humann.

# Please note since some pathways have reactions that overlap with other pathways
# more than just the selected pathways may appear for a specific species. 

# If running humann with the output from the script does not identify any species
# in the metaphlan step, recreate the humann input file with this script adding
# the option "--add-unintegrated". This option will add reads that map to a species
# in the humann nucleotide alignment step but the gene families associated with these
# alignments are not included in any pathways (and are counted as "unintegrated" in 
# the humann pathway abundance output file).

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
    parser.add_argument(
        "--percent",
        help="the percent of reads to output (subsampled by reaction set for requested pathways)\n[OPTIONAL]",
        type=int,
        default=100)
    parser.add_argument(
        "--add-markers",
        help="add in the reads that align to markers for the selected species using MetaPhlAn2\n[OPTIONAL]",
        action="store_true")
    parser.add_argument(
        "--min-markers",
        help="the minimum number of reads that map to markers for each species using MetaPhlAn2\n[OPTIONAL]",
        type=int,
        default=20)
    parser.add_argument(
        "--min-reads-per-marker",
        help="the minimum number of reads per marker for each species using MetaPhlAn2\n[OPTIONAL]",
        type=int,
        default=5)
    parser.add_argument(
        "--add-trimmable",
        help="add reads that would be trimmed based on length\n[OPTIONAL]",
        action="store_true")

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

def read_sam(file, gene_families_to_use):
    """ Read the sam file, only reading alignments """

    for line in open(file):
        data=line.rstrip().split("\t")
        if len(data) > config.sam_read_quality:
            read_name=data[config.sam_read_name_index]
            reference_name=data[config.sam_reference_index]
            sequence=data[config.sam_read_index]
            quality_scores=data[config.sam_read_quality]

            reference_data = reference_name.split(config.chocophlan_delimiter)
            try:
                uniref50=reference_data[-2]
                uniref90=reference_data[-3]
                species=reference_data[-4].split(".")[-1]
            except IndexError:
                # ignore reads that do not map to the reference
                continue
            
            gene_family=uniref90
            if gene_families_to_use == "UniRef50":
                gene_family=uniref50

            yield species, gene_family, read_name, sequence, quality_scores

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
    reactions={}
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
                reactions[species]=set()
            # get the reactions and then gene families for the pathway
            if pathway.upper() == "ALL":
                genefamilies[species]="ALL"
                reactions[species]="ALL"
            else:
                for reaction in pathways_database.find_reactions(pathway):
                    reactions[species].add(reaction)
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
    print("Reading sam file")
    selected_reads=set()
    reads_to_gene_familes={}
    reaction_totals={}
    
    # open a fasta file to write the species specific sequences for input to metaphlan
    if args.add_markers:
        species_fasta_file=args.output+".species_specific_reads.fasta"
        species_fasta_file_handle=open(species_fasta_file,"w")
        
    selected_reads_to_species={}
    for species, gene_family, read_name, sequence, quality_scores in read_sam(args.input_sam, args.gene_families):
        # record the reads that align to species/gene families requested
        requested_genes_for_species=genefamilies.get(species,[])
        if (gene_family in requested_genes_for_species) or ("ALL" in requested_genes_for_species) or (args.add_unintegrated and not(gene_family in all_genefamilies_in_pathways)):
            selected_reads.add(read_name)
            selected_reads_to_species[species]=selected_reads_to_species.get(species,0)+1
            if not read_name in reads_to_gene_familes:
                reads_to_gene_familes[read_name]=set()
            
            reads_to_gene_familes[read_name].add(gene_family)
            
            # add to the reaction totals count
            for reaction in reactions_database.find_reactions(gene_family):
                reaction_totals[reaction]=reaction_totals.get(reaction,0)+1.0
                
        # if this is a species in the list, and we are adding markers, write this to the fasta file for input to metaphlan
        if args.add_markers and species in genefamilies.keys():
            write_sequence(species_fasta_file_handle, read_name, sequence, quality_scores, "fasta")
           
    print("Total reads per species based on gene families")
    for species, count in selected_reads_to_species.items():
        print(species+"\t"+str(count))         

    # close the species fasta file handle
    if args.add_markers:
        species_fasta_file_handle.close()
            
    print("Total reads found: " + str(len(selected_reads)))
    
    # get the markers reads from the species fasta file
    marker_reads_to_add=set()
    read_to_species_marker={}
    read_to_marker_name={}
    all_reads_mapping_to_markers=set()
    if args.add_markers:
        print("Finding reads mapped to markers with MetaPhlAn2")
        
        print("Running MetaPhlAn2")
        metaphlan_marker_file=species_fasta_file+".marker_alignments.tsv"
        metaphlan_bowtie2_file=species_fasta_file+".bowtie2.tsv"
        
        try:
            # remove the bowtie2 output file if it exists to prevent metaphlan error
            os.remove(metaphlan_bowtie2_file)
        except EnvironmentError:
            pass
        
        output=subprocess.check_output(["metaphlan","--input_type","fasta",species_fasta_file,"-t","reads_map","-o",metaphlan_marker_file,"--bowtie2out",metaphlan_bowtie2_file])
        
        # read through the file to identify the maker reads to add
        for line in open(metaphlan_marker_file):
            if not line.startswith("#"):
                read_name, taxon = line.rstrip().split("\t")
                species = taxon.split("|")[-1]
                if species in genefamilies.keys():
                    marker_reads_to_add.add(read_name)
                    if not species in read_to_species_marker:
                        read_to_species_marker[species]=set()
                        
                    read_to_species_marker[species].add(read_name)
                    all_reads_mapping_to_markers.add(read_name)
                    
        # read through the file to identify the name of the marker the reads map to
        for line in open(metaphlan_bowtie2_file):
            try:
                read_name, marker_name = line.rstrip().split("\t")
            except IndexError:
                continue
            
            if read_name in all_reads_mapping_to_markers:
                read_to_marker_name[read_name]=marker_name
                    
        print("Found a total of "+str(len(marker_reads_to_add))+" reads to species markers")
    
    if args.percent < 100:
        print("Filtering reads by percent requested: " + str(args.percent))
        filtered_reads=set()
        current_reaction_totals={reaction:0 for reaction in reaction_totals.keys()}
        # go through the reads, adding until there is enough in the reaction list
        # to hit the percent requested
        for read_name in selected_reads:
            # check if this read is needed into increase the reaction percents
            add_read=False
            for gene in reads_to_gene_familes[read_name]:
                for reaction in reactions_database.find_reactions(gene):
                    if (current_reaction_totals[reaction]/reaction_totals[reaction])*100 < args.percent:
                        add_read=True
                        
            if add_read:
                filtered_reads.add(read_name)
                # update the current reaction counts
                for gene in reads_to_gene_familes[read_name]:
                    for reaction in reactions_database.find_reactions(gene):
                        current_reaction_totals[reaction]+=1
                        
        # update the selected reads to those that are filtered
        selected_reads=filtered_reads
        print("Total reads after filtering: "+str(len(selected_reads)))
    
    # check to make sure at least 200 marker reads (or command line setting) for each species are present
    # make sure there reads are spread over the markers for the sample
    if args.add_markers:
        # for each species, check the number of maker reads present
        print("Counting markers in set to determine if reads need to be added to meet min markers")
        for species, reads_for_species in read_to_species_marker.items():
            # count how many of the reads are in the selected reads list
            overlap=selected_reads.intersection(reads_for_species)
            total_overlap=len(list(overlap))
            if total_overlap < args.min_markers:
                # add in more reads to hit the min markers setting for this species
                max_reads_to_add=list(reads_for_species.difference(selected_reads))
                # group the reads based on the species markers
                to_add_by_markers={}
                for read_name in max_reads_to_add:
                    marker_for_read = read_to_marker_name[read_name]
                    if not marker_for_read in to_add_by_markers:
                        to_add_by_markers[marker_for_read]=set()
                    to_add_by_markers[marker_for_read].add(read_name)
                    
                selected_reads_by_markers={}
                for read_name in selected_reads:
                    try:
                        marker_for_read = read_to_marker_name[read_name]
                    except KeyError:
                        # ignore reads that do not map to markers
                        continue
                    if not marker_for_read in selected_reads_by_markers:
                        selected_reads_by_markers[marker_for_read]=set()
                    selected_reads_by_markers[marker_for_read].add(read_name)
                
                # count the total reads to add    
                total_to_add=args.min_markers - total_overlap
                total_added=0
                
                # get a list of all of the markers for the species
                all_species_markers=set(to_add_by_markers.keys())
                all_species_markers.update(selected_reads_by_markers.keys())
                
                # start with the markers not already included, adding one read to each
                for marker_name in all_species_markers.difference(selected_reads_by_markers.keys()):
                    total_markers_in_set=len(list(to_add_by_markers[marker_name]))
                    end_index = args.min_reads_per_marker if total_markers_in_set >= args.min_reads_per_marker else total_markers_in_set
                    
                    # add at most min reads per markers for this marker set to the set of selected reads
                    to_add=list(to_add_by_markers[marker_name])[:end_index]
                    selected_reads.update(to_add)
                    if not marker_name in selected_reads_by_markers:
                        selected_reads_by_markers[marker_name]=set()
                    selected_reads_by_markers[marker_name].update(to_add)
                    
                    # remove the added reads from the set of reads to add
                    to_add_by_markers[marker_name]=to_add_by_markers[marker_name].difference(to_add)
                    
                    # decrease the total amount to add
                    total_added+=end_index
                    
                    print("Adding " +str(end_index)+ " total reads to fill empty marker for species " + species)
                    
                # next add more reads for those sets that already have markers, starting with the smallest
                # set of markers to the largest
                marker_counts={marker:len(list(reads)) for marker, reads in selected_reads_by_markers.items()} 
                for marker_name in sorted(marker_counts, key=marker_counts.get):
                    try:
                        total_markers_in_set=len(list(to_add_by_markers[marker_name]))
                    except KeyError:
                        # ignore errors for markers that are not included in the to add list (only in the selected reads)
                        continue
                    
                    end_index = args.min_reads_per_marker if total_markers_in_set >= args.min_reads_per_marker else total_markers_in_set
                    
                    # add only so many reads to get to min reads per markers for this set
                    end_index=end_index-len(list(selected_reads_by_markers[marker_name]))
                    
                    # add at most min reads per marker reads for this marker set to the set of selected reads
                    if end_index > 0:
                        to_add=list(to_add_by_markers[marker_name])[:end_index]
                        selected_reads.update(to_add)
                        selected_reads_by_markers[marker_name].update(to_add)
                            
                        # remove the added reads from the set of reads to add
                        to_add_by_markers[marker_name]=to_add_by_markers[marker_name].difference(to_add)
                            
                        total_added+=end_index
                        print("Adding " +str(end_index)+ " total reads to fill slightly full marker for species " + species)
                    
                # add more reads to get to the min reads added value
                total_to_add=total_to_add-total_added if total_added < total_to_add else 0
                try:
                    selected_reads.update(max_reads_to_add[:total_to_add])
                    total_added+=total_to_add
                    print("Added "+str(total_to_add)+" reads for species "+species+" to meet min markers")
                except (TypeError, IndexError):
                    continue
                
            print("Total reads added for species "+species+" :"+str(total_added)) 
        print("Total reads after counting markers: " + str(len(selected_reads)))
   
    # determine how many trimmable reads to write
    total_trimmable = 0
    trimmable = []
    if args.add_trimmable:
        total_trimmable = int(random.uniform(50,100)) 
        print("Adding "+str(total_trimmable)+" total trimmable reads")
 
    print("Writing output file")
    for species, gene_family, read_name, sequence, quality_scores in read_sam(args.input_sam, args.gene_families):
        # write the read sequences requested once
        if read_name in selected_reads:
            write_sequence(file_handle, read_name, sequence, quality_scores, args.output_format)
            selected_reads.remove(read_name)
            if len(trimmable) < total_trimmable:
                trimmable.append([read_name, sequence, quality_scores])

    # write the trimmable reads
    for read_name, sequence, quality_scores in trimmable:
        new_length=int(len(sequence)/2.0)
        write_sequence(file_handle, "trimmable_"+read_name, sequence[0:new_length], quality_scores[0:new_length], args.output_format)

    print("Output file written: " + args.output)
    
if __name__ == "__main__":
    main()
