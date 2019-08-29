#!/usr/bin/env python

import sys
import os
import argparse

# This script will add the functional annotations from eggnog mapper and run_dbcan to a file of genes.

EGGNOG_KEGG_KO_INDEX = 8
EGGNOG_KEGG_PATHWAY_INDEX = 9
EGGNOG_KEGG_MODULE_INDEX = 10
EGGNOG_CAZY_INDEX = 15
DBCAN_DIAMOND_INDEX = 3

ANNOTATION_DELIMITER = ";"

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "Annotate Genome\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input-genome",
        help="the fasta file\n[REQUIRED]",
        metavar="<input.fasta>",
        required=True)
    parser.add_argument(
        "-e", "--input-eggnog",
        help="the eggnog annotations\n[REQUIRED]",
        metavar="<input.emapper.annotations>",
        required=True)
    parser.add_argument(
        "-d", "--input-dbcan",
        help="the dbcan annotations\n[REQUIRED]",
        metavar="<input.dbcan.overview.txt>",
        required=True)
    parser.add_argument(
        "-o", "--output",
        help="file to write the annotated genome\n[REQUIRED]",
        metavar="<output.fasta>",
        required=True)

    return parser.parse_args()

def format_cell(item, kegg_ko=False):
    # replace empty cells with none
    # remove ko: from each if kegg cell
    if item == "":
        item = "None"
    if kegg_ko:
        item = item.replace("ko:","")
    return item

def main():
    # read in run parameters
    args=parse_arguments(sys)

    # read in eggnog mapper annotations
    eggnog_kegg={}
    eggnog_cazy={}
    print("Reading eggnog annotations")
    with open(args.input_eggnog) as file_handle:
        for line in file_handle:
            if not line.startswith("#"):
                info = line.rstrip().split("\t")
                eggnog_kegg[info[0]]=(format_cell(info[EGGNOG_KEGG_KO_INDEX],kegg_ko=True),
                    format_cell(info[EGGNOG_KEGG_PATHWAY_INDEX]),
                    format_cell(info[EGGNOG_KEGG_MODULE_INDEX]))
                eggnog_cazy[info[0]]=format_cell(info[EGGNOG_CAZY_INDEX])

    # read in dbcan annotations
    run_dbcan={}
    print("Reading dbcan annotations")
    with open(args.input_dbcan) as file_handle:
        header = file_handle.readline()
        for line in file_handle:
            info = line.rstrip().split("\t")
            run_dbcan[info[0]]=format_cell(info[DBCAN_DIAMOND_INDEX])

    # write new fasta file with annotations
    print("Writing annotated fasta file")
    with open(args.output,"w") as file_handle_write:
        with open(args.input_genome) as file_handle_read:
            for line in file_handle_read:
                if line.startswith(">"):
                    info = line.rstrip().split(" ")
                    gene_id = info[0][1:]
                    name = "_".join(info[1:])
                    kegg_ko, kegg_pathway, kegg_module = eggnog_kegg.get(gene_id,["None","None","None"])
                    kegg = ANNOTATION_DELIMITER.join(["kegg_ko:"+kegg_ko, "kegg_pathway:"+kegg_pathway,"kegg_module:"+kegg_module])
                    cazy = "cazy:"+eggnog_cazy.get(gene_id,"None")
                    dbcan = "dbcan_diamond:"+run_dbcan.get(gene_id,"None")
                    file_handle_write.write(">"+";".join([gene_id,name,kegg,cazy,dbcan])+"\n")
                else:
                    file_handle_write.write(line)
            
if __name__ == "__main__":
    main()


