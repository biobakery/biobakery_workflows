#!/usr/bin/env python

"""
bioBakery Workflows: 16S workflow

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""


from anadama2 import Workflow
import os, sys, fnmatch

from biobakery_workflows.tasks import sixteen_s, dadatwo, general
from biobakery_workflows import utilities, config, files


# create a workflow instance, providing the version number and description
workflow = Workflow(version="3.1", description="A workflow for 16S sequencing data")

# add the custom arguments to the workflow
workflow_config = config.SixteenS()
workflow.add_argument("method", desc="method to process 16s workflow", default="vsearch", choices=["usearch","dada2","vsearch","its"])
workflow.add_argument("dada-db", desc="reference database for dada2 workflow", default="silva", choices=["gg","rdp","silva","unite"])
workflow.add_argument("usearch-db", desc="full paths for the reference databases (fna and taxonomy, comma delimited) for the usearch workflow",
    default=",".join([workflow_config.greengenes_fasta,workflow_config.greengenes_taxonomy]))
workflow.add_argument("bypass-functional-profiling", desc="bypass the functional profiling tasks", action="store_true")
workflow.add_argument("barcode-file", desc="the barcode file", default="")
workflow.add_argument("dual-barcode-file", desc="the string to identify the dual barcode file", default="")
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz", choices=["fastq.gz","fastq"])
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair, must proceed the file extension (ie R1_001.fastq.gz)", default="_R1_001")
workflow.add_argument("index-identifier", desc="the string to identify the index files", default="_I1_001")
workflow.add_argument("min-pred-qc-score", desc="the min phred quality score to use for demultiplexing", default=2)
workflow.add_argument("maxee", desc="the maxee value to use for quality control filtering", default=1)
workflow.add_argument("trunc-len-max", desc="the max length for truncating reads", default=200)
workflow.add_argument("trunc-len-rev-offset", desc="the offset (use negative number to subtract) to max length for truncating reads on the reverse strand (for dada2 only)", default=40)
workflow.add_argument("min-fold-parent-over-abundance", desc="the min fold difference between child and parent to call a sequence as chimeric (for dada2 only)", default=1)
workflow.add_argument("min-size", desc="the min size to use for clustering", default=2)
workflow.add_argument("bypass-primers-removal", desc="do not run remove primers tasks", action="store_true")
workflow.add_argument("fwd-primer", desc="forward primer, required for its workflow",default="")
workflow.add_argument("rev-primer", desc="reverse primer, required for its workflow",default="")
workflow.add_argument("amplicon-length", desc="length of the amplicon, required for figaro",default="")
workflow.add_argument("cutadapt-options", desc="additional options when running cutadapt",default="")
workflow.add_argument("minoverlap", desc="the min overlap required to merge pairs for the dada2 workflow", default=20)
workflow.add_argument("maxmismatch", desc="the max mismatch required to merge pairs for the dada2 workflow", default=0)
workflow.add_argument("tryRC", desc="try the reverse complement of the reads for the dada2 workflow", default="FALSE")
workflow.add_argument("percent-identity", desc="the percent identity to use for alignments", default=0.97)
workflow.add_argument("bypass-msa", desc="bypass running multiple sequence alignment and tree generation", action="store_true")
workflow.add_argument("picrust-version", desc="the picrust version to use", default="2")
workflow.add_argument("fastq-ascii", desc="the coding of Q scores", default="33", choices=["33","64"])
workflow.add_argument("min-len", desc="remove reads less then this length for DADA2", default="50")
workflow.add_argument("pooling", desc="how to pool samples for DADA2", default="FALSE", choices=["FALSE","TRUE","pseudo"])

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
# return an error if no files are found
input_files = utilities.find_files(args.input, extension=args.input_extension, exit_if_not_found=True)

# check for index files, do not error if they are not found
index_files = utilities.find_files(args.input, extension=args.index_identifier+"."+args.input_extension)

# remove the index files, if found, from the set of input files
input_files = list(filter(lambda file: not file in index_files, input_files))

# check for valid number for percent identity
if float(args.percent_identity) > 1:
    sys.exit("ERROR: Please set a percent identity that is less than 1")

# check for empty files in the input folder
utilities.check_for_empty_files(args.input, args.input_extension)

# if a dual index file is provided, then demultiplex dual indexing
if args.dual_barcode_file:
    barcode_files = fnmatch.filter(os.listdir(args.input), '*barcode*.fastq*')
    barcode_files = [os.path.join(args.input,file) for file in barcode_files]
    input_files = list(filter(lambda file: not file in barcode_files, input_files))

    demultiplexed_files, demultiplex_output_folder = general.demultiplex_dual(workflow,args.output, input_files,
             args.input_extension, barcode_files, args.dual_barcode_file, args.min_pred_qc_score, args.pair_identifier)

    # if the original files are gzipped, they will not be compressed after demultiplexing
    args.input_extension = args.input_extension.replace(".gz", "")

# if a barcode file is provided, then demultiplex
elif args.barcode_file:
    demultiplexed_files, demultiplex_output_folder=general.demultiplex(
            workflow, input_files, args.input_extension, args.output, args.barcode_file, index_files,
            args.min_pred_qc_score, args.pair_identifier)
    # if the original files are gzipped, they will not be compressed after demultiplexing
    args.input_extension = args.input_extension.replace(".gz","")
else:
    demultiplexed_files=input_files
    demultiplex_output_folder=args.input

    # check the max trunc len is not larger then the read length
    input_read_length_average = utilities.get_average_read_length_fastq(demultiplexed_files[0])
    if input_read_length_average < int(args.trunc_len_max):
        print("WARNING: The average input file read length ( {0} ) is less then the max trunc length provided ( {1} ). Please modify the max trunc length using the option '--trunc-len-max <200>'.".format(int(input_read_length_average), args.trunc_len_max))

if args.method == "dada2" or args.method == "its":

    # if its workflow remove primers first and set reference db to 'unite'
    primer_tasks=[]
    if args.method == "its":
        args.dada_db = "unite"
        args.trunc_len_max = 0

        if not args.bypass_primers_removal:
            if args.fwd_primer and args.rev_primer:
                cutadapt_folder=dadatwo.remove_primers(
                    workflow,args.fwd_primer,args.rev_primer,demultiplex_output_folder,args.output,args.pair_identifier,args.threads)
                demultiplex_output_folder=cutadapt_folder
            else:
                print("ITS workflow primers rmoval task requires fwd_primer and rev_primer arguments.")
                exit()
    elif args.fwd_primer and args.rev_primer:
        # check for pairs
        pair1, pair2=utilities.paired_files(demultiplex_output_folder, args.input_extension, args.pair_identifier)
        if pair1 and pair2:
            cutadapt_folder=dadatwo.remove_primers(
                workflow,args.fwd_primer,args.rev_primer,demultiplex_output_folder,args.output,args.pair_identifier,args.threads)
            demultiplex_output_folder=cutadapt_folder
        else:
            primer_tasks,cutadapt_files=general.remove_primers(
                workflow,args.fwd_primer,args.rev_primer,demultiplex_output_folder,args.output,args.pair_identifier,args.threads,args.input_extension,demultiplexed_files,args.cutadapt_options)
            args.input_extension=args.input_extension.replace(".gz","")
            demultiplex_output_folder=os.path.dirname(cutadapt_files[0])

    # run figaro to determine trunc length if needed
    figaro_csv=""
    if args.amplicon_length:
        figaro_csv = dadatwo.figaro(workflow,args.fwd_primer,args.rev_primer,args.amplicon_length,demultiplex_output_folder,args.output)

    # call dada2 workflow tasks
    # filter reads and trim
    read_counts_rds_path,  filtered_dir = dadatwo.filter_trim(
            workflow, demultiplex_output_folder,
            args.output,args.maxee,args.trunc_len_max,args.pair_identifier,args.threads,args.trunc_len_rev_offset,args.min_len,figaro_csv,primer_tasks)
    
    # learn error rates
    error_ratesF_path, error_ratesR_path = dadatwo.learn_error(
            workflow, args.output, filtered_dir, read_counts_rds_path, args.threads)
    
    # merge pairs
    mergers_file_path = dadatwo.merge_paired_ends(
            workflow, args.output, filtered_dir, error_ratesF_path, error_ratesR_path, args.threads, args.minoverlap, args.maxmismatch, args.pooling)

    # construct otu
    seqtab_file_path,read_counts_steps_path, seqs_fasta_path = dadatwo.const_seq_table(
            workflow, args.output, filtered_dir, mergers_file_path, args.threads, args.min_fold_parent_over_abundance)

    # centroid alignment
    centroid_fasta = files.SixteenS.path("msa_nonchimera", args.output)
    sixteen_s.centroid_alignment(workflow,
            seqs_fasta_path, centroid_fasta, args.threads, task_name="clustalo_nonchimera")

    # phylogenetic tree
    closed_tree = utilities.name_files("closed_reference.tre", args.output)
    sixteen_s.create_tree(workflow, centroid_fasta, closed_tree)

    # assign taxonomy
    closed_reference_tsv = dadatwo.assign_taxonomy(
            workflow, args.output, seqtab_file_path, args.dada_db, args.threads, args.tryRC)
    
    # functional profiling
    # check for picrust1 as not an option with this workflow
    if args.picrust_version == "1":
        print("WARNING: PICRUSt v1 is not compatible with ASV tables so will not be run for this workflow.")
    else:
        categorized_function = sixteen_s.functional_profile(workflow, closed_reference_tsv, seqs_fasta_path, 
                args.picrust_version, args.threads, args.output, otus=False, method=args.method)

else:
    cutadapt_files=demultiplexed_files
    if not args.bypass_primers_removal:
        if args.fwd_primer:
            primer_tasks,cutadapt_files=general.remove_primers(
                workflow,args.fwd_primer,args.rev_primer,demultiplex_output_folder,args.output,args.pair_identifier,args.threads,args.input_extension,demultiplexed_files,args.cutadapt_options)
            args.input_extension=args.input_extension.replace(".gz","")

    # call vsearch or usearch workflow tasks
    #  merge pairs, if paired-end, then rename so sequence id matches sample name then merge to single fastq file
    all_samples_fastq = sixteen_s.merge_samples_and_rename(
    	       workflow, args.method, cutadapt_files, args.input_extension, args.output, args.pair_identifier, args.threads, args.fastq_ascii)

	# add quality control tasks: generate qc report, filter by maxee, and truncate
    filtered_truncated_fasta, truncated_fasta, original_fasta = sixteen_s.quality_control(
            workflow, args.method, all_samples_fastq, args.output, args.threads, args.maxee, args.trunc_len_max, args.fastq_ascii)

    # taxonomic profiling (pick otus and then align creating otu tables, closed and open reference)
    try:
        usearch_fna_db, usearch_taxon_tsv = args.usearch_db.strip().split(",")
    except ValueError:
        sys.exit("ERROR: Please provide two comma seperated strings for the custom usearch databases")

    closed_reference_tsv, closed_ref_fasta = sixteen_s.taxonomic_profile(
            workflow, args.method, filtered_truncated_fasta, truncated_fasta, original_fasta, args.output,
            args.threads, args.percent_identity, usearch_fna_db, usearch_fna_db,
            usearch_taxon_tsv, args.min_size, args.bypass_msa)

    # functional profiling
    if not args.bypass_functional_profiling:
        categorized_function = sixteen_s.functional_profile(workflow, closed_reference_tsv, closed_ref_fasta, 
                args.picrust_version, args.threads, args.output, otus=True)

# start the workflow
workflow.go()
