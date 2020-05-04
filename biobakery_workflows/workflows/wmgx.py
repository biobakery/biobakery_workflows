#!/usr/bin/env python

"""
bioBakery Workflows: whole metagenome shotgun workflow

Copyright (c) 2016 Harvard School of Public Health

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
import sys
import os, fnmatch

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of biobakery_workflow tasks for shotgun sequences
from biobakery_workflows.tasks import shotgun, general

# import the utilities functions and config settings from biobakery_workflows
from biobakery_workflows import utilities, config

# create a workflow instance, providing the version number and description
# the version number will appear when running this script with the "--version" option
# the description will appear when running this script with the "--help" option
workflow = Workflow(version="0.1", description="A workflow for whole metagenome shotgun sequences")

# add the custom arguments to the workflow
workflow_config = config.ShotGun()
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz", choices=["fastq.gz","fastq","fq.gz","fq","fasta","fasta.gz","fastq.bz2","fq.bz2"])
workflow.add_argument("barcode-file", desc="the barcode file", default="")
workflow.add_argument("dual-barcode-file", desc="the string to identify the dual barcode file", default="")
workflow.add_argument("index-identifier", desc="the string to identify the index files", default="_I1_001")
workflow.add_argument("min-pred-qc-score", desc="the min phred quality score to use for demultiplexing", default=2)
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair", default=".R1")
workflow.add_argument("interleaved", desc="indicates whether or not sequence files are interleaved", default=False, action="store_true")
workflow.add_argument("bypass-quality-control", desc="do not run the quality control tasks", action="store_true")
workflow.add_argument("contaminate-databases", desc="the path (or comma-delimited paths) to the contaminate\nreference databases for QC", 
    default=",".join([workflow_config.kneaddata_db_human_genome]))
workflow.add_argument("qc-options", desc="additional options when running the QC step", default="")
workflow.add_argument("functional-profiling-options", desc="additional options when running the functional profiling step", default="")
workflow.add_argument("remove-intermediate-output", desc="remove intermediate output files", action="store_true")
workflow.add_argument("bypass-functional-profiling", desc="do not run the functional profiling tasks", action="store_true")
workflow.add_argument("bypass-strain-profiling", desc="do not run the strain profiling tasks (StrainPhlAn)", action="store_true")
workflow.add_argument("run-strain-gene-profiling", desc="run the gene-based strain profiling tasks (PanPhlAn)", action="store_true")
workflow.add_argument("bypass-taxonomic-profiling", desc="do not run the taxonomic profiling tasks (a tsv profile for each sequence file must be included in the input folder using the same sample name)", action="store_true")
workflow.add_argument("run-assembly", desc="run the assembly and annotation tasks", action="store_true")
workflow.add_argument("strain-profiling-options", desc="additional options when running the strain profiling step", default="")
workflow.add_argument("max-strains", desc="the max number of strains to profile", default=20, type=int)
workflow.add_argument("strain-list", desc="input file with list of strains to profile", default="")
workflow.add_argument("assembly-options", desc="additional options when running the assembly step", default="")

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
# return an error if no files are found
input_files = utilities.find_files(args.input, extension=args.input_extension, exit_if_not_found=True)

# check for index files, do not error if they are not found
index_files = utilities.find_files(args.input, extension=args.index_identifier+"."+args.input_extension)

# remove the index files, if found, from the set of input files
input_files = list(filter(lambda file: not file in index_files, input_files))

# if a dual index file is provided, then demultiplex dual indexing
if args.dual_barcode_file:
    if ".bz2" in args.input_extension:
        sys.exit("ERROR: Bz2 formatted files are not supported with demultiplexing")
    barcode_files = fnmatch.filter(os.listdir(args.input), '*barcode*.fastq*')
    barcode_files = [os.path.join(args.input,file) for file in barcode_files]
    input_files = list(filter(lambda file: not file in barcode_files, input_files))

    demultiplexed_files, demultiplex_output_folder = general.demultiplex_dual(workflow,args.output, input_files,
             args.input_extension, barcode_files, args.dual_barcode_file, args.min_pred_qc_score, args.pair_identifier)
    
    # if the original files are gzipped, they will not be compressed after demultiplexing
    args.input_extension = args.input_extension.replace(".gz", "")

# if a barcode file is provided, then demultiplex
elif args.barcode_file:
    if ".bz2" in args.input_extension:
        sys.exit("ERROR: Bz2 formatted files are not supported with demultiplexing")
    demultiplexed_files, demultiplex_output_folder=general.demultiplex(
            workflow, input_files, args.input_extension, args.output, args.barcode_file, index_files,
            args.min_pred_qc_score, args.pair_identifier)
    # if the original files are gzipped, they will not be compressed after demultiplexing
    args.input_extension = args.input_extension.replace(".gz","")
else:
    demultiplexed_files=input_files
    demultiplex_output_folder=args.input


### STEP #1: Run quality control on all input files ###
original_extension = args.input_extension
if args.bypass_quality_control:
    # merge files if they are paired
    qc_output_files, args.input_extension = shotgun.merge_pairs(workflow,
        demultiplexed_files, args.input_extension, args.pair_identifier, args.output)
    
elif not "fasta" in args.input_extension:
    qc_output_files, filtered_read_counts = shotgun.quality_control(workflow,
        demultiplexed_files, args.input_extension, args.output, args.threads, args.contaminate_databases,
        args.pair_identifier, args.qc_options, args.remove_intermediate_output)
    # get the new extension, if the original files were gzipped they will not be after quality control
    args.input_extension = args.input_extension.replace(".gz","")
    args.input_extension = args.input_extension.replace(".bz2","")
else:
    # if the input files are fasta, bypass quality control
    qc_output_files = demultiplexed_files

### STEP #2: Run taxonomic profiling on all of the filtered files ###
if not args.bypass_taxonomic_profiling:
    merged_taxonomic_profile, taxonomy_tsv_files, taxonomy_sam_files = shotgun.taxonomic_profile(workflow,
        qc_output_files,args.output,args.threads,args.input_extension)

elif not args.bypass_functional_profiling or not args.bypass_strain_profiling:
    # get the names of the taxonomic profiling files allowing for pairs
    input_pair1, input_pair2 = utilities.paired_files(demultiplexed_files, original_extension, args.pair_identifier)
    sample_names = utilities.sample_names(input_pair1 if input_pair1 else input_files,original_extension,args.pair_identifier)
    tsv_profiles = utilities.name_files(sample_names, demultiplex_output_folder, tag="taxonomic_profile", extension="tsv")
    # check all of the expected profiles are found
    if len(tsv_profiles) != len(list(filter(os.path.isfile,tsv_profiles))):
        sys.exit("ERROR: Bypassing taxonomic profiling but all of the tsv taxonomy profile files are not found in the input folder. Expecting the following input files:\n"+"\n".join(tsv_profiles))
    # run taxonomic profile steps bypassing metaphlan
    merged_taxonomic_profile, taxonomy_tsv_files, taxonomy_sam_files = shotgun.taxonomic_profile(workflow,
        tsv_profiles,args.output,args.threads,"tsv",already_profiled=True)
    # look for the sam profiles
    taxonomy_sam_files = utilities.name_files(sample_names, demultiplex_output_folder, tag="bowtie2", extension="sam")
    # if they do not all exist, then bypass strain profiling if not already set
    if not args.bypass_strain_profiling:
        if len(taxonomy_sam_files) != len(list(filter(os.path.isfile,taxonomy_sam_files))):
            print("Warning: Bypassing taxonomic profiling but not all taxonomy sam files are present in the input folder. Strain profiling will be bypassed. Expecting the following input files:\n"+"\n".join(taxonomy_sam_files))
            args.bypass_strain_profiling = True

### STEP #3: Run functional profiling on all of the filtered files ###
if not args.bypass_functional_profiling:
    genes_relab, ecs_relab, path_relab, genes, ecs, path = shotgun.functional_profile(workflow,
        qc_output_files,args.input_extension,args.output,args.threads,taxonomy_tsv_files,args.remove_intermediate_output,
        args.functional_profiling_options)

### STEP #4: Run strain profiling
# Provide taxonomic profiling output so top strains by abundance will be selected
if not args.bypass_strain_profiling:
    if args.bypass_taxonomic_profiling:
        sys.exit("ERROR: Taxonomic profiling must be run to also run strain profiling")
    shotgun.strain_profile(workflow,taxonomy_sam_files,args.output,args.threads,
        workflow_config.strainphlan_db_reference,workflow_config.strainphlan_db_markers,merged_taxonomic_profile,
        args.strain_profiling_options,args.max_strains,args.strain_list)

### STEP #5: Run gene-based strain profiling (optional)
if args.run_strain_gene_profiling:
    if args.bypass_taxonomic_profiling and not args.strain_list:
        sys.exit("ERROR: Taxonomic profiling must be run to also run gene-based strain profiling")
    if args.strain_list:
        merged_taxonomic_profile=args.strain_list
    shotgun.strain_gene_profile(workflow,qc_output_files,merged_taxonomic_profile,args.output,args.threads,workflow_config.panphlan_db,
    args.max_strains)

### STEP #6: Run assembly and annotation
if args.run_assembly:
    is_paired = args.interleaved or utilities.is_paired_end(input_files, original_extension, args.pair_identifier)
    assembled_contigs = shotgun.assemble(workflow, qc_output_files, args.input_extension, args.output, args.threads, args.pair_identifier, args.remove_intermediate_output, args.assembly_options, is_paired)
    shotgun.annotate(workflow, assembled_contigs, args.output, args.threads)
 
# start the workflow
workflow.go()
