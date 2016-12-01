#!/usr/bin/env python

import os
from anadama2 import Workflow

# create a workflow instance, providing the version number and description
# the version number will appear when running this script with the "--version" option
# the description will appear when running this script with the "--help" option
workflow = Workflow(version="0.1", description="A workflow for whole metagenome shotgun sequences")

# add the custom arguments to the workflow
workflow.add_argument("kneaddata-db", desc="the kneaddata database", required=True)
workflow.add_argument("input-extension", desc="the input file extension", default="fastq")
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
input_files = workflow.get_input_files(extension=args.input_extension)

""" STEP #1: Run Kneaddata on all input files """

# get a list of output files, one for each input file, with the kneaddata tag
kneaddata_output_files = workflow.name_output_files(name=input_files, tag="kneaddata", subfolder="kneaddata")
kneaddata_output_folder = os.path.dirname(kneaddata_output_files[0])

# create a task for each set of input and output files to run kneaddata
workflow.add_task_group(
    "kneaddata --input [depends[0]] --output [args[0]] --reference-db [args[1]] --threads [args[2]]",
    depends=input_files,
    targets=kneaddata_output_files,
    args=[kneaddata_output_folder, args.kneaddata_db, args.threads])

""" STEP #2: Run Metaphlan2 on all of the filtered fastq files (and merge tables)"""

# get a list of metaphlan2 output files, one for each input file
metaphlan2_profile_tag="taxonomic_profile"
metaphlan2_output_files_profile = workflow.name_output_files(name=input_files, subfolder="metaphlan2", tag=metaphlan2_profile_tag, extension="tsv")
metaphlan2_output_files_bowtie2 = workflow.name_output_files(name=input_files, subfolder="metaphlan2", tag="bowtie2", extension="tsv")
metaphlan2_output_files_sam = workflow.name_output_files(name=input_files, subfolder="metaphlan2", tag="bowtie2", extension="sam")
metaphlan2_output_folder = os.path.dirname(metaphlan2_output_files_profile[0])

# run metaphlan2 on each of the kneaddata output files
workflow.add_task_group(
    "metaphlan2.py [depends[0]] --input_type fastq --output_file [targets[0]] --bowtie2out [targets[1]] --samout [targets[2]] --nproc [args[0]]",
    depends=kneaddata_output_files,
    targets=zip(metaphlan2_output_files_profile, metaphlan2_output_files_bowtie2, metaphlan2_output_files_sam),
    args=[args.threads]) 

# merge all of the metaphlan taxonomy tables
metaphlan2_merged_output = workflow.name_output_files(name="taxonomic_profiles.tsv")

# run the humann2 join script to merge all of the metaphlan2 profiles
workflow.add_task(
    "humann2_join_tables --input [args[0]] --output [targets[0]] --file_name [args[1]]",
    depends=metaphlan2_output_files_profile,
    targets=metaphlan2_merged_output,
    args=[metaphlan2_output_folder, metaphlan2_profile_tag])

""" STEP #3: Run HUMAnN2 on the filtered fastq files (providing the MetaPhlAn2 taxonomic profiles) """

# get a list of output files, one for each input file, with the humann2 output file names
genefamiles = workflow.name_output_files(name=kneaddata_output_files, subfolder="humann2", tag="genefamilies", extension="tsv")
pathabundance = workflow.name_output_files(name=kneaddata_output_files, subfolder="humann2", tag="pathabundance", extension="tsv")
pathcoverage = workflow.name_output_files(name=kneaddata_output_files, subfolder="humann2", tag="pathcoverage", extension="tsv")
humann2_output_folder = os.path.dirname(genefamiles[0])

# create a task to run humann2 on each of the kneaddata output files
workflow.add_task_group(
    "humann2 --input [depends[0]] --output [args[0]] --taxonomic-profile [depends[1]] --threads [args[1]]",
    depends=zip(kneaddata_output_files,metaphlan2_output_files_profile),
    targets=zip(genefamiles, pathabundance, pathcoverage),
    args=[humann2_output_folder, args.threads])

""" STEP #4: Regroup UniRef90 gene families to ecs """

# get a list of all output ec files
ec_files = workflow.name_output_files(name=genefamiles, subfolder="humann2", tag="ecs")

# get ec files for all of the gene families files
workflow.add_task_group(
    "humann2_regroup_table --input [depends[0]] --output [targets[0]] --groups uniref90_level4ec",
    depends=genefamiles,
    targets=ec_files)

""" STEP #5: Normalize gene families, ecs, and pathway abundance to relative abundance (then merge files) """

# get a list of files for normalized ec, gene families, and pathway abundance
norm_genefamily_files = workflow.name_output_files(name=genefamiles, subfolder="genes", tag="relab")
norm_ec_files = workflow.name_output_files(name=ec_files, subfolder="ecs", tag="relab")
norm_pathabundance_files = workflow.name_output_files(name=pathabundance, subfolder="pathways", tag="relab")

# normalize the genefamily, ec, and pathabundance files
workflow.add_task_group(
    "humann2_renorm_table --input [depends[0]] --output [targets[0]] --units relab",
    depends=genefamiles + ec_files + pathabundance,
    targets=norm_genefamily_files + norm_ec_files + norm_pathabundance_files)

# get a list of merged files for ec, gene families, and pathway abundance
merged_genefamilies = workflow.name_output_files(name="genefamilies_relab.tsv")
merged_ecs = workflow.name_output_files(name="ecs_relab.tsv")
merged_pathabundance = workflow.name_output_files(name="pathabundance_relab.tsv")

# merge the ec, gene families, and pathway abundance files
all_depends=[norm_genefamily_files, norm_ec_files, norm_pathabundance_files]
all_targets=[merged_genefamilies, merged_ecs, merged_pathabundance]
for depends, targets in zip(all_depends, all_targets):
    workflow.add_task(
        "humann2_join_tables --input [args[0]] --output [targets[0]]",
        depends=depends,
        targets=targets,
        args=[os.path.dirname(depends[0])])

# start the workflow
workflow.go()
