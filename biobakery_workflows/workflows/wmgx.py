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

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of biobakery_workflow tasks for shotgun sequences
from biobakery_workflows.tasks import shotgun

# import the utilities functions and config settings from biobakery_workflows
from biobakery_workflows import utilities, config

# create a workflow instance, providing the version number and description
# the version number will appear when running this script with the "--version" option
# the description will appear when running this script with the "--help" option
workflow = Workflow(version="0.1", description="A workflow for whole metagenome shotgun sequences")

# add the custom arguments to the workflow
workflow_config = config.ShotGun()
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz", choices=["fastq.gz","fastq","fq.gz","fq","fasta","fasta.gz"])
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair", default=".R1")
workflow.add_argument("contaminate-databases", desc="the path (or comma-delimited paths) to the contaminate\nreference databases for QC", 
    default=",".join([workflow_config.kneaddata_db_human_genome,workflow_config.kneaddata_db_rrna]))
workflow.add_argument("qc-options", desc="additional options when running the QC step", default="")
workflow.add_argument("remove-intermediate-output", desc="remove intermediate output files", action="store_true")
workflow.add_argument("bypass-strain-profiling", desc="do not run the strain profiling tasks", action="store_true")
workflow.add_argument("strain-profiling-options", desc="additional options when running the strain profiling step", default="")

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
# return an error if no files are found
input_files = utilities.find_files(args.input, extension=args.input_extension, exit_if_not_found=True)

### STEP #1: Run quality control on all input files ###
if not "fasta" in args.input_extension:
    qc_output_files, filtered_read_counts = shotgun.quality_control(workflow, 
        input_files, args.input_extension, args.output, args.threads, args.contaminate_databases, 
        args.pair_identifier, args.qc_options, args.remove_intermediate_output)
    # get the new extension, if the original files were gzipped they will not be after quality control
    args.input_extension = args.input_extension.replace(".gz","")
else:
    # if the input files are fasta, bypass quality control
    qc_output_files = input_files

### STEP #2: Run taxonomic profiling on all of the filtered files ###
merged_taxonomic_profile, taxonomy_tsv_files, taxonomy_sam_files = shotgun.taxonomic_profile(workflow,
    qc_output_files,args.output,args.threads,args.input_extension)

### STEP #3: Run functional profiling on all of the filtered files ###
genes_relab, ecs_relab, path_relab, genes, ecs, path = shotgun.functional_profile(workflow,
    qc_output_files,args.input_extension,args.output,args.threads,taxonomy_tsv_files,args.remove_intermediate_output)

### STEP #4: Run strain profiling
if not args.bypass_strain_profiling:
    shotgun.strain_profile(workflow,taxonomy_sam_files,args.output,args.threads,
        workflow_config.strainphlan_db_reference,workflow_config.strainphlan_db_markers,args.strain_profiling_options)

# start the workflow
workflow.go()
