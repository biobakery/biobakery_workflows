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

import os

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of biobakery_workflow tasks for shotgun sequences
from biobakery_workflows.tasks import shotgun

# import the utilities functions from biobakery_workflows
from biobakery_workflows import utilities

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with two input folder options
workflow = Workflow(version="0.1", remove_options=["input"],
                    description="A workflow for whole metagenome and metatranscriptome shotgun sequences")

# add the custom arguments to the workflow
workflow.add_argument("input-metagenome",desc="the input folder of whole metagenome shotgun sequences", required=True)
workflow.add_argument("input-metatranscriptome",desc="the input folder of whole metatranscriptome shotgun sequences", required=True)
workflow.add_argument("input-mapping",desc="the mapping file of metatranscriptome samples to metagenome samples")
workflow.add_argument("kneaddata-db", desc="the kneaddata database", required=True)
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair", default=".R1.")

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
input_files_metagenome = utilities.find_files(args.input_metagenome, extension=args.input_extension, exit_if_not_found=True)
input_files_metatranscriptome = utilities.find_files(args.input_metatranscriptome, extension=args.input_extension, exit_if_not_found=True)

### STEP #1: Run quality control on all input files ###
wms_output_folder = os.path.join(args.output,"whole_metagenome_shotgun")
wts_output_folder = os.path.join(args.output,"whole_metatranscriptome_shotgun")
wms_qc_output_files, wms_filtered_read_count = shotgun.quality_control(workflow, input_files_metagenome, wms_output_folder, args.threads, args.kneaddata_db, args.pair_identifier)
wts_qc_output_files, wts_filtered_read_count = shotgun.quality_control(workflow, input_files_metatranscriptome, wts_output_folder, args.threads, args.kneaddata_db, args.pair_identifier)

### STEP #2: Run taxonomic profiling on all of the metagenome filtered files (and metatranscriptome if mapping not provided)###
wms_taxonomic_profile, wms_taxonomy_tsv_files, wms_taxonomy_sam_files = shotgun.taxonomic_profile(workflow,wms_qc_output_files,wms_output_folder,args.threads)

if not args.input_mapping:
    wts_taxonomic_profile, wts_taxonomy_tsv_files, wts_taxonomy_sam_files = shotgun.taxonomic_profile(workflow,wts_qc_output_files,wts_output_folder,args.threads)   

### STEP #3: Run functional profiling on all of the filtered files ###

# run the wms samples with their taxonomy files
wms_genefamilies, wms_ecs, wms_pathabundance = shotgun.functional_profile(workflow,wms_qc_output_files,wms_output_folder,args.threads,wms_taxonomy_tsv_files)

# provide the taxonomy files for the wms samples to the wts samples, if mapping file provided
if args.input_mapping:
    # get the mapped taxonomy files for the quality controlled input files
    filtered_fastq, matched_taxonomic_profile = utilities.match_files(wts_qc_output_files,wms_taxonomy_tsv_files,args.input_mapping)
    wts_genefamilies, wts_ecs, wts_pathabundance = shotgun.functional_profile(workflow,filtered_fastq,wts_output_folder,args.threads,matched_taxonomic_profile)
else:
    # if no mapping file is provided then run to get a taxonomic profile from the wts samples
    wts_genefamilies, wts_ecs, wts_pathabundance = shotgun.functional_profile(workflow,wts_qc_output_files,wts_output_folder,args.threads,wts_taxonomy_tsv_files)

# start the workflow
workflow.go()
