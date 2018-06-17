#!/usr/bin/env python

"""
bioBakery Workflows: DADA2 workflow

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

import os

from anadama2 import Workflow

from biobakery_workflows.tasks import dadatwo
from biobakery_workflows import utilities, config



# create a workflow instance, providing the version number and description
workflow = Workflow(version="0.1", description="A workflow for DADA2 sequencing data")

# add the custom arguments to the workflow
workflow_config = config.DADA2()
workflow.add_argument("barcode-file", desc="the barcode file", default="")
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz", choices=["fastq.gz","fastq"])
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair", default="_R1_001")
workflow.add_argument("index-identifier", desc="the string to identify the index files", default="_I1_001")
workflow.add_argument("min-pred-qc-score", desc="the min phred quality score to use for demultiplexing", default=2)
workflow.add_argument("maxee", desc="the maxee value to use for quality control filtering", default=1)
workflow.add_argument("trunc-len-max", desc="the max length for truncating reads", default=200)
workflow.add_argument("min-size", desc="the min size to use for clustering", default=2)
workflow.add_argument("percent-identity", desc="the percent identity to use for alignments", default=0.97)
workflow.add_argument("pool", desc="if provided then samples are pooled for analysis, default setting is analysis of each sample independently", default="FALSE")

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
# return an error if no files are found
input_files = utilities.find_files(args.input, extension=args.input_extension, exit_if_not_found=True)

# check for index files, do not error if they are not found
index_files = utilities.find_files(args.input, extension=args.index_identifier+"."+args.input_extension)

# remove the index files, if found, from the set of input files
input_files = list(filter(lambda file: not file in index_files, input_files))

# if a barcode file is provided, then demultiplex
if args.barcode_file:
    demultiplexed_files=sixteen_s.demultiplex(
        workflow, input_files, args.input_extension, args.output, args.barcode_file, index_files,
        args.min_pred_qc_score, args.pair_identifier)
    # if the original files are gzipped, they will not be compressed after demultiplexing
    args.input_extension = args.input_extension.replace(".gz","")
else:
    demultiplexed_files=input_files


        
# call  workflow tasks
dadatwo.filter_trim(workflow,args.input, args.output, args.pool)
dadatwo.learn_error(workflow, args.output, args.pool)
dadatwo.merge_paired_ends(workflow, args.input, args.output, args.pool)
# start the workflow
workflow.go()
