#!/usr/bin/env python

"""
bioBakery Workflows: isolate assembly

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
workflow = Workflow(version="0.1", description="A workflow for isolate assembly")

# add the custom arguments to the workflow
workflow_config = config.ShotGun()
workflow.add_argument("species-name", desc="the species name", required=True)
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz", choices=["fastq.gz","fastq","fq.gz","fq"])
workflow.add_argument("threads", desc="number of threads/cores for each task to use", default=1)
workflow.add_argument("pair-identifier", desc="the string to identify the first file in a pair", default="_R1_001")
workflow.add_argument("reference-database", desc="the path to the reference database for quality assessment", default="")

# get the arguments from the command line
args = workflow.parse_args()

# get all input files with the input extension provided on the command line
# return an error if no files are found
input_files = utilities.find_files(args.input, extension=args.input_extension, exit_if_not_found=True)

### STEP #1: Run quality control on all input files ###
sample_names=utilities.sample_names(input_files,args.input_extension)
input_pair1, input_pair2 = utilities.paired_files(input_files, args.input_extension, args.pair_identifier)
paired = False
if input_pair1:
    qc_targets=utilities.name_files([sample_names[0]+".trimmed.1.fastq",sample_names[0]+".trimmed.2.fastq",sample_names[0]+".trimmed.single.1.fastq",sample_names[0]+".trimmed.single.2.fastq",sample_names[0]+".trimmed.single.12.fastq"], args.output, subfolder="kneaddata", create_folder=True)
    paired = True
    workflow.add_task(
        "kneaddata --run-fastqc-start --input [depends[0]] --input [depends[1]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]] && cat [args[3]] [args[4]] > [targets[2]]",
        depends=[input_pair1[0], input_pair2[0]],
        targets=[qc_targets[0],qc_targets[1],qc_targets[4]],
        args=[os.path.dirname(qc_targets[0]),args.threads,sample_names[0],qc_targets[2],qc_targets[3]])
else:
    qc_targets=[os.path.join(qc_folder, sample_names[0]+".fastq")]
    workflow.add_task(
        "kneaddata --run-fastqc-start --input [depends[0]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]]",
        depends=[input_pair1[0]],
        targets=[qc_targets[0]],
        args=[os.path.dirname(qc_targets[0]),args.threads,sample_names[0]])

### STEP #2: Run assembly ###
assembly_targets = utilities.name_files("contigs.fasta", args.output, subfolder="spades", create_folder=True)
if paired:
    workflow.add_task(
        "spades.py --pe1-1 [depends[0]] --pe1-2  [depends[1]] --pe1-s [depends[2]] --careful --cov-cutoff auto -o [args[0]] --threads [args[1]]",
        depends=[qc_targets[0],qc_targets[1],qc_targets[4]],
        targets=assembly_targets,
        args=[os.path.dirname(assembly_targets),args.threads])
else:
    workflow.add_task(
        "spades.py -s [depends[0]] --careful --cov-cutoff auto -o [args[0]] --threads [args[1]]",
        depends=qc_targets[0],
        targets=assembly_targets,
        args=[os.path.dirname(assembly_targets),args.threads])

### STEP #3: Annotate assembly ###
annotation_targets = utilities.name_files(args.species_name +".faa", args.output, subfolder="prokka")
workflow.add_task(
    "prokka --outdir [args[0]] --prefix [args[1]] [depends[0]] --cpus [args[2]]",
    depends=assembly_targets,
    targets=annotation_targets,
    args=[os.path.dirname(annotation_targets),args.species_name,args.threads])

### STEP #4: Run quality assessment with quast ###
quast_targets = utilities.name_files("quast_stdout.txt", args.output, subfolder="quast", create_folder=True)
optional_database = ""
if args.reference_database:
    optional_database = " -r " + args.reference_database
workflow.add_task(
    "quast [depends[0]] --output-dir [args[0]] --threads [args[1]] "+optional_database+" > [targets[0]]",
    depends=assembly_targets,
    targets=quast_targets,
    args=[os.path.dirname(quast_targets), args.threads])

### STEP #5: Run quality assessment with checkm ###
checkm_targets = utilities.name_files("checkm_stdout.log", args.output, subfolder="checkm", create_folder=True)
workflow.add_task(
    "checkm lineage_wf [args[0]] [args[1]] > [targets[0]]",
    depends=annotation_targets,
    targets=checkm_targets,
    args=[os.path.dirname(annotation_targets), os.path.dirname(checkm_targets)])

### STEP #6: Functional annotations ###
functional_targets = utilities.name_files(["eggnog_mapper.log", "eggnog_mapper.emapper.annotations"], args.output, subfolder="eggnog_mapper", create_folder=True)
workflow.add_task(
    "emapper.py -o [args[0]] --output_dir [args[1]] -i [depends[0]] -m diamond --cpu [args[2]] > [targets[0]]",
    depends=annotation_targets,
    targets=functional_targets,
    args=[os.path.basename(functional_targets[0]).split(".")[0],os.path.dirname(functional_targets[0]),args.threads])

# start the workflow
workflow.go()
