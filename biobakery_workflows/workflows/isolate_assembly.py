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
from anadama2.tracked import TrackedExecutable

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
workflow.add_argument("dbcan-path", desc="the path to the run_dbcan.py script", default="/app/")

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
    sample_names=utilities.sample_names(input_pair1,args.input_extension,args.pair_identifier)
    qc_targets=[utilities.name_files([name+".trimmed.1.fastq",name+".trimmed.2.fastq",name+".trimmed.single.1.fastq",name+".trimmed.single.2.fastq",name+".trimmed.single.12.fastq"], args.output, subfolder="kneaddata", create_folder=True) for name in sample_names]
    paired = True
    for target_set,input_R1,input_R2,name in zip(qc_targets,input_pair1,input_pair2,sample_names):
        workflow.add_task(
            "kneaddata --run-fastqc-start --input [depends[0]] --input [depends[1]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]] && cat [args[3]] [args[4]] > [targets[2]]",
            depends=[input_R1, input_R2, TrackedExecutable("kneaddata")],
            targets=[target_set[0],target_set[1],target_set[4]],
            args=[os.path.dirname(target_set[0]),args.threads,name,target_set[2],target_set[3]])
else:
    qc_targets=utilities.name_files(sample_names, args.output, subfolder="kneaddata", create_folder=True, extension="trimmed.fastq")
    for target_file,input_file,name in zip(qc_targets,input_files,sample_names):
        workflow.add_task(
            "kneaddata --run-fastqc-start --input [depends[0]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]]",
            depends=[input_file, TrackedExecutable("kneaddata")],
            targets=[target_file],
            args=[os.path.dirname(target_file),args.threads,name])

### STEP #2: Run assembly ###
assembly_inputs = ""
assembly_depends=[]
assembly_targets = utilities.name_files("contigs.fasta", args.output, subfolder="spades", create_folder=True)
if paired:
    for i, inputs in enumerate(qc_targets):
        assembly_inputs += "--pe{0}-1 {1} --pe{0}-2  {2} --pe{0}-s {3} ".format(i+1, inputs[0], inputs[1], inputs[4])
        assembly_depends+=[inputs[0],inputs[1],inputs[4]]
    workflow.add_task(
        "spades.py "+assembly_inputs+" --careful --cov-cutoff auto -o [args[0]] --threads [args[1]]",
        depends=assembly_depends+[TrackedExecutable("spades.py")],
        targets=assembly_targets,
        args=[os.path.dirname(assembly_targets),args.threads])
else:
    for input_file in qc_targets:
        assembly_inputs += "-s {0} ".format(input_file)
    workflow.add_task(
        "spades.py "+assembly_inputs+" --careful --cov-cutoff auto -o [args[0]] --threads [args[1]]",
        depends=qc_targets+[TrackedExecutable("spades.py")],
        targets=assembly_targets,
        args=[os.path.dirname(assembly_targets),args.threads])

### STEP #3: Annotate assembly ###
annotation_targets = utilities.name_files([args.species_name +".faa", args.species_name +".ffn"], args.output, subfolder="prokka")
workflow.add_task(
    "prokka --outdir [args[0]] --prefix [args[1]] [depends[0]] --cpus [args[2]]",
    depends=[assembly_targets,TrackedExecutable("prokka")],
    targets=annotation_targets,
    args=[os.path.dirname(annotation_targets[0]),args.species_name,args.threads])

### STEP #4: Run quality assessment with quast ###
quast_targets = utilities.name_files("quast_stdout.txt", args.output, subfolder="quast", create_folder=True)
optional_database = ""
if args.reference_database:
    optional_database = " -r " + args.reference_database
workflow.add_task(
    "quast [depends[0]] --output-dir [args[0]] --threads [args[1]] "+optional_database+" > [targets[0]]",
    depends=[assembly_targets,TrackedExecutable("quast")],
    targets=quast_targets,
    args=[os.path.dirname(quast_targets), args.threads])

### STEP #5: Run quality assessment with checkm ###
checkm_targets = utilities.name_files("checkm_stdout.log", args.output, subfolder="checkm", create_folder=True)
workflow.add_task(
    "checkm lineage_wf [args[0]] [args[1]] > [targets[0]]",
    depends=[annotation_targets[0],TrackedExecutable("checkm",version_command="{}")],
    targets=checkm_targets,
    args=[os.path.dirname(annotation_targets[0]), os.path.dirname(checkm_targets)])

### STEP #6: Functional annotations with emapper ###
functional_targets = utilities.name_files(["eggnog_mapper.log", "eggnog_mapper.emapper.annotations"], args.output, subfolder="eggnog_mapper", create_folder=True)
workflow.add_task(
    "emapper.py -o [args[0]] --output_dir [args[1]] -i [depends[0]] -m diamond --cpu [args[2]] > [targets[0]]",
    depends=[annotation_targets[0],TrackedExecutable("emapper.py")],
    targets=functional_targets,
    args=[os.path.basename(functional_targets[0]).split(".")[0],os.path.dirname(functional_targets[0]),args.threads])

### STEP #7: Functional annotations with dbcan ###
dbcan_targets = utilities.name_files("overview.txt", args.output, subfolder="run_dbcan")
workflow.add_task(
    "cd [args[0]] && python run_dbcan.py [depends[0]] protein --out_dir [args[1]]",
    depends=annotation_targets[0],
    targets=dbcan_targets,
    args=[args.dbcan_path,os.path.dirname(dbcan_targets)],
    name="run_dbcan")

### STEP #8: Add functional annotations to genome files ###
final_genomes = utilities.name_files([args.species_name +".faa", args.species_name +".ffn"], args.output, subfolder="annotated_genomes", create_folder=True)
for annotation_file, genome_file in zip(annotation_targets, final_genomes):
    workflow.add_task(
        "annotate_genome.py --input-genome [depends[0]] --input-dbcan [depends[1]] --input-eggnog [depends[2]] --output [targets[0]]",
        depends=[annotation_file,dbcan_targets,functional_targets[1]],
        targets=genome_file)

# start the workflow
workflow.go()
