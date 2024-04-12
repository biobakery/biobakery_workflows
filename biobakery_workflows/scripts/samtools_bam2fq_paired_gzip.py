#!/usr/bin/env python

from anadama2 import Workflow

workflow = Workflow()

input_files = workflow.get_input_files(extension="bam")
# rename input files to replace period with underscore (so not to look like extension)
input_files_renamed=[file.replace(".","_").replace("_bam",".fastq") for file in input_files]

output_files = workflow.name_output_files(name=input_files_renamed, tag="R1", extension="fastq")

for output, input in zip(output_files,input_files):
    all_file=output.replace("R1","ALL")
    output_r2=output.replace("R1","R2")
    workflow.add_task_gridable(
        "samtools bam2fq [depends[0]] > [args[0]] && "+\
        "grep '^@.*/1$' [args[0]] -A 3 --no-group-separator > [args[1]] && "+\
        "grep '^@.*/2$' [args[0]] -A 3 --no-group-separator > [args[2]] && "+\
        "gzip -f [args[1]] && "+\
        "gzip -f [args[2]] && "+\
        "rm [args[0]] ",
        depends=input,
        args=[all_file,output,output_r2],
        targets=[output+".gz",output_r2+".gz"],
        time=60*3, # 3 hours
        mem=4*1024, # 4 GB
        cores=1)

workflow.go()
