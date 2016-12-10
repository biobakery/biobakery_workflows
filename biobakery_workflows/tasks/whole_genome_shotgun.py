"""
bioBakery Workflows: tasks.whole_genome_shotgun module
A collection of tasks for wgs workflows

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

from biobakery_workflows import utilities

def kneaddata(workflow, input_files, threads, paired=None, databases=None):
    """Run kneaddata
    
    This set of tasks will run kneaddata on the input files provided. It will run with
    single-end or paired-end input files.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to kneaddata. This
        is a list of lists if the input files are paired.
        threads (int): The number of threads/cores for kneaddata to use.
        paired (bool): This indicates if the input files are paired.
        databases (string/list): The databases to use with kneaddata (optional).
        
    Requires:
        kneaddata v0.5.3+: A tool to perform quality control on metagenomic and
            metatranscriptomic sequencing data
        
    Returns:
        list: A list of the filtered fastq files created by kneaddata.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import whole_genome_shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = whole_genome_shotgun.kneaddata(workflow,
            ["demo.fastq","demo2.fastq"], 1)
            
        # run the workflow
        workflow.go()
    """

    # get the sample basenames from the input files
    if paired:
        sample_names=utilities.sample_names(input_files[0])
    else:
        sample_names=utilities.sample_names(input_files)
        

    if paired:
        # reorder the input files so they are a set of paired files
        input_files=zip(input_files[0],input_files[1])
        # get a list of output files, one for each pair of input files
        # get the names of the output files that are always written
        kneaddata_output_files = workflow.name_output_files(name=sample_names, tag="paired_1", subfolder="kneaddata", extension="fastq")
        kneaddata_output_files = zip(kneaddata_output_files, 
            workflow.name_output_files(name=sample_names, tag="paired_2", subfolder="kneaddata", extension="fastq"))
        
        # get the names for the output files that are sometimes written based on mapping
        kneaddata_optional_output_files = workflow.name_output_files(name=sample_names, tag="unmatched_1", subfolder="kneaddata", extension="fastq")
        kneaddata_optional_output_files = zip(kneaddata_optional_output_files,
            workflow.name_output_files(name=sample_names, tag="unmatched_2", subfolder="kneaddata", extension="fastq"))
        
        # add the second input file to the kneaddata arguments
        second_input_option=" --input [depends[1]] "
        # get the output folder
        kneaddata_output_folder = os.path.dirname(kneaddata_output_files[0][0])
    else:
        # get a list of output files one for each single-end input file
        kneaddata_output_files = workflow.name_output_files(name=sample_names, subfolder="kneaddata", extension="fastq")
        # the second input option is not used since these are single-end input files
        second_input_option=" "
        # get the output folder
        kneaddata_output_folder = os.path.dirname(kneaddata_output_files[0])
        
    # create the database command option string to provide zero or more databases to kneaddata
    if databases is None:
        optional_arguments=""
    elif isinstance(databases,list):
        # start the string with the kneaddata option and add an option for each database
        optional_arguments=" --reference-db "+" --reference-db ".join(databases)
    else:
        optional_arguments=" --reference-db " + databases
    
    # create a task for each set of input and output files to run kneaddata
    for sample, depends, targets in zip(sample_names, input_files, kneaddata_output_files):
        workflow.add_task_gridable(
            "kneaddata --input [depends[0]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]] "+second_input_option+optional_arguments,
            depends=depends,
            targets=targets,
            args=[kneaddata_output_folder, threads, sample],
            time=6*60, # 6 hours
            mem=12*1024, # 12 GB
            cores=threads) # time/mem based on 8 cores
    
    if paired:
        # if the inputs are paired, merge the final fastq output files into a single file
        # two output files are always written while the second two are sometimes printed
        # those that are optionally printed are not included in the depends
        kneaddata_merged_files = workflow.name_output_files(name=sample_names, subfolder="kneaddata", extension="fastq")
        for depends, args, target in zip(kneaddata_output_files, kneaddata_optional_output_files, kneaddata_merged_files):
            workflow.add_task_gridable(
                "cat [depends[0]] [depends[1]] [args[0]] [args[1]] > [targets[0]]",
                depends=depends,
                targets=target,
                args=args,
                time=10, # 10 minutes
                mem=1*1024, # 1 GB
                cores=1)
            
        # set the final output files to be the list of the merged files
        kneaddata_output_files = kneaddata_merged_files
    
    return kneaddata_output_files


def quality_control(workflow, input_files, threads, databases=None, pair_identifier=None):
    """Quality control tasks for whole genome shotgun sequences
    
    This set of tasks performs quality control on whole genome shotgun
    input files of single-end fastq format. It runs kneaddata using all of the 
    databases provided. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to kneaddata.
        threads (int): The number of threads/cores for kneaddata to use.
        databases (string/list): The databases to use with kneaddata (optional).
        pair_identifer (string): The string in the file basename to identify
            the first pair in the set (optional).
        
    Requires:
        kneaddata v0.5.3+: A tool to perform quality control on metagenomic and
            metatranscriptomic sequencing data
        
    Returns:
        list: A list of the filtered fastq files created by kneaddata.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import whole_genome_shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = whole_genome_shotgun.quality_control(workflow,
            ["demo.fastq","demo2.fastq"], 1)
            
        # run the workflow
        workflow.go()
        
    Todo:
        * Add option parameter to allow for setting options like fastqc.
    """
    
    # check for paired input files
    input_pair1, input_pair2 = utilities.paired_files(input_files, pair_identifier)

    
    paired = False
    if input_pair1:
        paired = True
        input_files = [input_pair1, input_pair2]
    
    # create a task for each set of input and output files to run kneaddata
    kneaddata_output_files=kneaddata(workflow, input_files, threads, paired, databases)
    
    return kneaddata_output_files


def taxonomic_profile(workflow,input_files,threads):
    """Taxonomic profile for whole genome shotgun sequences
    
    This set of tasks performs taxonomic profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files already run through quality control.
        threads (int): The number of threads/cores for metaphlan2 to use.
        
    Requires:
        metaphlan2 v2.5.0+: A tool to profile the composition of microbial communities.
        humann2 v0.9.6+: A tool for functional profiling (only humann2_join_tables is required).
        
    Returns:
        string: A file of the merged taxonomic profiles from all samples.
        list: A list of the sam files generated by metaphlan2.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import whole_genome_shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = whole_genome_shotgun.quality_control(workflow,
            ["demo.fastq","demo2.fastq"], 1)
        
        # run taxonomic profile
        taxonomic_profile, sam_outputs = whole_genome_shotgun.taxonomic_profile(
            workflow, filtered_fastq, 1) 
            
        # run the workflow
        workflow.go()
        
    Todo:
        * Add option for fasta input files.
        * Add option to provide paired input files which are merged then run.
    """
    
    # get a list of metaphlan2 output files, one for each input file
    metaphlan2_profile_tag="taxonomic_profile"
    metaphlan2_output_files_profile = workflow.name_output_files(name=input_files, subfolder="metaphlan2", tag=metaphlan2_profile_tag, extension="tsv")
    metaphlan2_output_files_sam = workflow.name_output_files(name=input_files, subfolder="metaphlan2", tag="bowtie2", extension="sam")
    metaphlan2_output_folder = os.path.dirname(metaphlan2_output_files_profile[0])
    
    # run metaphlan2 on each of the kneaddata output files
    workflow.add_task_group_gridable(
        "metaphlan2.py [depends[0]] --input_type fastq --output_file [targets[0]] --samout [targets[1]] --nproc [args[0]] --no_map --tmp_dir [args[1]]",
        depends=input_files,
        targets=zip(metaphlan2_output_files_profile, metaphlan2_output_files_sam),
        args=[threads,metaphlan2_output_folder],
        time=3*60, # 3 hours
        mem=12*1024, # 12 GB
        cores=threads) # time/mem based on 8 cores
    
    # merge all of the metaphlan taxonomy tables
    metaphlan2_merged_output = workflow.name_output_files(name="taxonomic_profiles.tsv")
    
    # run the humann2 join script to merge all of the metaphlan2 profiles
    workflow.add_task(
        "humann2_join_tables --input [args[0]] --output [targets[0]] --file_name [args[1]]",
        depends=metaphlan2_output_files_profile,
        targets=metaphlan2_merged_output,
        args=[metaphlan2_output_folder, metaphlan2_profile_tag])
    
    return metaphlan2_merged_output, metaphlan2_output_files_profile, metaphlan2_output_files_sam

def functional_profile(workflow,input_files,threads,taxonomic_profiles=None):
    """Functional profile for whole genome shotgun sequences
    
    This set of tasks performs functional profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. Optionally the taxonomic
    profiles can be provided for the samples.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq (or fasta) files already run through quality control.
        threads (int): The number of threads/cores for kneaddata to use.
        taxonomic_profiles (list): A set of taxonomic profiles, one per sample (optional).
        
    Requires:
        humann2 v0.9.6+: A tool for functional profiling.
        
    Returns:
        string: A file of the merged gene families (relative abundance) for all samples.
        string: A file of the merged ecs (relative abundance) for all samples.
        string: A file of the merged pathway abundances (relative abundance) for all samples.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import whole_genome_shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = whole_genome_shotgun.quality_control(workflow,
            ["demo.fastq","demo2.fastq"], 1)
        
        # run functional profiling
        genefamilies_file, ecs_file, pathabundance_file = whole_genome_shotgun.functional_profile(
            workflow, filtered_fastq, 1) 
            
        # run the workflow
        workflow.go()
    """
    
    ### Step 1: Run humann2 on all input files ###

    # get a list of output files, one for each input file, with the humann2 output file names
    genefamiles = workflow.name_output_files(name=input_files, subfolder="humann2", tag="genefamilies", extension="tsv")
    pathabundance = workflow.name_output_files(name=input_files, subfolder="humann2", tag="pathabundance", extension="tsv")
    pathcoverage = workflow.name_output_files(name=input_files, subfolder="humann2", tag="pathcoverage", extension="tsv")
    
    humann2_output_folder = os.path.dirname(genefamiles[0])
    
    # if taxonomic profiles are provided, add these to the targets and the command option
    if taxonomic_profiles:
        optional_profile_args=" --taxonomic-profile [depends[1]] "
        depends=zip(input_files,taxonomic_profiles)
    else:
        optional_profile_args=""
        depends=input_files
    
    # create a task to run humann2 on each of the kneaddata output files
    workflow.add_task_group_gridable(
        "humann2 --input [depends[0]] --output [args[0]] --threads [args[1]]"+optional_profile_args,
        depends=depends,
        targets=zip(genefamiles, pathabundance, pathcoverage),
        args=[humann2_output_folder, threads],
        time=24*60, # 24 hours
        mem=36*1024, # 36 GB
        cores=threads)
    
    ### STEP #2: Regroup UniRef90 gene families to ecs ###
    
    # get a list of all output ec files
    ec_files = workflow.name_output_files(name=genefamiles, subfolder="humann2", tag="ecs")
    
    # get ec files for all of the gene families files
    workflow.add_task_group_gridable(
        "humann2_regroup_table --input [depends[0]] --output [targets[0]] --groups uniref90_level4ec",
        depends=genefamiles,
        targets=ec_files,
        time=10*60, # 10 minutes
        mem=5*1024, # 5 GB
        cores=1)
    
    ### STEP #3: Normalize gene families, ecs, and pathway abundance to relative abundance (then merge files) ###
    
    # get a list of files for normalized ec, gene families, and pathway abundance
    norm_genefamily_files = workflow.name_output_files(name=genefamiles, subfolder="genes", tag="relab")
    norm_ec_files = workflow.name_output_files(name=ec_files, subfolder="ecs", tag="relab")
    norm_pathabundance_files = workflow.name_output_files(name=pathabundance, subfolder="pathways", tag="relab")
    
    # normalize the genefamily, ec, and pathabundance files
    workflow.add_task_group_gridable(
        "humann2_renorm_table --input [depends[0]] --output [targets[0]] --units relab",
        depends=genefamiles + ec_files + pathabundance,
        targets=norm_genefamily_files + norm_ec_files + norm_pathabundance_files,
        time=5*60, # 5 minutes
        mem=5*1024, # 5 GB
        cores=1)
    
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
        
    return merged_genefamilies, merged_ecs, merged_pathabundance
