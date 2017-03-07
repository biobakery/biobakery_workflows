"""
bioBakery Workflows: tasks.shotgun module
A collection of tasks for workflows with shotgun sequences

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

def kneaddata(workflow, input_files, output_folder, threads, paired=None, databases=None, pair_identifier=None):
    """Run kneaddata
    
    This set of tasks will run kneaddata on the input files provided. It will run with
    single-end or paired-end input files.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to kneaddata. This
        is a list of lists if the input files are paired.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for kneaddata to use.
        paired (bool): This indicates if the input files are paired.
        databases (string/list): The databases to use with kneaddata (optional).
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        
    Requires:
        kneaddata v0.5.4+: A tool to perform quality control on metagenomic and
            metatranscriptomic sequencing data
        
    Returns:
        list: A list of the filtered fastq files created by kneaddata.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = shotgun.kneaddata(workflow,
            ["demo.fastq","demo2.fastq"], 1)
            
        # run the workflow
        workflow.go()
    """

    # get the sample basenames from the input files
    if paired:
        sample_names=utilities.sample_names(input_files[0],pair_identifier)
    else:
        sample_names=utilities.sample_names(input_files)
        
    # get the kneaddata final output files
    kneaddata_output_fastq = utilities.name_files(sample_names, output_folder, subfolder="kneaddata", extension="fastq", create_folder=True)
    kneaddata_output_logs = utilities.name_files(sample_names, output_folder, subfolder="kneaddata", extension="log")
    kneaddata_output_files = zip(kneaddata_output_fastq, kneaddata_output_logs)

    # get the output folder
    kneaddata_output_folder = os.path.dirname(kneaddata_output_files[0][0])
        
    if paired:
        # reorder the input files so they are a set of paired files
        input_files=zip(input_files[0],input_files[1])
        # add the second input file to the kneaddata arguments
        # also add the option to cat the final output files into a single file
        second_input_option=" --input [depends[1]] --cat-final-output"
    else:
        # the second input option is not used since these are single-end input files
        second_input_option=" "

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
    
    return kneaddata_output_fastq, kneaddata_output_logs

def kneaddata_read_count_table(workflow, input_files, output_folder):
    """Run kneaddata_read_count_table
    
    This task will run kneaddata_read_count_table on the input files provided.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to log files for input to kneaddata. 
            All input files are expected to be in the same folder and all should
            have the log extension. 
        output_folder (string): The path of the output folder.
        
    Requires:
        kneaddata v0.5.4+: A tool to perform quality control on metagenomic and
            metatranscriptomic sequencing data
        
    Returns:
        string: The path to the read count table written.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add a task to run kneaddata_read_count_table
        read_count_file = shotgun.kneaddata_read_count_table(workflow,
            ["demo.log","demo2.log"])
            
        # run the workflow
        workflow.go()
    """
    
    # get the folder of the log files
    input_folder=os.path.dirname(input_files[0])
    
    # get the name for the output file
    kneaddata_read_count_file = utilities.name_files("kneaddata_read_count_table.tsv",output_folder,subfolder="counts",create_folder=True)
    
    # add the task (which is not gridable as this task should take under 5 minutes)
    workflow.add_task("kneaddata_read_count_table --input [args[0]] --output [targets[0]]",
        depends=input_files,
        targets=kneaddata_read_count_file,
        args=[input_folder])
    
    return kneaddata_read_count_file


def quality_control(workflow, input_files, output_folder, threads, databases=None, pair_identifier=None):
    """Quality control tasks for whole genome shotgun sequences
    
    This set of tasks performs quality control on whole genome shotgun
    input files of single-end fastq format. It runs kneaddata using all of the 
    databases provided. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to kneaddata.
        threads (int): The number of threads/cores for kneaddata to use.
        output_folder (string): The path of the output folder.
        databases (string/list): The databases to use with kneaddata (optional).
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        
    Requires:
        kneaddata v0.5.4+: A tool to perform quality control on metagenomic and
            metatranscriptomic sequencing data
        
    Returns:
        list: A list of the filtered fastq files created by kneaddata.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = shotgun.quality_control(workflow,
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
    kneaddata_output_fastq, kneaddata_output_logs=kneaddata(workflow, input_files, output_folder, threads, paired, databases, pair_identifier)
    
    # create the read count table
    kneaddata_read_count_file=kneaddata_read_count_table(workflow, kneaddata_output_logs, output_folder)
    
    return kneaddata_output_fastq, kneaddata_read_count_file


def taxonomic_profile(workflow,input_files,output_folder,threads):
    """Taxonomic profile for whole genome shotgun sequences
    
    This set of tasks performs taxonomic profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files already run through quality control.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for metaphlan2 to use.
        
    Requires:
        metaphlan2 v2.5.0+: A tool to profile the composition of microbial communities.
        humann2 v0.9.6+: A tool for functional profiling (only humann2_join_tables is required).
        
    Returns:
        string: A file of the merged taxonomic profiles from all samples.
        list: A list of the sam files generated by metaphlan2.
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = shotgun.quality_control(workflow,
            ["demo.fastq","demo2.fastq"], 1)
        
        # run taxonomic profile
        taxonomic_profile, sam_outputs = shotgun.taxonomic_profile(
            workflow, filtered_fastq, 1) 
            
        # run the workflow
        workflow.go()
        
    Todo:
        * Add option for fasta input files.
        * Add option to provide paired input files which are merged then run.
    """
    
    # get the sample names from the input files
    sample_names=utilities.sample_names(input_files)
    
    # get a list of metaphlan2 output files, one for each input file
    metaphlan2_profile_tag="taxonomic_profile"
    metaphlan2_output_files_profile = utilities.name_files(sample_names, output_folder, subfolder="metaphlan2", tag=metaphlan2_profile_tag, extension="tsv", create_folder=True)
    metaphlan2_output_files_sam = utilities.name_files(sample_names, output_folder, subfolder="metaphlan2", tag="bowtie2", extension="sam")
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
    metaphlan2_merged_output = utilities.name_files("taxonomic_profiles.tsv", output_folder)
    
    # run the humann2 join script to merge all of the metaphlan2 profiles
    workflow.add_task(
        "humann2_join_tables --input [args[0]] --output [targets[0]] --file_name [args[1]]",
        depends=metaphlan2_output_files_profile,
        targets=metaphlan2_merged_output,
        args=[metaphlan2_output_folder, metaphlan2_profile_tag])
   
    # get the name for the file to write the species counts
    metaphlan2_species_counts_file = utilities.name_files("metaphlan2_species_counts_table.tsv",output_folder,subfolder="counts",create_folder=True)

    # create a file of species counts
    workflow.add_task(
    "count_features.py --input [depends[0]] --output [targets[0]] --include s__ --filter t__ --reduce-sample-name"
    depends=metaphlan2_merged_output,
    targets=metaphlan2_species_counts_file) 

    return metaphlan2_merged_output, metaphlan2_output_files_profile, metaphlan2_output_files_sam

def functional_profile(workflow,input_files,output_folder,threads,taxonomic_profiles=None):
    """Functional profile for whole genome shotgun sequences
    
    This set of tasks performs functional profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. Optionally the taxonomic
    profiles can be provided for the samples.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq (or fasta) files already run through quality control.
        output_folder (string): The path of the output folder.
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
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        filtered_fastq = shotgun.quality_control(workflow,
            ["demo.fastq","demo2.fastq"], 1)
        
        # run functional profiling
        genefamilies_file, ecs_file, pathabundance_file = shotgun.functional_profile(
            workflow, filtered_fastq, 1) 
            
        # run the workflow
        workflow.go()
    """
    
    # get the sample names from the input files
    sample_names=utilities.sample_names(input_files)
    
    ### Step 1: Run humann2 on all input files ###

    # get a list of output files, one for each input file, with the humann2 output file names
    genefamiles = utilities.name_files(sample_names, output_folder, subfolder="humann2", tag="genefamilies", extension="tsv", create_folder=True)
    pathabundance = utilities.name_files(sample_names, output_folder, subfolder="humann2", tag="pathabundance", extension="tsv")
    pathcoverage = utilities.name_files(sample_names, output_folder, subfolder="humann2", tag="pathcoverage", extension="tsv")

    # get the list of the log files that will be created, one for each input file, each in a temp humann2 folder
    log_files = [os.path.join(output_folder, "humann2", sample+"_humann2_temp", sample)+".log" for sample in sample_names]
    # get the name for the file of read and species counts created from the humann2 log outputs
    log_counts = utilities.name_files("humann2_read_and_species_count_table.tsv",output_folder,subfolder="counts",create_folder=True)

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
        targets=zip(genefamiles, pathabundance, pathcoverage, log_files),
        args=[humann2_output_folder, threads],
        time=24*60, # 24 hours
        mem=36*1024, # 36 GB
        cores=threads)

    # create a task to get the read and species counts for each humann2 run from the log files
    workflow.add_task(
        "get_counts_from_humann2_logs.py --input [args[0]] --output [targets[0]]",
        depends=log_files,
        targets=log_counts,
        args=humann2_output_folder)
    
    ### STEP #2: Regroup UniRef90 gene families to ecs ###
    
    # get a list of all output ec files
    ec_files = utilities.name_files(genefamiles, output_folder, subfolder="humann2", tag="ecs")
    
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
    norm_genefamily_files = utilities.name_files(genefamiles, output_folder, subfolder="genes", tag="relab", create_folder=True)
    norm_ec_files = utilities.name_files(ec_files, output_folder,  subfolder="ecs", tag="relab", create_folder=True)
    norm_pathabundance_files = utilities.name_files(pathabundance, output_folder,  subfolder="pathways", tag="relab", create_folder=True)
    
    # normalize the genefamily, ec, and pathabundance files
    workflow.add_task_group_gridable(
        "humann2_renorm_table --input [depends[0]] --output [targets[0]] --units relab",
        depends=genefamiles + ec_files + pathabundance,
        targets=norm_genefamily_files + norm_ec_files + norm_pathabundance_files,
        time=5*60, # 5 minutes
        mem=5*1024, # 5 GB
        cores=1)
    
    # get a list of merged files for ec, gene families, and pathway abundance
    merged_genefamilies = utilities.name_files("genefamilies_relab.tsv", output_folder)
    merged_ecs = utilities.name_files("ecs_relab.tsv", output_folder)
    merged_pathabundance = utilities.name_files("pathabundance_relab.tsv", output_folder)
    
    # merge the ec, gene families, and pathway abundance files
    all_depends=[norm_genefamily_files, norm_ec_files, norm_pathabundance_files]
    all_targets=[merged_genefamilies, merged_ecs, merged_pathabundance]
    for depends, targets in zip(all_depends, all_targets):
        workflow.add_task(
            "humann2_join_tables --input [args[0]] --output [targets[0]]",
            depends=depends,
            targets=targets,
            args=[os.path.dirname(depends[0])])

    # get feature counts for the ec, gene families, and pathways
    genefamilies_counts = utilities.name_files("humann2_genefamilies_relab_counts.tsv", output_folder, subfolder="counts")
    ecs_counts = utilities.name_files("humann2_ecs_relab_counts.tsv", output_folder, subfolder="counts")
    pathabundance_counts = utilities.name_files("humann2_pathabundance_relab_counts.tsv", output_folder, subfolder="counts")
    workflow.add_task_group(
        "count_features.py --input [depends[0]] --output [targets[0]] --reduce-sample-name --ignore-un-features --ignore-stratification",
        depends=[merged_genefamilies, merged_ecs, merged_pathabundance],
        targets=[genefamilies_counts, ecs_counts, pathabundance_counts])

        
    return merged_genefamilies, merged_ecs, merged_pathabundance
