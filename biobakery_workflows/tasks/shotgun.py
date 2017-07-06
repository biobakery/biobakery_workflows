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
import subprocess

from anadama2.tracked import TrackedExecutable

from biobakery_workflows import utilities
from biobakery_workflows import files
from biobakery_workflows import data

def kneaddata(workflow, input_files, output_folder, threads, paired=None, 
    databases=None, pair_identifier=None, additional_options=None, remove_intermediate_output=None):
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
            Allow for a single path or multiple paths in one string comma-delimited.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        additional_options (string): Additional options when running kneaddata (optional).
        remove_intermediate_output (bool): Remove intermediate output files.
        
    Requires:
        kneaddata v0.6.1+: A tool to perform quality control on metagenomic and
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
        second_input_option=" --input [depends[1]] --cat-final-output --serial "
        # determine time/memory equations based on the two input files
        time_equation="6*60 if ( file_size('[depends[0]]') + file_size('[depends[1]]') ) < 25 else 4*6*60"
        mem_equation="12*1024 if ( file_size('[depends[0]]') + file_size('[depends[1]]') ) < 25 else 2*12*1024"
    else:
        # the second input option is not used since these are single-end input files
        second_input_option=" "
        # determine time/memory equations based on the single input file
        time_equation="6*60 if file_size('[depends[0]]') < 25 else 4*6*60"
        mem_equation="12*1024 if file_size('[depends[0]]') < 25 else 2*12*1024"
        
    # set additional options to empty string if not provided
    if additional_options is None:
        additional_options=""

    # create the database command option string to provide zero or more databases to kneaddata
    if databases is None:
        optional_arguments=""
    elif isinstance(databases,list):
        # start the string with the kneaddata option and add an option for each database
        optional_arguments=" --reference-db "+" --reference-db ".join(databases)
    elif isinstance(databases,basestring) and "," in databases:
        # split the paths by comma
        database_list=list(filter(lambda x: x, databases.split(",")))
        # start the string with the kneaddata option and add an option for each database
        optional_arguments=" --reference-db "+" --reference-db ".join(database_list)        
    else:
        optional_arguments=" --reference-db " + databases
        
    # add option to remove intermediate output, if set
    if remove_intermediate_output:
        additional_options+=" --remove-intermediate-output "
    
    # create a task for each set of input and output files to run kneaddata
    for sample, depends, targets in zip(sample_names, input_files, kneaddata_output_files):
        workflow.add_task_gridable(
            "kneaddata --input [depends[0]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]] "+second_input_option+optional_arguments+" "+additional_options,
            depends=utilities.add_to_list(depends,TrackedExecutable("kneaddata")),
            targets=targets,
            args=[kneaddata_output_folder, threads, sample],
            time=time_equation, # 6 hours or more depending on file size
            mem=mem_equation, # 12 GB or more depending on file size
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
        kneaddata v0.6.1+: A tool to perform quality control on metagenomic and
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
    kneaddata_read_count_file = files.ShotGun.path("kneaddata_read_counts",output_folder,create_folder=True)
    
    # add the task (which is not gridable as this task should take under 5 minutes)
    workflow.add_task("kneaddata_read_count_table --input [args[0]] --output [targets[0]]",
        depends=input_files,
        targets=kneaddata_read_count_file,
        args=[input_folder],
        name="kneaddata_read_count_table")
    
    return kneaddata_read_count_file


def quality_control(workflow, input_files, output_folder, threads, databases=None, 
    pair_identifier=None, additional_options=None, remove_intermediate_output=None):
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
        additional_options (string): Additional options when running kneaddata (optional).
        remove_intermediate_output (bool): Remove intermediate output files.
        
    Requires:
        kneaddata v0.6.1+: A tool to perform quality control on metagenomic and
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
    if pair_identifier:
        input_pair1, input_pair2 = utilities.paired_files(input_files, pair_identifier)
    else:
        input_pair1 = []
    
    paired = False
    if input_pair1:
        paired = True
        input_files = [input_pair1, input_pair2]
    
    # create a task for each set of input and output files to run kneaddata
    kneaddata_output_fastq, kneaddata_output_logs=kneaddata(workflow, input_files, 
        output_folder, threads, paired, databases, pair_identifier, additional_options,
        remove_intermediate_output)
    
    # create the read count table
    kneaddata_read_count_file=kneaddata_read_count_table(workflow, kneaddata_output_logs, output_folder)
    
    return kneaddata_output_fastq, kneaddata_read_count_file


def taxonomic_profile(workflow,input_files,output_folder,threads,input_extension):
    """Taxonomic profile for whole genome shotgun sequences
    
    This set of tasks performs taxonomic profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files already run through quality control.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for metaphlan2 to use.
        input_extension (string): The extension for the input files.
        
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
    
    # determine the input file type based on the extension
    if input_extension in ["fasta","fasta.gz","fa","fa.gz"]:
        input_type="fasta"
    else:
        input_type="fastq"
    
    # run metaphlan2 on each of the kneaddata output files
    for depend_fastq, target_profile, target_sam in zip(input_files, metaphlan2_output_files_profile, metaphlan2_output_files_sam):
        workflow.add_task_gridable(
            "metaphlan2.py [depends[0]] --input_type [args[2]] --output_file [targets[0]] --samout [targets[1]] --nproc [args[0]] --no_map --tmp_dir [args[1]]",
            depends=[depend_fastq,TrackedExecutable("metaphlan2.py")],
            targets=[target_profile,target_sam],
            args=[threads,metaphlan2_output_folder,input_type],
            time=3*60, # 3 hours
            mem=12*1024, # 12 GB
            cores=threads) # time/mem based on 8 cores
    
    # merge all of the metaphlan taxonomy tables
    metaphlan2_merged_output = files.ShotGun.path("taxonomic_profile", output_folder)
    
    # run the humann2 join script to merge all of the metaphlan2 profiles
    workflow.add_task(
        "humann2_join_tables --input [args[0]] --output [targets[0]] --file_name [args[1]]",
        depends=metaphlan2_output_files_profile,
        targets=metaphlan2_merged_output,
        args=[metaphlan2_output_folder, metaphlan2_profile_tag])
   
    # get the name for the file to write the species counts
    metaphlan2_species_counts_file = files.ShotGun.path("species_counts",output_folder,create_folder=True)

    # create a file of species counts
    workflow.add_task(
    "count_features.py --input [depends[0]] --output [targets[0]] --include s__ --filter t__ --reduce-sample-name",
    depends=metaphlan2_merged_output,
    targets=metaphlan2_species_counts_file) 

    return metaphlan2_merged_output, metaphlan2_output_files_profile, metaphlan2_output_files_sam

def functional_profile(workflow,input_files,output_folder,threads,taxonomic_profiles=None, remove_intermediate_output=None):
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
        remove_intermediate_output (bool): Remove intermediate output files.
        
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

    # get the list of the log files that will be created, one for each input file, set to the same folder as the main output files
    log_files = utilities.name_files(sample_names, output_folder, subfolder="humann2", extension="log")
    # get the name for the file of read and species counts created from the humann2 log outputs
    log_counts = files.ShotGun.path("humann2_read_counts",output_folder,create_folder=True)

    humann2_output_folder = os.path.dirname(genefamiles[0])
    
    # if taxonomic profiles are provided, add these to the targets and the command option
    if taxonomic_profiles:
        optional_profile_args=" --taxonomic-profile [depends[1]] "
        depends=zip(input_files,taxonomic_profiles)
    else:
        optional_profile_args=""
        depends=input_files
        
    # remove intermediate temp files, if set
    if remove_intermediate_output:
        optional_profile_args+=" --remove-temp-output "
    
    # create a task to run humann2 on each of the kneaddata output files
    for depend_fastq, target_gene, target_path, target_coverage, target_log in zip(depends, genefamiles, pathabundance, pathcoverage, log_files):
        workflow.add_task_gridable(
            "humann2 --input [depends[0]] --output [args[0]] --o-log [targets[3]] --threads [args[1]]"+optional_profile_args,
            depends=utilities.add_to_list(depend_fastq,TrackedExecutable("humann2")),
            targets=[target_gene, target_path, target_coverage, target_log],
            args=[humann2_output_folder, threads],
            time="24*60 if file_size('[depends[0]]') < 25 else 4*24*60", # 24 hours or more depending on file size
            mem="32*1024 if file_size('[depends[0]]') < 25 else 2*32*1024", # 32 GB or more depending on file size
            cores=threads)

    # create a task to get the read and species counts for each humann2 run from the log files
    workflow.add_task(
        "get_counts_from_humann2_logs.py --input [args[0]] --output [targets[0]]",
        depends=log_files,
        targets=log_counts,
        args=humann2_output_folder)
    
    ### STEP #2: Regroup UniRef90 gene families to ecs ###
    
    # get a list of all output ec files
    ec_files = utilities.name_files(sample_names, output_folder, subfolder="humann2", tag="ecs", extension="tsv")
    
    # get ec files for all of the gene families files
    workflow.add_task_group_gridable(
        "humann2_regroup_table --input [depends[0]] --output [targets[0]] --groups uniref90_level4ec",
        depends=genefamiles,
        targets=ec_files,
        time=10*60, # 10 minutes
        mem=5*1024, # 5 GB
        cores=1)
    
    ### STEP #3: Merge gene families, ecs, and pathway abundance files

    # get a list of merged files for ec, gene families, and pathway abundance
    merged_genefamilies = files.ShotGun.path("genefamilies", output_folder)
    merged_ecs = files.ShotGun.path("ecs", output_folder)
    merged_pathabundance = files.ShotGun.path("pathabundance", output_folder)
    
    # merge the ec, gene families, and pathway abundance files
    all_depends=[genefamiles, ec_files, pathabundance]
    all_targets=[merged_genefamilies, merged_ecs, merged_pathabundance]
    file_basenames=["genefamilies","ecs","pathabundance"]
    for depends, targets, basename in zip(all_depends, all_targets, file_basenames):
        workflow.add_task(
            "humann2_join_tables --input [args[0]] --output [targets[0]] --file_name [args[1]]",
            depends=depends,
            targets=targets,
            args=[os.path.dirname(depends[0]),basename])
    
    ### STEP #4: Normalize gene families, ecs, and pathway abundance to relative abundance (then merge files) ###
    
    # get a list of files for normalized ec, gene families, and pathway abundance
    norm_genefamily_files = utilities.name_files(genefamiles, output_folder, subfolder="genes", tag="relab", create_folder=True)
    norm_ec_files = utilities.name_files(ec_files, output_folder,  subfolder="ecs", tag="relab", create_folder=True)
    norm_pathabundance_files = utilities.name_files(pathabundance, output_folder,  subfolder="pathways", tag="relab", create_folder=True)
    
    # normalize the genefamily, ec, and pathabundance files
    # do not include special features (ie UNMAPPED, UNINTEGRATED, UNGROUPED) in norm computation
    workflow.add_task_group_gridable(
        "humann2_renorm_table --input [depends[0]] --output [targets[0]] --units relab --special n",
        depends=genefamiles + ec_files + pathabundance,
        targets=norm_genefamily_files + norm_ec_files + norm_pathabundance_files,
        time=5*60, # 5 minutes
        mem=5*1024, # 5 GB
        cores=1)
    
    # get a list of merged files for ec, gene families, and pathway abundance
    merged_genefamilies_relab = files.ShotGun.path("genefamilies_relab", output_folder)
    merged_ecs_relab = files.ShotGun.path("ecs_relab", output_folder)
    merged_pathabundance_relab = files.ShotGun.path("pathabundance_relab", output_folder)
    
    # merge the ec, gene families, and pathway abundance files
    all_depends=[norm_genefamily_files, norm_ec_files, norm_pathabundance_files]
    all_targets=[merged_genefamilies_relab, merged_ecs_relab, merged_pathabundance_relab]
    for depends, targets in zip(all_depends, all_targets):
        workflow.add_task(
            "humann2_join_tables --input [args[0]] --output [targets[0]]",
            depends=depends,
            targets=targets,
            args=[os.path.dirname(depends[0])])

    # get feature counts for the ec, gene families, and pathways
    genefamilies_counts = files.ShotGun.path("genefamilies_relab_counts", output_folder)
    ecs_counts = files.ShotGun.path("ecs_relab_counts", output_folder)
    pathabundance_counts = files.ShotGun.path("pathabundance_relab_counts", output_folder)
    workflow.add_task_group(
        "count_features.py --input [depends[0]] --output [targets[0]] --reduce-sample-name --ignore-un-features --ignore-stratification",
        depends=[merged_genefamilies_relab, merged_ecs_relab, merged_pathabundance_relab],
        targets=[genefamilies_counts, ecs_counts, pathabundance_counts])
    
    # merge the feature counts into a single file
    all_feature_counts = files.ShotGun.path("feature_counts", output_folder)
    workflow.add_task(
        "humann2_join_tables --input [args[0]] --output [targets[0]] --file_name _relab_counts.tsv",
        depends=[genefamilies_counts, ecs_counts, pathabundance_counts],
        targets=all_feature_counts,
        args=[os.path.dirname(genefamilies_counts)])

        
    return merged_genefamilies_relab, merged_ecs_relab, merged_pathabundance_relab, merged_genefamilies, merged_ecs, merged_pathabundance

def norm_ratio(workflow, wms_genes, wms_ecs, wms_paths, wts_genes, wts_ecs, wts_paths, output_folder, mapping=None):
    """Compute a rna/dna normalized ratio
    
    Compute the ratio for the genes, ecs, and pathway abundances using the original
    output files generated by HUMAnN2 (merged for all samples). These files have not
    had any normalization applied.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        wms_genes (string): A path to the gene file for the metagenome samples.
        wms_ecs (string): A path to the ecs file for the metagenome samples.
        wms_paths (string): A path to the pathway abundance file for the metagenome samples.
        wts_genes (string): A path to the gene file for the metatranscriptome samples.
        wts_ecs (string): A path to the ecs file for the metatranscriptome  samples.
        wts_paths (string): A path to the pathway abundance file for the metatranscriptome  samples..
        output_folder (string): The path of the output folder.
        mapping (string): A file of names mapping rna samples to dna samples (optional).
        
    Requires:
        None.
        
    Returns:
        string: A file of the normed gene families for all samples.
        string: A file of the normed ecs for all samples.
        string: A file of the normed pathway abundances for all samples.
        
    """
    
    # get the names of the output files
    norm_ratio_genes = files.ShotGun.path("genefamilies_norm_ratio",output_folder,create_folder=True)
    norm_ratio_ecs = files.ShotGun.path("ecs_norm_ratio",output_folder)
    norm_ratio_pathway = files.ShotGun.path("paths_norm_ratio",output_folder)
    
    # if a mapping file is provided, set the mapping option
    mapping_option=""
    if mapping:
        mapping_option=" --mapping [args[1]]"
    
    # add tasks to renorm each set of files (genes, ecs, pathway abundance)
    all_depends=[(wms_genes,wts_genes),(wms_ecs,wts_ecs),(wms_paths,wts_paths)]
    all_targets=[norm_ratio_genes,norm_ratio_ecs,norm_ratio_pathway]
    for depend, target in zip(all_depends, all_targets):
        workflow.add_task(
            "rna_dna_norm.py --input-dna [depends[0]] --input-rna [depends[1]] --output [args[0]] --reduce-sample-name" + mapping_option,
            depends=depend,
            targets=target,
            args=[os.path.dirname(target), mapping])
        
    return norm_ratio_genes, norm_ratio_ecs, norm_ratio_pathway

def strainphlan(task,threads,clade_number,clade_list,reference_folder,marker_folder):
    """ Run strainphlan for the specific clade
    
    Args:
        task: (anadama2.task): An instance of the task class.
        threads: (int): Total threads to use for task.
        clade_number: (int): The number of clade to run.
        clade_list: (string): The path to the clade list.
        reference_folder (string): The folder containing the reference files.
        marker_folder (string): The folder containing the marker files.

    Requires:
        StrainPhlAn: A tool for strain profiling.    
        
    """        
    
    # find the name of the clade in the list
    with open(clade_list) as file_handle:
        clades=[line.strip() for line in filter(lambda line: line.startswith("s__"), file_handle.readlines())]
        try:
            profile_clade=clades[clade_number]
        except IndexError:
            profile_clade=None
            
    if profile_clade:
        command = "strainphlan.py --ifn_samples [args[0]]/*.markers --output_dir [args[1]] "+\
            "--clades [args[2]] --nprocs_main [args[3]] --keep_alignment_files"
            
        # add the marker files to the command
        all_marker_file=os.path.join(marker_folder,"all_markers.fasta")
        marker_file=os.path.join(marker_folder,profile_clade+".markers.fasta")
        command += " --ifn_markers "+marker_file
        
        # generate the marker file if it does not already exist
        if not os.path.isfile(marker_file):
            # get the pkl file relative to the strainphlan install
            try:
                strainphlan_pkl=os.path.join(os.path.dirname(subprocess.check_output(["which","strainphlan.py"])),"db_v20","mpa_v20_m200.pkl")
            except subprocess.CalledProcessError:
                raise EnvironmentError("Unable to find strainphlan install.")
            
            marker_command="extract_markers.py --mpa_pkl [depends[0]] "+\
                "--ifn_markers [depends[1]] --clade [args[0]] "+\
                "--ofn_markers [targets[0]]"
                
            # create the marker file in the output folder
            marker_file=os.path.join(os.path.dirname(task.targets[0].name),profile_clade+".markers.fasta")
                
            return_code = utilities.run_task(marker_command,depends=[strainphlan_pkl,all_marker_file], 
                targets=[marker_file],args=[profile_clade])
        
        # check that the marker file exists
        if not os.path.isfile(marker_file):
            raise EnvironmentError("Unable to find StrainPhlAn markers for clade "+ profile_clade)
        
        # get the list of reference genomes
        genomes=set()
        with open(data.get_file("strainphlan_species_gcf.tsv")) as file_handle:
            for line in file_handle:
                if line.startswith(profile_clade):
                    genomes.add(os.path.join(reference_folder,line.rstrip().split("\t")[-1]+".fna.bz2"))
        
        # only use those references found in the folder
        genomes = filter(os.path.isfile,genomes)
        
        
        # add the reference genome files to the command, if any are found
        if len(list(genomes)):
            command += " --ifn_ref_genomes " + " --ifn_ref_genomes ".join(genomes)
        
        # write the output to the log
        command += " > [targets[0]]"
        
    else:
        # there is not a clade of this number, create an empty output file
        command = "touch [targets[0]]"
        
    # run the task
    return_code = utilities.run_task(command, depends=task.depends, targets=task.targets, 
        args=[os.path.dirname(task.depends[0].name),os.path.dirname(task.targets[0].name),profile_clade,threads])

def strain_profile(workflow,sam_files,output_folder,threads,reference_folder,marker_folder,max_species=10):
    """Strain profile for whole genome shotgun sequences
    
    This set of tasks performs strain profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        sam_files (list): A list of paths to sam files generated by MetaPhlAn2.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores to use.
        reference_folder (string): The folder containing the reference files.
        marker_folder (string): The folder containing the marker files.
        max_species (int): The maximum number of species to profile.
        
    Requires:
        StrainPhlAn: A tool for strain profiling.
        
    Returns:
        None
        
    """    

    ### STEP #1: Identify markers for each of the samples
    # name the marker files based on the sam files
    strainphlan_markers = utilities.name_files(sam_files, output_folder, subfolder="strainphlan", extension="markers", create_folder=True)
     
    # create a marker file from each sam file
    for sam, markers in zip(sam_files, strainphlan_markers):
        workflow.add_task_gridable(
            "sample2markers.py --ifn_samples [depends[0]] --input_type sam --output_dir [args[0]] --nprocs [args[1]]",
            depends=sam,
            targets=markers,
            args=[os.path.dirname(markers),threads],
            time=15*60, # 15 minutes
            mem=5*1024, # 5 GB
            cores=threads)

    ### STEP #2: Find the top species by average abundance
    clade_list = utilities.name_files("clades_list.txt", output_folder, subfolder="strainphlan")
    
    workflow.add_task(
        "strainphlan.py --ifn_samples [args[0]]/*.markers --output_dir [args[0]] --print_clades_only > [targets[0]]",
        depends=strainphlan_markers,
        targets=clade_list,
        args=os.path.dirname(strainphlan_markers[0]),
        name="strainphlan_print_clades")
    
    ### STEP #3: Run strainphlan on the top set of clades identified
    clade_logs = utilities.name_files(map(str,range(10)), output_folder, tag="clade", subfolder="strainphlan", extension="log")
    for clade_number in range(max_species):
        workflow.add_task_gridable(
            utilities.partial_function(strainphlan,threads=threads,clade_number=clade_number,
                clade_list=clade_list,reference_folder=os.path.abspath(reference_folder),marker_folder=os.path.abspath(marker_folder)),
            depends=strainphlan_markers,
            targets=clade_logs[clade_number-1],                           
            time=6*60, # 6 hours
            mem=10*1024, # 10 GB
            cores=threads)

     
