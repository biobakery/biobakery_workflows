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

import sys
import os
import subprocess
import itertools

from anadama2.tracked import TrackedExecutable, TrackedDirectory

from biobakery_workflows import utilities
from biobakery_workflows import files
from biobakery_workflows import data

# constants
BOWTIE2_EXTENSION=".1.bt2"

def kneaddata(workflow, input_files, extension, output_folder, threads, paired=None, 
    databases=None, pair_identifier=None, additional_options=None, remove_intermediate_output=None):
    """Run kneaddata
    
    This set of tasks will run kneaddata on the input files provided. It will run with
    single-end or paired-end input files.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to kneaddata. This
        is a list of lists if the input files are paired.
        extension (string): The extension for all files.
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
        sample_names=utilities.sample_names(input_files[0],extension,pair_identifier)
    else:
        sample_names=utilities.sample_names(input_files,extension)
        
    # get the kneaddata final output files
    main_folder=os.path.join("kneaddata","main")
    kneaddata_output_repeats_removed_fastq = utilities.name_files(sample_names, output_folder, subfolder=main_folder, extension="repeats.removed.fastq", create_folder=True)
    kneaddata_output_fastq = utilities.name_files(sample_names, output_folder, subfolder=main_folder, extension="fastq", create_folder=True)
    kneaddata_output_logs = utilities.name_files(sample_names, output_folder, subfolder=main_folder, extension="log")
    kneaddata_output_files = list(zip(kneaddata_output_fastq, kneaddata_output_logs))

    # get the output folder
    kneaddata_output_folder = os.path.dirname(kneaddata_output_files[0][0])

    rename_final_output = ""        
    if paired:
        # reorder the input files so they are a set of paired files
        input_files=list(zip(input_files[0],input_files[1]))
        # add the second input file to the kneaddata arguments
        # also add the option to cat the final output files into a single file
        second_input_option=" --input [depends[1]] --cat-final-output "
        # determine time/memory equations based on the two input files
        time_equation="3*6*60 if ( file_size('[depends[0]]') + file_size('[depends[1]]') ) < 10 else 5*6*60"
        mem_equation="3*12*1024 if ( file_size('[depends[0]]') + file_size('[depends[1]]') ) < 10 else 6*12*1024"
    else:
        # the second input option is not used since these are single-end input files
        second_input_option=" "
        # determine time/memory equations based on the single input file
        time_equation="3*6*60 if file_size('[depends[0]]') < 10 else 5*6*60"
        mem_equation="3*12*1024 if file_size('[depends[0]]') < 10 else 6*12*1024"
        # need to rename the final output file here to the sample name
        rename_final_output = " && mv [args[3]] [targets[0]]"
        
    # set additional options to empty string if not provided
    if additional_options is None:
        additional_options=""

    # always run with the serial option, which presents read counts in log in the manner expected
    # by the visualization workflows (in serial filtering in the order of the databases provided)
    # also run tandem repeat filtering
    additional_options+=" --serial --run-trf "

    # create the database command option string to provide zero or more databases to kneaddata
    if databases is None or databases=="":
        optional_arguments=""
    elif isinstance(databases,list):
        # start the string with the kneaddata option and add an option for each database
        optional_arguments=" --reference-db "+" --reference-db ".join(databases)
    elif isinstance(databases,str) and "," in databases:
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
    # rename file with repeats in name to only sample name
    for sample, depends, targets, intermediate_file in zip(sample_names, input_files, kneaddata_output_files, kneaddata_output_repeats_removed_fastq):
        workflow.add_task_gridable(
            "kneaddata --input [depends[0]] --output [args[0]] --threads [args[1]] --output-prefix [args[2]] "+second_input_option+optional_arguments+" "+additional_options+rename_final_output,
            depends=utilities.add_to_list(depends,TrackedExecutable("kneaddata")),
            targets=targets,
            args=[kneaddata_output_folder, threads, sample, intermediate_file],
            time=time_equation, # 6 hours or more depending on file size
            mem=mem_equation, # 12 GB or more depending on file size
            cores=threads, # time/mem based on 8 cores
            name=utilities.name_task(sample,"kneaddata")) # name task based on sample name
    
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


def quality_control(workflow, input_files, extension, output_folder, threads, databases=None, 
    pair_identifier=None, additional_options=None, remove_intermediate_output=None):
    """Quality control tasks for whole genome shotgun sequences
    
    This set of tasks performs quality control on whole genome shotgun
    input files of single-end fastq format. It runs kneaddata using all of the 
    databases provided. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to kneaddata.
        extension (string): The extension for all files.
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
        input_pair1, input_pair2 = utilities.paired_files(input_files, extension, pair_identifier)
    else:
        input_pair1 = []
    
    paired = False
    if input_pair1:
        paired = True
        input_files = [input_pair1, input_pair2]
    
    # create a task for each set of input and output files to run kneaddata
    kneaddata_output_fastq, kneaddata_output_logs=kneaddata(workflow, input_files, extension,
        output_folder, threads, paired, databases, pair_identifier, additional_options,
        remove_intermediate_output)
    
    # create the read count table
    kneaddata_read_count_file=kneaddata_read_count_table(workflow, kneaddata_output_logs, output_folder)
    
    return kneaddata_output_fastq, kneaddata_read_count_file


def taxonomic_profile(workflow,input_files,output_folder,threads,input_extension,already_profiled=False):
    """Taxonomic profile for whole genome shotgun sequences
    
    This set of tasks performs taxonomic profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files already run through quality control.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for metaphlan to use.
        input_extension (string): The extension for the input files.
        already_profiled (bool): Indicates if the input files need to be run through metaphlan.
            If not, just join profiles and count species.
        
    Requires:
        metaphlan v3.0.0+: A tool to profile the composition of microbial communities.
        humann v0.9.6+: A tool for functional profiling (only humann_join_tables is required).
        
    Returns:
        string: A file of the merged taxonomic profiles from all samples.
        list: A list of the sam files generated by metaphlan.
        
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
    sample_names=utilities.sample_names(input_files,input_extension)
    
    # get a list of metaphlan output files, one for each input file
    metaphlan_profile_tag="taxonomic_profile"
    main_folder=os.path.join("metaphlan","main")
    metaphlan_output_files_profile = utilities.name_files(sample_names, output_folder, subfolder=main_folder, tag=metaphlan_profile_tag, extension="tsv", create_folder=True)
    metaphlan_output_files_sam = utilities.name_files(sample_names, output_folder, subfolder=main_folder, tag="bowtie2", extension="sam")
    metaphlan_output_folder = os.path.dirname(metaphlan_output_files_profile[0])
    
    # determine the input file type based on the extension
    if input_extension in ["fasta","fasta.gz","fa","fa.gz"]:
        input_type="fasta"
    else:
        input_type="fastq"
    
    # run metaphlan on each of the kneaddata output files
    if not already_profiled:
        for sample, depend_fastq, target_profile, target_sam in zip(sample_names, input_files, metaphlan_output_files_profile, metaphlan_output_files_sam):
            workflow.add_task_gridable(
                "metaphlan [depends[0]] --input_type [args[2]] --output_file [targets[0]] --samout [targets[1]] --nproc [args[0]] --no_map --tmp_dir [args[1]]",
                depends=[depend_fastq,TrackedExecutable("metaphlan")],
                targets=[target_profile,target_sam],
                args=[threads,metaphlan_output_folder,input_type],
                time="3*60 if file_size('[depends[0]]') < 25 else 4*3*60", # 3 hours or more depending on input file size
                mem="12*1024 if file_size('[depends[0]]') < 25 else 4*12*1024", # 12 GB or more depending on input file size
                cores=threads, # time/mem based on 8 cores
                name=utilities.name_task(sample,"metaphlan"))
    else:
        # set the names of the already profiled outputs
        metaphlan_output_files_profile = input_files
        metaphlan_output_folder = os.path.dirname(input_files[0])
    
    # merge all of the metaphlan taxonomy tables
    metaphlan_merged_output = files.ShotGun.path("taxonomic_profile", output_folder)
    
    # run the join taxonomic profiles script to merge all of the metaphlan profiles
    workflow.add_task(
        "join_taxonomic_profiles.py --input [args[0]] --output [targets[0]] --file_name [args[1]]",
        depends=metaphlan_output_files_profile,
        targets=metaphlan_merged_output,
        args=[metaphlan_output_folder, metaphlan_profile_tag],
        name="metaphlan_join_taxonomic_profiles")
   
    # get the name for the file to write the species counts
    metaphlan_species_counts_file = files.ShotGun.path("species_counts",output_folder,create_folder=True)

    # create a file of species counts
    workflow.add_task(
    "count_features.py --input [depends[0]] --output [targets[0]] --include s__ --filter t__ --reduce-sample-name",
    depends=metaphlan_merged_output,
    targets=metaphlan_species_counts_file,
    name="metaphlan_count_species") 

    return metaphlan_merged_output, metaphlan_output_files_profile, metaphlan_output_files_sam

def merge_pairs(workflow,input_files,extension,pair_identifier,output_folder):
    """ Merge the paired files into a single file 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq (or fasta) files already run through quality control.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.  
        
    Requires:
        none
        
    Returns:
        list: A list of all of the merged files (the same as the input files if files are not paired)
        string: The extension for the merged files
    """
    
    # check for paired input files
    input_pair1, input_pair2 = utilities.paired_files(input_files, extension, pair_identifier)
    
    # set the default output extension to be the same as the input
    output_extension = extension
    
    if input_pair1:
        sample_names=utilities.sample_names(input_pair1,extension,pair_identifier)
        # determine the command based on the pair identifier
        if extension.endswith(".gz"):
            command="gunzip -c "
            output_extension=extension.replace(".gz","")
        else:
            command="cat "
        # determine the targets based on the output folder and extension
        merged_files = utilities.name_files(sample_names, output_folder, subfolder="input_merged", extension=output_extension, create_folder=True)  
        
        # add a task to merge each pair set
        for pair1, pair2, target in zip(input_pair1, input_pair2, merged_files):
            workflow.add_task_gridable(
                command+" [depends[0]] [depends[1]] > [targets[0]]",
                depends=[pair1,pair2],
                time=10,
                mem=1*1024,
                cores=1,
                targets=target)
    else:
        merged_files=input_files
        
    return merged_files, output_extension
          

def functional_profile(workflow,input_files,extension,output_folder,threads,taxonomic_profiles=None, remove_intermediate_output=None,
    options=None):
    """Functional profile for whole genome shotgun sequences
    
    This set of tasks performs functional profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control. Optionally the taxonomic
    profiles can be provided for the samples.
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq (or fasta) files already run through quality control.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for kneaddata to use.
        taxonomic_profiles (list): A set of taxonomic profiles, one per sample (optional).
        remove_intermediate_output (bool): Remove intermediate output files.
        
    Requires:
        humann v0.9.6+: A tool for functional profiling.
        
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
    sample_names=utilities.sample_names(input_files,extension)
    
    ### Step 1: Run humann on all input files ###

    # get a list of output files, one for each input file, with the humann output file names
    main_folder=os.path.join("humann","main")
    genefamiles = utilities.name_files(sample_names, output_folder, subfolder=main_folder, tag="genefamilies", extension="tsv", create_folder=True)
    pathabundance = utilities.name_files(sample_names, output_folder, subfolder=main_folder, tag="pathabundance", extension="tsv")
    pathcoverage = utilities.name_files(sample_names, output_folder, subfolder=main_folder, tag="pathcoverage", extension="tsv")

    # get the list of the log files that will be created, one for each input file, set to the same folder as the main output files
    log_files = utilities.name_files(sample_names, output_folder, subfolder=main_folder, extension="log")
    # get the name for the file of read and species counts created from the humann log outputs
    log_counts = files.ShotGun.path("humann_read_counts",output_folder,create_folder=True)

    humann_output_folder = os.path.dirname(genefamiles[0])
    
    # if taxonomic profiles are provided, add these to the targets and the command option
    if taxonomic_profiles:
        optional_profile_args=" --taxonomic-profile [depends[1]] "
        depends=list(zip(input_files,taxonomic_profiles))
    else:
        optional_profile_args=""
        depends=input_files
        
    # remove intermediate temp files, if set
    if remove_intermediate_output:
        optional_profile_args+=" --remove-temp-output "
   
    if not options is None:
        optional_profile_args+=" "+options+" "
 
    # create a task to run humann on each of the kneaddata output files
    for sample, depend_fastq, target_gene, target_path, target_coverage, target_log in zip(sample_names, depends, genefamiles, pathabundance, pathcoverage, log_files):
        workflow.add_task_gridable(
            "humann --input [depends[0]] --output [args[0]] --o-log [targets[3]] --threads [args[1]]"+optional_profile_args,
            depends=utilities.add_to_list(depend_fastq,TrackedExecutable("humann")),
            targets=[target_gene, target_path, target_coverage, target_log],
            args=[humann_output_folder, threads],
            time="2*24*60 if file_size('[depends[0]]') < 25 else 8*24*60", # 24 hours or more depending on file size
            mem="2*32*1024 if file_size('[depends[0]]') < 25 else 3*48*1024", # 32 GB or more depending on file size
            cores=threads,
            name=utilities.name_task(sample,"humann"))

    # create a task to get the read and species counts for each humann run from the log files
    workflow.add_task(
        "get_counts_from_humann_logs.py --input [args[0]] --output [targets[0]]",
        depends=log_files,
        targets=log_counts,
        args=humann_output_folder,
        name="humann_count_alignments_species")
    
    ### STEP #2: Regroup UniRef90 gene families to ecs ###
    
    # get a list of all output ec files
    ec_files = utilities.name_files(sample_names, output_folder, 
        subfolder=os.path.join("humann","regrouped"), tag="ecs", 
        extension="tsv", create_folder=True)
    
    # get ec files for all of the gene families files
    workflow.add_task_group_gridable(
        "humann_regroup_table --input [depends[0]] --output [targets[0]] --groups uniref90_level4ec",
        depends=genefamiles,
        targets=ec_files,
        time=10, # 10 minutes
        mem=5*1024, # 5 GB
        cores=1,
        name=map(lambda sample: utilities.name_task(sample,"humann_regroup_UniRef2EC"), sample_names))

    
    ### STEP #3: Merge gene families, ecs, and pathway abundance files

    # get a list of merged files for ec, gene families, and pathway abundance
    merged_genefamilies = files.ShotGun.path("genefamilies", output_folder, create_folder=True)
    merged_ecs = files.ShotGun.path("ecs", output_folder)
    merged_pathabundance = files.ShotGun.path("pathabundance", output_folder)
    
    # merge the ec, gene families, and pathway abundance files
    all_depends=[genefamiles, ec_files, pathabundance]
    all_targets=[merged_genefamilies, merged_ecs, merged_pathabundance]
    file_basenames=["genefamilies","ecs","pathabundance"]
    for depends, targets, basename in zip(all_depends, all_targets, file_basenames):
        workflow.add_task(
            "humann_join_tables --input [args[0]] --output [targets[0]] --file_name [args[1]]",
            depends=depends,
            targets=targets,
            args=[os.path.dirname(depends[0]),basename],
            name="humann_join_tables_"+basename)
    
    ### STEP #4: Normalize gene families, ecs, and pathway abundance to relative abundance (then merge files) ###
    
    # get a list of files for normalized ec, gene families, and pathway abundance
    relab_folder=os.path.join("humann","relab")
    norm_genefamily_files = utilities.name_files(genefamiles, output_folder, subfolder=os.path.join(relab_folder,"genes"), tag="relab", create_folder=True)
    norm_ec_files = utilities.name_files(ec_files, output_folder,  subfolder=os.path.join(relab_folder,"ecs"), tag="relab", create_folder=True)
    norm_pathabundance_files = utilities.name_files(pathabundance, output_folder,  subfolder=os.path.join(relab_folder,"pathways"), tag="relab", create_folder=True)
    
    # normalize the genefamily, ec, and pathabundance files
    # do not include special features (ie UNMAPPED, UNINTEGRATED, UNGROUPED) in norm computation
    renorm_task_names=len(norm_genefamily_files)*["genes"] + len(norm_ec_files)*["ecs"] + len(norm_pathabundance_files)*["pathways"]
    renorm_task_names=[utilities.name_task(sample,"humann_renorm_"+type+"_relab") for sample, type in zip(itertools.cycle(sample_names),renorm_task_names)]
    workflow.add_task_group_gridable(
        "humann_renorm_table --input [depends[0]] --output [targets[0]] --units relab --special n",
        depends=genefamiles + ec_files + pathabundance,
        targets=norm_genefamily_files + norm_ec_files + norm_pathabundance_files,
        time=15, # 15 minutes
        mem=5*1024, # 5 GB
        cores=1,
        name=renorm_task_names)

    
    # get a list of merged files for ec, gene families, and pathway abundance
    merged_genefamilies_relab = files.ShotGun.path("genefamilies_relab", output_folder)
    merged_ecs_relab = files.ShotGun.path("ecs_relab", output_folder)
    merged_pathabundance_relab = files.ShotGun.path("pathabundance_relab", output_folder)
    
    # merge the ec, gene families, and pathway abundance files
    all_depends=[norm_genefamily_files, norm_ec_files, norm_pathabundance_files]
    all_targets=[merged_genefamilies_relab, merged_ecs_relab, merged_pathabundance_relab]
    all_types=["genes_relab","ecs_relab","pathways_relab"]
    for depends, targets, input_type in zip(all_depends, all_targets, all_types):
        workflow.add_task(
            "humann_join_tables --input [args[0]] --output [targets[0]]",
            depends=depends,
            targets=targets,
            args=[os.path.dirname(depends[0])],
            name="humann_join_tables_"+input_type)

    # get feature counts for the ec, gene families, and pathways
    genefamilies_counts = files.ShotGun.path("genefamilies_relab_counts", output_folder)
    ecs_counts = files.ShotGun.path("ecs_relab_counts", output_folder)
    pathabundance_counts = files.ShotGun.path("pathabundance_relab_counts", output_folder)
    workflow.add_task_group(
        "count_features.py --input [depends[0]] --output [targets[0]] --reduce-sample-name --ignore-un-features --ignore-stratification",
        depends=[merged_genefamilies_relab, merged_ecs_relab, merged_pathabundance_relab],
        targets=[genefamilies_counts, ecs_counts, pathabundance_counts],
        name=["humann_count_features_genes","humann_count_features_ecs","humann_count_features_pathways"])
    
    # merge the feature counts into a single file
    all_feature_counts = files.ShotGun.path("feature_counts", output_folder)
    workflow.add_task(
        "humann_join_tables --input [args[0]] --output [targets[0]] --file_name _relab_counts.tsv",
        depends=[genefamilies_counts, ecs_counts, pathabundance_counts],
        targets=all_feature_counts,
        args=[os.path.dirname(genefamilies_counts)],
        name="humann_merge_feature_counts")

        
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

def strainphlan(task,threads,clade_number,clade_list,reference_folder,marker_folder,options):
    """ Run strainphlan for the specific clade
    
    Args:
        task: (anadama2.task): An instance of the task class.
        threads: (int): Total threads to use for task.
        clade_number: (int): The number of clade to run.
        clade_list: (string): The path to the clade list.
        reference_folder (string): The folder containing the reference files.
        marker_folder (string): The folder containing the marker files.
        options (string): Options to apply when running strainphlan.

    Requires:
        StrainPhlAn: A tool for strain profiling.    
        
    """        
    
    # find the name of the clade in the list
    with open(clade_list) as file_handle:
        clades=[line.strip().split(" ")[0] for line in filter(lambda line: line.startswith("s__") or line.startswith("g__"), file_handle.readlines())]
        try:
            profile_clade=clades[clade_number]
        except IndexError:
            profile_clade=None
            
    if profile_clade:
        command = "strainphlan --samples [args[0]]/*/*.pkl --output_dir [args[1]] "+\
            "--clade [args[2]] --nprocs [args[3]] "+options
            
        # add the marker files to the command
        all_marker_file=os.path.join(marker_folder,"all_markers.fasta")
        marker_file=os.path.join(marker_folder,profile_clade+".fna")
        
        # generate the marker file if it does not already exist
        if not os.path.isfile(marker_file):
            # get the pkl file relative to the strainphlan install
            try:
                import metaphlan
                strainphlan_pkl=os.path.join(os.path.dirname(metaphlan.__file__),"metaphlan_databases","mpa_v30_CHOCOPhlAn_201901.pkl")
            except subprocess.CalledProcessError:
                raise EnvironmentError("Unable to find strainphlan install.")
            
            marker_command="extract_markers.py --database [depends[0]] "+\
                "--clade [args[0]] --output_dir [args[1]]"
                
            # create the marker file in the output folder
            marker_file=os.path.join(os.path.dirname(task.targets[0].name),profile_clade+".fna")
                
            return_code = utilities.run_task(marker_command,depends=[strainphlan_pkl], 
                targets=[marker_file],args=[profile_clade,os.path.dirname(marker_file)])
        
        # check that the marker file exists
        if not os.path.isfile(marker_file):
            raise EnvironmentError("Unable to find StrainPhlAn markers for clade "+ profile_clade)

        # add the marker file location (default or per workflow) to the command
        command += " --clade_markers "+marker_file       
 
        # get the list of reference genomes
        genomes=set()
        with open(data.get_file("strainphlan_species_gcf.tsv")) as file_handle:
            for line in file_handle:
                if line.startswith(profile_clade):
                    genomes.add(os.path.join(reference_folder,line.rstrip().split("\t")[-1]+".fna.bz2"))
        
        # only use those references found in the folder
        genomes = list(filter(os.path.isfile,genomes))
    
        # if more than six references are found, limit total used
        # this resolves e. coli which has so many references it can exceed command line length
        if len(genomes) > 6:
            genomes = genomes[:6]    
        
        # add the reference genome files to the command, if any are found
        if len(genomes):
            command += " --references " + " --references ".join(genomes)
        
        # write the output to the log
        command += " > [targets[0]] && touch [targets[1]] && if [ -f [args[1]]/RAxML_bestTree.[args[2]].tree ]; then cp [args[1]]/RAxML_bestTree.[args[2]].tree [targets[1]]; fi"
        
    else:
        # there is not a clade of this number, create an empty output file
        command = "touch [targets[0]] && touch [targets[1]]"
        
    # run the task
    return_code = utilities.run_task(command, depends=task.depends, targets=task.targets, 
        args=[os.path.abspath(os.path.join(os.path.dirname(task.depends[0].name),"..")),os.path.dirname(task.targets[0].name),profile_clade,threads])
    

def strain_profile(workflow,sam_files,output_folder,threads,reference_folder,marker_folder,abundance_file,options="",max_species=20,strain_list=""):
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
        abundance_file (string): The merged taxonomic abundance file to be used to select
            the top strains by average abundance.
        options (string): Options to apply when running strainphlan.
        max_species (int): The maximum number of species to profile.
        strain_list (string): The path to a file with the list of strains to profile
    Requires:
        StrainPhlAn: A tool for strain profiling.
        
    Returns:
        None
        
    """    

    ### STEP #1: Identify markers for each of the samples
    # name the marker files based on the sam files
    strainphlan_markers_temp = utilities.name_files(sam_files, output_folder, subfolder="strainphlan", extension="pkl", create_folder=True)
    # place each in its own output folder to allow for unique temp output folders for each run
    strainphlan_markers=[]
    for filename in strainphlan_markers_temp:
        sample_name = os.path.basename(filename).split(".")[0]
        strainphlan_markers.append(os.path.join(os.path.dirname(filename), sample_name, os.path.basename(filename)))

    # create a marker file from each sam file, require min depth
    for sam, markers in zip(sam_files, strainphlan_markers):
        sample_name=os.path.basename(sam).replace("_bowtie2.sam","")
        workflow.add_task_gridable(
            "mkdir -p [args[0]] && sample2markers.py --input [depends[0]] --input_format sam --output_dir [args[0]] --nprocs [args[1]]",
            depends=sam,
            targets=markers,
            args=[os.path.dirname(markers),threads],
            time=30, # 30 minutes
            mem=5*1024, # 5 GB
            cores=threads,
            name=utilities.name_task(sample_name,"strainphlan_sample2markers"))

    ### STEP #2: Find the top species by average abundance
    clade_list = utilities.name_files("clades_list.txt", output_folder, subfolder="strainphlan")
    
    workflow.add_task(
        "strainphlan --samples [args[0]]/*/*.pkl --output_dir [args[0]] --print_clades_only > [targets[0]] "+options,
        depends=strainphlan_markers,
        targets=clade_list,
        args=os.path.abspath(os.path.join(os.path.dirname(strainphlan_markers[0]),"..")),
        name="strainphlan_print_clades")
    
    # order the clades list by average abundance if list not provided
    if not strain_list:
        ordered_clade_list = utilities.name_files("clades_list_order_by_average_abundance.txt", output_folder, subfolder="strainphlan")
        workflow.add_task(
            utilities.partial_function(utilities.order_clade_list,clade_list=clade_list,abundance_file=abundance_file,output_file=ordered_clade_list),
            depends=[clade_list, abundance_file],
            targets=ordered_clade_list)
    else:
        # use the list provided for the strains to profile
        ordered_clade_list = strain_list
    
    ### STEP #3: Run strainphlan on the top set of clades identified
    clade_logs = utilities.name_files(map(str,range(max_species)), output_folder, tag="clade", subfolder="strainphlan", extension="log")
    clade_tree = utilities.name_files(map(str,range(max_species)), output_folder, tag="clade", subfolder="strainphlan", extension="tree")
    for clade_number in range(max_species):
        workflow.add_task_gridable(
            utilities.partial_function(strainphlan,threads=threads,clade_number=clade_number,
                clade_list=ordered_clade_list,reference_folder=os.path.abspath(reference_folder),marker_folder=os.path.abspath(marker_folder),
                options=options),
            depends=strainphlan_markers+[ordered_clade_list],
            targets=[clade_logs[clade_number-1],clade_tree[clade_number-1]], 
            time=6*60, # 6 hours
            mem=50*1024, # 50 GB
            cores=threads,
            name="strainphlan_clade_"+str(clade_number))

def get_panphlan_species_name(abundance_file, species_number, panphlan_db):
    """ Get the panphlan species name based on the abundance file """

    species_ranked = utilities.rank_species_average_abundance(abundance_file)
    
    try:
        selected_species_info = species_ranked[species_number].split("_")
        selected_species = selected_species_info[-2].lower()[0]+selected_species_info[-1].lower()
    except IndexError:
        selected_species = None
        selected_species_with_version = None

    if selected_species:
        # get the version number based on the files found
        possible_dbs=[]
        for file in os.listdir(panphlan_db):
            if file.startswith("panphlan_"+selected_species) and file.endswith(BOWTIE2_EXTENSION):
                possible_dbs.append(file)

        if not possible_dbs:
            raise EnvironmentError("Unable to find panphlan database "+selected_species+" in folder "+panphlan_db)

        selected_species_with_version = sorted(possible_dbs)[-1].split(".")[0]

    return selected_species_with_version

def panphlan_map(task,species_number,threads,panphlan_db,output_file):
    """ Run the panphlan map step for the input file and the database selected """

    # get the species for this task
    selected_species = get_panphlan_species_name(task.depends[0].name, species_number, panphlan_db)

    if selected_species:
        # get a single folder for this species database
        species_db = os.path.join(panphlan_db, "panphlan_"+selected_species+BOWTIE2_EXTENSION)
    
        # run the task
        return_code = utilities.run_task(
            "panphlan_map.py -c [args[0]] -i [depends[1]] -o [args[1]] --i_bowtie2_indexes [args[2]] --tmp [args[3]] --nproc [args[4]] --verbose > [targets[0]]", 
            depends=task.depends+[species_db], targets=task.targets, 
            args=[selected_species,output_file,panphlan_db,os.path.dirname(task.targets[0].name),threads])
    else:
        # there is not a clade of this number, create an empty output file
        utilities.run_task("touch [targets[0]]", targets=task.targets)


def panphlan_profile(task,species_number,panphlan_db):
    """ Run panphlan profile on the set of input files and species """

    # get the species for this task
    selected_species = get_panphlan_species_name(task.depends[0].name, species_number, panphlan_db)

    if selected_species:
        # get a single folder for this species database
        species_db = os.path.join(panphlan_db, "panphlan_"+selected_species+BOWTIE2_EXTENSION)

        # get the name of the species gene file target
        output_folder = os.path.dirname(task.depends[1].name)
        gene_target=os.path.join(output_folder,selected_species+"_gene_presence_absence.tsv")
        # run the task
        return_code = utilities.run_task(
            "panphlan_profile.py -c [args[0]] -i [args[1]] --o_dna [targets[0]] --add_strains --i_bowtie2_indexes [args[2]] --verbose > [targets[1]]",
            depends=task.depends+[species_db], targets=[gene_target]+task.targets,
            args=[selected_species,output_folder,panphlan_db])
    else:
        # there is not a clade of this number, create an empty output file
        utilities.run_task("touch [targets[0]]", targets=task.targets)


def strain_gene_profile(workflow,qc_files,abundance_file,output,threads,panphlan_db,max_species):
    """Strain profile, gene-based for whole genome shotgun sequences
   
    This set of tasks performs strain profiling on whole genome shotgun
    input files. For paired-end files, first merge and provide a single file per sample.
    Input files should first be run through quality control.
   
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        qc_files (list): A list of paths to qc files generated by Kneaddata.
        abundance_file (string): The merged taxonomic abundance file to be used to select
            the top strains by average abundance.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores to use.
        panphlan_db (string): The folder containing the database files.
        max_species (int): The maximum number of species to profile.
    Requires:
        PanPhlAn: A tool for strain profiling.

    Returns:
        None
    """

    ### STEP #1: Run panphlan map on all of the samples for each of the top species
    out_files = workflow.name_output_files(name=qc_files, tag="panphlan_map", extension="csv.bz2")
    log_files = workflow.name_output_files(name=qc_files, tag="panphlan_map", extension="log")
    for clade_number in range(max_species):
        map_out = [os.path.join(output,"panphlan",str(clade_number),os.path.split(file)[-1].replace("panphlan_map","panphlan_map_clade"+str(clade_number))) for file in out_files]
        log_out = [os.path.join(output,"panphlan",str(clade_number),os.path.split(file)[-1].replace("panphlan_map","panphlan_map_clade"+str(clade_number))) for file in log_files]

        # create the output folder
        subfolder=os.path.dirname(map_out[0])
        if not os.path.isdir(subfolder):
            os.makedirs(subfolder)

        for input_file, output_file, sample_map_out, sample_log_out in zip(qc_files, out_files, map_out, log_out):
            file_sample = os.path.split(output_file)[-1].split(".")[0]
            workflow.add_task_gridable(
                utilities.partial_function(panphlan_map,species_number=clade_number,
                    threads=threads, panphlan_db=os.path.abspath(panphlan_db),output_file=sample_map_out),
                depends=[abundance_file, input_file],
                targets=[sample_log_out],
                time=2*60, # 2 hours
                mem=5*1024, # 5 GB
                cores=threads,
                name=file_sample+"_clade_"+str(clade_number))

        ### STEP #2: Run panphlan profile on each clade
        profile_output_log = os.path.join(output,"panphlan","clade_"+str(clade_number)+".log")
        workflow.add_task_gridable(
            utilities.partial_function(panphlan_profile,species_number=clade_number,
                panphlan_db=os.path.abspath(panphlan_db)),
            depends=[abundance_file]+log_out,
            targets=profile_output_log,
            mem=5*1024, # 5 GB
            time=60, # 1 hour
            cores=1,
            name="panphlan_profile_clade_"+str(clade_number))

def megahit(workflow, input_files, extension, output_folder, threads, remove_intermediate_output=True,
    additional_options=None, interleaved=False):
    """Run MEGAHIT.
    
    This set of tasks will run MEGAHIT on the input files provided to produce contigs via metagenomic
    assembly. 
   
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to megahit. A list of lists is 
            provided here containing the reads to assemble and any singleton reads post-QC.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for kneaddata to use.
        paired (boolean): Whether or not the input files are paired-end
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        remove_intermediate_output (bool): Remove intermediate output files.
        additional_options (string): Any additional options to pass onto megahit.
        interleaved (bool): Whether or not sequence files passed in are interleaved
        
    Requires:
        MEGAHIT v1.1.3: A tool to perform metagenomic assembly on 
            metagenomic sequencing data.
        
    Returns:
        list: A list of the assembled contigs created by MEGAHIT.
    """
    megahit_contigs = []
    sample_names = utilities.sample_names(input_files[0], extension)

    time_equation="8*60 if file_size('[depends[0]]') < 10 else 10*60"
    mem_equation="2*12*1024 if file_size('[depends[0]]') < 10 else 4*12*1024"
        
    assembly_dir = os.path.join(output_folder, "assembly", "main")
    depends = []
    megahit_template = "megahit %s -t [args[0]] -m 0.99 -f -o [args[3]] --out-prefix [args[1]] [args[2]]"

    workflow.add_task('mkdir -p [targets[0]]',
                      depends=[output_folder],
                      targets=[assembly_dir])

    for (sample_name, input_reads, orphan_reads) in zip(sample_names, input_files[0], input_files[1]):
        sample_name = sample_name.replace("_sorted_final","")

        megahit_contig_dir = os.path.join(assembly_dir, sample_name)
        intermediate_dir = os.path.join(megahit_contig_dir, 'intermediate_contigs')
        megahit_contig = os.path.join(megahit_contig_dir, '%s.contigs.fa' % sample_name)
        completed_file = os.path.join(megahit_contig_dir, 'done')

        if interleaved:
            megahit_cmd = megahit_template % "--12 [depends[0]] -r [depends[1]]"
            depends = [input_reads, orphan_reads]
        else:
            megahit_cmd = megahit_template % "-r [depends[0]]"
            depends = [input_reads]

        workflow.add_task_gridable(megahit_cmd,
                                   depends=depends,
                                   targets=[megahit_contig, completed_file],
                                   args=[threads, sample_name, additional_options, megahit_contig_dir],
                                   cores=threads,
                                   mem=mem_equation,
                                   time=time_equation)

        if remove_intermediate_output:
            workflow.add_task('rm -rf [depends[0]]',
                              depends=[intermediate_dir])

        megahit_contigs.append(megahit_contig)
    
    return megahit_contigs


def assemble(workflow, input_files, extension, output_folder, threads, pair_identifier=None,
    remove_intermediate_output=None, additional_options=None, interleaved=False):
    """Metagenomic assembly for whole genome shotgun sequences.

    This set of tasks performs metagenomic assembly on whole genome shotgun input files in either 
    paired-end or single-end fastq format. 
        
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list) A list of lists of fastq files that have been run through quality control and 
            produced cleaned sequences and any orphan reads.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores to use during assembly.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        remove_intermediate_output (bool): Remove intermediate output files.
        additional_options (string): Any additional options to pass onto megahit.
        interleaved (bool): Whether or not sequence files passed in are interleaved
        
    Requires:
        MEGAHIT v1.1.3: A tool to perform metagenomic assembly on 
            metagenomic sequencing data.
        
    Returns:
        list: A list of contig files for all samples
        list: A list of contigs that pass a length filter (DEFAULT: 500bp)
        
    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        (assembled_contigs, assembled_contigs_filtered) = shotgun.assemble(["cleaned_seqsA.fastq"])
        
        # run the workflow
        workflow.go()
    """
    # Start out by sorting our sequences to make sure they are in the proper order (if interleaved)
    sample_names = utilities.sample_names(input_files, extension)
    product_extension = "fastq"

    if interleaved:
        sort_dir = os.path.join(output_folder, "sort", "main")
        sorted_sequences = utilities.name_files(sample_names, sort_dir, tag="sorted", extension=product_extension, create_folder=True)

        workflow.add_task_group_gridable(utilities.sort_fastq_file,
                                        depends = input_files,
                                        targets = sorted_sequences,
                                        time = "2*60 if file_size('[depends[0]]') < 6 else 4*60",
                                        mem="12*1024 if file_size('[depends[0]]') < 6 else 2*12*1024",
                                        cores = threads)

        orphans_dir = os.path.join(output_folder, "extract_orphans", "main")
        orphan_seqs_files = utilities.name_files(sorted_sequences, orphans_dir, tag="orphans", extension=product_extension, create_folder=True)
        balanced_seqs_files = utilities.name_files(sorted_sequences, orphans_dir, tag="final", extension=product_extension)

        for (sorted_seq, balanced_seq, orphan_seq) in zip(sorted_sequences, balanced_seqs_files, orphan_seqs_files):
            workflow.add_task_gridable(utilities.extract_orphan_reads,
                                       depends = sorted_seq,
                                       targets = [balanced_seq, orphan_seq],
                                       cores = threads,
                                       mem = "4*1024 if file_size('[depends[0]]') < 10 else 2*4*1024",
                                       time = "1*60 if file_size('[depends[0]]') < 10 else 2*60")

        input_files = [balanced_seqs_files, orphan_seqs_files]
    else:
        input_files = [input_files, [None] * len(input_files)]

    assembled_contig_files = megahit(workflow, input_files, product_extension, output_folder, threads,
                                     remove_intermediate_output, additional_options, interleaved)

    return assembled_contig_files


def prodigal(workflow, contigs, output_folder, threads):
    """Open reading frame prediction for genes identified in the supplied contigs.

    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        contigs (list): A list of FASTQ files containing contigs to predict genes from.
        output_folder (string): The path of the output folder.
        threads (int) The number of threads/cores to use during sorting.

    Requires:
        prodigal v2.6.3+: Fast, reliable protein-coding gene prediction for prokaryotic genomes.

    Returns: 
        list: A list of GFF3 files containing predicted genes and associated annotations.
        list: A list of predicted gene nucleotide coding sequences.
        list: A list of predicted gene amino acid coding sequences.

    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        (assembled_contigs, assembled_contigs_filtered) = shotgun.assemble(["cleaned_seqsA.fastq"])
        shotgun.annoate(workflow, assembled_contigs, "/home/output_folder", 8)

        # run the workflow
        workflow.go()
    """
    sample_names = utilities.sample_names(contigs, ".fa")

    time_equation="2*60 if file_size('[depends[0]]') < 10 else 2*2*60"
    mem_equation="2*12*1024 if file_size('[depends[0]]') < 10 else 4*12*1024"

    annotation_dir = os.path.join(output_folder, "annotation", "main")
    gff3_files = utilities.name_files(sample_names, annotation_dir, create_folder=True, extension="gff")
    nuc_cds_files = utilities.name_files(sample_names, annotation_dir, extension="fna")
    aa_cds_files = utilities.name_files(sample_names, annotation_dir, extension="faa")

    for (input_contig, gff3_file, nuc_cds_file, aa_cds_file) in zip(contigs, gff3_files, nuc_cds_files, aa_cds_files):
        workflow.add_task_gridable("prodigal -q -p meta -i [depends[0]] -f gff -o [targets[0]] -d [targets[1]] -a [targets[2]]",
                                   depends=[input_contig, TrackedDirectory(annotation_dir)],
                                   targets=[gff3_file, nuc_cds_file, aa_cds_file],
                                   cores=threads,
                                   mem=mem_equation,
                                   time=time_equation)

    return (gff3_files, nuc_cds_files, aa_cds_files)


def prokka(workflow, contigs, output_folder, threads):
    """Whole genome annotation for the supplied contigs.

    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        contigs (list): A list of FASTQ files containing contigs to predict genes from.
        output_folder (string): The path of the output folder.
        threads (int) The number of threads/cores to use during sorting.

    Requires:
        PROKKA v1.13.3+: Rapid prokaryotic genome annotation

    Returns: 
        list: A list of GFF3 files containing predicted genes and associated annotations.
        list: A list of predicted gene nucleotide coding sequences.
        list: A list of predicted gene amino acid coding sequences.
        list: A list of table files containing all features predicted [locus_tag,ftype,len_bp,gene,EC_number,COG,product]

    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun
        
        # create an anadama2 workflow instance
        workflow=Workflow()
        
        # add quality control tasks for the fastq files
        (assembled_contigs, assembled_contigs_filtered) = shotgun.assemble(["cleaned_seqsA.fastq"])
        shotgun.annoate(workflow, assembled_contigs, "/home/output_folder", 8)

        # run the workflow
        workflow.go()
    """
    sample_names = utilities.sample_names(contigs, ".contigs.fa")

    time_equation="20*60 if file_size('[depends[0]]') < 10 else 30*60"
    mem_equation="2*12*1024 if file_size('[depends[0]]') < 10 else 4*12*1024"

    annotation_dir = os.path.join(output_folder, "annotation", "main")
    gff3_files = utilities.name_files(sample_names, annotation_dir, create_folder=True, extension="gff")
    nuc_cds_files = utilities.name_files(sample_names, annotation_dir, extension="fna")
    aa_cds_files = utilities.name_files(sample_names, annotation_dir, extension="faa")
    feature_tables = utilities.name_files(sample_names, annotation_dir, extension="tsv")

    for (sample_name, input_contig, gff3_file, nuc_cds_file, aa_cds_file, feature_table) in zip(sample_names, contigs, 
                                                                                                gff3_files, nuc_cds_files, aa_cds_files, 
                                                                                                feature_tables):
        workflow.add_task_gridable("prokka --force --outdir [args[2]] --prefix [args[0]] [depends[0]] --cpus [args[1]]",
                                   depends=[input_contig],
                                   targets=[gff3_file, nuc_cds_file, aa_cds_file, feature_table],
                                   args=[sample_name, threads, annotation_dir],
                                   cores=threads,
                                   mem=mem_equation,
                                   time=time_equation)

    return (gff3_files, nuc_cds_files, aa_cds_files, feature_tables)


def annotate(workflow, contigs, output_folder, threads):
    """Annotates the provided contig files.

    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        contigs (list): A list of paths to fastq files to predict genes from.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads/cores for kneaddata to use.

    Requires:
        prodigal v2.6.3+: Fast, reliable protein-coding gene prediction for prokaryotic genomes.

    Returns:
        list: A list of GFF3 files containing predicted genes and associated annotations.
        list: A list of predicted gene nucleotide coding sequences.
        list: A list of predicted gene amino acid coding sequenecs.

    Example:
        from anadama2 import Workflow
        from biobakery_workflows.tasks import shotgun

        # create an anadama2 workflow instance
        workflow=Workflow()

        out_folder = "/home/gene_predictions"
        contigs = ['/home/contigsA.fasta', '/home/contigsB.fasta']

        # predict genes
        (gff3_files, nuc_cds_seqs, aa_cds_seqs) = shotgun.predict_genes(workflow, contigs, out_folder, 4)
            
        # run the workflow
        workflow.go()
    """
    (gff3_files, nuc_cds_seqs, aa_cds_seqs, feature_tables) = prokka(workflow, contigs, output_folder, threads)

    return (gff3_files, nuc_cds_seqs, aa_cds_seqs, feature_tables)
