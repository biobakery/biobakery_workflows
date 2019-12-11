"""
bioBakery Workflows: tasks.sixteen_s module
A collection of tasks for workflows with shotgun sequences

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

from anadama2.tracked import TrackedExecutable

from biobakery_workflows import utilities
from biobakery_workflows import files

    
def quality_control(workflow, method, fastq_file, output_folder, threads, maxee, trunc_len):
    """ Create a quality report, filter fastq, and then truncate fasta files
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        fastq_file (string): The path to the fastq file.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        maxee (int): The maxee value to use for filtering.
        trunc_len (int): The value to use for max length.
        
    Requires:
        usearch: tools for sequence analysis.
        
    Returns:
        list: Paths to the truncated files (full fasta, and filtered fasta) and full file
        
    """
    
    # generate a quality report that can be used again for filtering
    qc_report = quality_report(workflow, method, fastq_file, output_folder, threads)
        
    # filter the fastq file with the maxee scores
    filtered_truncated_fasta, fasta = filter_fastq(workflow, method, fastq_file, output_folder, threads, maxee, trunc_len)
    
    # truncate reads based on length
    truncated_fastas = truncate(workflow, method, [fasta], output_folder, threads, trunc_len)
    
    return filtered_truncated_fasta, truncated_fastas[0], fasta

def merge_samples_and_rename(workflow, method, input_files, extension, output_folder, pair_identifier, threads):
    """ Merge the files, first if pairs, then rename sequence ids to match sample id
         Then merge all files into a single fastq file

    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set.
        threads (int): The number of threads for each task.
        
    Requires:
        usearch: tools for sequence analysis.
        
    Returns:
        string: The path to the merged fastq file for all samples

    """
    
    # merge the files, if pairs, and then rename sequence ids to match sample ids
    renamed_files = merge_pairs_and_rename(workflow, method, input_files, extension, output_folder, pair_identifier, threads)
    
    # merge the renamed files into a single fastq file
    all_samples_fastq = merge_fastq(workflow, renamed_files, output_folder)
    
    return all_samples_fastq


def merge_pairs(task, method, threads=1):
    """ Merge the pair files, allowing for empty input files 
    
        Args:
            task: (anadama2.task): An instance of the task class.
            method (string): tools for sequence analysis - vsearch or usearch(default)
            threads: (int): Total threads to use for usearch task.
        
        Requires:
            usearch or vsearch
        """
    
    if os.path.getsize(task.depends[0].name) > 0:
        # since the input files are not empty, run usearch merge pairs
        # using vsearch
        if method == 'vsearch':
            command = "export OMP_NUM_THREADS=[args[0]]; " + \
                      "vsearch -fastq_mergepairs [depends[0]] -reverse [depends[1]]  -fastqout [targets[0]] -fastqout_notmerged_fwd [targets[1]] -threads [args[0]]"
        # using usearch
        else:
            command="export OMP_NUM_THREADS=[args[0]]; " +\
            "usearch -fastq_mergepairs [depends[0]] -reverse [depends[1]]  -fastqout [targets[0]] -fastqout_notmerged_fwd [targets[1]] -threads [args[0]]"
        
    else:
        # the input files are empty, create empty output files
        command="touch [targets[0]] [targets[1]]"
        
    # run the task
    return_code = utilities.run_task(command, depends=task.depends, targets=task.targets, args=threads)
    

def merge_pairs_and_rename(workflow, method, input_files, extension, output_folder, pair_identifier, threads):
    """ Merge the files if pairs and rename sequence ids to match sample id
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis, usearch default or vsearch
        input_files (list): A list of paths to fastq files.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set.
        threads (int): The number of threads for each task.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: A list of the renamed files.
        
    """   
    
    pair1, pair2=utilities.paired_files(input_files, extension, pair_identifier)
    
    if pair1 and pair2:
        # paired input files were found

        # if the files are gzipped, first decompress as fastq_mergepairs will take in fastq.gz but the output will not be correctly formatted
        if pair1[0].endswith(".gz"):
            # get the names of the decompressed output files
            decompressed_pair1=utilities.name_files([os.path.basename(file).replace(".gz","") for file in pair1],output_folder,subfolder="merged_renamed")
            # get the names of the decompressed output files
            decompressed_pair2=utilities.name_files([os.path.basename(file).replace(".gz","") for file in pair2],output_folder,subfolder="merged_renamed")
            
            # add tasks to decompress the files
            workflow.add_task_group(
                "gunzip -c [depends[0]] > [targets[0]]",
                depends=pair1+pair2,
                targets=decompressed_pair1+decompressed_pair2)
            
            # the pair files to be used for the remaining tasks are those that are decompressed
            pair1=decompressed_pair1
            pair2=decompressed_pair2


        # get the sample names from the input file names
        sample_names=[os.path.basename(file).replace(pair_identifier+".fastq","") for file in pair1]
        
        # get the names of the output files
        stitched_files=utilities.name_files(sample_names,output_folder,subfolder="merged_renamed",tag="stitched",extension="fastq",create_folder=True)
        unjoined_files=utilities.name_files(sample_names,output_folder,subfolder="merged_renamed",tag="unjoined",extension="fastq")
        
        # run usearch to merge pairs, if input files are non-empty
        for read1, read2, stitched_output, unjoined_output in zip(pair1,pair2,stitched_files,unjoined_files):
            if method == 'vsearch':
                workflow.add_task(
                    utilities.partial_function(merge_pairs, method="vsearch", threads=threads),
                    depends=[read1, read2, TrackedExecutable("vsearch")],
                    targets=[stitched_output, unjoined_output],
                    name="vsearch_fastq_mergepairs")
            else:
                workflow.add_task(
                    utilities.partial_function(merge_pairs,method="userach", threads=threads),
                    depends=[read1, read2, TrackedExecutable("usearch")],
                    targets=[stitched_output, unjoined_output],
                    name="usearch_fastq_mergepairs")
        
        # merge the stitched and unjoined from the prior step
        renamed_files=utilities.name_files(sample_names,output_folder,subfolder="merged_renamed",tag="renamed",extension="fastq")
        workflow.add_task_group(
            "merge_and_rename_fastq.py [depends[0]] [depends[1]] _stitched [targets[0]]",
            depends=zip(stitched_files,unjoined_files),
            targets=renamed_files)
    
    else:
        # these files are not pairs and do not need to be merged
        # rename the files
        renamed_files=utilities.name_files(input_files,output_folder,subfolder="merged_renamed",tag="renamed",extension="fastq",create_folder=True)
        workflow.add_task_group(
            "merge_and_rename_fastq.py [depends[0]] '' '' [targets[0]]",
            depends=input_files,
            targets=renamed_files)
        
    return renamed_files

def merge_fastq(workflow, input_files, output_folder):
    """ Merge all of the fastq files into a single fastq file
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files.
        output_folder (string): The path of the output folder.
        
    Requires:
        None
        
    Returns:
        string: A path to the merged file
        
    """       
        
    # get the name of the final merged fastq file
    all_samples_fastq = utilities.name_files("all_samples_concatenated.fastq",output_folder)

    workflow.add_task(
        "merge_fastq.py [args[0]] _renamed.fastq [targets[0]]",
        depends=input_files,
        targets=all_samples_fastq,
        args=os.path.dirname(input_files[0]))
    
    return all_samples_fastq


def quality_report(workflow, method, fastq_file, output_folder, threads, qmax=45):

    """ Generate a qc report from the fastq file of all samples
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearch(default) or vsearch
        fastq_file (string): The path to the fastq file.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        qmax (int): Max qvalue increased from the default of 43 to allow for Ion Torrent data
    Requires:
        usearch or vsearch
        
    Returns:
        string: A path to the qc report file
        
    """       
        
    # get the name of the final merged fastq file
    qc_file = files.SixteenS.path("eestats2", output_folder)
    if method == 'vsearch':
        workflow.add_task(
            "export OMP_NUM_THREADS=[args[0]]; " + \
            "vsearch -fastq_eestats2 [depends[0]] -output [targets[0]] -threads [args[0]]",
            depends=[fastq_file, TrackedExecutable("vsearch")],
            targets=qc_file,
            args=threads,
            name="vsearch_fastq_eestats2")
    else:
        workflow.add_task(
            "export OMP_NUM_THREADS=[args[0]]; "+\
            "usearch -fastq_eestats2 [depends[0]] -output [targets[0]] -threads [args[0]] -fastq_qmax [args[1]]",
            depends=[fastq_file,TrackedExecutable("usearch")],
            targets=qc_file,
            args=[threads,qmax],
            name="usearch_fastq_eestats2")
    
    return qc_file
  
  
def filter_fastq(workflow, method, fastq_file, output_folder, threads, maxee, trunc_len, qmax=45):
    """ Filter the fastq files using the maxee value
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        fastq_file (string): The path to the fastq file.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        maxee (int): The maxee value to use for filtering.
        trunc_len (int): The value to use for max length.
        qmax (int): Max qvalue increased from the default of 43 to allow for Ion Torrent data
    Requires:
        usearch or vsearch
        
    Returns:
        string: A path to the filtered fasta file
        string: A path to the full fasta file
        
    """  
  
    # get the name of the final merged fastq file
    fasta_filtered_file = utilities.name_files("all_samples_concatenated_filtered.fasta",output_folder)
    fasta_discarded_file = utilities.name_files("all_samples_concatenated_discarded.fasta",output_folder)
    if method == "vsearch":
        workflow.add_task(
            "export OMP_NUM_THREADS=[args[0]]; " + \
            "vsearch -fastq_filter [depends[0]] -fastq_maxee [args[1]] -fastaout [targets[0]] -threads [args[0]] -fastaout_discarded [targets[1]] -fastq_trunclen [args[2]]",
            depends=[fastq_file, TrackedExecutable("vsearch")],
            targets=[fasta_filtered_file, fasta_discarded_file],
            args=[threads, maxee, trunc_len],
            name="vsearch_fastq_filter")
    else:
        workflow.add_task(
            "export OMP_NUM_THREADS=[args[0]]; "+\
            "usearch -fastq_filter [depends[0]] -fastq_maxee [args[1]] -fastaout [targets[0]] -threads [args[0]] -fastaout_discarded [targets[1]] -fastq_trunclen [args[2]] -fastq_qmax [args[3]]",
            depends=[fastq_file,TrackedExecutable("usearch")],
            targets=[fasta_filtered_file, fasta_discarded_file],
            args=[threads, maxee, trunc_len, qmax],
            name="usearch_fastq_filter")
    
    # create a fasta file of all reads (included the discarded
    fasta_file = utilities.name_files("all_samples_concatenated.fasta",output_folder)
    workflow.add_task(
        "cat [depends[0]] [depends[1]] > [targets[0]]",
        depends=[fasta_filtered_file, fasta_discarded_file],
        targets=fasta_file)
    
    return fasta_filtered_file, fasta_file


def truncate(workflow, method, input_files, output_folder, threads, trunc_len):
    """ Truncate the fasta sequences by length
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        input_files (list): A list of paths to fastq files.
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        trunc_len (int): The value to use for max length.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: Paths to the truncated files
        
    """  

    # get the name of the output files
    output_files = utilities.name_files(input_files, output_folder, tag="truncated")
    if method == "vsearch":
        workflow.add_task_group(
            "vsearch --fastx_filter [depends[0]]  --fastq_trunclen [args[0]]  -fastaout [targets[0]]",
            depends=input_files,
            targets=output_files,
            args=trunc_len,
            name="vsearch_fastx_truncate")
    else:
        workflow.add_task_group(
            "usearch -fastx_truncate [depends[0]] -trunclen [args[0]] -fastaout [targets[0]]",
            depends=input_files,
            targets=output_files,
            args=trunc_len,
            name="usearch_fastx_truncate")

    return output_files
  
def pick_otus(workflow, method, fasta_file, reference_fasta, output_folder, threads, min_size):
    """ Dereplicate, sort by size, and then cluster otus
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        fasta_file (string): The path to the fasta file (filtered and dereplicated).
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        min_size (int): Min size of the reads to filter.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: Path to the otu fasta file
    """
    
    dereplicated_fasta = dereplicate(workflow, method, fasta_file, output_folder, threads)
    
    sorted_fasta = sort_by_size(workflow, method, dereplicated_fasta, output_folder, min_size)
    
    otu_fasta = cluster_otus(workflow, method, sorted_fasta, reference_fasta, output_folder)
    
    return otu_fasta
    
    
def dereplicate(workflow, method, fasta_file, output_folder, threads):
    """ Dereplicate reads
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        fasta_file (string): The path to the fasta file (filtered and dereplicated).
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: Path to the dereplicated fasta file
    """    
    
    # get the name of the output files
    output_file = utilities.name_files("all_samples_dereplicated.fasta", output_folder)

    if method == "vsearch":
        workflow.add_task(
            "export OMP_NUM_THREADS=[args[0]]; " + \
            "vsearch --derep_fulllength [depends[0]] --output [targets[0]] --sizein --sizeout --threads [args[0]]",
            depends=[fasta_file, TrackedExecutable("vsearch")],
            targets=output_file,
            args=threads,
            name="vsearch_derep_fulllength")
    else:
        workflow.add_task(
            "export OMP_NUM_THREADS=[args[0]]; "+\
            "usearch -derep_fulllength [depends[0]] -fastaout [targets[0]] -sizeout -threads [args[0]]",
            depends=[fasta_file,TrackedExecutable("usearch")],
            targets=output_file,
            args=threads,
            name="usearch_derep_fulllength")

    return output_file

def sort_by_size(workflow, method, fasta_file, output_folder, min_size):
    """ Sort reads by size, removing those that are not of min size
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        fasta_file (string): The path to the fasta file (filtered and dereplicated).
        output_folder (string): The path of the output folder.
        min_size (int): Min size of the reads to filter.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: Path to the fasta file sorted by size
    """    
    
    # get the name of the output files
    output_file = utilities.name_files("all_samples_sorted.fasta", output_folder)
    if method == "vsearch":
        workflow.add_task(
            "vsearch --sortbysize [depends[0]] --output [targets[0]] --minsize [args[0]]",
            depends=[fasta_file, TrackedExecutable("vsearch")],
            targets=output_file,
            args=min_size,
            name="vsearch_sortbysize")
    else:
        workflow.add_task(
            "usearch -sortbysize [depends[0]] -fastaout [targets[0]] -minsize [args[0]]",
            depends=[fasta_file,TrackedExecutable("usearch")],
            targets=output_file,
            args=min_size,
            name="usearch_sortbysize")

    return output_file

def cluster_otus(workflow, method, fasta_file, reference_fasta, output_folder):
    """ Cluster the otus with usearch
    
    Args:
        workflow (anadama2.workflow): an instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        fasta_file (string): the path to the fasta file (filtered and dereplicated).
        reference_fasta (string): the path to reference fasta db
        output_folder (string): the path of the output folder.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: Path to the fasta file sorted by size

    """    
    
    # get the name of the output files
    output_fasta = utilities.name_files("all_samples_otus_nonchimeras.fasta", output_folder)

    if method == "vsearch":
        output_txt = utilities.name_files("all_samples_vsearch_otus.txt", output_folder)
        all_otus = utilities.name_files("all_otus.fasta", output_folder)
        workflow.add_task(
            "vsearch --cluster_size [depends[0]] --consout [targets[0]] --id 0.97 --relabel 'OTU' --uc [targets[1]]",
            depends=[fasta_file, TrackedExecutable("vsearch")],
            targets=[all_otus, output_txt],
            name="vsearch_cluster_otus")

        workflow.add_task(
            "vsearch --uchime_ref [depends[0]] --nonchimeras [targets[0]] --strand plus --db [args[0]]",
            depends=[all_otus, TrackedExecutable("vsearch")],
            targets=[output_fasta],
            args=[reference_fasta],
            name="vsearch_nochim")
    else:
        output_txt = utilities.name_files("all_samples_uparse_otus.txt", output_folder)
        workflow.add_task(
            "usearch -cluster_otus [depends[0]] -otus [targets[0]] -relabel 'OTU' -uparseout [targets[1]]",
            depends=[fasta_file,TrackedExecutable("usearch")],
            targets=[output_fasta, output_txt],
            name="usearch_cluster_otus")

    return output_fasta 
  
  
def taxonomic_profile(workflow, method, filtered_fasta_file, truncated_fasta_file, original_fasta_file, output_folder, threads, percent_identity,
    reference_usearch, reference_fasta, reference_taxonomy, min_size, bypass_msa=False):
    """ Pick otus, cluster centroids, otu mapping, reference mapping to create open/closed reference taxonomy files
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearhc(default) or vsearch
        filtered_fasta_file (string): The path to the fasta file (filtered and dereplicated).
        truncated_fasta_file (string): The path to the fasta file (truncated not qced).
        original_fasta_file (string): The path to the fasta file (not qc or truncated).
        output_folder (string): The path of the output folder.
        threads (int): The number of threads for each task.
        percent_identity (float): The percent identity to use for alignments.
        reference_usearch (string): The path to the reference usearch formatted database.
        reference_fasta (string): The path to the reference fasta file.
        reference_taxonomy (string): The path to the reference taxonomy file.
        min_size (int): Min size of the reads to filter.
        bypass_msa (bool): Bypass msa clustering and tree generation.        

    Requires:
        usearch(as of usearch v9: has built-in de novo chimera filtering) or vsearch
        clustal omega: multiple sequence alignment for proteins 
        
    Returns:
        list: Path to the fasta file sorted by size

    """      
    
    # first pick otus
    otu_fasta = pick_otus(workflow, method, filtered_fasta_file, reference_fasta, output_folder, threads, min_size)
    
    # centroid OTU sequence alignment
    # get the name of the output files
    if not bypass_msa:
        centroid_fasta = files.SixteenS.path("msa_nonchimera", output_folder)
        centroid_alignment(workflow, otu_fasta, centroid_fasta, threads, task_name="clustalo_nonchimera")
    
    # align the reads to the otus
    otu_alignment_uc = utilities.name_files("all_samples_otu_mapping_results.uc", output_folder)
    otu_alignment_tsv = utilities.name_files("all_samples_otu_mapping_results.tsv", output_folder)
    global_alignment(workflow, method, truncated_fasta_file, otu_fasta, percent_identity, threads, otu_alignment_uc, otu_alignment_tsv)
    
    # align the otus to the reference database
    reference_alignment_uc = utilities.name_files("all_samples_green_genes_mapping_results.uc", output_folder)
    reference_alignment_tsv = utilities.name_files("all_samples_green_genes_mapping_results.tsv", output_folder)
    global_alignment(workflow, method, otu_fasta, reference_usearch, percent_identity, threads, reference_alignment_uc, reference_alignment_tsv, top_hit_only=True)
    
    # create the open/cosed reference tables
    closed_reference_tsv, closed_ref_fasta = build_otu_tables(workflow, reference_taxonomy, reference_fasta, reference_alignment_uc, otu_alignment_uc, otu_fasta, original_fasta_file, output_folder)
    
    # cluster the closed reference otu fasta sequences and generate a tree
    if not bypass_msa:
        centroid_closed_fasta = files.SixteenS.path("msa_closed_reference", output_folder)
        closed_tree = utilities.name_files("closed_reference.tre", output_folder)
        centroid_alignment(workflow, closed_ref_fasta, centroid_closed_fasta, threads, task_name="clustalo_closed_reference")
        create_tree(workflow, centroid_closed_fasta, closed_tree)    

    return closed_reference_tsv, closed_ref_fasta

def centroid_alignment(workflow, fasta_file, output_fasta, threads, task_name=None):
    """ Run clustalo for centroid alignment
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        fasta_file (string): The path to the fasta file (otu sequences).
        output_fasta (string): The path of the output file.
        threads (int): The number of threads/cores for each task.
        task_name (string): The custom name of the task.
        
    Requires:
        clustal omega: multiple sequence alignment for proteins 
        
    Returns:
        string: Path to the clustered otu file

    """     
    
    # remove existing output file if already exists as clustalo will not overwrite
    workflow.add_task(
        "remove_if_exists.py [targets[0]] ; "
        "clustalo -i [depends[0]] -o [targets[0]] --threads [args[0]]",
        depends=[fasta_file,TrackedExecutable("clustalo",version_command="echo 'clustalo' `clustalo --version`")],
        targets=output_fasta,
        args=threads,
        name=task_name if task_name else "clustalo")

def create_tree(workflow, msa_file, tree_file):
    """ Run fasttree to generate phylogenetic tree
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        msa_file (string): The path to the multiple sequence alignment file.
        tree_file (string): The path of the output tree file.

    Requires:
        fasttree: phylogenetic trees from alignments of nucleotide or protein sequences

    Returns:
        None
    """

    # run fasttree on msa file with default settings
    workflow.add_task(
        "FastTree -gtr -nt [depends[0]] > [targets[0]]",
        depends=msa_file,
        targets=tree_file,
        name="fasttree")

def global_alignment(workflow, method, fasta_file, database_file, id, threads, output_file_uc, output_file_tsv, top_hit_only=None):
    """ Run global alignment with the database provided 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        method (string): tools for sequence analysis - usearch(default) or vsearch
        fasta_file (string): The path to the fasta file (filtered and dereplicated).
        database_file (string): Path to the database file (fasta or usearch format)
        id (float): The percent identity for alignment
        threads (int): The number of threads/cores for each task
        output_file_uc (string): The name for the uc output file 
        output_file_tsv (string): The name for the tsv output file 
        top_hit_only (bool): If set, only get the top hits.
        
    Requires:
        usearch or vsearch
        
    Returns:
        list: Path to the mapping results files

    """    
    
    optional_flags=""
    
    # remove existing output file if already exists as clustalo will not overwrite
    if method == "vsearch":
        if top_hit_only:
            optional_flags = " -top_hits_only"
        workflow.add_task_gridable(
            "export OMP_NUM_THREADS=[args[0]]; " + \
            "vsearch -usearch_global [depends[0]] -db [depends[1]] -strand plus -id [args[1]] -uc [targets[0]] -otutabout [targets[1]] -threads [args[0]]" + optional_flags,
            depends=[fasta_file, database_file, TrackedExecutable("vsearch")],
            targets=[output_file_uc, output_file_tsv],
            args=[threads, id],
            name="vsearch_global",
            time=60,  # 60 minutes
            mem=2 * 1024,  # 2 GB
            cores=threads)  # time/mem based on 8 cores
    else:
        if top_hit_only:
            optional_flags = " -top_hit_only"
        workflow.add_task_gridable(
            "export OMP_NUM_THREADS=[args[0]]; "+\
            "usearch -usearch_global [depends[0]] -db [depends[1]] -strand 'both' -id [args[1]] -uc [targets[0]] -otutabout [targets[1]] -threads [args[0]]"+optional_flags,
            depends=[fasta_file, database_file, TrackedExecutable("usearch")],
            targets=[output_file_uc, output_file_tsv],
            args=[threads, id],
            name="usearch_global",
            time=60, # 60 minutes
            mem=2*1024, # 2 GB
            cores=threads) # time/mem based on 8 cores
   
   
def build_otu_tables(workflow, reference_taxonomy, reference_fasta, reference_mapping_results_uc, otu_mapping_results_uc, otu_fasta, original_fasta, output_folder):
    """ Build the open/closed reference otu tables, denovo table, and corresponding fasta files

    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        reference_taxonomy (string): The path to the reference taxonomy file.
        reference_fasta (string): The path to the reference fasta file.
        reference_mapping_results_uc (string): The path to the reference mapping uc results file.
        otu_mapping_results_uc (string): The path to the otu mapping uc results file.
        otu_fasta (string): The path to the fasta file of otu sequences.
        original_fasta (string): The path to the fasta file (not qc or truncated).
        output_folder (string): The path of the output folder.
        
    Requires:
        None
        
    Returns:
        list: The path to the closed reference otu files (tsv and fasta)

    """
    
    # name the output files
    open_ref_tsv = files.SixteenS.path("otu_table_open_reference", output_folder) 
    open_ref_fasta = utilities.name_files("all_samples_open_reference.fasta", output_folder)
    closed_ref_tsv = files.SixteenS.path("otu_table_closed_reference", output_folder)
    closed_ref_fasta = utilities.name_files("all_samples_closed_reference.fasta", output_folder)
    denovo_tsv = utilities.name_files("all_samples_denovo_otu_table.tsv", output_folder)
    read_counts = files.SixteenS.path("read_count_table", output_folder)
  

    workflow.add_task(
        "create_otu_tables_from_alignments.py [depends[0]] [depends[1]] [depends[2]] [depends[3]] [depends[4]] [depends[5]] "+\
        "[targets[0]] [targets[1]] [targets[2]] [targets[3]] [targets[4]] [targets[5]]",
        depends=[reference_taxonomy, reference_fasta, reference_mapping_results_uc, otu_mapping_results_uc, otu_fasta, original_fasta],
        targets=[open_ref_tsv,open_ref_fasta,closed_ref_tsv,closed_ref_fasta,denovo_tsv,read_counts])
    
    return closed_ref_tsv, closed_ref_fasta

def run_picrust2(task, threads, otus=False):
    """ Run picrust2, first changing sequence ids to avoid all numeric (as per picrust2 tutorial) """

    picrust2_input_dir = os.path.dirname(task.depends[0].name)
    picrust2_output_dir = os.path.dirname(task.targets[0].name)

    if otus:
        reformat_input_fasta = utilities.name_files(task.depends[0].name, picrust2_input_dir, tag="picrust_reformatted_input", create_folder=True)
        with open(task.depends[0].name) as file_handle:
            with open(reformat_input_fasta, "w") as file_handle_write:
                for line in file_handle:
                    if line.startswith(">"):
                        line=line.replace(">",">seq")
                    file_handle_write.write(line)
    else:
        reformat_input_fasta = task.depends[0].name

    reformat_input_tsv = utilities.name_files(task.depends[1].name, picrust2_input_dir, tag="picrust_reformatted_input")
    with open(task.depends[1].name) as file_handle:
        with open(reformat_input_tsv, "w") as file_handle_write:
            header = file_handle.readline()
            file_handle_write.write(header)
            for line in file_handle:
                if otus:
                    line="seq"+line
                line="\t".join(line.split("\t")[:-1])+"\n"
                file_handle_write.write(line)

    utilities.run_task("remove_if_exists.py [args[0]] --is-folder ; picrust2_pipeline.py -s [args[1]] -i [args[2]] -o [args[0]] -p [args[3]]",
        depends=task.depends,
        targets=task.depends,
        args=[picrust2_output_dir, reformat_input_fasta, reformat_input_tsv, threads])


def functional_profile(workflow, closed_reference_tsv, closed_reference_fasta, picrust_version, threads, output_folder, otus):
    """ Run picrust for functional profiling
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        closed_reference_tsv (string): The path to the closed reference tsv file.
        closed_reference_fasta (string): The path to the closed reference fasta file.
        picrust_version (str): The version of picrust to use.
        threads (int): The number of threads/cores for each task.
        output_folder (string): The path of the output folder.
        otus (bool): Are the inputs from OTUs (so all numerical ids).
        
    Requires:
        Picrust v1.1 or v2: Software to predict metagenome function.
        Biom v2: A tool for general use formatting of biological data.
        
    Returns:
        string: The path to the functional data file in tsv format.

    """
    
    if picrust_version == "1":
        # convert the tsv file to biom format
        closed_reference_biom_file = utilities.name_files(closed_reference_tsv,output_folder,extension="biom")
        convert_to_biom_from_tsv(workflow,closed_reference_tsv,closed_reference_biom_file,options="--process-obs-metadata=taxonomy --output-metadata-id=taxonomy")

        # run picrust to get functional data
        functional_data_categorized_biom, functional_data_predicted_biom = picrust(workflow, closed_reference_biom_file,
                                                                                   output_folder)

        # convert the predited biom file to tsv
        functional_data_predicted_tsv = utilities.name_files(functional_data_predicted_biom, output_folder, extension="tsv")
        convert_from_biom_to_tsv(workflow, functional_data_predicted_biom, functional_data_predicted_tsv)

        # convert the categorized biom file to tsv
        functional_data_categorized_tsv = utilities.name_files(functional_data_categorized_biom, output_folder,
                                                           extension="tsv")
        convert_from_biom_to_tsv(workflow, functional_data_categorized_biom, functional_data_categorized_tsv)

        return functional_data_categorized_tsv

    else:
        # run the v2 pipeline
        functional_data_predicted_tre = utilities.name_files("out.tre", output_folder, subfolder="picrust2")
        workflow.add_task(
            utilities.partial_function(run_picrust2,threads=threads,otus=otus),
            depends=[closed_reference_fasta, closed_reference_tsv],
            targets=functional_data_predicted_tre)
        return functional_data_predicted_tre

def picrust(workflow,otu_table_biom,output_folder):
    """ Runs picrust normalize, then predict
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        out_table_biom (string): The path to the biom file (closed reference otu table).
        output_folder (string): The path of the output folder.
        
    Requires:
        Picrust v1.1: Software to predict metagenome function.
        
    Returns:
        string: The path to the functional data file in biom format.
    
    """

    # normalize the otu table
    normalized_otu_table=utilities.name_files("all_samples_normalize_by_copy_number.biom", output_folder)
    # first remove target file as picrust will not overwrite
    # expects biom file is json (not hdf5) format
    workflow.add_task(
        "remove_if_exists.py [targets[0]] ; "+\
        "normalize_by_copy_number.py -i [depends[0]] -o [targets[0]]",
        depends=[otu_table_biom,TrackedExecutable("normalize_by_copy_number.py")],
        targets=normalized_otu_table,
        name="normalize_by_copy_number.py")
    
    # predict metagenomes
    predict_metagenomes_table=utilities.name_files("all_samples_predict_metagenomes.biom", output_folder)
    # first remove target file as picrust will not overwrite
    workflow.add_task(
        "remove_if_exists.py [targets[0]] ; "+\
        "predict_metagenomes.py -i [depends[0]] -o [targets[0]]",
        depends=[normalized_otu_table,TrackedExecutable("predict_metagenomes.py")],
        targets=predict_metagenomes_table,
        name="predict_metagenomes.py")

    # categorize by function
    categorized_function_table = utilities.name_files("all_samples_categorize_by_function.biom", output_folder)
    # first remove target file as picrust will not overwrite
    workflow.add_task(
        "remove_if_exists.py [targets[0]] ; " + \
        "categorize_by_function.py -i [depends[0]] -o [targets[0]] --level 3 -c KEGG_Pathways",
        depends=[predict_metagenomes_table, TrackedExecutable("categorize_by_function.py")],
        targets=categorized_function_table,
        name="categorize_by_function.py")

    return categorized_function_table, predict_metagenomes_table
    

def convert_to_biom_from_tsv(workflow, tsv_file, biom_file, table_type="OTU table", options=""):
    """ Convert a tsv file to a biom file 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        tsv_file (string): The path to a tsv file (otu table format) with taxonomy.
        biom_file (string): The path to write the new biom file
        table_type (string): The type of table to convert
        options (string): Additional options to provide to the convert function
        
    Requires:
        Biom v2: A tool for general use formatting of biological data.
        
    """
    
    # first remove biom file if exists as biom will not overwrite
    workflow.add_task(
        "remove_if_exists.py [targets[0]] ; "+\
        "biom convert -i [depends[0]] -o [targets[0]] --table-type='"+table_type+"' --to-json "+options,
        depends=[tsv_file,TrackedExecutable("biom")],
        targets=biom_file,
        name="biom")
    
def convert_from_biom_to_tsv(workflow, biom_file, tsv_file, table_type="OTU table", options=""):
    """ Convert from a biom file to a tsv file 
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        biom_file (string): The path to the biom file.
        tsv_file (string): The path to write the new tsv file.
        table_type (string): The type of table to convert
        options (string): Additional options to provide to the convert function
        
    Requires:
        Biom v2: A tool for general use formatting of biological data.
        
    """
    
    # first remove biom file if exists as biom will not overwrite
    workflow.add_task(
        "remove_if_exists.py [targets[0]] ; "+\
        "biom convert -i [depends[0]] -o [targets[0]] --table-type='"+table_type+"' --to-tsv "+options,
        depends=[biom_file,TrackedExecutable("biom")],
        targets=tsv_file,
        name="biom")

      
