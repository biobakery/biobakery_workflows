"""
bioBakery Workflows: tasks.dada module
A collection of tasks for DADA2 workflow with 16s amplicon sequences

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
from biobakery_workflows import files, config
from anadama2.tracked import TrackedDirectory

import os.path


def filter_trim(workflow, input_folder, output_folder, maxee, trunc_len_max, pair_id, threads):
    """ Filters samples by threshould and trims them

       Args:
           workflow (anadama2.workflow): An instance of the workflow class.
           input_folder (string): The path of the input folder.
           output_folder (string): The path of the output folder.

       Requires:
           None

       Returns:
           Creates folder with filtered and trimed sample files
           Creates file that contains read counts before and after filtering
    """
    reads_plotF_png = files.SixteenS.path("readF_qc", output_folder)
    reads_plotR_png = files.SixteenS.path("readR_qc", output_folder)

    readcounts_tsv_path = output_folder + "/Read_counts_after_filtering.tsv"
    readcounts_rds_path = output_folder + "/Read_counts_after_filtering.rds"
    dir_path = os.path.dirname(os.path.realpath(__file__))
    main_folder = dir_path + "/.."
    print(main_folder)

    workflow.add_task(
        "[vars[0]]/scripts/filter_and_trim.py \
          --input=[args[0]]\
          --output=[args[1]]\
          --maxee=[args[2]]\
          --trunc_len_max=[args[3]]\
          --readcounts_tsv_path=[targets[0]]\
          --readcounts_rds_path=[targets[1]]\
          --plotF=[targets[2]]\
          --plotR=[targets[3]]\
          --pair_id=[args[4]]\
          --threads=[args[5]]",
        depends=TrackedDirectory(input_folder),
        targets=[readcounts_tsv_path, readcounts_rds_path, reads_plotF_png, reads_plotR_png],
        args=[input_folder, output_folder, maxee, trunc_len_max, pair_id, threads],
        vars=main_folder,
        name="filter_and_trim"
    )
    return readcounts_tsv_path,readcounts_rds_path


def learn_errors(workflow, output_folder, readcounts_tsv_path, threads):
    """ Learns error rates for each sample

       Args:
           workflow (anadama2.workflow): An instance of the workflow class.
           output_folder (string): The path of the output folder.

       Requires:
           Path to read counts file

       Returns:
           Creates files that contain forward and reverse error rates for each sample
    """

    error_ratesF_png = files.SixteenS.path("error_ratesF", output_folder)
    error_ratesR_png = files.SixteenS.path("error_ratesR", output_folder)

    error_ratesF_rds = output_folder + "/error_ratesFWD.rds"
    error_ratesR_rds = output_folder + "/error_ratesREV.rds"
    filt_path = output_folder + "/filtered_input"

    dir_path = os.path.dirname(os.path.realpath(__file__))
    main_folder = dir_path + "/.."

    workflow.add_task(
        "[vars[0]]/scripts/learn_error_rates.py\
          --filt_path=[vars[1]]\
          --error_ratesF_png=[targets[0]]\
          --error_ratesR_png=[targets[1]]\
          --error_ratesF_rds=[targets[2]]\
          --error_ratesR_rds=[targets[3]]\
          --threads=[args[0]]",
        depends=[readcounts_tsv_path],
        targets=[error_ratesF_png, error_ratesR_png, error_ratesF_rds, error_ratesR_rds],
        args=[threads],
        vars=[main_folder,filt_path],
        name="learn_error_rates"
    )
    return error_ratesF_rds, error_ratesR_rds


def merge_paired_ends(workflow, output_folder, error_ratesF_rds, error_ratesR_rds, threads):
    """ Dereplicates and merges paired reads

        Args:
            workflow (anadama2.workflow): An instance of the workflow class.
            input_folder (string): The path of the input folder.
            output_folder (string): The path of the output folder.

        Requires:
            Path to error rates files

        Returns:
            Creates files that contain merged and dereplicated reads
     """

    mergers_rds = output_folder + "/mergers.rds"
    sample_names_rds = output_folder + "/sample_names.rds"
    filt_path = output_folder + "/filtered_input"
    dir_path = os.path.dirname(os.path.realpath(__file__))
    main_folder = dir_path + "/.."

    workflow.add_task(
        "[vars[0]]/scripts/merge_paired_ends.py\
          --filt_path=[vars[1]]\
          --error_ratesF_rds=[depends[0]]\
          --error_ratesR_rds=[depends[1]]\
          --mergers_rds=[targets[0]]\
          --sample_names_rds=[targets[1]]\
          --threads=[args[0]]",
        depends=[error_ratesF_rds, error_ratesR_rds],
        targets=[mergers_rds, sample_names_rds],
        args=threads,
        vars=[main_folder, filt_path],
        name="dereplicate_and_merge"
    )



def const_seq_table(workflow, output_folder, threads):
    """ Builds sequence table, removes chimeras

       Args:
           workflow (anadama2.workflow): An instance of the workflow class.
           input_folder (string): The path of the input folder.
           output_folder (string): The path of the output folder.

       Requires:
           Path to merged reads

       Returns:
           Creates file that contains sequence table
           Creates file that contains read counts on each step
    """

    read_counts_steps_path = files.SixteenS.path("counts_each_step", output_folder)

    mergers_rds = output_folder + "/mergers.rds"
    seqtab_rds = output_folder + "/seqtab_final.rds"
    sample_names_rds = output_folder + "/sample_names.rds"
    all_counts_tsv = output_folder + "/all_samples_SV_counts.tsv"
    readcounts_rds = output_folder + "/Read_counts_after_filtering.rds"

    dir_path = os.path.dirname(os.path.realpath(__file__))
    main_folder = dir_path + "/.."

    workflow.add_task(
        "[vars[0]]/scripts/const_seq_table.py\
          --mergers_rds=[depends[0]]\
          --sample_names_rds=[depends[1]]\
          --readcounts_steps_path=[targets[0]]\
          --seqtab_rds=[targets[1]]\
          --readcounts_rds=[depends[2]]\
          --all_counts_tsv=[targets[2]]\
          --threads=[args[0]]",
        depends=[mergers_rds, sample_names_rds, readcounts_rds],
        targets=[read_counts_steps_path, seqtab_rds, all_counts_tsv],
        args=threads,
        vars=main_folder,
        name="construct_sequence_table"
    )
    return  read_counts_steps_path


def phylogeny(workflow,output_folder):
    """ Aligns sequences and reconstructs phylogeny

       Args:
           workflow (anadama2.workflow): An instance of the workflow class.
           output_folder (string): The path of the output folder.

       Requires:
           Path to sequence table data

       Returns:
           Creates file that contains pylogeny
    """

    msa_fasta_path = files.SixteenS.path("msa_nonchimera", output_folder)
    seqtab_rds = output_folder + "/seqtab_final.rds"

    dir_path = os.path.dirname(os.path.realpath(__file__))
    main_folder = dir_path + "/.."

    workflow.add_task(
        "[vars[0]]/scripts/phylogeny.py\
          --seqtab_rds=[depends[0]]\
          --msa_fasta_path=[targets[0]]",
        depends=[seqtab_rds],
        targets=[msa_fasta_path],
        args=[output_folder],
        vars=main_folder,
        name="phylogeny"
    )
    return msa_fasta_path


def fasttree(workflow, output_folder, msa_fasta_path):
    """ Generates a phylogenetic tree

       Args:
           workflow (anadama2.workflow): An instance of the workflow class.
           output_folder (string): The path of the output folder.

       Requires:
           Path to phylogeny file

       Returns:
           Creates file that contains pylogenic tree
    """

    fasttree_file_path = output_folder + "/closed_reference.tre"

    workflow.add_task(
        "FastTree -gtr -nt  [depends[0]] >  [targets[0]]",
        depends=[msa_fasta_path],
        targets=[fasttree_file_path],
        args=[output_folder],
        name="fasttree"
    )
    return fasttree_file_path


def assign_taxonomy(workflow, output_folder, ref_path, threads):
    """ Assigns taxonomy using gg, silva, or rdp database

       Args:
           workflow (anadama2.workflow): An instance of the workflow class.
           output_folder (string): The path of the output folder.

       Requires:
           Sequence table data

       Returns:
           Creates closed reference file
    """
    # check what reference db to use for taxonomy assignment

    seqtab_rds = output_folder + "/seqtab_final.rds"

    otu_closed_ref_path = files.SixteenS.path("otu_table_closed_reference", output_folder)

    if ref_path == "silva":
        refdb_path = config.SixteenS().silva_dada2
        refdb_species_path = config.SixteenS().silva_species_dada2
    elif ref_path == "rdp":
        refdb_path = config.SixteenS().rdp_dada2
        refdb_species_path = config.SixteenS().rdp_species_dada2
    else:
        refdb_path = config.SixteenS().greengenes_dada2
        refdb_species_path = "None"

    dir_path = os.path.dirname(os.path.realpath(__file__))
    main_folder = dir_path + "/.."


    workflow.add_task(
        "[vars[2]]/scripts/assign_taxonomy.py\
          --refdb_path=[vars[0]]\
          --refdb_species_path=[vars[1]]\
          --seqtab_rds=[depends[0]]\
          --otu_closed_ref_path=[targets[0]]\
          --threads=[args[0]]",
        depends=[seqtab_rds],
        targets=[otu_closed_ref_path],
        args=threads,
        vars=[refdb_path, refdb_species_path, main_folder],
        name="assign_taxonomy"
    )

    return otu_closed_ref_path
