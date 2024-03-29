"""
bioBakery Workflows: tasks.files module
A collection of file names used by tasks

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
import copy
import sys

from anadama2 import reporters

from .utilities import name_files

class FileInfo(object):
    def __init__(self, name=None, subfolder=None, tag=None, extension=None, description=None):
        # set a list of non-path keywords
        self.non_path_keywords=["description"]
        
        # concat multiple strings if present in description
        if description and isinstance(description,tuple):
            description="\n".join(description)
        
        keywords={"names":name, "subfolder":subfolder, "tag":tag, "extension":extension, "description":description}
        self.keywords = {key:value for key, value in keywords.items() if value}
        
    def get_path_keywords(self):
        info = copy.copy(self.keywords)
        
        # remove the non-path keywords
        for key in self.non_path_keywords:
            try:
                del info[key]
            except KeyError:
                pass
        
        return info
    
    def __getitem__(self, key):
        """ Return the file info """
        try:
            value=self.keywords[key]
        except KeyError:
            value=""
        return value

class Workflow(object):
    file_info = {}
    
    file_info["log"]=FileInfo(reporters.LOG_FILE_NAME,
        description="The AnADAMA2 workflow log.")
    
    @classmethod
    def path(cls, name, main_folder="", none_if_not_found=None, error_if_not_found=None, search_for_file=False, **keywords):
        merged_keywords = copy.copy(keywords)
        merged_keywords.update(cls.file_info[name].get_path_keywords())
        file_path=name_files(folder=main_folder, **merged_keywords)

        # if the file is not found, then look in the input folder
        if not os.path.isfile(file_path) and search_for_file:
            file_name = cls.filename(name)
            file_path=os.path.join(main_folder, file_name)

        # if set, error if the file does not exist
        if error_if_not_found and not os.path.isfile(file_path):
            message="\nERROR: Unable to find file: "+file_path
            desc=cls.description(name)
            if desc:
                message+="\n\nFile description:\n"+desc
            sys.exit(message)
        
        # if set, check if the file exists, if not return None
        if none_if_not_found and not os.path.isfile(file_path):
                file_path = None
           
        return file_path
    
    @classmethod
    def description(cls,name):
        try:
            desc=cls.file_info[name].keywords["description"]
        except (KeyError, AttributeError):
            desc=""
            
        return desc

    @classmethod
    def filename(cls, name):
        try: 
            fname=cls.file_info[name].keywords["names"]
        except (KeyError, AttributeError):
            fname=""

        return fname
    
    @classmethod
    def list_file_path_description(cls,folder,input_files):
        """ List the file names and descriptions in a format to be used in an argument help description """
        
        desc=""
        for required in input_files:
            desc+="\n\n".join(["* "+cls.path(name,folder)+ " ( " + required + " )\n-- "+cls.description(name) for name in input_files[required]])+"\n"
            
        return desc

    @classmethod
    def list_file_description(cls,input_files):
        """ List the file names and descriptions in a format to be used in an argument help description """
        
        desc=""
        for required in input_files:
            desc+="\n".join(["* "+cls.filename(name)+ " ( " + required + " )\n-- "+cls.description(name)+"\n\n" for name in input_files[required]])+"\n"
            
        return desc

class ShotGun(Workflow):
    """ A collection of information of folders/files created by the shotgun tasks """
    
    file_info=copy.copy(Workflow.file_info)
    
    # set the folder names for wmgx_wmtx data workflows
    wmgx_folder_name="whole_metagenome_shotgun"
    wmtx_folder_name="whole_metatranscriptome_shotgun"
    
    # set the kneaddata file name
    file_info["kneaddata_read_counts"]=FileInfo("kneaddata_read_count_table.tsv",subfolder=os.path.join("kneaddata","merged"),
        description=("A tab-delimited file with samples as rows and read counts as columns.",
            "This file is generated by compiling information from the KneadData log files. ",
            "The file contains read counts after each step of the Kneaddata workflow. ",
            "Depending on the input to Kneaddata this file will contain counts for ",
            "single or paired/orphan reads. The read counts from the filtering steps ",
            "are for filtering applied in a serial manner with the reads that do not ",
            "map to the first reference used as input in filtering with the next ",
            "reference database continuing until all databases have been run."))

    # set the taxonomy file names
    file_info["taxonomic_profile"]=FileInfo("metaphlan_taxonomic_profiles.tsv",subfolder=os.path.join("metaphlan","merged"),
        description=("A tab-delimited file with samples as columns and relative abundance as rows.",
            "This file contains the merged taxonomic profiles computed by MetaPhlAn for all samples."))
    file_info["species_counts"]=FileInfo("metaphlan_species_counts_table.tsv",subfolder=os.path.join("metaphlan","merged"),
        description=("A tab-delimited file with samples as rows and counts as columns.",
            "This file contains the counts of total species for each sample using the ",
            "species identified by MetaPhlAn."))
    
    # set the merged feature file names
    file_info["genefamilies"]=FileInfo("genefamilies.tsv",subfolder=os.path.join("humann","merged"),
        description=("A tab-delimited file with samples as columns and gene families ",
            "as rows. This file is a merged set of gene families for all samples ",
            "computed by HUMAnN. This file contains stratified counts as RPKs."))
    file_info["ecs"]=FileInfo("ecs.tsv", subfolder=os.path.join("humann","merged"),
        description=("A tab-delimited file with samples as columns and ecs as rows. ",
            "This file is a merged set of ecs for all samples generated from the gene ",
            "families computed by HUMAnN. This file contains stratified counts as RPKs."))
    file_info["pathabundance"]=FileInfo("pathabundance.tsv", subfolder=os.path.join("humann","merged"),
        description=("A tab-delimited file with samples as columns and pathways ",
            "as rows. This file is a merged set of pathway abundances for all ",
            "samples computed by HUMAnN. This file contains stratified counts ",
            "of non-normalized abundances."))
    
    # set the normed feature file names
    file_info["genefamilies_relab"]=FileInfo("genefamilies_relab.tsv", subfolder=os.path.join("humann","merged"))
    file_info["ecs_relab"]=FileInfo("ecs_relab.tsv", subfolder=os.path.join("humann","merged"),
        description=("A tab-delimited file with samples as columns and ECs as rows.",
                "This file is a merged set of EC abundances for all samples computed ",
                "by HUMAnN. This file contains relative abundances."))
    file_info["pathabundance_relab"]=FileInfo("pathabundance_relab.tsv", subfolder=os.path.join("humann","merged"),
        description=("A tab-delimited file with samples as columns and pathways ",
                "as rows. This file is a merged set of pathway abundances for all ",
                "samples computed by HUMAnN. This file contains stratified counts ",
                "of relative abundances."))
    
    # set the feature count file names
    file_info["genefamilies_relab_counts"]=FileInfo("humann_genefamilies_relab_counts.tsv", subfolder=os.path.join("humann","counts"))
    file_info["ecs_relab_counts"]=FileInfo("humann_ecs_relab_counts.tsv", subfolder=os.path.join("humann","counts"))
    file_info["pathabundance_relab_counts"]=FileInfo("humann_pathabundance_relab_counts.tsv", subfolder=os.path.join("humann","counts"))
    
    # set the all feature counts file names
    file_info["feature_counts"]=FileInfo("humann_feature_counts.tsv", subfolder=os.path.join("humann","counts"),
        description=("A tab-delimited file with samples as rows and features ",
            "as columns. This file includes the total feature counts (non-stratified)",
            "for the features computed by HUMAnN (genes, ecs, and pathways)."))
    file_info["humann_read_counts"]=FileInfo("humann_read_and_species_count_table.tsv", subfolder=os.path.join("humann","counts"),
        description=("A tab-delimited file with samples as rows and counts as columns.",
            "This file was created using the HUMAnN logs. It includes the total number ",
            "of species used to generate the custom database, the total number of initial",
            "reads, and the total reads aligning for both search steps."))
    
    # set the names for the rna/dna normed files
    file_info["genefamilies_norm_ratio"]=FileInfo("rna_dna_relative_expression_unstratified.tsv",subfolder=os.path.join("humann","rna_dna_norm","genes"),
        description=("A tab-delimited file with samples as columns and genes as rows. ",
            "This file includes the normalized RNA abundances as a ratio to DNA abundance. ",
            "This file does not include stratified features."))
    file_info["ecs_norm_ratio"]=FileInfo("rna_dna_relative_expression_unstratified.tsv",subfolder=os.path.join("humann","rna_dna_norm","ecs"),
        description=("A tab-delimited file with samples as columns and ecs as rows. ",
            "This file includes the normalized RNA abundances as a ratio to DNA abundance. ",
            "This file does not include stratified features."))
    file_info["paths_norm_ratio"]=FileInfo("rna_dna_relative_expression_unstratified.tsv",subfolder=os.path.join("humann","rna_dna_norm","paths"),
        description=("A tab-delimited file with samples as columns and pathways as rows. ",
            "This file includes the normalized RNA abundances as a ratio to DNA abundance. ",
            "This file does not include stratified features."))
    
class ShotGunVis(Workflow):
    """ A collection of information of folders/files created by the shotgun vis templates """
    
    file_info=copy.copy(Workflow.file_info)
    
    file_info["microbial_counts"]=FileInfo("microbial_counts_table.tsv",
        description="A tab-delimited file with samples as rows and ratios as "+\
            "columns. Includes the proportion of reads remaining after "+\
            "trimming and filtering in the quality control workflow.")
    file_info["rna_microbial_counts"]=FileInfo("rna_microbial_counts_table.tsv",
        description="A tab-delimited file with RNA samples as rows and ratios as "+\
            "columns. Includes the proportion of reads remaining after "+\
            "trimming and filtering in the quality control workflow.")
    file_info["qc_counts"]=FileInfo("qc_counts_table.tsv",
        description="A tab-delimited file with samples as rows and read counts "+\
            "as columns. Includes the read counts for trimming and filtering steps "+\
            "in the quality control workflow. The reads are single end.")
    file_info["qc_counts_paired"]=FileInfo("qc_counts_pairs_table.tsv",
        description="A tab-delimited file with samples as rows and read counts "+\
            "as columns. Includes the read counts for trimming and filtering steps "+\
            "in the quality control workflow. The reads are paired end and these "+\
            "counts are only for pairs.")
    file_info["qc_counts_orphan"]=FileInfo("qc_counts_orphans_table.tsv",
        description="A tab-delimited file with samples as rows and read counts "+\
            "as columns. Includes the read counts for trimming and filtering steps "+\
            "in the quality control workflow. The reads are paired end and these "+\
            "counts are only for orphans.")
    file_info["rna_qc_counts_paired"]=FileInfo("rna_qc_counts_pairs_table.tsv",
        description="A tab-delimited file with RNA samples as rows and read counts "+\
            "as columns. Includes the read counts for trimming and filtering steps "+\
            "in the quality control workflow. The reads are paired end and these "+\
            "counts are only for pairs.")
    file_info["rna_qc_counts_orphan"]=FileInfo("rna_qc_counts_orphans_table.tsv",
        description="A tab-delimited file with RNA samples as rows and read counts "+\
            "as columns. Includes the read counts for trimming and filtering steps "+\
            "in the quality control workflow. The reads are paired end and these "+\
            "counts are only for orphans.")
    file_info["taxa_counts"]=FileInfo("taxa_counts_table.tsv",
        description="A tab-delimited file with samples as rows and counts as "+\
            "columns. These are the total number of species/genera identified for each "+\
            "sample before and after filtering.")

class SixteenS(Workflow):
    """ A collection of information of folders/files created by the 16s tasks """
    
    file_info=copy.copy(Workflow.file_info)
    
    # set the names for the otu table and read count files
    file_info["otu_table_closed_reference"]=FileInfo("all_samples_taxonomy_closed_reference.tsv",
        description=("A tab-delimited file with samples/taxonomy as columns and taxonomy as rows. ",
            "First column is the OTU id and the last column is the taxonomy. The remaining",
            "columns are sample names. Values are counts."))
    file_info["otu_table_open_reference"]=FileInfo("all_samples_taxonomy_open_reference.tsv",
        description=("A tab-delimited file with samples/taxonomy as columns and taxonomy as rows. ",
            "First column is the OTU id and the last column is the taxonomy. The remaining",
            "columns are sample names. Values are counts. All OTUs without taxonomy are labeled Unclassified."))
    file_info["read_count_table"]=FileInfo("all_samples_read_counts.tsv",
        description=("A tab-delimited file with samples as rows and counts as columns. ",
            "The counts included are the original read count, total number of reads ",
            "mapping to an OTU with known taxonomy, and total reads mapping to an ",
            "unclassified OTU."))
    file_info["eestats2"]=FileInfo("all_samples_eestats2.txt",
        description=("A file with maxee as columns and read lengths as rows."))
    file_info["msa_nonchimera"]=FileInfo("all_samples_clustalo_aligned_nonchimera.fasta",
        description=("A multiple sequence alignment file generated from the nonchimera sequences ",
            "using Clustalo."))
    file_info["msa_closed_reference"]=FileInfo("all_samples_clustalo_aligned_closed_reference.fasta",
        description=("A multiple sequence alignment file generated from the closed reference sequences ",
            "using Clustalo."))
    # set the names for the otu tables and read count files
    file_info["error_ratesF"]=FileInfo("Error_rates_per_sample_FWD.png",
        description=("Plots of forward read error rates in dada2 workflow for each sample"))
    file_info["error_ratesR"]=FileInfo("Error_rates_per_sample_REV.png",
        description=("Plots of reverse read error rates in dada2 workflow for each sample"))
    file_info["readF_qc"]=FileInfo("FWD_read_plot.png",
        description=("Plots of quality of forward reads for each sample"))
    file_info["readR_qc"]=FileInfo("REV_read_plot.png",
        description=("Plots of quality of reverse reads for each sample"))
    file_info["counts_each_step"]=FileInfo("Read_counts_at_each_step.tsv",
        description=("A tab-delimited file with samples as rows and counts as columns. ",
            "The counts included in each step of workflow process"))
    file_info["filtN"]=FileInfo("filtN",
        description=("Folder with N filtered files"))

        
