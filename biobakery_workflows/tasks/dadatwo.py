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

def filter_trim(workflow, input_folder, output_folder):
    
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
         readcounts_rds_path = output_folder + "/Read_counts_filt.rds"
       
         workflow.add_task(
             "biobakery_workflows/scripts/filter_and_trim.R \
               --input_dir=[args[0]]\
               --output_dir=[args[1]]\
               --readcounts_tsv_path=[targets[0]]\
               --readcounts_rds_path=[targets[1]]\
               --reads_plotF=[targets[2]]\
               --reads_plotR=[targets[3]]",
             depends = input_folder,
             targets = [readcounts_tsv_path, readcounts_rds_path, reads_plotF_png, reads_plotR_png],
             args = [input_folder, output_folder],
             name ="filter_and_trim"
             )
         return readcounts_tsv_path
     

def learn_error(workflow,output_folder,readcounts_tsv_path):
    
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
         
         error_ratesF_path= output_folder + "/error_ratesFWD.rds"
         error_ratesR_path = output_folder + "/error_ratesREV.rds"
       
         workflow.add_task(
             "biobakery_workflows/scripts/learn_error_rates.R\
               --output_dir=[args[0]]\
               --error_ratesF_png=[targets[0]]\
               --error_ratesR_png=[targets[1]]\
               --error_ratesF_path=[targets[2]]\
               --error_ratesR_path=[targets[3]]",
             depends = [readcounts_tsv_path],
             targets = [error_ratesF_png, error_ratesR_png, error_ratesF_path, error_ratesR_path],  
             args = [output_folder],
             name = "learn_error_rates"
             )
         return error_ratesF_path, error_ratesR_path
     

def merge_paired_ends(workflow, input_folder, output_folder, error_ratesF_path, error_ratesR_path):
    
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

        mergers_file_path = output_folder + "/mergers.rds"
        
        workflow.add_task(
            "biobakery_workflows/scripts/merge_paired_ends.R\
              --input_dir=[args[0]]\
              --output_dir=[args[1]]\
              --error_ratesF_path=[depends[0]]\
              --error_ratesR_path=[depends[1]]\
              --mergers_file_path=[targets[0]]",
            depends = [error_ratesF_path, error_ratesR_path],
            targets = [mergers_file_path],                       
            args = [input_folder, output_folder],
            name = "dereplicate_and_merge"
            )
        return mergers_file_path
    

def const_seq_table(workflow, input_folder, output_folder, mergers_file_path):
    
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
         
         seqtab_file_path = output_folder + "/seqtab_final.rds"

         workflow.add_task(
            "biobakery_workflows/scripts/const_seq_table.R\
              --input_dir=[args[0]]\
              --output_dir=[args[1]]\
              --merged_file=[depends[0]]\
              --read_counts_steps_path=[targets[0]]\
              --seqtab_file_path=[targets[1]]",
            depends = [mergers_file_path],
            targets = [read_counts_steps_path, seqtab_file_path],                              
            args = [input_folder, output_folder],
            name = "construct_sequence_table"
            )
         return seqtab_file_path, read_counts_steps_path


def phylogeny(workflow, output_folder, seqtab_file_path ):
    
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

         workflow.add_task(
            "biobakery_workflows/scripts/phylogeny.R\
              --output_dir=[args[0]]\
              --seqtab_file_path=[depends[0]]\
              --msa_fasta_path=[targets[0]]",
            depends = [seqtab_file_path],
            targets = [msa_fasta_path],                              
            args = [output_folder],
            name = "phylogeny"
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
            depends = [msa_fasta_path],
            targets = [fasttree_file_path],                              
            args = [output_folder],
            name = "fasttree"
            )
         return fasttree_file_path


def assign_taxonomy(workflow, output_folder, seqtab_file_path, ref_path):
    
         """ Assigns taxonomy using gg, silva, or rdp database
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Sequence table data
                
            Returns:
                Creates closed reference file
         """
 #check what reference db to use for taxonomy assignment

         otu_closed_ref_path  = files.SixteenS.path("otu_table_closed_reference", output_folder)
         
         if ref_path == "silva": 
             refdb_path = config.SixteenS().silva_dada2  
             workflow.add_task(
                "biobakery_workflows/scripts/assign_taxonomy_silva_rdp.R\
                  --output_dir=[args[0]]\
                  --refdb_path=[vars[0]]\
                  --seqtab_file_path=[depends[0]]\
                  --otu_closed_ref_path=[targets[0]]",
                depends = [seqtab_file_path],
                targets = [otu_closed_ref_path],                              
                args = [output_folder],
                vars =[refdb_path],
                name = "assign_silva"
                
                )
         elif ref_path == "rdp":
             refdb_path = config.SixteenS().rdp_dada2  
             workflow.add_task(
                "biobakery_workflows/scripts/assign_taxonomy_silva_rdp.R\
                  --output_dir=[args[0]]\
                  --refdb_path=[vars[0]]\
                  --seqtab_file_path=[depends[0]]\
                  --otu_closed_ref_path=[targets[0]]",
                depends = [seqtab_file_path],
                targets = [otu_closed_ref_path],                              
                args = [output_folder],
                vars = [refdb_path],
                name = "assign_rdp"
                
                )
         else:    
             refdb_path = config.SixteenS().greengenes_dada2
             all_samples_sv_gg13_path = output_folder + "/all_samples_GG13-8-taxonomy.tsv"
     
             workflow.add_task(
                "biobakery_workflows/scripts/assign_taxonomy.R\
                  --output_dir=[args[0]]\
                  --refdb_path=[vars[0]]\
                  --seqtab_file_path=[depends[0]]\
                  --otu_closed_ref_path=[targets[0]]\
                  --allsamples_gg_tax_path=[targets[1]]",
                depends = [seqtab_file_path],
                targets = [otu_closed_ref_path, all_samples_sv_gg13_path],                              
                args = [output_folder],
                vars = [refdb_path],
                name = "assign_gg"
            )
             
         return otu_closed_ref_path
     
 
def remove_tmp_files(workflow, output_folder, otu_closed_ref_path,
                      msa_fasta_path, fasttree_file_path):
    
         """ Removes temporary .rds files
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Path to output folder
                
            Returns:
                None
          """
         
         rm_out_file = output_folder + "/out.txt"
         rm_error_file = output_folder + "/error.txt"
             
         workflow.add_task(
                     "rm  [args[0]]/*.rds  > [targets[0]] 2 > [targets[1]] ",
                     depends = [otu_closed_ref_path, msa_fasta_path, fasttree_file_path],
                     args = [output_folder],
                     targets = [rm_out_file, rm_error_file],
                     name = "rm_tmp_files"
                     )
         