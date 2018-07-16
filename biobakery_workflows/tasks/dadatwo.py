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

         read_counts_file = "/Read_counts_after_filtering.tsv"
         read_counts_file_path = output_folder + read_counts_file
       
         workflow.add_task(
             "biobakery_workflows/scripts/filter_and_trim.R \
               --input_dir=[args[0]]\
               --output_dir=[args[1]]\
               --readcounts_file=[args[2]]",
             depends = input_folder,
             targets = [read_counts_file_path],
             args = [input_folder, output_folder, read_counts_file],
             name ="filter_and_trim"
             )
         return read_counts_file

def learn_error(workflow,output_folder,read_counts_file):
    
         """ Learns error rates for each sample
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Path to read counts file
                
            Returns:
                Creates files that contain forward and reverse error rates for each sample
         """

         read_counts_file_path = output_folder + read_counts_file
         error_rates_F_file = "/error_rates_F.rds"
         error_rates_R_file = "/error_rates_R.rds"
         error_ratesF_path = output_folder + error_rates_F_file
         error_ratesR_path = output_folder + error_rates_R_file
       
         workflow.add_task(
             "biobakery_workflows/scripts/learn_error_rates.R\
               --output_dir=[args[0]]\
               --error_ratesF=" + error_rates_F_file + " --error_ratesR=" + error_rates_R_file,
             depends = [read_counts_file_path],
             targets = [error_ratesF_path, error_ratesR_path],  
             args = [output_folder],
             name = "learn_error_rates"
             )
         return error_rates_F_file, error_rates_R_file

def merge_paired_ends(workflow, input_folder, output_folder, error_ratesF, error_ratesR):
    
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

        mergers_file = "/mergers.rds"
        mergers_file_path = output_folder + mergers_file
        error_ratesF_path = output_folder + error_ratesF
        error_ratesR_path = output_folder + error_ratesR
        workflow.add_task(
            "biobakery_workflows/scripts/merge_paired_ends.R\
              --input_dir=[args[0]]\
              --output_dir=[args[1]]\
              --error_ratesF=" + error_ratesF + " --error_ratesR=" + error_ratesR,
            depends = [error_ratesF_path, error_ratesR_path],
            targets = [mergers_file_path],                       
            args = [input_folder, output_folder],
            name = "dereplicate_and_merge"
            )
        return mergers_file

def const_seq_table(workflow, input_folder, output_folder, merged_file):
    
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

         merged_file_path = output_folder + merged_file
         read_counts_steps_file = "/Read_counts_at_each_step.tsv"
         read_counts_steps_path = output_folder + read_counts_steps_file
         seqtab_data_file = "/seqtab_final.rds"
         seqtab_data_file_path = output_folder + seqtab_data_file

         workflow.add_task(
            "biobakery_workflows/scripts/const_seq_table.R\
              --input_dir=[args[0]]\
              --output_dir=[args[1]]\
              --merged_file= " + merged_file + " --read_counts_steps=" + read_counts_steps_file,
            depends = [merged_file_path],
            targets = [read_counts_steps_path, seqtab_data_file_path],                              
            args = [input_folder, output_folder],
            name = "construct_sequence_table"
            )
         return seqtab_data_file, read_counts_steps_file


def phylogeny(workflow, output_folder, seqtab_data_file ):
    
         """ Aligns sequences and reconstructs phylogeny
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Path to sequence table data
                
            Returns:
                Creates file that contains pylogeny
         """

         seqtab_data_path = output_folder + seqtab_data_file
         msa_fasta_file = "/all_samples_clustalo_aligned_nonchimera.fasta"
         msa_fasta_path = output_folder + msa_fasta_file

         workflow.add_task(
            "biobakery_workflows/scripts/phylogeny.R\
              --output_dir=[args[0]]\
              --seqtab_file=" + seqtab_data_file,
            depends = [seqtab_data_path],
            targets = [msa_fasta_path],                              
            args = [output_folder],
            name = "phylogeny"
            )
         return msa_fasta_file

def fasttree(workflow, output_folder):
    
         """ Generates a phylogenetic tree
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Path to phylogeny file
                
            Returns:
                Creates file that contains pylogenic tree
         """

         msa_fasta_file = output_folder + "/all_samples_clustalo_aligned_nonchimera.fasta"
         fasttree_file = output_folder + "/closed_reference.tre"

         workflow.add_task(
            "FastTree -gtr -nt  [depends[0]] >  [targets[0]]",
            depends = [msa_fasta_file],
            targets = [fasttree_file],                              
            args = [output_folder],
            name = "fasttree"
            )


def assign_taxonomy(workflow, output_folder, gg_path):
    
         """ Assigns taxonomy using gg database
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Sequence table data
                
            Returns:
                Creates closed reference file
         """

         seqtab_data_path = output_folder + "/seqtab_final.rds"
         all_samples_taxonomy = output_folder + "/all_samples_GG13-8-taxonomy.tsv"
         all_samples_sv_gg13 = output_folder + "/all_samples_taxonomy_closed_reference.tsv"
 
         workflow.add_task(
            "biobakery_workflows/scripts/assign_taxonomy.R --output_dir=[args[0]] --gg_path=[args[1]]",
            depends = [seqtab_data_path],
            targets = [all_samples_taxonomy, all_samples_sv_gg13],                              
            args = [output_folder, gg_path],
            name = "assign_taxonomy"
            )
               

def assign_silva_rdp(workflow, output_folder, rdp_path, silva_path):
    
         """ Assigns taxonomy using silva and rdp database
            
            Args:
                workflow (anadama2.workflow): An instance of the workflow class.
                output_folder (string): The path of the output folder.
                
            Requires:
                Sequence table data
                
            Returns:
                Creates closed reference silva and rdp files
         """

         seqtab_data_path = output_folder + "/seqtab_final.rds"
         all_samples_silva = output_folder + "/all_samples_taxonomy_closed_reference_silva.tsv"
         all_samples_silva_rdp = output_folder + "/all_samples_taxonomy_closed_reference_rdp.tsv"
 
         workflow.add_task(
            "biobakery_workflows/scripts/assign_silva_rdp.R --output_dir=[args[0]] --rdp_path=[args[1]] --silva_path=[args[2]]",
            depends = [seqtab_data_path],
            targets = [all_samples_silva, all_samples_silva_rdp],                              
            args = [output_folder, rdp_path, silva_path],
            name = "assign_silva_rdp"
            
            )

   
