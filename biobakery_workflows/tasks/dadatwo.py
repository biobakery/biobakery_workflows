"""
bioBakery Workflows: tasks.dadatwo module
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
from anadama2.tracked import TrackedDirectory, TrackedExecutable
from biobakery_workflows import files, config, utilities
import os,fnmatch

def identify_primers(workflow,fwd_primer,rev_primer,input_folder,output_folder,pair_id):
    """ Identifies primers and N filters samples
       Args:
           workflow (anadama2.workflow): an instance of the workflow class
           input_folder (string): path to input folder
           output_folder (string):  path to output folder
           fwd_primer (string): forward primer
           rev_primer (string): reverse primer
           pair_id (string): pair identifier

       Requires:
          dada2, ShortRead, Biostrings, tools r packages

       Returns:
           string: path to folder with N filtered files
           string: path to file with fwd_primer and its reverse
           string: path to file with rev_primer and its reverse
    """
    script_path = utilities.get_package_file("identify_primers", "Rscript")
    filtN_folder = os.path.join(output_folder,"filtN")
    primers_folder = os.path.join(output_folder,"primers")
    fwd_primer_file = os.path.join(primers_folder,"fwd_primer_file.txt")
    rev_primer_file = os.path.join(primers_folder,"rev_primer_file.txt")

    workflow.add_task(
        "[vars[0]]  \
          --input_dir=[args[3]] \
          --output_dir=[vars[1]] \
          --primers_dir=[vars[2]] \
          --fwd_primer_file=[targets[0]] \
          --rev_primer_file=[targets[1]] \
          --fwd_primer=[args[0]] \
          --rev_primer=[args[1]] \
          --pair_id=[args[2]]",
        depends=[TrackedDirectory(input_folder)],
        targets=[fwd_primer_file, rev_primer_file,TrackedDirectory(filtN_folder)],
        args=[fwd_primer, rev_primer, pair_id,input_folder],
        vars=[script_path,filtN_folder,primers_folder],
        name="identify_primers"
    )

    return filtN_folder,fwd_primer_file,rev_primer_file

def remove_primers(workflow,input_folder,fwd_primer_file,rev_primer_file, pair_id):
    """ Removes primers from samples
           Args:
               workflow (anadama2.workflow): an instance of the workflow class
               input_folder (string):  path to folder with input files
               fwd_primer_file (string): path to file with fwd primer and its reverse
               rev_primer_file (string): path to file with rev primer and its reverse
               pair_id (string): pair identifier

           Requires:
              cutadapt

           Returns:
               string: path to folder with cleaned from primers files
    """
    with open(fwd_primer_file) as f:
        FWD = f.read().splitlines()
    with open(rev_primer_file) as f:
        REV = f.read().splitlines()
    pair_id2 = pair_id.replace("1", "2")
    fwd_reads=fnmatch.filter(os.listdir(input_folder), "*"+pair_id+"*.fastq*")
    rev_reads=fnmatch.filter(os.listdir(input_folder),"*"+ pair_id2+"*.fastq*")
    fwd_reads_inp = sorted([os.path.join(input_folder,f) for f in fwd_reads])
    rev_reads_inp = sorted([os.path.join(input_folder, f) for f in rev_reads])
    cutadapt_folder = os.path.join(input_folder, "cutadapt")
    if not os.path.exists(cutadapt_folder):
        os.mkdir(cutadapt_folder)
    fwd_reads_out = sorted([os.path.join(cutadapt_folder,f) for f in fwd_reads])
    rev_reads_out = sorted([os.path.join(cutadapt_folder, f) for f in rev_reads])

    for i in range(0,len(rev_reads_inp)):
        workflow.add_task(
            "cutadapt -g [vars[0]] -a [vars[3]] -G [vars[2]] -A [vars[1]] -n 2 \
             -o [targets[0]] -p [targets[1]] [vars[4]] [vars[5]]",
            depends=[fwd_primer_file,rev_primer_file],
            targets=[fwd_reads_out[i], rev_reads_out[i]],
            vars=[FWD[0],FWD[1], REV[0],REV[1],fwd_reads_inp[i], rev_reads_inp[i]],
            name="remove_primers"
        )


    return cutadapt_folder

def filter_trim(workflow, input_folder, output_folder, maxee, trunc_len_max, pair_id,threads):
    
         """ Filters samples by maxee and trims them, renders quality control plots
         of forward and reverse reads for each sample, creates read counts tsv and rds files.
            
            Args:
                workflow (anadama2.workflow): an instance of the workflow class
                input_folder (string): path to input folder
                output_folder (string):  path to output folder
                maxee (string): maxee value to use for filtering
                trunc_len_max (string): max length for truncating
                pair_id (string): pair identifier
                threads (int): number of threads
                
            Requires:
               dada2, gridExtra,tools r packages
                
            Returns:
                string: path to file that contains read counts before and after filtering
                string: path to folder with filtered and trimmed sample files
         """
         reads_plotF_png = files.SixteenS.path("readF_qc", output_folder)
         reads_plotR_png = files.SixteenS.path("readR_qc", output_folder)

         readcounts_tsv_path = os.path.join(output_folder, "Read_counts_after_filtering.tsv")
         readcounts_rds_path = os.path.join(output_folder, "Read_counts_filt.rds")
         filtered_dir = "filtered_input"
         script_path = utilities.get_package_file("filter_and_trim", "Rscript")
         workflow.add_task(
             "[vars[0]] \
               --input_dir=[args[0]]\
               --output_dir=[args[1]]\
               --filtered_dir=[vars[1]]\
               --maxee=[args[2]]\
               --trunc_len_max=[args[3]]\
               --readcounts_tsv_path=[targets[0]]\
               --readcounts_rds_path=[targets[1]]\
               --reads_plotF=[targets[2]]\
               --reads_plotR=[targets[3]]\
               --pair_id=[args[4]]\
               --threads=[args[5]]",
             depends =[TrackedDirectory(input_folder)],
             targets = [readcounts_tsv_path, readcounts_rds_path, reads_plotF_png, reads_plotR_png],
             args = [input_folder, output_folder, maxee, trunc_len_max, pair_id, threads],
             vars = [script_path,filtered_dir],
             name ="filter_and_trim"
             )
         return readcounts_tsv_path, filtered_dir
     

def learn_error(workflow, output_folder, filtered_dir, readcounts_tsv_path, threads):
    
         """ Learns error rates for each sample, renders error rates plots for forward and reverse reads
            
            Args:
                workflow (anadama2.workflow): an instance of the workflow class
                output_folder (string): path to output folder
                filtered_dir (string): path to directory with filtered files
                readcounts_tsv_path (string): path to read counts after filtering tsv file
                threads (int): number of threads

            Requires:
                dada2, ggplot2 r packages

            Returns:
                string: path to file that contains error rates of forward reads
                string: path to file that contains error rates of reverse reads
         """

         error_ratesF_png = files.SixteenS.path("error_ratesF", output_folder)
         error_ratesR_png = files.SixteenS.path("error_ratesR", output_folder)
         
         error_ratesF_path= os.path.join(output_folder, "error_ratesFWD.rds")
         error_ratesR_path =os.path.join(output_folder, "error_ratesREV.rds")

         script_path = utilities.get_package_file("learn_error_rates", "Rscript")
       
         workflow.add_task(
             "[vars[0]] \
               --output_dir=[args[0]]\
               --filtered_dir=[args[1]]\
               --error_ratesF_png=[targets[0]]\
               --error_ratesR_png=[targets[1]]\
               --error_ratesF_path=[targets[2]]\
               --error_ratesR_path=[targets[3]]\
               --threads=[vars[1]]",
             depends = [readcounts_tsv_path],
             targets = [error_ratesF_png, error_ratesR_png, error_ratesF_path, error_ratesR_path],  
             args = [output_folder, filtered_dir],
             vars = [script_path, threads],
             name = "learn_error_rates"
             )
         return error_ratesF_path, error_ratesR_path
     

def merge_paired_ends(workflow, output_dir, filtered_dir, error_ratesF_path, error_ratesR_path, threads):
    
        """ Dereplicates and merges paired reads
            
            Args:
                workflow (anadama2.workflow): an instance of the workflow class
                output_folder (string): path to output folder
                filtered_dir (string): path to directory with filtered files
                error_ratesF_path (string): path to rds file that contains error rates of forward reads
                error_ratesR_path (string): path to rds file that contains error rates of reverse reads
                threads (int): number fo threads

            Requires:
                dada2, tools r packages
                
            Returns:
                string: path to rds file that contains merged and dereplicated reads
         """

        mergers_file_path = os.path.join(output_dir, "mergers.rds")
        script_path = utilities.get_package_file("merge_paired_ends", "Rscript")
        
        workflow.add_task(
            "[vars[0]] \
              --output_dir=[args[0]]\
              --filtered_dir=[args[1]]\
              --error_ratesF_path=[depends[0]]\
              --error_ratesR_path=[depends[1]]\
              --mergers_file_path=[targets[0]]\
              --threads=[vars[1]]",
            depends = [error_ratesF_path, error_ratesR_path],
            targets = [mergers_file_path],                       
            args = [output_dir, filtered_dir],
            vars = [script_path, threads],
            name = "dereplicate_and_merge"
            )
        return mergers_file_path
    

def const_seq_table(workflow, output_folder, filtered_dir,  mergers_file_path, threads):
    
         """ Builds ASV table, removes chimeras, creates read counts at each step, and fasta file with all sequences
            
            Args:
                workflow (anadama2.workflow): an instance of the workflow class
                output_folder (string):  path to output folder
                filtered_dir (string): path to directory with filtered files
                mergers_file_path (string): path to rds file that contains merged reads
                threads (int): number of threads
                
            Requires:
                dada2, tools, seqinr r packages

            Returns:
                string: path to rds file that contains ASV data
                string: path to read counts at each step tsv file
                string: path to fasta file with all sequences
         """
         
         read_counts_steps_path = files.SixteenS.path("counts_each_step", output_folder)
         
         seqtab_file_path = os.path.join(output_folder, "seqtab_final.rds")
         seqs_fasta_path = os.path.join(output_folder, "sequences.fasta")
         readcounts_rds = "Read_counts_filt.rds"
         asv_tsv = "all_samples_SV_counts.tsv"

         script_path = utilities.get_package_file("const_seq_table", "Rscript")
         version_script = utilities.get_package_file("dada2_version", "Rscript")

         version_command = """echo 'r' `r -e 'packageVersion("dada2")' | grep -C 1 dada2`"""

         workflow.add_task(
            "[vars[0]] \
              --output_dir=[args[0]]\
              --filtered_dir=[args[1]]\
              --merged_file_path=[depends[0]]\
              --read_counts_steps_path=[targets[0]]\
              --readcounts_rds=[vars[2]]\
              --asv_tsv=[vars[3]]\
              --seqtab_file_path=[targets[1]]\
              --seqs_fasta_path=[targets[2]]\
              --threads=[vars[1]]",
            depends = [mergers_file_path,TrackedExecutable("R", version_command="echo '" +  version_script + "' `" + version_script + "`")],
            targets = [read_counts_steps_path, seqtab_file_path, seqs_fasta_path],
            args = [output_folder, filtered_dir],
            vars = [script_path, threads, readcounts_rds, asv_tsv ],
            name = "construct_sequence_table"
            )
         return seqtab_file_path, read_counts_steps_path, seqs_fasta_path


def assign_taxonomy(workflow, output_folder, seqtab_file_path, ref_path, threads):
    
         """ Assigns taxonomy using green genes, silva, or rdp database, creates closed reference file
            
            Args:
                workflow (anadama2.workflow): an instance of the workflow class
                output_folder (string): path to output folder
                seqtab_file_path (string): path to rds file that contains ASV data
                ref_path (string): reference database name
                threads (int):
                
            Requires:
                dada2 r package
                
            Returns:
                string: path to closed reference file
         """

         otu_closed_ref_path  = files.SixteenS.path("otu_table_closed_reference", output_folder)

         # check what reference db to use for taxonomy assignment
         if ref_path == "unite":
             refdb_path = config.SixteenS().unite
             refdb_species_path = "None"
         elif ref_path == "silva":
             refdb_path = config.SixteenS().silva_dada2 
             refdb_species_path = config.SixteenS().silva_species_dada2 
         elif ref_path == "rdp":
             refdb_path = config.SixteenS().rdp_dada2 
             refdb_species_path = config.SixteenS().rdp_species_dada2  
         else:    
             refdb_path = config.SixteenS().greengenes_dada2
             refdb_species_path = "None"

         script_path = utilities.get_package_file("assign_taxonomy", "Rscript")
             
         workflow.add_task(
            "[vars[2]] \
              --output_dir=[args[0]]\
              --refdb_path=[vars[0]]\
              --refdb_species_path=[vars[1]]\
              --seqtab_file_path=[depends[0]]\
              --otu_closed_ref_path=[targets[0]]\
              --threads=[vars[3]]",
            depends = [seqtab_file_path],
            targets = [otu_closed_ref_path],                              
            args = [output_folder],
            vars =[refdb_path, refdb_species_path, script_path, threads],
            name = "assign_taxonomy"
            )         
     
         return otu_closed_ref_path
 
 
def remove_tmp_files(workflow, output_folder, otu_closed_ref_path, msa_fasta_path, fasttree_file_path):
    
         """ Removes temporary rds files
            
            Args:
                workflow (anadama2.workflow): an instance of the workflow class
                output_folder (string): path to output folder.
                otu_closed_ref_path (string): path to closed reference file
                msa_fasta_path (string): path to msa file
                fasttree_file_path (string): path to phylogenetic tree file

            Requires:
               None
                
            Returns:
               None
          """
         
         rm_out_file = os.path.join(output_folder, "tmp_rm.txt")
             
         workflow.add_task(
                     "rm  [args[0]]/*.rds  &>[targets[0]] ",
                     depends = [otu_closed_ref_path, msa_fasta_path, fasttree_file_path],
                     args = [output_folder],
                     targets = [rm_out_file],
                     name = "rm_tmp_files"
                     )
         
