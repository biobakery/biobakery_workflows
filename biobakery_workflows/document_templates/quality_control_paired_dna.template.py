
#+ echo=False
import os
import numpy

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

# set the max number of table rows to be shown
# and set the messages to print for each table size
max_table_rows=20
table_message="A data file exists of this table: "
large_table_message="The table is too large to include the full table in this document."+\
    " A partial table is shown which includes only "+str(max_table_rows)+" samples."+\
    " Please see the data file for the full table: "

# get the name of the contaminate database used
db_name = vars["contaminate_database"]

#' # Quality Control

#' This report section contains information about the quality control processing
#' for all samples. These samples were
#' run through [KneadData](http://huttenhower.sph.harvard.edu/kneaddata).
#' Samples were first trimmed then filtered against a contaminate reference database.

#' Data is organized by paired and orphan reads. When 
#' one read in a pair passes a filtering step and the other does not the surviving
#' read is an orphan. The tables and plots are annotated as follows:

#' * raw : Untouched fastq reads.
#' * trim : Number of reads remaining after trimming bases with Phred score < 20. If the 
#' trimmed reads is <70% of original length then it is removed altogether.
#' * <% print(db_name) %> : Number of reads remaining after depleting reads against a contaminate reference database.

#+ echo=False

# read in the DNA samples
columns, dna_samples, dna_paired_data = document.read_table(vars["dna_read_counts"], only_data_columns=(0,2,6), format_data=int)
dna_paired_columns = ["Raw","Trim",db_name]

columns, dna_samples, dna_orphan_data = document.read_table(vars["dna_read_counts"], only_data_columns=(4,5,8,9), format_data=int)
dna_orphan_columns = ["Trim orphan1", "Trim orphan2", db_name+" orphan1", db_name+" orphan2"]

#' ## DNA Samples Quality Control

#' ### DNA Samples Tables of Filtered Reads

#+ echo=False

# create a table of the paired counts
paired_counts_file = os.path.join(document.data_folder,"paired_counts_table.tsv")
document.write_table(["# Sample"]+dna_paired_columns, dna_samples, dna_paired_data, paired_counts_file)

if len(dna_samples) <= max_table_rows:
    document.show_table(dna_paired_data, dna_samples, dna_paired_columns, "DNA Paired end reads", 
        format_data_comma=True)
else:
    document.show_table(dna_paired_data[:max_table_rows], dna_samples[:max_table_rows], 
        dna_paired_columns, "DNA Paired end reads (partial table)", format_data_comma=True)
        
#' <% print(large_table_message) if len(dna_samples) > max_table_rows else print(table_message) %>
#' [paired_counts_table.tsv](data/paired_counts_table.tsv)    

#+ echo=False

# create a table of the orphan counts
orphan_counts_file = os.path.join(document.data_folder,"orphan_counts_table.tsv")
document.write_table(["# Sample"]+dna_orphan_columns, dna_samples, dna_orphan_data, orphan_counts_file)

if len(dna_samples) <= max_table_rows:
    document.show_table(dna_orphan_data, dna_samples, dna_orphan_columns, "DNA Orphan reads", 
        format_data_comma=True)
else:
    document.show_table(dna_orphan_data[:max_table_rows], dna_samples[:max_table_rows], 
        dna_orphan_columns, "DNA Orphan reads (partial table)", format_data_comma=True)    

#' <% print(large_table_message) if len(dna_samples) > max_table_rows else print(table_message) %>
#' [orphan_counts_table.tsv](data/orphan_counts_table.tsv) 

#+ echo=False
# plot the microbial reads ratios    
dna_microbial_reads, dna_microbial_labels = utilities.microbial_read_proportion(dna_paired_data, dna_orphan_data, database_name=db_name)

# create a table of the microbial reads
microbial_counts_file = os.path.join(document.data_folder,"microbial_counts_table.tsv")
document.write_table(["# Sample"]+dna_microbial_labels, dna_samples, dna_microbial_reads, microbial_counts_file)

if len(dna_samples) <= max_table_rows:
    document.show_table(dna_microbial_reads, dna_samples, dna_microbial_labels,
        "DNA microbial read proportion")
else:
    document.show_table(dna_microbial_reads[:max_table_rows], dna_samples[:max_table_rows], 
        dna_microbial_labels, "DNA microbial read proportion (partial table)")    
        
#' Proportion of reads remaining after removing host reads relative to the number of: i) quality-trimmed reads, and ii) raw unfiltered reads.  
        
#' <% print(large_table_message) if len(dna_samples) > max_table_rows else print(table_message) %>
#' [microbial_counts_table.tsv](data/microbial_counts_table.tsv) 
        
#' ### DNA Samples Plots of Filtered Reads

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(dna_paired_data), row_labels=dna_paired_columns, 
    column_labels=dna_samples, title="DNA Paired end reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(dna_orphan_data), row_labels=dna_orphan_columns, 
    column_labels=dna_samples, title="DNA Orphan reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)





