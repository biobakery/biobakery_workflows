
#+ echo=False
import os
import numpy

from biobakery_workflows import utilities
from biobakery_workflows import visualizations

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

# read in the DNA samples
columns, dna_samples, dna_data = document.read_table(vars["dna_read_counts"], only_data_columns=(0,1,3), format_data=int)
dna_columns = ["Raw","Trim",db_name]


#' # Quality Control

#' 
#' <%= visualizations.ShotGun.format_caption("qc_intro",total_samples=len(dna_samples),seq_type="single-end") %>
#' Reads were first trimmed then filtered against a contaminate reference database.

#' * raw : Untouched fastq reads.
#' * trim : Number of reads remaining after trimming bases with Phred score < 20. If the 
#' trimmed reads is <70% of original length then it is removed altogether.
#' * <% print(db_name) %> : Number of reads remaining after depleting reads against the contaminate reference database.

#' ## DNA Samples Quality Control

#' ### DNA Samples Tables of Filtered Reads

#+ echo=False

# create a table of the paired counts
qc_counts_file = os.path.join(document.data_folder,"qc_counts_table.tsv")
document.write_table(["# Sample"]+dna_columns, dna_samples, dna_data, qc_counts_file)

if len(dna_samples) <= max_table_rows:
    document.show_table(dna_data, dna_samples, dna_columns, "DNA reads", 
        format_data_comma=True)
else:
    document.show_table(dna_data[:max_table_rows], dna_samples[:max_table_rows], 
        dna_columns, "DNA reads (partial table)", format_data_comma=True)
        
#' <% print(large_table_message) if len(dna_samples) > max_table_rows else print(table_message) %>
#' [qc_counts_table.tsv](data/qc_counts_table.tsv)

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
# plot the microbial reads ratios    
dna_microbial_reads, dna_microbial_labels = utilities.microbial_read_proportion(dna_data, database_name=db_name)

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
document.plot_grouped_barchart(numpy.transpose(dna_data), row_labels=dna_columns, 
    column_labels=dna_samples, title="DNA reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)





