
#+ echo=False
import os
import numpy

from biobakery_workflows import utilities, visualizations, files

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False
    
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
document.write_table(["# Sample"]+dna_columns, dna_samples, dna_data,
    files.ShotGunVis.path("qc_counts",document.data_folder))

table_message=visualizations.show_table_max_rows(document, dna_data, dna_samples,
    dna_columns, "DNA reads", files.ShotGunVis.path("qc_counts"), format_data_comma=True)
        
#' <%= table_message %>

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
# plot the microbial reads ratios    
dna_microbial_reads, dna_microbial_labels = utilities.microbial_read_proportion(dna_data, database_name=db_name)

# create a table of the microbial reads
document.write_table(["# Sample"]+dna_microbial_labels, dna_samples, 
    dna_microbial_reads, files.ShotGunVis.path("microbial_counts",document.data_folder))

table_message=visualizations.show_table_max_rows(document, dna_microbial_reads, dna_samples,
    dna_microbial_labels, "DNA microbial read proportion",
    files.ShotGunVis.path("microbial_counts"))  
    
#' <%= visualizations.ShotGun.captions["microbial_ratios"] %>  
        
#' <%= table_message %>
        
#' ### DNA Samples Plots of Filtered Reads

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(dna_data), row_labels=dna_columns, 
    column_labels=dna_samples, title="DNA reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)





