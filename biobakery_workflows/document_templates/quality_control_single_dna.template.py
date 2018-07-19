
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

# read in the DNA samples
dna_columns, dna_samples, dna_data = visualizations.qc_read_counts(document, vars["dna_read_counts"])

#' # Quality Control

#' <% visualizations.ShotGun.print_qc_intro_caption(len(dna_samples), dna_columns[2:]) %>

#' ## DNA Samples Quality Control

#' ### DNA Samples Tables of Filtered Reads

#+ echo=False

# create a table of the paired counts
document.write_table(["# Sample"]+dna_columns, dna_samples, dna_data,
    files.ShotGunVis.path("qc_counts",document.data_folder))

table_message=visualizations.show_table_max_rows(document, dna_data, dna_samples,
    dna_columns, "DNA reads", files.ShotGunVis.path("qc_counts"), format_data_comma=True)
        
#' <%= table_message %>

#+ echo=False
# compute and plot the microbial reads ratios
# compute ratios for each database used for qc
dna_microbial_reads, dna_microbial_labels = utilities.microbial_read_proportion_multiple_databases(dna_data, dna_columns)

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





