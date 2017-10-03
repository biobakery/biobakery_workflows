
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
(dna_paired_columns, dna_orphan_columns), dna_samples, (dna_paired_data, dna_orphan_data) = visualizations.qc_read_counts(document, vars["dna_read_counts"])


#' # Quality Control

#' <% visualizations.ShotGun.print_qc_intro_caption(len(dna_samples), dna_paired_columns[2:], paired=True) %>

#' ## DNA Samples Quality Control

#' ### DNA Samples Tables of Filtered Reads

#+ echo=False

# create a table of the paired counts
document.write_table(["# Sample"]+dna_paired_columns, dna_samples, dna_paired_data,
    files.ShotGunVis.path("qc_counts_paired",document.data_folder))

table_message=visualizations.show_table_max_rows(document, dna_paired_data, dna_samples,
    dna_paired_columns, "DNA Paired end reads", files.ShotGunVis.path("qc_counts_paired"),
    format_data_comma=True)
        
#' <%= table_message %>

#+ echo=False

# create a table of the orphan counts
document.write_table(["# Sample"]+dna_orphan_columns, dna_samples, dna_orphan_data,
    files.ShotGunVis.path("qc_counts_orphan",document.data_folder))

table_message=visualizations.show_table_max_rows(document, dna_orphan_data, dna_samples,
    dna_orphan_columns, "DNA Orphan reads", files.ShotGunVis.path("qc_counts_orphan"),
    format_data_comma=True)
        
#' <%= table_message %>  

#+ echo=False
# plot the microbial reads ratios    
dna_microbial_reads, dna_microbial_labels = utilities.microbial_read_proportion_multiple_databases(
    dna_paired_data, dna_paired_columns, dna_orphan_data)

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
# sort the samples/data by read count with the largest original read count first
def sort_samples_reads_decreasing(read_data, read_samples):
    """ Sort the reads from largest to smallest total read count """
        
    sorted_samples, sorted_total_reads = utilities.sort_data(read_data[0], read_samples)
    sorted_all_read_data = []
    for data_set in read_data:
        sorted_all_read_data.append([data_set[read_samples.index(sample)] for sample in sorted_samples])
    
    return sorted_samples, sorted_all_read_data

sorted_samples, sorted_all_read_data = sort_samples_reads_decreasing(numpy.transpose(dna_paired_data), dna_samples)

document.plot_grouped_barchart(sorted_all_read_data, row_labels=dna_paired_columns, 
    column_labels=sorted_samples, title="DNA Paired end reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)

#+ echo=False
sorted_samples, sorted_all_read_data = sort_samples_reads_decreasing(numpy.transpose(dna_orphan_data), dna_samples)

document.plot_grouped_barchart(sorted_all_read_data, row_labels=dna_orphan_columns, 
    column_labels=sorted_samples, title="DNA Orphan reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)





