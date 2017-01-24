
#+ echo=False
import numpy

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

#' # Quality Control

#' This report section contains information about the quality control processing
#' for all samples. These samples were
#' run through [KneadData](http://huttenhower.sph.harvard.edu/kneaddata).
#' Samples were first trimmed then filtered using the human genome (hg38).

#' Data is organized by paired and orphan reads. When 
#' one read in a pair passes a filtering step and the other does not the surviving
#' read is an orphan. The tables and plots are annotated as follows:

#' * raw: Untouched fastq reads.
#' * trim: Number of reads remaining after trimming bases with Phred score < 20. If the 
#' trimmed reads is <70% of original length then it is removed altogether.
#' * hg38: Number of reads remaining after depleting reads against the human genome (hg38).

#+ echo=False

# read in the DNA samples
columns, dna_samples, dna_paired_data = document.read_table(vars["dna_read_counts"], only_data_columns=(0,2,6))
dna_paired_columns = ["Raw","Trim","hg38"]

columns, dna_samples, dna_orphan_data = document.read_table(vars["dna_read_counts"], only_data_columns=(4,5,8,9))
dna_orphan_columns = ["Trim orphan1", "Trim orphan2", "hg38 orphan1", "hg38 orphan2"]

#' ## DNA Samples Quality Control

#' ### DNA Samples Tables of Filtered Reads

#+ echo=False
document.show_table(dna_paired_data, dna_samples, dna_paired_columns, "DNA Paired end reads", 
    format_data_comma=True)

#+ echo=False
document.show_table(dna_orphan_data, dna_samples, dna_orphan_columns, "DNA Orphan reads", 
    format_data_comma=True)

#+ echo=False
# plot the microbial reads ratios    
dna_microbial_reads, dna_microbial_labels = utilities.microbial_read_proportion(dna_paired_data, dna_orphan_data)
document.show_table(dna_microbial_reads, dna_samples, dna_microbial_labels,
                    "DNA microbial read proportion", column_width=0.25)

#' ### DNA Samples Plots of Filtered Reads

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(dna_paired_data), row_labels=dna_paired_columns, 
    column_labels=dna_samples, title="DNA Paired end reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(dna_orphan_data), row_labels=dna_orphan_columns, 
    column_labels=dna_samples, title="DNA Orphan reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)





