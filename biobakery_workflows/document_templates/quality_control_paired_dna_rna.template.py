
#' # Quality Control

#' This report section contains information about the quality control processing
#' for all samples. These samples were
#' run through [KneadData](http://huttenhower.sph.harvard.edu/kneaddata).
#' Samples were first trimmed then filtered using the human genome (hg38). RNA
#' samples where than filtered using the human transcriptome (hg38 mRNA).

#' Data is organized by paired and orphan reads. When 
#' one read in a pair passes a filtering step and the other does not the surviving
#' read is an orphan. The tables and plots are annotated as follows:

#' * raw: Untouched fastq reads.
#' * trim: Number of reads remaining after trimming bases with Phred score < 20. If the 
#' trimmed reads is <70% of original length then it is removed altogether.
#' * hg38: Number of reads remaining after depleting reads against the human genome (hg38).
#' * mRNA: Number of reads remaining after depleting reads against the human genome (hg38)
#' and the human transcriptome (hg38 mRNA). (RNA samples only)

#+ echo=False

# read in the DNA samples
columns, dna_samples, dna_paired_data = document.read_table(vars["dna_read_counts"], only_data_columns=(0,2,6))
dna_paired_columns = ["Raw","Trim","hg38"]

columns, dna_samples, dna_orphan_data = document.read_table(vars["dna_read_counts"], only_data_columns=(4,5,8,9))
dna_orphan_columns = ["Trim orphan1", "Trim orphan2", "hg38 orphan1", "hg38 orphan2"]

# read in the RNA samples
columns, rna_samples, rna_paired_data = document.read_table(vars["rna_read_counts"], only_data_columns=(0,2,6,14))
rna_paired_columns = ["Raw","Trim","hg38","hg38 mRNA"]

columns, rna_samples, rna_orphan_data = document.read_table(vars["rna_read_counts"], only_data_columns=(4,5,10,11,16,17))
rna_orphan_columns = ["Trim orphan1", "Trim orphan2", "hg38 orphan1", "hg38 orphan2", 
    "mRNA orphan1", "mRNA orphan2"]

#' ## DNA Samples Quality Control

#' ### DNA Samples Tables of Filtered Reads

#+ echo=False
document.show_table(dna_paired_data, dna_samples, dna_paired_columns, "DNA Paired end reads", 
    format_data_comma=True)

#+ echo=False
document.show_table(dna_orphan_data, dna_samples, dna_orphan_columns, "DNA Orphan reads", 
    format_data_comma=True)

#+ echo=False
def microbial_read_proportion(paired_data, orphan_data):
    # compute the proportion of reads relative to the trim/raw reads
    proportion_decontaminated = []
    for paired_row, orphan_row in zip(paired_data, orphan_data):
        decontaminated_sum = 2.0 * paired_row[-1] + orphan_row[-1] + orphan_row[-2]
        decon_trim = decontaminated_sum / (2.0 * paired_row[1] + orphan_row[0] + orphan_row[1])
        decon_raw = decontaminated_sum / (2.0 * paired_row[0])
        proportion_decontaminated.append(["{0:.5f}".format(i) for i in [decon_trim, decon_raw]])
    return proportion_decontaminated
    
dna_microbial_read_proportion = microbial_read_proportion(dna_paired_data, dna_orphan_data)
document.show_table(dna_microbial_read_proportion, dna_samples, ["hg38 / Trim","hg38 / Raw"],
                    "DNA microbial read proportion")

#' ### DNA Samples Plots of Filtered Reads

#+ echo=False
import numpy

document.plot_grouped_barchart(numpy.transpose(dna_paired_data), row_labels=dna_paired_columns, 
    column_labels=dna_samples, title="DNA Paired end reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(dna_orphan_data), row_labels=dna_orphan_columns, 
    column_labels=dna_samples, title="DNA Orphan reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)

#' ## RNA Samples Quality Control

#' ### RNA Samples Tables of Filtered Reads

#+ echo=False
document.show_table(rna_paired_data, rna_samples, rna_paired_columns, "RNA Paired end reads", 
    format_data_comma=True)

#+ echo=False
document.show_table(rna_orphan_data, rna_samples, rna_orphan_columns, "RNA Orphan reads", 
    format_data_comma=True, column_width=0.10)

#+ echo=False
rna_microbial_read_proportion = microbial_read_proportion(rna_paired_data, rna_orphan_data)
document.show_table(rna_microbial_read_proportion, rna_samples, ["hg38 mRNA / Trim","hg38 mRNA / Raw"],
                    "RNA microbial read proportion")

#' ### RNA Samples Plots of Filtered Reads

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(rna_paired_data), row_labels=rna_paired_columns, 
    column_labels=rna_samples, title="RNA Paired end reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)

#+ echo=False
document.plot_grouped_barchart(numpy.transpose(rna_orphan_data), row_labels=rna_orphan_columns, 
    column_labels=rna_samples, title="RNA Orphan reads", ylabel="Read count (in millions)",
    legend_title="Filter", yaxis_in_millions=True)




