
#+ echo=False
import numpy

min_abundance=0.01
min_samples=10

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

# read in the read count table
# columns expected are total reads, reads that map to OTUs with taxonomy,
# and reads that map to OTUs without taxonomy
columns, samples, data = document.read_table(vars["read_count_table"])

#' There were a total of <% print(len(samples))%> samples processed with this workflow for this project.

#' <% if pdf_format: print("\clearpage") %>

#' # Quality Control

#+ echo=False
# read the eestats2 table and add table to document
eestats_rows, eestats_columns, eestats_data, overall_stats = utilities.read_eestats2(vars["eestats_table"])

document.show_table(eestats_data, eestats_rows, eestats_columns, 
    "Expected error filter by read length",font="10")

#' The general stats for this data set are: <%= overall_stats %> .
#' This table shows the number of reads based on length for different error filters.

#' # Read Count

#+ echo=False

# sort the samples/data by read count with the largest original read count first
total_reads=[row[0] for row in data]
sorted_samples, sorted_total_reads = utilities.sort_data(total_reads, samples)
sorted_all_read_data = [data[samples.index(sample)] for sample in sorted_samples]

known_reads = [row[1] for row in sorted_all_read_data]
unknown_reads = [row[2] for row in sorted_all_read_data]
unmapped_reads = [row[0]-(row[1]+row[2]) for row in sorted_all_read_data]

# plot the read counts
document.plot_stacked_barchart([known_reads,unknown_reads,unmapped_reads], ["classified","unclassified","unmapped"], sorted_samples, 
    title="Read counts by Sample", ylabel="Total Reads", xlabel="Samples")

#' This figure shows counts of reads in three categories: 1) classified: reads that align to OTUs with known taxonomy,
#' 2) reads that align to OTUs of unknown taxonomy, 3) reads that do not align to any OTUs. The sum of these
#' three read counts for each sample is the total original read count not including filtering prior to OTU clustering.

#' <% if pdf_format: print("\clearpage") %>

#' # Taxonomy

#+ echo=False
import numpy

# read in the otu table data
samples, ids, taxonomy, data = utilities.read_otu_table(vars["otu_table"])

# plot the top taxa by genus level, plotting the relative abundance values
max_taxa=15
# get the relative abundance values for the samples
relab_data = utilities.relative_abundance(data)
# get the taxa summarized by genus level
genus_level_taxa, genus_level_data = utilities.taxa_by_level(taxonomy, relab_data, level=5)
# get the top rows of the relative abundance data
top_taxa, top_data = utilities.top_rows(genus_level_taxa, genus_level_data, max_taxa, function="average")
# shorten the top taxa names to just the genus level for plotting
top_taxa_short_names = utilities.taxa_shorten_name(top_taxa, level=5, remove_identifier=True)

# sort the data so those with the top genera are shown first
sorted_samples, sorted_data = utilities.sort_data(top_data[0], samples)
transpose_top_data = numpy.transpose(top_data)
sorted_top_data = numpy.transpose([transpose_top_data[samples.index(sample)] for sample in sorted_samples])

document.plot_stacked_barchart(sorted_top_data, row_labels=top_taxa_short_names, 
    column_labels=sorted_samples, title="Top "+str(max_taxa)+" genera by average abundance",
    ylabel="Relative abundance", legend_title="Genera")

#+ echo=False

# plot the relative abundance of the top terminal taxa
# get the terminal taxa
terminal_taxa_relab, terminal_data_relab = utilities.terminal_taxa(taxonomy, relab_data)
# get the top rows of terminal taxa
top_terminal_taxa, top_terminal_data = utilities.top_rows(terminal_taxa_relab, terminal_data_relab, max_taxa, function="average")

# reduce the taxa names to just the most specific identifier
shorted_names=utilities.taxonomy_trim(top_terminal_taxa)

# sort the data with the samples with the top terminal taxa first
sorted_samples_terminal, sorted_data_terminal = utilities.sort_data(top_terminal_data[0], samples)
transpose_top_terminal_data = numpy.transpose(top_terminal_data)
sorted_top_terminal_data = numpy.transpose([transpose_top_terminal_data[samples.index(sample)] for sample in sorted_samples_terminal])

document.plot_stacked_barchart(sorted_top_terminal_data, row_labels=shorted_names, 
    column_labels=sorted_samples_terminal, title="Top "+str(max_taxa)+" terminal taxa by average abundance",
    ylabel="Relative abundance", legend_title="Terminal taxa")

#' # Ordination

#+ echo=False

# plot the top terminal node taxa in a PCOA
# provide data as values [0-1] organized as samples as columns and features as rows

# filter the data by min abundance and min samples
filtered_taxonomy, filtered_data = utilities.filter_taxa(terminal_taxa_relab, terminal_data_relab, min_abundance, min_samples)

document.show_pcoa(samples, filtered_taxonomy, filtered_data, title="PCOA Ordination of terminal taxa using Bray-Curtis similarity")

#' For the PCoA plot, relative abundances are passed through a basic filter requiring each terminal taxa
#' to have at least <% print(min_abundance)%> % abundance in at least 
#' <% print(min_samples) %> % of all samples.





