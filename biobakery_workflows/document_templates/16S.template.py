
#+ echo=False
import numpy

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

#' # Read Count

#+ echo=False

# read in the read count table
columns, samples, data = document.read_table(vars["read_count_table"])

# sort the samples/data by read count with the largest first
sorted_samples, sorted_data = utilities.sort_data(data, samples)

# plot the read counts
document.plot_barchart(sorted_data, sorted_samples, title="Read counts by Sample",
    ylabel="Total Reads", xlabel="Samples")

#' # Taxonomy

#+ echo=False

# read in the otu table data
samples, ids, taxonomy, data = utilities.read_otu_table(vars["otu_table"])

# create a plot of total counts from the otu table
# get the total number of counts for each sample
counts=[sum(column) for column in numpy.transpose(data)]

# order the counts to put the samples with the largest counts first in the plot
sorted_samples, sorted_counts = utilities.sort_data(counts,samples)

document.plot_barchart(sorted_counts, sorted_samples, title="OTU Counts by Sample",
    ylabel="Counts", xlabel="Samples")

#+ echo=False

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

document.plot_stacked_barchart(top_data, row_labels=top_taxa_short_names, 
    column_labels=samples, title="Top "+str(max_taxa)+" genera by average abundance",
    ylabel="Relative abundance", legend_title="Genera")

#+ echo=False

# plot the relative abundance of the top terminal taxa
# get the terminal taxa
terminal_taxa_relab, terminal_data_relab = utilities.terminal_taxa(taxonomy, relab_data)
# get the top rows of terminal taxa
top_terminal_taxa, top_terminal_data = utilities.top_rows(terminal_taxa_relab, terminal_data_relab, max_taxa, function="average")

# reduce the taxa names to just the most specific identifier
shorted_names=[taxon.split(";")[-1] for taxon in top_terminal_taxa]

document.plot_stacked_barchart(top_terminal_data, row_labels=shorted_names, 
    column_labels=samples, title="Top "+str(max_taxa)+" terminal taxa by average abundance",
    ylabel="Relative abundance", legend_title="Terminal taxa")

#+ echo=False

# plot the top terminal node taxa in a PCOA
# provide data as values [0-1] organized as samples as columns and features as rows
# filter out any zero rows
taxa_nonzero, data_nonzero = utilities.filter_zero_rows(terminal_taxa_relab, terminal_data_relab)

document.show_pcoa(samples, taxa_nonzero, data_nonzero, title="Ordination of terminal taxa")





