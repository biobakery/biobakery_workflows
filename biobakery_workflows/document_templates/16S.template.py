
#+ echo=False

min_abundance=0.01
min_samples=10

#' # Taxonomy
    
#' ## Genera
 
#+ echo=False

# read in the otu table data
samples, ids, taxonomy, data = utilities.read_otu_table(vars["otu_table"])

# plot the top taxa by genus level, plotting the relative abundance values
max_taxa=15

# get the relative abundance values for the samples
relab_data = utilities.relative_abundance(data, percent=True)

# get the top taxa by genus level
max_taxa = 15
sorted_samples, sorted_top_data, top_data, top_taxa_short_names, legend_size = visualizations.get_top_taxonomy_by_level(taxonomy, samples, relab_data, max_taxa)

# add other to the taxonomy data, other represents total genera not shown on plot
top_taxa_short_names_plus_other, sorted_top_data_plus_other = visualizations.fill_taxonomy_other(top_taxa_short_names, sorted_top_data)

document.plot_stacked_barchart(sorted_top_data_plus_other, row_labels=top_taxa_short_names_plus_other,
    column_labels=sorted_samples, title="Top "+str(max_taxa)+" genera by average abundance",
    ylabel="Relative abundance", legend_title="Genera", legend_style="italic", legend_size=legend_size)

#+ echo=False
# plot grouped and average barplots for metadata if provided
visualizations.plot_grouped_and_average_barplots_taxonomy(document, vars, sorted_samples, sorted_top_data_plus_other, top_taxa_short_names_plus_other, max_taxa, feature="genera")

#' <% if pdf_format: print("\clearpage") %>

#' ## Terminal Taxa

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

# add the remaining terminal taxa as "other" to the data
shorted_names_plus_other, sorted_top_terminal_data_plus_other = visualizations.fill_taxonomy_other(shorted_names, sorted_top_terminal_data)

document.plot_stacked_barchart(sorted_top_terminal_data_plus_other, row_labels=shorted_names_plus_other,
    column_labels=sorted_samples_terminal, title="Top "+str(max_taxa)+" terminal taxa by average abundance",
    ylabel="Relative abundance", legend_title="Terminal taxa")

visualizations.plot_grouped_and_average_barplots_taxonomy(document, vars, sorted_samples_terminal, sorted_top_terminal_data_plus_other, shorted_names_plus_other, max_taxa, feature="terminal taxa")

#' # Heatmaps

#+ echo=False
max_sets_heatmap=25

#' <%= visualizations.ShotGun.format_caption("heatmap_intro",max_sets=max_sets_heatmap,type="genera and terminal taxa",method="Spearman and Bray-Curtis", data_type="taxa") %>

#' ## Genera

#+ echo=False
# update the figure size based on output format for the heatmaps
utilities.change_pweave_figure_size_heatmap(pdf_format)

visualizations.plot_heatmap(document,vars,samples,top_taxa_short_names,top_data,
    pdf_format,"Top {} genera by average abundance (Spearman)".format(max_sets_heatmap),max_sets_heatmap)
utilities.reset_pweave_figure_size()
#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
utilities.change_pweave_figure_size_heatmap(pdf_format)
visualizations.plot_heatmap(document,vars,samples,top_taxa_short_names,top_data,
    pdf_format,"Top {} genera by average abundance (Bray-Curtis)".format(max_sets_heatmap),max_sets_heatmap,method="lbraycurtis")
utilities.reset_pweave_figure_size()



#' ## Terminal Taxa

#+ echo=False
utilities.change_pweave_figure_size_heatmap(pdf_format)
visualizations.plot_heatmap(document,vars,samples,shorted_names,top_terminal_data,
    pdf_format,"Top {} terminal taxa by average abundance (Spearman)".format(max_sets_heatmap),max_sets_heatmap)
utilities.reset_pweave_figure_size()

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
utilities.change_pweave_figure_size_heatmap(pdf_format)
visualizations.plot_heatmap(document,vars,samples,shorted_names,top_terminal_data,
    pdf_format,"Top {} terminal taxa by average abundance (Bray-Curtis)".format(max_sets_heatmap),max_sets_heatmap,method="lbraycurtis")
utilities.reset_pweave_figure_size()


#' # Ordination

#' ## Genera

#+ echo=False
import numpy

# filter then get top genera
filtered_taxonomy_all, filtered_relab_data_all = utilities.filter_taxa_abundance(taxonomy, relab_data, min_abundance, min_samples)
samples_genera, sorted_top_data_genera, top_data_genera, top_taxa_genera, legend_size = visualizations.get_top_taxonomy_by_level(filtered_taxonomy_all, samples, filtered_relab_data_all, max_sets_heatmap)

# provide data as values [0-1] organized as samples as columns and features as rows
top_filtered_data_pcoa=numpy.array(sorted_top_data_genera)/100.0
document.show_pcoa(samples_genera, top_taxa_genera, top_filtered_data_pcoa, title="PCoA Ordination of top {} genera using Bray-Curtis similarity".format(max_sets_heatmap))

#' For the PCoA plots, relative abundances are passed through a basic filter requiring each taxon
#' have at least <% print(min_abundance)%> % abundance in at least
#' <% print(min_samples) %> % of all samples.

#+ echo=False
visualizations.show_pcoa_metadata(document, vars, samples_genera, top_taxa_genera, top_filtered_data_pcoa,
    title="PCoA Ordination of top {} genera".format(max_sets_heatmap))

#' ## Terminal taxa

#+ echo=False
# plot the top terminal node taxa in a PCoA

# filter the data by min abundance and min samples
filtered_taxonomy, filtered_data = utilities.filter_taxa_abundance(terminal_taxa_relab, terminal_data_relab, min_abundance, min_samples)

# get the top set of terminal taxa
top_filtered_taxonomy, top_filtered_data = utilities.top_rows(filtered_taxonomy, filtered_data, max_sets_heatmap, function="average")

# provide data as values [0-1] organized as samples as columns and features as rows
top_filtered_data_pcoa=numpy.array(top_filtered_data)/100.0
document.show_pcoa(samples, top_filtered_taxonomy, top_filtered_data_pcoa, title="PCoA Ordination of top {} terminal taxa using Bray-Curtis similarity".format(max_sets_heatmap))

visualizations.show_pcoa_metadata(document, vars, samples, top_filtered_taxonomy, top_filtered_data_pcoa,
    title="PCoA Ordination of top {} terminal taxa".format(max_sets_heatmap))

#' <% if pdf_format: print("\clearpage") %>
