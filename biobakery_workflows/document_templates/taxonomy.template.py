
#+ echo=False
import os
import math

min_abundance=0.01
min_samples=10
max_sets_heatmap=25
max_sets_barplot=15

from biobakery_workflows import utilities, visualizations, files

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' # Taxonomic Profiling of Metagenomic Reads

#' This report section contains information about the taxonomy
#' for all DNA samples. These samples were
#' run through [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2).

#' Species abundances are passed through a basic filter requiring each species
#' to have at least <% print(min_abundance)%> % abundance in at least 
#' <% print(min_samples) %> % of all samples.

#+ echo=False

# read in the taxonomy data
samples, taxonomy, data = document.read_table(vars["taxonomic_profile"])

# remove extra information from sample name if included from workflow join
samples=[s.replace("_taxonomic_profile","") for s in samples]

# filter to only include data for the species level
# get the rows with species but not strain information
species_taxonomy, species_data = utilities.filter_species(taxonomy,data)

# now filter species also applying min abundance and min samples
filtered_species_taxonomy, filtered_species_data = utilities.filter_species(taxonomy,
    data, min_abundance=min_abundance, min_samples=min_samples)

#' A total of <% print(len(species_taxonomy)) %> species were identified. After basic
#' filtering <% print(len(filtered_species_taxonomy)) %> species remained. 

#' ## Species Count Table

#+ echo=False
import numpy

# count the number of species that pass filters for each sample
def count_filtered_columns(data, min):
    # first transpose the data
    data=numpy.transpose(data)
    # get the items from a row that pass the min filter
    # then count the total items per row
    return [len(list(filter(None,filter(lambda x: x>min,row)))) for row in data]

species_counts=count_filtered_columns(species_data, min=0)
species_counts_after_filter=count_filtered_columns(filtered_species_data, min=0)

all_species_counts=[[a,b] for a,b in zip(species_counts, species_counts_after_filter)]

# create a table of the data in the output folder
document.write_table(["# Sample","Total","After filter"],samples, all_species_counts,
    files.ShotGunVis.path("species_counts",document.data_folder))

# show the table, reducing the rows if there are lots of samples
table_message=visualizations.show_table_max_rows(document, all_species_counts, samples,
    ["Total","After filter"],"Total species per sample",files.ShotGunVis.path("species_counts"))

#' <%= table_message %>

#' ## Ordination

#+ echo=False
# get the top species by average abundance
top_taxonomy, top_data = utilities.top_rows(species_taxonomy, species_data, max_sets_heatmap,
    function="average") 

# compute the pcoa and plot
# provide data as range of [0-1] organised as samples as rows and features as columns
pcoa_data=numpy.array(top_data)/100.0
caption=document.show_pcoa(samples,top_taxonomy,pcoa_data,"Ordination of species abundances")

#' <% print(caption) %>

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
if 'metadata' in vars and vars['metadata']:
    # organize metadata for plot if available
    sample_metadata=vars["metadata"][0]
    for category in vars["metadata"][1:]:
        name=category[0]
        metadata_mapping=dict((x,y) for x,y in zip(sample_metadata[1:],category[1:]))

        document.show_pcoa(samples, top_taxonomy, pcoa_data, title="PCOA Ordination of terminal taxa using Bray-Curtis similarity - "+name,
            metadata=metadata_mapping)

#' ## Heatmap

#' <%= visualizations.ShotGun.format_caption("heatmap_intro",max_sets=max_sets_heatmap,type="species and genera",method="Spearman and Bray-Curtis") %>

#+ echo=False
# update the figure size based on output format for the heatmaps
utilities.change_pweave_figure_size_heatmap(pdf_format)

visualizations.plot_heatmap(document,vars,samples,top_taxonomy,top_data,pdf_format,max_sets_heatmap=max_sets_heatmap)

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
# reset the figure size to the defaults
utilities.reset_pweave_figure_size()

#' ## Barplot

#+ echo=False
top_taxonomy, top_data = utilities.top_rows(species_taxonomy, species_data, max_sets_barplot,
    function="average") 

sorted_data, sorted_samples = visualizations.sort_data(document, top_data, samples)

# add other to the taxonomy data
# other represents the total abundance of all species not included in the top set
top_taxonomy.append("other")
other_abundances=[]
for column in numpy.transpose(sorted_data):
    other_abundances.append(100-sum(column))
sorted_data.append(other_abundances)

#+ echo=False
# plot all samples taxonomy
document.plot_stacked_barchart(sorted_data, row_labels=top_taxonomy, 
    column_labels=sorted_samples, title="Top "+str(max_sets_barplot)+" species by average abundance",
    ylabel="Relative abundance", legend_title="Species", legend_style="italic")

#' Stacked barplot of <% print(max_sets_barplot) %> most abundant species among samples.
#' Samples in the plot were sorted on the species with the highest mean abundances among samples, in decreasing order. 

#+ echo=False
# plot a barchart for a set of categorical data
if visualizations.metadata_provided(vars):
    categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples, 
        sorted_data)
    for cat_metadata in categorical_metadata:
        visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, top_taxonomy,
            samples_found,title="Top {} species by average abundance".format(max_sets_barplot))
    # plot average for all samples grouped by categorical metadata
    for cat_metadata in categorical_metadata:
        visualizations.plot_average_taxonomy(document, ordered_sorted_data, samples_found, top_taxonomy, cat_metadata, max_sets_barplot, legend_title="species")
    
#' <% if categorical_metadata: print("Stacked barplot of average abundance among samples grouped by metadata.") %>

