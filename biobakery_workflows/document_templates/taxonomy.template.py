
#+ echo=False
import math

min_abundance=0.01
min_samples=10
max_sets_heatmap=25
max_sets_barplot=15

from biobakery_workflows import utilities

#' # Taxonomy

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
species_counts_after_filter_min_ab=count_filtered_columns(filtered_species_data, min=min_abundance)

all_species_counts=[[a,b,c] for a,b,c in zip(species_counts, species_counts_after_filter, species_counts_after_filter_min_ab)]

document.show_table(all_species_counts,samples,["Total","After filter","After filter min abun "+str(min_abundance)+"%"],
    title="Total species per sample", column_width=0.4)

#' ## Ordination

#+ echo=False
# get the top species by average abundance
top_taxonomy, top_data = utilities.top_rows(species_taxonomy, species_data, max_sets_heatmap,
    function="average") 

# compute the pcoa and plot
# provide data as range of [0-1] organised as samples as rows and features as columns
document.show_pcoa(samples,top_taxonomy,numpy.array(top_data)/100.0,"Ordination of species abundances")

#' ## Heatmap

#+ echo=False
document.show_hclust2(samples,top_taxonomy,top_data,
                      title="Top "+str(max_sets_heatmap)+" species by average abundance")

#' ## Barplot

#+ echo=False
top_taxonomy, top_data = utilities.top_rows(species_taxonomy, species_data, max_sets_barplot,
    function="average") 

document.plot_stacked_barchart(top_data, row_labels=top_taxonomy, 
    column_labels=samples, title="Top "+str(max_sets_barplot)+" species by average abundance",
    ylabel="Predicted community composition (% of total)", legend_title="Species")

