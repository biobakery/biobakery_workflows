#' % <% from anadama2 import PweaveDocument; document=PweaveDocument(); vars = document.get_vars(); print(vars["title"]) %>
#' % Project: <% print(vars["project"]) %>
#' % Date: <% import time; print(time.strftime("%m/%d/%Y")) %>

#' # Introduction

#+ echo=False
# get the variable settings from the data processing workflow
from anadama2.reporters import LoggerReporter
try:
    workflow_settings = LoggerReporter.read_log(vars["log"],"variables")
except AttributeError:
    workflow_settings = []

# print a warning if the variables could not be read
if isinstance(workflow_settings, list):
    print("WARNING: Unable to read workflow settings from log file.")
    workflow_settings={}

maxee = workflow_settings.get("maxee","UNK")
trunc_len_max = workflow_settings.get("trunc_len_max","UNK")
percent_identity = workflow_settings.get("percent_identity","UNK")
min_cluster_size = workflow_settings.get("min_size","UNK")

# read in the read count table
# columns expected are total reads, reads that map to OTUs with taxonomy,
# and reads that map to OTUs without taxonomy
columns, samples, data = document.read_table(vars["read_count_table"])

#' The <% print(len(samples)) %> samples from this project were run through the standard 16S workflow. The workflow
#' follows the UPARSE OTU analysis pipeline for OTU calling and taxonomy prediction with percent identity 
#' of <%= percent_identity %> and minimum cluster size of <%= min_cluster_size %>. 
#' The GreenGenes 16S RNA Gene Database version 13_8 was used for taxonomy prediction.
#' Reads were filtered for quality control using a MAXEE score of <%= maxee %>. Filtered reads were
#' used to generate the OTUs. Reads not passing quality control were kept and used in the step
#' assigning reads to OTUs. First these reads were truncated to a max length of <%= trunc_len_max %> bases.

#+ echo=False
import os
import numpy

min_abundance=0.01
min_samples=10

from biobakery_workflows import utilities, visualizations

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' <% if pdf_format: print("\clearpage") %>

#' # Quality Control

#+ echo=False
# read the eestats2 table and add table to document
eestats_rows, eestats_columns, eestats_data, overall_stats = utilities.read_eestats2(vars["eestats_table"])

document.show_table(eestats_data, eestats_rows, eestats_columns, 
    "Expected error filter by read length",font="10")

#' The general stats for this data set are: <%= overall_stats %> .
#' This table shows the number of reads based on length for different error filters.

#+ echo=False
# if picard files are provided, then plot those that do not meet a threshold
picard_text = ""
if vars["picard"]:
    picard_min_threshold=20
    # read through all picard files, capturing those based on the threshold
    data_below_threshold={}
    samples_above_threshold=[]
    for file in utilities.get_files(vars["picard"], extension=vars["picard_ext"]):
        picard_sample_name = os.path.basename(file).replace("."+vars["picard_ext"],"")
        picard_data, below_threshold = utilities.read_picard(file, threshold=picard_min_threshold)
        if picard_data:
            if below_threshold:
                data_below_threshold[picard_sample_name]=[a[1] for a in picard_data]
            else:
                samples_above_threshold.append(picard_sample_name)

    # plot each sample that is below the threshold
    for sample, picard_data in data_below_threshold.items():
        document.plot_barchart(picard_data,
            title="Picard quality scores for "+sample, xlabel="Cycle", ylabel="Quality Score")
    
    
    if len(data_below_threshold.keys()) > 0:
        picard_text="Only the samples with at least one quality score below the threshold "+\
            " ( "+str(picard_min_threshold)+" ) are shown as quality scores for each base. These samples are: "+\
            ",".join(list(data_below_threshold.keys()))+"."
    else:
        picard_text="All samples had all quality scores above the threshold "+\
            " ( "+str(picard_min_threshold)+" )."
    
    # list those that were above the threshold
    above_threshold_list=",".join(samples_above_threshold)
    if above_threshold_list:
        picard_text+=" The following samples did not have any quality scores below the threshold: " + above_threshold_list + "."
    else:
        picard_text+=" None of the samples had all quality scores above the threshold."
                    

#' <% if picard_text: print(picard_text) %>

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

#+ echo=False
# plot grouped taxonomy for all categorical data provided
if visualizations.metadata_provided(vars):
    categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples, 
        [known_reads,unknown_reads,unmapped_reads])
    for cat_metadata in categorical_metadata:
        visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, ["classified","unclassified","unmapped"], 
            samples_found, title="Read counts by Sample", ylabel="Total Reads", legend_title="")

#' # Taxonomy

#' ## Top Genera

#+ echo=False
import numpy

# read in the otu table data
samples, ids, taxonomy, data = utilities.read_otu_table(vars["otu_table"])

# get the relative abundance values for the samples
relab_data = utilities.relative_abundance(data)

# get the top taxa by genus level
max_taxa = 15
sorted_samples, sorted_top_data, top_data, top_taxa_short_names, legend_size = visualizations.get_top_taxonomy_by_level(taxonomy, samples, relab_data, max_taxa)

document.plot_stacked_barchart(sorted_top_data, row_labels=top_taxa_short_names, 
    column_labels=sorted_samples, title="Top "+str(max_taxa)+" genera by average abundance",
    ylabel="Relative abundance", legend_title="Genera", legend_style="italic", legend_size=legend_size)

#+ echo=False
if visualizations.metadata_provided(vars):
    categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples, sorted_top_data)
    # plot taxonomy grouped by feature
    for cat_metadata in categorical_metadata:
        visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, top_taxa_short_names, 
            samples_found, title="Top {} genera by average abundance".format(max_taxa), ylabel="Relative abundance", 
            legend_title="Genera", legend_size=legend_size)
    # plot taxonomy average for all samples, grouped by feature
    for cat_metadata in categorical_metadata:
        visualizations.plot_average_taxonomy(document, ordered_sorted_data, samples_found, top_taxa_short_names,
            cat_metadata, max_taxa, legend_title="Genera")

#' <% if pdf_format: print("\clearpage") %>

#' ## Top Terminal Taxa

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

if visualizations.metadata_provided(vars):
    categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples_terminal, 
        sorted_top_terminal_data)
    for cat_metadata in categorical_metadata:
        visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, shorted_names, samples_found,
            title="Top {} terminal taxa by average abundance".format(max_taxa), ylabel="Relative abundance", legend_title="Terminal taxa")
    for cat_metadata in categorical_metadata:
        visualizations.plot_average_taxonomy(document, ordered_sorted_data, samples_found, shorted_names, 
            cat_metadata, max_taxa, legend_title="Terminal taxa")

#' # Heatmaps

#+ echo=False
max_sets_heatmap=25

#' <%= visualizations.ShotGun.format_caption("heatmap_intro",max_sets=max_sets_heatmap,type="genera and terminal taxa",method="Spearman and Bray-Curtis") %>

#' ## Top Genera

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

#' ## Top Terminal Taxa

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

#' <% if pdf_format: print("\clearpage") %>

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


#+ echo=False
if 'metadata' in vars and vars['metadata']:
    # organize metadata for plot if available
    sample_metadata=vars["metadata"][0]
    for category in vars["metadata"][1:]:
        name=category[0]
        metadata_mapping=dict((x,y) for x,y in zip(sample_metadata[1:],category[1:]))

        document.show_pcoa(samples, filtered_taxonomy, filtered_data, title="PCOA Ordination of terminal taxa using Bray-Curtis similarity "+name,
            metadata=metadata_mapping)

