#' % <% from anadama2 import PweaveDocument; document=PweaveDocument(); vars = document.get_vars(); print(vars["title"]) %>
#' % Project: <% print(vars["project"]) %>
#' % Date: <% import time; print(time.strftime("%m/%d/%Y")) %>
#+ echo=False

# get the variable settings from the data processing workflow
from anadama2.reporters import LoggerReporter
from biobakery_workflows import visualizations

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
dada_db = workflow_settings.get("dada_db", "UNK").upper()

method=vars["method"]

if method == "dada2" or method == "its":
    if method == "its":
        dada_db = "UNITE"
    dadadb_info="\nThe " + str(dada_db) + " database was used for taxonomy prediction."
    columns, samples, data = document.read_table(vars["counts_each_step"])

else:
    columns, samples, data = document.read_table(vars["read_count_table"])
    
    usearchintro="The " + str(len(samples)) + "  samples from this project were run through the standard 16S workflow.  \
        follows the UPARSE OTU analysis pipeline for OTU calling and taxonomy prediction with percent identity " \
        + str(percent_identity) + " and minimum cluster size of " + str(min_cluster_size) + "." \
        + "\n\nThe GreenGenes 16S RNA Gene Database version 13_8 was used for taxonomy prediction.\
        \n\nReads were filtered for quality control using a MAXEE score of " + str(maxee) + ". Filtered reads were \
        used to generate the OTUs. Reads not passing quality control were kept and used in the step \
        assigning reads to OTUs. First these reads were truncated to a max length of " + str(trunc_len_max) + " bases.\n"
        

import os
import numpy

# determine the document format
pdf_format=True if vars["format"] == "pdf" else False

#' <% if pdf_format: print("\clearpage") %>

#' # Introduction
#+ echo=False
#' <% if method== "its": print(visualizations.Sixteen_S.captions["itsintro"]) %>
#' <% if method == "dada2": print(visualizations.Sixteen_S.captions["dada2intro"]) %>
#' <% if method == "dada2" or method == "its": print(visualizations.Sixteen_S.captions["dada2stepsinfo"]) %>
#' <% if method == "dada2" or method == "its": print(dadadb_info) %>
#' <% if method != "dada2" and method != "its": print(usearchintro) %>

#+ echo=False


min_abundance=0.01
min_samples=10

from biobakery_workflows import utilities

# read in the read count table
# columns expected are total reads, reads that map to OTUs with taxonomy,
# and reads that map to OTUs without taxonomy

#' <% if pdf_format: print("\clearpage") %>


#' # Quality Control


#' <% if method == "dada2" or method == "its": print("## Forward Reads") %>   \
#' <% if method == "dada2" or method == "its": print("![FWD Read](" + vars["readF_qc"] + ")") %>
#' <% if method == "dada2" or method == "its": print("\clearpage")  %>

#' <% if method == "dada2" or method == "its": print("## Reverse Reads") %>   \
#' <% if method == "dada2" or method == "its": print("![REV Read](" + vars["readR_qc"] + ")") %>
#' <% if method == "dada2" or method == "its": print("\clearpage") %>
#+ echo=False

if method != "dada2" and method != "its":
    eestats_rows, eestats_columns, eestats_data, overall_stats = utilities.read_eestats2(vars["eestats_table"])
    document.show_table(eestats_data, eestats_rows, eestats_columns,"Expected error filter by read length",font="10")
    
#' <% if method != "dada2" and method != "its": print("The general stats for this data set are:" + str(overall_stats)) %>
#' <% if method != "dada2" and method != "its": print("This table shows the number of reads based on length for different error filters.") %>

#' <% if method == "dada2" or method == "its": print("# Error rates") %>
    
#' <% if method == "dada2" or method == "its": print(visualizations.Sixteen_S.captions["dada2errorintro"]) %>
    
#' <% if method == "dada2" or method == "its": print("## Forward Reads") %>  \
#' <% if method == "dada2" or method == "its": print("![FWD Error Rates](" + vars["error_ratesF"] +")") %>
#' <% if method == "dada2" or method == "its": print("\clearpage") %>

#' <% if method == "dada2" or method == "its": print("## Reverse Reads") %>   \
#' <% if method == "dada2" or method == "its": print("![REV Error Rates](" + vars["error_ratesR"] + ")") %>

#' <% if method == "dada2" or method == "its": print(visualizations.Sixteen_S.captions["dada2errorinfo"]) %>
#' <% if pdf_format: print("\clearpage") %>
      
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

if method == "dada2" or method == "its":

    filtered_reads = [row[1] for row in sorted_all_read_data]
    merged_reads = [row[3] for row in sorted_all_read_data]
    tabled_reads = [row[4] for row in sorted_all_read_data]
    nochim_reads = [row[5] for row in sorted_all_read_data]
    # plot the read counts
    document.plot_grouped_barchart([sorted_total_reads,filtered_reads,merged_reads,tabled_reads,nochim_reads],
           ["Original","Filtered","Merged","Tabled","Nochimera"], sorted_samples,
           title="Read counts by Sample", xlabel="Samples", ylabel="Total Reads")
    # plot grouped taxonomy for all categorical data provided
    if visualizations.metadata_provided(vars):
        categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples,
            [total_reads,filtered_reads,merged_reads,tabled_reads,nochim_reads])
        for cat_metadata in categorical_metadata:
            visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, 
                ["total","filtered","merged","tabled","nochimera"], samples_found, 
                title="Read counts by sample", ylabel="Total Reads", legend_title="")

else:
    known_reads = [row[1] for row in sorted_all_read_data]
    unknown_reads = [row[2] for row in sorted_all_read_data]
    unmapped_reads = [row[0]-(row[1]+row[2]) for row in sorted_all_read_data]
    # plot the read counts
    document.plot_stacked_barchart([known_reads,unknown_reads,unmapped_reads], ["classified","unclassified","unmapped"], sorted_samples, 
        title="Read counts by Sample", ylabel="Total Reads", xlabel="Samples")
    # plot grouped taxonomy for all categorical data provided
    if visualizations.metadata_provided(vars):
        categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = visualizations.merge_categorical_metadata(vars, sorted_samples,
            [known_reads,unknown_reads,unmapped_reads])
        for cat_metadata in categorical_metadata:
            visualizations.plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, ["classified","unclassified","unmapped"],
                samples_found, title="Read counts by sample", ylabel="Total Reads", legend_title="")

#' <% if method == "dada2" or method == "its": print(visualizations.Sixteen_S.captions["dada2countsinfo"]) %>
#' <% if method != "dada2" and method != "its": print(visualizations.Sixteen_S.captions["usearchcountsinfo"]) %>

#' <% if pdf_format: print("\clearpage") %>

#' # Taxonomy
    
#' <% if method == "dada2" or method == "its": print(visualizations.Sixteen_S.captions["dada2taxinfo"]) %>
#' <% if method == "dada2" or method == "its": print(dadadb_info) %>

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

#' <%= visualizations.ShotGun.format_caption("heatmap_intro",max_sets=max_sets_heatmap,type="genera and terminal taxa",method="Spearman and Bray-Curtis") %>

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
