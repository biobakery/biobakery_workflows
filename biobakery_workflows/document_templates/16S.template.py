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

method=vars["method"]

if method == "dada2":
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
#' <% if vars["method"] == "dada2": print(visualizations.Sixteen_S.captions["dada2intro"]) %>
#' <% if vars["method"] == "dada2": print(visualizations.Sixteen_S.captions["dada2stepsinfo"]) %>
#' <% if vars["method"] != "dada2": print(usearchintro) %>

#+ echo=False


min_abundance=0.01
min_samples=10

from biobakery_workflows import utilities

# read in the read count table
# columns expected are total reads, reads that map to OTUs with taxonomy,
# and reads that map to OTUs without taxonomy

#' <% if pdf_format: print("\clearpage") %>


#' # Quality Control


#' <% if method == "dada2": print("## Forward Read Quality Plot by Sample") %>   \
#' <% if method == "dada2": print("![FWD Read](" + vars["readF_qc"] + ")") %> 
#' <% if method == "dada2": print("\clearpage")  %> 

#' <% if method == "dada2": print("## Reverse Read Quality Plot  by Sample") %>   \
#' <% if method == "dada2": print("![REV Read](" + vars["readR_qc"] + ")") %>
#' <% if method == "dada2": print("\clearpage") %>
#+ echo=False

if method != "dada2":
    eestats_rows, eestats_columns, eestats_data, overall_stats = utilities.read_eestats2(vars["eestats_table"])
    document.show_table(eestats_data, eestats_rows, eestats_columns,"Expected error filter by read length",font="10")
    
#' <% if method != "dada2": print("The general stats for this data set are:" + str(overall_stats)) %>
#' <% if method != "dada2": print("This table shows the number of reads based on length for different error filters.") %>    

#' <% if method == "dada2": print("# Error rates") %>
    
#' <% if method == "dada2": print(visualizations.Sixteen_S.captions["dada2errorintro"]) %>  
    
#' <% if method == "dada2": print("## Forward Read Error Rates by Sample") %>  \
#' <% if method == "dada2": print("![FWD Error Rates](" + vars["error_ratesF"] +")") %>
#' <% if method == "dada2": print("\clearpage") %>

#' <% if method == "dada2": print("## Reverse Read Error Rates by Sample") %>   \
#' <% if method == "dada2": print("![REV Error Rates](" + vars["error_ratesR"] + ")") %>

#' <% if method == "dada2": print(visualizations.Sixteen_S.captions["dada2errorinfo"]) %> 
#' <% if pdf_format: print("\clearpage") %>
      
#+ echo=False
def sort_data(top_data, samples):
    # sort the top data so it is ordered with the top sample/abundance first
    sorted_sample_indexes=sorted(range(len(samples)),key=lambda i: top_data[0][i],reverse=True)
    sorted_samples=[samples[i] for i in sorted_sample_indexes]
    sorted_data=[]
    for row in top_data:
        sorted_data.append([row[i] for i in sorted_sample_indexes])
    return sorted_data, sorted_samples

def plot_grouped_taxonomy_subsets(sorted_data, cat_metadata, taxa, title, samples_found, ylabel, legend_title, legend_size):
    """ Plot the grouped taxonomy sorted by species abundance for a single feature """
    # group the samples by metadata
    sorted_data_grouped, sorted_samples_grouped = utilities.group_samples_by_metadata(cat_metadata, sorted_data, samples_found)
    # sort the data by abundance
    for metadata_type in sorted_data_grouped:
        sorted_data_grouped[metadata_type], sorted_samples_grouped[metadata_type] = sort_data(sorted_data_grouped[metadata_type], sorted_samples_grouped[metadata_type])

    # print out a plot for each group of metadata if there are lots of categories
    sorted_metadata_subsets=sorted(sorted_data_grouped.keys())
    
    try:
        sorted_metadata_subsets=sorted(sorted_data_grouped.keys(), key=float)
    except ValueError:
        pass
    
    max_subsets=8
    
    # split into subsets
    split_sorted_metadata_subsets = [sorted_metadata_subsets[x:x+max_subsets] for x in range(0, len(sorted_metadata_subsets), max_subsets)]
    
    # make sure the last group is not just a single data set
    if len(split_sorted_metadata_subsets[-1]) == 1:
        last_set = split_sorted_metadata_subsets.pop()
        split_sorted_metadata_subsets[-1].append(last_set[0])
    
    for metadata_subset in split_sorted_metadata_subsets:
        subset_sorted_data_grouped=dict((key, sorted_data_grouped[key]) for key in metadata_subset)
        subset_sorted_samples_grouped=dict((key, sorted_samples_grouped[key]) for key in metadata_subset)
        
        # get title addition for subset
        title_add=""
        if len(sorted_metadata_subsets) > max_subsets:
            title_add=" "+metadata_subset[0]+" to "+metadata_subset[-1]

        document.plot_stacked_barchart_grouped(subset_sorted_data_grouped, row_labels=taxa, 
            column_labels_grouped=subset_sorted_samples_grouped, title=title+" - "+str(cat_metadata[0])+title_add,
            ylabel=ylabel, legend_title=legend_title, legend_style="italic", legend_size=legend_size)
        

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

def plot_all_categorical_metadata(sorted_samples, sorted_data, labels, title, ylabel, legend_title="", legend_size=7):
    """ Generate a plot of each set of categorical metadata """
    if 'metadata' in vars and vars['metadata'] and 'metadata_labels' in vars and vars['metadata_labels']:
        # get the metadata organized into the same sample columns as the data
        new_data, samples_found = utilities.merge_metadata(vars['metadata'], sorted_samples, sorted_data, values_without_names=True)
        # split the data and metadata 
        ordered_metadata=new_data[0:len(vars['metadata'])-1]
        ordered_sorted_data=new_data[len(vars['metadata'])-1:]
        # get the categorical metadata
        categorical_metadata=utilities.filter_metadata_categorical(ordered_metadata, vars['metadata_labels'])
        # plot a bar chart for a set of categorical data
        for cat_metadata in categorical_metadata:
            plot_grouped_taxonomy_subsets(ordered_sorted_data, cat_metadata, labels, title, samples_found, ylabel, legend_title, legend_size)

if method == "dada2":

    filtered_reads = [row[1] for row in sorted_all_read_data]
    merged_reads = [row[3] for row in sorted_all_read_data]
    tabled_reads = [row[4] for row in sorted_all_read_data]
    nochim_reads = [row[5] for row in sorted_all_read_data]
    # plot the read counts
    document.plot_grouped_barchart([sorted_total_reads,filtered_reads,merged_reads,tabled_reads,nochim_reads],
           ["Original","Filtered","Merged","Tabled","Nochimera"], sorted_samples,
           title="Read counts by Sample", xlabel="Samples", ylabel="Total Reads")
    plot_all_categorical_metadata(sorted_samples, [total_reads,filtered_reads,merged_reads,tabled_reads,nochim_reads], 
    ["total","filtered","merged","tabled","nochimera"], title="Read counts in each step by sample", ylabel="Total Reads")


else:
    known_reads = [row[1] for row in sorted_all_read_data]
    unknown_reads = [row[2] for row in sorted_all_read_data]
    unmapped_reads = [row[0]-(row[1]+row[2]) for row in sorted_all_read_data]
    # plot the read counts
    document.plot_stacked_barchart([known_reads,unknown_reads,unmapped_reads], ["classified","unclassified","unmapped"], sorted_samples, 
        title="Read counts by Sample", ylabel="Total Reads", xlabel="Samples")
    plot_all_categorical_metadata(sorted_samples, [known_reads,unknown_reads,unmapped_reads], 
        ["classified","unclassified","unmapped"], title="Read counts by Sample", ylabel="Total Reads")
 
#' <% if vars["method"] == "dada2": print(visualizations.Sixteen_S.captions["dada2countsinfo"]) %>  
#' <% if vars["method"] != "dada2": print(visualizations.Sixteen_S.captions["usearchcountsinfo"]) %> 

#' <% if pdf_format: print("\clearpage") %>

#' # Taxonomy
#' ## Average Abundance
    
#' <% if vars["method"] == "dada2": print(visualizations.Sixteen_S.captions["dada2taxinfo"]) %>      
    
#+ echo=False

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

# check for duplicate genera in list
legend_size = 7
if len(top_taxa_short_names) != len(list(set(top_taxa_short_names))):
    # if duplicate names, then add family to the taxonomy
    top_taxa_short_names = [family+"."+genus for family, genus in zip(utilities.taxa_shorten_name(top_taxa, level=4),utilities.taxa_shorten_name(top_taxa, level=5))]
    # reduce legend size to fit names
    legend_size = 5

# sort the data so those with the top genera are shown first
sorted_samples, sorted_data = utilities.sort_data(top_data[0], samples)

transpose_top_data = numpy.transpose(top_data)
sorted_top_data = numpy.transpose([transpose_top_data[samples.index(sample)] for sample in sorted_samples])

document.plot_stacked_barchart(sorted_top_data, row_labels=top_taxa_short_names, 
    column_labels=sorted_samples, title="Top "+str(max_taxa)+" genera by average abundance",
    ylabel="Relative abundance", legend_title="Genera", legend_style="italic", legend_size=legend_size)


plot_all_categorical_metadata(sorted_samples, sorted_top_data, top_taxa_short_names,
    title="Top "+str(max_taxa)+" genera by average abundance", ylabel="Relative abundance", legend_title="Genera", legend_size=legend_size)


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

document.plot_stacked_barchart(sorted_top_terminal_data, row_labels=shorted_names, 
    column_labels=sorted_samples_terminal, title="Top "+str(max_taxa)+" terminal taxa by average abundance",
    ylabel="Relative abundance", legend_title="Terminal taxa")

plot_all_categorical_metadata(sorted_samples_terminal, sorted_top_terminal_data, shorted_names,
    title="Top "+str(max_taxa)+" terminal taxa by average abundance", ylabel="Relative abundance", legend_title="Terminal taxa")



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
    pdf_format,"Top {} terminal taxa by average abundance (Bray-Curtis)".format(max_sets_heatmap),max_sets_heatmap)
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