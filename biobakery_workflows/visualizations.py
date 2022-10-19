"""
bioBakery Workflows: visualizations module
A collection of utilities for the document visualizations

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os
import re
import copy
import sys
import subprocess
import tempfile
import shutil

from . import utilities

# if there are not names with the ecs, add them to the file for the plots generated
def add_ec_names(workflow, ecsabundance, output, template_depends):

    # allow for optional ec file
    if not ecsabundance:
        return ecsabundance

    # check for names in the ec file
    with open(ecsabundance) as file_handle:
        header = file_handle.readline()
        top_ec = file_handle.readline()
        if ":" in top_ec:
            return ecsabundance

    ecsabundance_withnames=os.path.join(output,"ecs",os.path.basename(ecsabundance))
    rename_task=workflow.add_task(
        "mkdir -p [args[0]] && humann_rename_table --input [depends[0]] --output [targets[0]] --names ec",
        depends=ecsabundance,
        targets=ecsabundance_withnames,
        args=os.path.dirname(ecsabundance_withnames),
        name="humann_rename_table_ecs")

    # add the task as a dependency of the template task
    template_depends+=[rename_task]

    return ecsabundance_withnames

# sort the samples/data by read count with the largest original read count first
def sort_samples_reads_decreasing(read_data, read_samples):
    """ Sort the reads from largest to smallest total read count """

    sorted_samples, sorted_total_reads = utilities.sort_data(read_data[0], read_samples)
    sorted_all_read_data = []
    for data_set in read_data:
        sorted_all_read_data.append([data_set[read_samples.index(sample)] for sample in sorted_samples])

    return sorted_samples, sorted_all_read_data

def plot_grouped_and_average_barplots_taxonomy(document, vars, sorted_samples, sorted_data, top_taxonomy,
    max_sets_barplot, max_groups_barplot, feature="species", sort_by_name=False, sort_by_name_inverse=False, ylabel="Relative abundance"):
    """ Plot grouped barplots and average barplots for all of the features provided.

        Args:
            document (anadama2.document): Document to use to add plots
            vars (dict): The dictionary of input variables provided for the visualization run
            sorted_samples (list): The sample names organized to match the data
            sorted_data (list): The data organized to match the sample/taxonomy
            top_taxonomy (list): The list of full taxonomy names
            max_sets_barplot (int): The max number of features (default species) to plot
            max_groups_barplot (int): The max number of grouped barplots to show for a metadata variable
            feature (string): The data being plotted (default species)
            sort_by_name (bool): If set, sort data by sample name instead of abundance
            sort_by_name_inverse (bool): If true, sort by the inverse of the name (so the reverse of the string)
                this is useful for samples with sample name plus features
            ylabel (string): The label for the y-axis for the average plots

        Returns: 
            list: A list of metadata strings (one for each categorical feature plotted)
    """

    # plot a barchart for a set of categorical data
    categorical_metadata = []
    if metadata_provided(vars):
        categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found = merge_categorical_metadata(vars, sorted_samples,
            sorted_data)
        for cat_metadata in categorical_metadata:
            plot_grouped_taxonomy_subsets(document, ordered_sorted_data, cat_metadata, top_taxonomy,
                samples_found,title="Top {} {} by average abundance".format(max_sets_barplot,feature),ylabel=ylabel,sort_by_name=sort_by_name,
                sort_by_name_inverse=sort_by_name_inverse,feature=feature,max_groups_barplot=max_groups_barplot,legend_title=feature[0].upper()+feature[1:])
        # plot average for all samples grouped by categorical metadata
        for cat_metadata in categorical_metadata:
            plot_average_taxonomy(document, ordered_sorted_data, samples_found, top_taxonomy, cat_metadata, max_sets_barplot, legend_title=feature, ylabel=ylabel)

    return categorical_metadata

def fill_taxonomy_other(top_taxonomy, sorted_data):
    """ Fill the taxonomy abundances provided with 'other' so each totals 100.

        Args:
            top_taxonomy (list): The list of full taxonomy names
            sorted_data (list): The data sorted by abundance matched to the taxonomy list

        Returns:
            list : The list of full taxonomy names plus 'other'
            list : The list of data filled with the 'other' abundance
    """

    import numpy

    # add other to the taxonomy data
    # other represents the total abundance of all species not included in the top set
    new_top_taxonomy = top_taxonomy + ["Other"]
    new_sorted_data = sorted_data
    other_abundances=[]
    for column in numpy.transpose(sorted_data):
        other_abundances.append(100-sum(column))

    if isinstance(new_sorted_data,numpy.ndarray):
        new_sorted_data = new_sorted_data.tolist()
    new_sorted_data.append(other_abundances)

    return new_top_taxonomy, new_sorted_data


def show_pcoa_metadata(document, vars, samples, top_taxonomy, pcoa_data, title, data_type):
    """ Plot pcoa for each feature if metadata has been provided

        Args:
            document (anadama2.document): Document to use to add plots
            vars (dict): The dictionary of input variables provided for the visualization run
            samples (list): The sample names organized to match the data
            top_taxonomy (list): The full taxonomic names organized to match the data
            pcoa_data (list): The data (range 0 to 1) organized to match the samples and taxonomy
            title (string): The base string for the title for the plots
            data_type (string): The type of data (taxonomic level)
        Return: None

    """

    if 'metadata' in vars and vars['metadata']:
        # organize metadata for plot if available
        sample_metadata=vars["metadata"][0]
        for category in vars["metadata"][1:]:
            name=category[0]
            metadata_mapping=dict((x,y) for x,y in zip(sample_metadata[1:],category[1:]))
            metadata_dict = vars['metadata_labels']
            metadata_type = metadata_dict[name]
            
            document.show_pcoa(samples, top_taxonomy, pcoa_data,title=title+" - "+name,
                               metadata=metadata_mapping, metadata_type=metadata_type, outfilename=os.path.join(document.figures_folder,"pcoa_"+name+"_"+data_type+".png"))

def get_top_taxonomy_by_level(taxonomy, samples, relab_data, max_taxa, taxa_level=5):
    """ Get the top, based on average abundance, set of taxa at the genus level

        Args:
            taxonomy (list): The full taxonomic names organized to match the data
            samples (list): The sample names organized to match the data
            relab_data (list): The relative abundance data
            max_taxa (int): The max number of genera to identify
            taxa_level (int): The taxonomy level to select (default of 5 is genus, using 4 would be family)

        Return:
            list : The sorted sample names (sorted to list the top genera first)
	    list : The data filtered to only the top genera (sorted by decreasing abundance)
            list : The data filtered to only the top genera 
            list : The taxonomic names (short for plotting) of the top taxa
            int : The recommended font size for the ploting legend (based on the taxonomic names)
    """ 

    import numpy

    # get the taxa summarized by genus level
    genus_level_taxa, genus_level_data = utilities.taxa_by_level(taxonomy, relab_data, level=taxa_level)
    # get the top rows of the relative abundance data
    top_taxa, top_data = utilities.top_rows(genus_level_taxa, genus_level_data, max_taxa, function="average")
    # shorten the top taxa names to just the genus level for plotting
    top_taxa_short_names = utilities.taxa_shorten_name(top_taxa, level=5, remove_identifier=True)
    # check for duplicate genera in list
    legend_size = 7
    if len(top_taxa_short_names) != len(list(set(top_taxa_short_names))):
        # if duplicate names, then add family to the taxonomy
        top_taxa_short_names = [family+"."+genus for family, genus in zip(utilities.taxa_shorten_name(top_taxa, level=taxa_level-1),utilities.taxa_shorten_name(top_taxa, level=taxa_level))]
        # reduce legend size to fit names
        legend_size = 5

    # sort the data so those with the top genera are shown first
    sorted_samples, sorted_data = utilities.sort_data(top_data[0], samples)
    transpose_top_data = numpy.transpose(top_data)
    sorted_top_data = numpy.transpose([transpose_top_data[samples.index(sample)] for sample in sorted_samples])

    return sorted_samples, sorted_top_data, top_data, top_taxa_short_names, legend_size


def metadata_provided(vars):
    """ Check if metadata was provided by the user for this visualization run 

        Args:
            vars (dict): The dictionary of input variables provided for the visualization run 

        Returns:
            bool : True if metadata was provided from the user for this run

    """

    if 'metadata' in vars and vars['metadata'] and 'metadata_labels' in vars and vars['metadata_labels']:
        return True
    else:
        return False

def merge_categorical_metadata(vars, sorted_samples, sorted_data):
    """ Merge data with metadata and filter continuous data 

        Args:
            vars (dict): The dictionary of input variables provided for the visualization run
            sorted_samples (list): A list of sample names sorted to match the data provided
            sorted_data (list): The full set of data (taxonomy or other)

        Returns: 
            list: A list of metadata strings (one for each categorical feature)
            list: The data modified to only include those samples with metadata
            list: The metadata sorted to match the order of the data (categorical and continuous)
            list: The metadata sorted to match the order of the data (only categorical)
    """

    # get the metadata organized into the same sample columns as the data
    new_data, samples_found = utilities.merge_metadata(vars['metadata'], sorted_samples, sorted_data, values_without_names=True)
    # split the data and metadata
    ordered_metadata=new_data[0:len(vars['metadata'])-1]
    ordered_sorted_data=new_data[len(vars['metadata'])-1:]
    # get the categorical metadata
    categorical_metadata=utilities.filter_metadata_categorical(ordered_metadata, vars['metadata_labels'])

    return categorical_metadata, ordered_sorted_data, ordered_metadata, samples_found


def plot_average_taxonomy(document, ordered_sorted_data, samples_found, top_taxonomy, cat_metadata, max_sets_barplot, legend_title, ylabel="Relative Abundance"):
    """ Plot the average taxonomy sorted by species abundance for a single feature.

        Args:
            document (anadama2.document): Document to use to add plots
            ordered_sorted_data (list): A list of sorted/ordered data (to match the sample names and metadata)
            samples_found (list): An ordered list of sample names (to match the data/metadata)
            top_taxonomy (list): A list of the taxonomy names ordered to match the data
            cat_metadata (list): The ordered categorical metadata (for a single feature)
            max_sets_barplot (int): The max number of taxonomic values (to use for the title)
            legend_title (string): The legend title for the plots
            ylabel (string): The label for the y-axis

        Returns: None
    """

    # group the samples by metadata
    sorted_data_grouped, sorted_samples_grouped = utilities.group_samples_by_metadata(cat_metadata, ordered_sorted_data, samples_found)
    metadata_names=[]
    average_data=[]
    for name, row in sorted_data_grouped.items():
        metadata_names.append(name)
        average_data.append([sum(group)/(1.0*len(group)) for group in row])

    # reorder average data so it is grouped by taxonomy
    average_data = [list(a) for a in zip(*average_data)]

    # sort the data by abundance or name (depending on metadata type, use name if numeric)
    sort_by_name=False
    if document.sorted_data_numerical_or_alphabetical(metadata_names) != sorted(metadata_names):
        sort_by_name = True

    sorted_data, sorted_names = sort_data(document, average_data, metadata_names, sort_by_name=sort_by_name)

    document.plot_stacked_barchart(sorted_data, row_labels=top_taxonomy,
        column_labels=sorted_names, 
        title="Top {} {} group average - {}".format(max_sets_barplot, legend_title, cat_metadata[0]),
        ylabel=ylabel, legend_title=legend_title[0].upper()+legend_title[1:], legend_style="italic", outfilename=os.path.join(document.figures_folder,legend_title+"_"+cat_metadata[0]+"_average_taxonomy.png"),legend_reverse=True)

def plot_stacked_barchart_taxonomy(document, samples, taxonomy, data, max_sets_barplot, taxonomy_level):
    # for the taxonomy data, organize and then plot the stacked barchart

    top_taxonomy, top_data = utilities.top_rows(taxonomy, data, max_sets_barplot, function="average")
    sorted_data, sorted_samples = sort_data(document, top_data, samples)

    top_taxonomy, sorted_data = fill_taxonomy_other(top_taxonomy, sorted_data)

    document.plot_stacked_barchart(sorted_data, row_labels=top_taxonomy,
        column_labels=sorted_samples, title="Top "+str(max_sets_barplot)+" "+taxonomy_level+" by average abundance",
        ylabel="Relative abundance", legend_title=taxonomy_level[0].upper()+taxonomy_level[1:], legend_style="italic", 
        outfilename=os.path.join(document.figures_folder,taxonomy_level+"_average_abundance.png"),legend_reverse=True)

    return sorted_samples, sorted_data, top_taxonomy


def sort_data(document, top_data, samples, sort_by_name=False, sort_by_name_inverse=False):
    # sort the top data so it is ordered with the top sample/abundance first
    if sort_by_name:
        sorted_sample_indexes=[samples.index(a) for a in document.sorted_data_numerical_or_alphabetical(samples)]
    elif sort_by_name_inverse:
        inverse_samples = [sample[::-1] for sample in samples]
        sorted_sample_indexes=[inverse_samples.index(a) for a in document.sorted_data_numerical_or_alphabetical(inverse_samples)]
    else:
        sorted_sample_indexes=sorted(range(len(samples)),key=lambda i: top_data[0][i],reverse=True)

    sorted_samples=[samples[i] for i in sorted_sample_indexes]
    sorted_data=[]
    for row in top_data:
        sorted_data.append([row[i] for i in sorted_sample_indexes])
    return sorted_data, sorted_samples

def plot_grouped_taxonomy_subsets(document, sorted_data, cat_metadata, top_taxonomy, samples_found, title, 
    ylabel="Relative abundance", legend_title="Species", legend_size=7, max_subsets=2, sort_by_name=False,
    sort_by_name_inverse=False, feature="", max_groups_barplot=5):
    """ Plot the grouped taxonomy with samples sorted by species abundance for each feature.

        Args:
            document (anadama2.document): Document to use to add plots
            sorted_data (list): A list of sorted/ordered data (to match the sample names and metadata)
            cat_metadata (list): A list of metadata (feature) names
            top_taxonomy (list): A list of the taxonomy names ordered to match the data
            samples_found (list): An ordered list of sample names (to match the data/metadata)
            title (string): The base string for the title for all plots
            ylabel (string): The label for the y-axis
            legend_title (string): The legend title for the plots
            legend_size (int): The font size for the legend
            max_subsets (int): The max number of subsets to use for each plot
            sort_by_name (bool): If set, then sort by sample name instead of abundance
            sort_by_name_inverse (bool): If true, sort by the inverse of the name (so the reverse of the string)
                this is useful for samples with sample name plus features

        Return: None
    """

    # group the samples by metadata
    sorted_data_grouped, sorted_samples_grouped = utilities.group_samples_by_metadata(cat_metadata, sorted_data, samples_found)
    # sort the data by abundance or sample name
    for metadata_type in sorted_data_grouped:
        sorted_data_grouped[metadata_type], sorted_samples_grouped[metadata_type] = sort_data(document, sorted_data_grouped[metadata_type], sorted_samples_grouped[metadata_type], sort_by_name=sort_by_name, sort_by_name_inverse=sort_by_name_inverse)

    # print out a plot for each group of metadata
    sorted_metadata_subsets=document.sorted_data_numerical_or_alphabetical(list(sorted_data_grouped.keys()))
   
    # split into subsets
    split_sorted_metadata_subsets = [sorted_metadata_subsets[x:x+max_subsets] for x in range(0, len(sorted_metadata_subsets), max_subsets)]
   
    # make sure the last group is not just a single data set
    if len(split_sorted_metadata_subsets[-1]) == 1:
        last_set = split_sorted_metadata_subsets.pop()
        split_sorted_metadata_subsets[-1].append(last_set[0])

    # set a max number of subsets
    index=0
    for metadata_subset in split_sorted_metadata_subsets:
        subset_sorted_data_grouped=dict((key, sorted_data_grouped[key]) for key in metadata_subset)
        subset_sorted_samples_grouped=dict((key, sorted_samples_grouped[key]) for key in metadata_subset)

        # get title addition for subset
        title_add=""
        if len(sorted_metadata_subsets) > max_subsets:
            title_add=" "+metadata_subset[0]+" to "+metadata_subset[-1]

        document.plot_stacked_barchart_grouped(subset_sorted_data_grouped, row_labels=top_taxonomy,
            column_labels_grouped=subset_sorted_samples_grouped, title=title+" - "+str(cat_metadata[0])+title_add,
            ylabel=ylabel, legend_title=legend_title, legend_style="italic", legend_size=legend_size, outfilename=os.path.join(document.figures_folder,"grouped_taxonomy_"+feature+"_"+str(cat_metadata[0])+"_"+str(metadata_subset[0].replace("-","_"))+".png"),legend_reverse=True)
        index+=1

        if index>= max_groups_barplot:
            break

def zscore_heatmap(document, dna_samples, dna_top_average_pathways, dna_top_average_data, merged_data, metadata_pathways, metadata_samples, data_type="pathways"):
    # if there is metadata, add it to the heatmap
    if 'metadata' in vars and vars['metadata'] and 'metadata_labels' in vars and vars['metadata_labels']:
        # get the total number of features
        total_features=len(vars['metadata'])-1
        # for the zscore to be applied first all of the categorical data must be removed
        filtered_metadata_pathways=[]
        filtered_merged_data=[]
        filtered_row_count=0
        for i in range(total_features):
            label = vars['metadata_labels'].get(metadata_pathways[i],"cat")
            if label != "cat":
                filtered_metadata_pathways.append(metadata_pathways[i])
                filtered_merged_data.append(merged_data[i])
            else:
                filtered_row_count+=1
        filtered_metadata_rows=range(1,total_features-filtered_row_count+1)

        # add abundance data to filtered metadata
        filtered_metadata_pathways+=metadata_pathways[total_features:]
        filtered_merged_data+=merged_data[total_features:]

        document.show_hclust2(metadata_samples, filtered_metadata_pathways, filtered_merged_data,
            title="Top "+str(max_sets_heatmap)+" "+data_type+" by average abundance",
            log_scale=False,zscore=True,
            metadata_rows=filtered_metadata_rows)
    else:
        document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
            title="Top "+str(max_sets_heatmap)+" "+data_type+" by average abundance",
            log_scale=False,zscore=True)


def log10_heatmap(document, dna_samples, dna_top_average_pathways, dna_top_average_data, data_type="pathways"):
    merged_data=[]
    metadata_pathways=[]
    metadata_samples=[]
    if 'metadata' in vars and vars['metadata']:
        merged_data, metadata_samples=merge_metadata(vars['metadata'], dna_samples,
            [[dna_top_average_pathways[i]]+dna_top_average_data[i] for i in range(len(dna_top_average_pathways))])
        metadata_pathways=[row.pop(0) for row in merged_data]
        # get the metadata row numbers
        metadata_rows=range(1,len(vars['metadata']))
        document.show_hclust2(metadata_samples, metadata_pathways, merged_data,
            title="Top "+str(max_sets_heatmap)+" "+data_type+" by average abundance",
            metadata_rows=metadata_rows)
    else:
        document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
            title="Top "+str(max_sets_heatmap)+" "+data_type+" by average abundance")

    return merged_data, metadata_pathways, metadata_samples

def remove_large_metadata_levels(merged_data, row_names, total_metadata_variables, max_levels):
    new_data=[]
    new_names=[]
    total_metadata=0
    for index in range(0, total_metadata_variables):
        grouped_values = list(set(merged_data[index]))
        if len(grouped_values) <= max_levels:
            new_data.append(merged_data[index])
            new_names.append(row_names[index])
            total_metadata+=1

    for index in range(total_metadata_variables, len(merged_data)):
        new_data.append(merged_data[index])
        new_names.append(row_names[index])

    return total_metadata, new_data, new_names

def plot_heatmap(document,vars,samples,top_taxonomy,top_data,pdf_format,filename,title=None,max_sets_heatmap=25,method="correlation",zscore=False):
    """ Generate a heatmap using the doc function. Include metadata if available. """

    # set the default title if not provided
    if not title:
        title = "Top {} species by average abundance".format(max_sets_heatmap)

    # if there is metadata, add it to the top taxonomy data
    if 'metadata' in vars and vars['metadata']:
        merged_data, metadata_samples=utilities.merge_metadata(vars['metadata'], samples,
            [[top_taxonomy[i]]+top_data[i] for i in range(len(top_taxonomy))])
        metadata_taxonomy=[row.pop(0) for row in merged_data]
       

        total_metadata, merged_data, metadata_taxonomy=remove_large_metadata_levels(merged_data, metadata_taxonomy, len(vars['metadata']), max_sets_heatmap)
        # if the total metadata are reduced to remove large metadata levels add one to the range to allow for the sample ids
        if len(vars['metadata']) > total_metadata:
            total_metadata+=1
        metadata_rows=range(1,total_metadata)

        document.show_hclust2(metadata_samples, metadata_taxonomy, merged_data,
            title=title, metadata_rows=metadata_rows, method=method,outfilename=os.path.join(document.figures_folder,filename),zscore=zscore)
    else:
        document.show_hclust2(samples,top_taxonomy,top_data,title=title,method=method,outfilename=os.path.join(document.figures_folder,filename),zscore=zscore)


def plot_pcoa_top_average_abundance(document, samples, feature_names, feature_data, feature_type, scale_data=None, legend_title="% Abundance", max_sets=6):
    """ Plot multiple pcoa in a single figure for the top abundances for the feature set """

    # if function is provided, scale the abundance data
    if scale_data:
        new_data = []
        for row in feature_data:
            new_data.append(map(scale_data,row))
        feature_data = new_data

    # get the top features by average abundance
    top_names, top_data = utilities.top_rows(feature_names, feature_data, max_sets, function="average")

    top_abundances=dict((x,y) for x,y in zip(["#"+str(i+1)+" "+x for i,x in enumerate(top_names[:max_sets])], top_data[:max_sets]))
    document.show_pcoa_multiple_plots(samples, feature_names, feature_data,
        "PCoA Ordination of "+feature_type+", top "+str(max_sets)+" "+feature_type+" by average abundance", top_abundances, legend_title)
    
def qc_read_counts(document, file):
    """ Read in the file of read counts compiled from kneaddata logs with the utility script """
    
    columns, samples, data = document.read_table(file, format_data=int)
    
    # shorten known database names
    columns=[name.replace("SILVA_128_LSUParc_SSUParc_ribosomal_RNA","rRNA") for name in columns]
    columns=[name.replace("Homo_sapiens_hg38","hg38") for name in columns]
    columns=[name.replace("human_hg38_refMrna","mRNA") for name in columns]
    columns=[name.replace("hg37_and_human_contamination","hg37") for name in columns]
    
    # change the names of the raw and trimmed columns (to include case)
    columns=[name.replace("raw","Raw").replace("trimmed","Trim").replace("decontaminated","") for name in columns]
    
    # check if this is single or paired end data
    if list(filter(lambda x: "single" in x, columns)):
        # remove single from the column name
        columns=[name.replace("single","").strip() for name in columns[:-1]]
        
        # return all but the final filtered column as this is not needed
        data=[row[:-1] for row in data]
        
    else:
        # remove the final columns from the list as these are not needed
        columns=list(filter(lambda x: not "final" in x, columns))
        
        # organize the data into pairs and orphans
        pairs_index = [index for index, name in enumerate(columns) if "pair1" in name]
        pairs_columns = [name.replace("pair1","").strip() for index, name in enumerate(columns) if index in pairs_index]
        orphans_index = [index for index, name in enumerate(columns) if "orphan" in name]
        orphans_columns = [name for index, name in enumerate(columns) if index in orphans_index]
        
        pairs_data=[]
        orphans_data=[]
        for row in data:
            pairs_data.append([])
            orphans_data.append([])
            for index in range(len(row)):
                if index in pairs_index:
                    pairs_data[-1].append(row[index])
                if index in orphans_index:
                    orphans_data[-1].append(row[index])
                
        columns = pairs_columns, orphans_columns
        data = pairs_data, orphans_data
    
    return columns, samples, data
    

def feature_counts(document, read_counts_file, feature_counts_file):
    """ Compute feature counts from the humann log read counts file and the feature counts file """
    
    # read in the read count and feature count files
    read_type, read_samples, read_count_data = document.read_table(read_counts_file)
    feature_type, feature_samples, feature_count_data = document.read_table(feature_counts_file)
    
    # remove any samples for which the prescreen did not find any species so nucleotide search was bypassed
    # these samples will have NA read counts but could have non-zero species count (species is the last column)
    read_samples, read_count_data = utilities.filter_zero_rows(read_samples, read_count_data, ignore_index=-1)
    
    # get the total reads for samples along with those for nucleotide alignment and translated
    # convert values to log10

    total_reads=[utilities.try_log10(row[read_type.index("total reads")]) for row in read_count_data]
    nucleotide_reads=[utilities.try_log10(row[read_type.index("total nucleotide aligned")]) for row in read_count_data]
    translated_reads=[utilities.try_log10(row[read_type.index("total translated aligned")]) for row in read_count_data]
    
    # sort the feature counts so they are in the same sample order as the read counts
    all_feature_counts={sample:row for sample, row in zip(feature_samples, feature_count_data)}
    sorted_feature_count_data=[all_feature_counts[sample] for sample in read_samples]
    
    # get the counts by each feature type
    # convert values to log10
    genefamilies_counts=[utilities.try_log10(row[feature_type.index("humann_genefamilies_relab_counts")]) for row in sorted_feature_count_data]
    ecs_counts=[utilities.try_log10(row[feature_type.index("humann_ecs_relab_counts")]) for row in sorted_feature_count_data]
    pathabundance_counts=[utilities.try_log10(row[feature_type.index("humann_pathabundance_relab_counts")]) for row in sorted_feature_count_data]
    
    return total_reads, nucleotide_reads, translated_reads, genefamilies_counts, ecs_counts, pathabundance_counts

def write_pathway_average_variance_table(document, file_name, data, names_and_descriptions, format_table_decimal="{:.3}"):
    """ Write a file of the pathways including averages and variance """

    # get the average abundances and descriptions for the pathways
    top_average_pathways_file = os.path.join(document.data_folder, file_name)
    
    # get the average abundances, formatting as a single value per row 
    average_abundance_variance=[]
    for average, variance in zip(utilities.row_average(data), utilities.row_variance(data)):
        average_abundance_variance.append([format_table_decimal.format(average),format_table_decimal.format(variance)])
    
    document.write_table(["# Pathway","Average abundance", "Variance"], names_and_descriptions, 
                         average_abundance_variance, top_average_pathways_file)
    
    return average_abundance_variance

def top_average_pathways(document, file, max_sets, get_all=False, filter_correlation=False, correlation_threshold=0.7):
    """ Read the pathways file and get the top average pathways """
    
    # read in the samples and get the data with out the stratification by bug
    samples, pathways, data = document.read_table(file)
    pathway_names = utilities.pathway_names(pathways)
    pathways, data = utilities.remove_stratified_pathways(pathways, 
        data, remove_description=True)
    
    # remove extra identifier from sample name if included in workflow
    samples = [sample.replace("_Abundance","").replace("-RPKs","") for sample in samples]
    
    # get the average abundance for the pathways
    if get_all:
        top_pathways, top_data = pathways, data
    else:
        if filter_correlation:
            pathways, data = utilities.filter_correlation(pathways, data, correlation_threshold, max_return=max_sets*2)

        top_pathways, top_data = utilities.top_rows(pathways,
            data, max_sets, function="average")
    
    # get the top names with descriptions
    top_names_and_descriptions = [name+":"+pathway_names[name] for name in top_pathways]
    
    return samples, top_pathways, top_data, top_names_and_descriptions

def remove_unexpected_chars(value_list):
    return [re.sub(r"[^0-9a-zA-Z.\-\s/]+", "_",val) for val in value_list]

def show_table_max_rows(document, data, row_labels, column_labels, title, table_file,
    max_rows=20, format_data_comma=None, location="center", font=None, max_columns=7, outfilename=None):
    """ For large numbers of samples, only show a reduced table """
    
    table_message="A data file exists of this table: "
    large_table_message="The table is too large to include the full table in this document."+\
        " A partial table is shown which includes only {max} {item}."+\
        " Please see the data file for the full table: "
        
    # check if there are too many rows
    partial_table_rows=False
    if len(row_labels) > max_rows:
        data=data[:max_rows]
        row_labels=row_labels[:max_rows]
        partial_table_rows=True
        
    # check if there are too many columns
    partial_table_columns=False
    if len(column_labels) > max_columns-1:
        data=[row[:max_columns-1] for row in data]
        column_labels=column_labels[:max_columns-1]
        partial_table_columns=True
        
    # determine the message and title if this is a full or partial table
    if partial_table_rows or partial_table_columns:
        if partial_table_rows:
            message=large_table_message.format(max=max_rows,item="rows")
        else:
            message=large_table_message.format(max=max_columns,item="columns")
        title+=" (partial table)"
    else:
        message=table_message
        
    # render the table
    document.show_table(data, remove_unexpected_chars(row_labels), remove_unexpected_chars(column_labels), 
        title, format_data_comma=format_data_comma, location=location, font=font, outfilename=outfilename)
    
    message+="[{file}](data/{file})".format(file=os.path.basename(table_file))
        
    return message

def print_pathways_urls(names, descriptions, total):
    """ List pathways with urls, including descriptions """
    
    print("Detailed functions of the top {} pathways can be found on the following MetaCyc pages:  ".format(total))
    
    print("")
    for pathway, desc in zip(remove_unexpected_chars(names[:total]), remove_unexpected_chars(descriptions[:total])):
        print(" * ["+desc+"]("+utilities.metacyc_url(pathway)+")  ")
        
    print("")
    print("To learn more about other pathways, search for the pathway by name on the [MetaCyc](https://metacyc.org/) website.")

class Workflow(object):
    
    @classmethod
    def format_caption(cls,name,**keywords):
        return cls.captions[name].format(**keywords)


class ShotGun(Workflow):
    captions={}
    
    # add captions for functional data section
    captions["functional_intro"]="This report section contains preliminary "+\
        "exploratory figures that summarize HUMAnN functional profiling of "+\
        "all samples. HUMAnN performs species-specific and species-agnostic "+\
        " quantification of gene families, EC enzyme modules, and pathways, "+\
        "using the UniRef and MetaCyc databases. For more information on "+\
        "functional profiling and the databases used, see websites for "+\
        "[HUMAnN](http://huttenhower.sph.harvard.edu/humann), "+\
        "[UniRef](http://www.uniprot.org/help/uniref), "+\
        "and [MetaCyc](https://metacyc.org/)."
        
    captions["heatmap_intro"]="Hierarchical clustering of samples "+\
        "and {type}, using top {max_sets} {type} "+\
        "with highest mean relative abundance among samples. "+\
        "The 'average linkage' clustering on the Euclidean "+\
        "distance metric was used to cluster samples. The {type} "+\
        "dendrogram is based on pairwise ( {method} ) correlation between {data_type}. "+\
        "Samples are columns and {type} are rows. The heatmaps were generated "+\
        "with [Hclust2](https://bitbucket.org/nsegata/hclust2)."
        
    captions["feature_detection"]="Feature detection as a function of sequencing "+\
        "depth. Effect of sample sequencing depth on the ability to detect "+\
        "microbiome functional features in {seq_type} sequence data. HUMAnN "+\
        "functional profiling of {seq_short_type} quality filtered reads was performed on "+\
        "individual samples in species-specific mode (blue), i.e. nucleotide "+\
        "alignment against pangenomes of species identified in the sample "+\
        "with MetaPhlAn, and in combined species-specific and -agnostic "+\
        "(orange) mode, in which reads not matching any pangenome reference "+\
        "sequences were subjected to translated searching against the "+\
        "UniRef90 database. Each profiled sample is represented by a "+\
        "orange and blue point in each plot. Linear regression fit is "+\
        "represented by straight lines in each plot."
        
    captions["pathway_abundance_heatmap"]="Abundances were {norm} transformed "+\
        "prior to clustering. The color bar represents relative abundances on a {norm} scale."  
        
    captions["scatter_reads_aligned"]="Number of aligned reads in species-specific "+\
        "(nucleotide search) and species-agnostic (translated search) HUMAnN mode "+\
        "as a function of input reads."
        
    captions["scatter_features"]="Detection of UniRef90 gene families, enzyme modules,"+\
        " and pathways as a function of aligned reads."
        
    captions["microbial_ratios"]="Proportion of reads remaining after removing host"+\
        " reads relative to the number of: i) quality-trimmed reads, and ii) raw "+\
        "unfiltered reads."
        
    captions["qc_intro"]="This report section contains information about the "+\
        "quality control processing for all {total_samples} {seq_type} fastq input "+\
        "files. These files were run through the "+\
        "[KneadData](http://huttenhower.sph.harvard.edu/kneaddata) QC pipeline. "+\
        "Reads were first trimmed then filtered against contaminate reference database{dbs}. "
    
    captions["qc_intro_multiple_db"]="Reads were filtered sequentially "+\
        "with those reads passing the first filtering step used as input to the next "+\
        "filtering step. This chain of filtering removes reads from all references in serial."
        
    captions["qc_intro_paired"]="\n \nData is organized by paired "+\
        "and orphan reads. When one read in a pair passes a filtering step and "+\
        "the other does not the surviving read is an orphan."
        
    captions["qc_intro_table"]="\nThe tables and plots are annotated as follows:\n \n"+\
        " * raw : Untouched fastq reads.\n"+\
        " * trim : Number of reads remaining after trimming bases with Phred score < 20. If the "+\
        "trimmed reads is < 50% of original length then it is removed altogether.\n"
        
    # set descriptions for command qc databases
    captions["qc_databases"]={}
    captions["qc_databases"]["hg38"]="The human genome database is used to remove "+\
        "reads originating from the host DNA."
    captions["qc_databases"]["mRNA"]="The human transcriptome (hg38 mRNA) database "+\
        "is used to remove reads originating from host gene isoforms."
    captions["qc_databases"]["rRNA"]="The SILVA (rRNA) database is used to remove "+\
        " small and large subunit ribosomal RNA."
    
    @classmethod
    def print_qc_intro_caption(cls, total_samples, databases, paired=None):
        """ Generate the qc intro caption based on the samples and databases """
        
        caption=cls.captions["qc_intro"]
        
        # if there are multiple databases, add the description about sequential filtering
        if len(databases) > 1:
            caption+=cls.captions["qc_intro_multiple_db"]
        
        if paired:
            caption+=cls.captions["qc_intro_paired"]
        caption+=cls.captions["qc_intro_table"]
        
        # add each database to the list
        dbs_list=[]
        for db in databases:
            dbs_list.append(db)
            desc=cls.captions["qc_databases"].get(db,"")
            caption+=" * {db} : Number of reads after depleting against reference database {list}. {desc}\n".format(
                db=db,list=" and ".join(dbs_list), desc=desc)
        caption+=" \n"
        
        # get the sequence type string
        seq_type="single-end"
        if paired:
            seq_type="paired-end"
        
        # format the caption to include the specific details for this data set
        caption=caption.format(total_samples=total_samples, seq_type=seq_type,
            dbs="s: "+", ".join(databases[:-1])+" and "+databases[-1] if len(databases) > 1 else " "+databases[0])
        
        for line in caption.split("\n"):
            print(line)
        

class Sixteen_S(Workflow):
    @classmethod
    def compile_default_intro(cls,vars):
        # read through the workflow log to gather information and compile the default intro for the report

        # get the variable settings from the data processing workflow
        from anadama2.reporters import LoggerReporter
        from anadama2 import PweaveDocument
        document=PweaveDocument()

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
        usearch_db = workflow_settings.get("usearch_db", "UNK").lower()

        method = vars["method"]
        if method == "dada2" or method == "its":
            if method == "its":
                dada_db = "UNITE"
            dadadb_info="\nThe " + str(dada_db) + " database was used for taxonomy prediction."
    
        else:
            try:
                columns, samples, data = document.read_table(vars["read_count_table"])
            except IOError:
                columns, samples, data = "", [], []

            if "silva" in usearch_db:
                db_info = "SILVA database"
            else:
                db_info = "GreenGenes 16S RNA Gene Database version 13_8"

            samples_text=str(len(samples))+" " if samples else ""
            usearchintro="The " + samples_text + "  samples from this project were run through the standard 16S workflow.  \
                follows the UPARSE OTU analysis pipeline for OTU calling and taxonomy prediction with percent identity " \
                + str(percent_identity) + " and minimum cluster size of " + str(min_cluster_size) + "." \
                + "\n\nThe " + db_info + " was used for taxonomy prediction.\
                \n\nReads were filtered for quality control using a MAXEE score of " + str(maxee) + ". Filtered reads were \
                used to generate the OTUs. Reads not passing quality control were kept and used in the step \
                assigning reads to OTUs. First these reads were truncated to a max length of " + str(trunc_len_max) + " bases.\n"

        if method == "its":
            flowchart="![]({0})\n\n".format(utilities.get_package_file("dada2","image"))
            return flowchart+cls.captions["itsintro"]+"\n\n"+dadadb_info
        elif method == "dada2":
            flowchart="![]({0})\n\n".format(utilities.get_package_file("dada2","image"))
            return flowchart+cls.captions["dada2intro"]+"\n\n"+dadadb_info
        else:
            flowchart="![]({0})\n\n".format(utilities.get_package_file("16s_workflow","image"))
            return flowchart+usearchintro


    captions={}

    captions["itsintro"]="The [DADA2 ITS pipeline](https://benjjneb.github.io/dada2/ITS_workflow.html) is the DADA2 pipeline, Callahan BJ, et all (2016). “DADA2: High-resolution sample inference from Illumina amplicon data.” Nature Methods, 13, 581-583., with small changes to adapt to sequencing of the ITS region.\n\n\n" 

    captions["dada2intro"]="The [DADA2 pipeline](https://benjjneb.github.io/dada2/tutorial.html), Callahan BJ, et all (2016). “DADA2: High-resolution sample inference from Illumina amplicon data.” Nature Methods, 13, 581-583., resolves sequence variants from 16S rRNA to generate a amplicon sequence variant (ASV) table.\n\n\n"
    
    captions["dada2errorintro"]="The read quality profiles show the quality scores at each base for the forward and reverse reads. The green line shows the mean and the quartiles are shown with the orange lines.\n"
    
    captions["dada2countsinfo"]="This figure shows stacked counts of reads in four categories: \n \
        \n1) Original: total count of raw reads,\n2) Filtered: number of reads after filtering for length and quality,\n3) Merged: number of reads where the pairs merge,\n4) NonChimera: total remaining after chimera removal.\n\n"


    captions["usearchcountsinfo"]="This figure shows stacked counts of reads in three categories: \n \
        \n1) classified: reads that align to OTUs with known taxonomy,\n2) unclassified: reads that align to OTUs of unknown taxonomy,\n3) unmapped: reads that do not align to any OTUs.\n\n The sum of these\
        three read counts for each sample is the total original read count not including filtering prior to OTU clustering.\n"
    

class Stats(Workflow): 
    captions = {}

    captions["sections"]=['\n\n1. All against all : A mantel test, which computes the correlation between two matrices of the same dimension, is run to compare all points of each data set provided with all points of all other data sets. This computation is only run if there are multiple data sets (eg taxonomy and functional data).',
        '\n\n2. One against all : A PERMANOVA, used to compare groups of objects, is run with two different approaches to the statistical analysis for each of the data types, eg taxonomy, provided for the study. First the PERMANOVA is run to compare a single metadata variable with the data set. Next, in a multivariable analysis, all of the metadata variables are used for the analysis against the data set. For a longitudinal study, where there are multiple time points for each individual, the metadata co-variates that do not vary for each individual are factored into the analysis.',
        '\n\n3. Each metadata variable against all features individually: MaAsLin 2 filters, transforms, and then performs a linear model to fit metadata variables to feature data (e.g. taxonomy, pathways), one at a time. If both taxonomy and pathways are provided, plots for the significant pathways, stratified by species are generated. ',
        '\n\n4. Each data type against all metadata variables : HAllA tests all possible associations of each feature in a data set against all metadata variables. It is run to compare each of the data sets provided for the study against all variables. \n\n\n ']

    captions["intro"]='The data for this project was run through the standard stats workflow. The workflow is composed of four sections.\n\n'+"".join(captions["sections"])
    captions["intro_bypass_halla"]='The data for this project was run through the standard stats workflow. The workflow is composed of three sections.\n\n'+"".join(captions["sections"][0:3])
