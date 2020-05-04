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
import copy
import sys
import subprocess
import tempfile
import shutil

from . import utilities

def plot_grouped_and_average_barplots_taxonomy(document, vars, sorted_samples, sorted_data, top_taxonomy,
    max_sets_barplot, feature="species", sort_by_name=False, sort_by_name_inverse=False, ylabel="Relative abundance"):
    """ Plot grouped barplots and average barplots for all of the features provided.

        Args:
            document (anadama2.document): Document to use to add plots
            vars (dict): The dictionary of input variables provided for the visualization run
            sorted_samples (list): The sample names organized to match the data
            sorted_data (list): The data organized to match the sample/taxonomy
            top_taxonomy (list): The list of full taxonomy names
            max_sets_barplot (int): The max number of features (default species) to plot
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
                sort_by_name_inverse=sort_by_name_inverse)
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
    new_top_taxonomy = top_taxonomy + ["other"]
    new_sorted_data = sorted_data
    other_abundances=[]
    for column in numpy.transpose(sorted_data):
        other_abundances.append(100-sum(column))

    if isinstance(new_sorted_data,numpy.ndarray):
        new_sorted_data = new_sorted_data.tolist()
    new_sorted_data.append(other_abundances)

    return new_top_taxonomy, new_sorted_data


def show_pcoa_metadata(document, vars, samples, top_taxonomy, pcoa_data, title):
    """ Plot pcoa for each feature if metadata has been provided

        Args:
            document (anadama2.document): Document to use to add plots
            vars (dict): The dictionary of input variables provided for the visualization run
            samples (list): The sample names organized to match the data
            top_taxonomy (list): The full taxonomic names organized to match the data
            pcoa_data (list): The data (range 0 to 1) organized to match the samples and taxonomy
            title (string): The base string for the title for the plots
 
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
                               metadata=metadata_mapping, metadata_type=metadata_type)

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
    average_data = zip(*average_data)

    # sort the data by abundance or name (depending on metadata type, use name if numeric)
    sort_by_name=False
    if document.sorted_data_numerical_or_alphabetical(metadata_names) != sorted(metadata_names):
        sort_by_name = True

    sorted_data, sorted_names = sort_data(document, average_data, metadata_names, sort_by_name=sort_by_name)

    document.plot_stacked_barchart(sorted_data, row_labels=top_taxonomy,
        column_labels=sorted_names, 
        title="Top {} {} group average - {}".format(max_sets_barplot, legend_title, cat_metadata[0]),
        ylabel=ylabel, legend_title=legend_title, legend_style="italic")

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
    sort_by_name_inverse=False):
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
    sorted_metadata_subsets=document.sorted_data_numerical_or_alphabetical(sorted_data_grouped.keys())
   
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

        document.plot_stacked_barchart_grouped(subset_sorted_data_grouped, row_labels=top_taxonomy,
            column_labels_grouped=subset_sorted_samples_grouped, title=title+" - "+str(cat_metadata[0])+title_add,
            ylabel=ylabel, legend_title=legend_title, legend_style="italic", legend_size=legend_size)

def plot_heatmap(document,vars,samples,top_taxonomy,top_data,pdf_format,title=None,max_sets_heatmap=25,method="correlation"):
    """ Generate a heatmap using the doc function. Include metadata if available. """

    # set the default title if not provided
    if not title:
        title = "Top {} species by average abundance".format(max_sets_heatmap)

    # if there is metadata, add it to the top taxonomy data
    if 'metadata' in vars and vars['metadata']:
        merged_data, metadata_samples=utilities.merge_metadata(vars['metadata'], samples,
            [[top_taxonomy[i]]+top_data[i] for i in range(len(top_taxonomy))])
        metadata_taxonomy=[row.pop(0) for row in merged_data]
        # get the metadata row numbers
        metadata_rows=range(1,len(vars['metadata']))
        document.show_hclust2(metadata_samples, metadata_taxonomy, merged_data,
            title=title, metadata_rows=metadata_rows, method=method)
    else:
        document.show_hclust2(samples,top_taxonomy,top_data,title=title,method=method)


def plot_hallagram(document, feature_set_1_data, feature_set_2_data, axis1_label, axis2_label, model=None, title=None,
    output_folder=None, strongest=10, show_table=False,q_value=0.1):
    """ Run halla on the two feature sets and plot the hallagram with the strongest associations 
        The first data line for each feature set should be a header of sample ids
        Requires halla v0.8.7
        If show_table is set, then instead of including heatmap include a table of associations
    """

    from matplotlib._png import read_png
    import matplotlib.pyplot as pyplot

    # create a temp output folder, if not provided
    if output_folder:
        outfolder = output_folder
        # if folder does not exist, then create
        if not os.path.isdir(outfolder):
            os.makedirs(outfolder)
    else:
        outfolder=tempfile.mkdtemp(suffix="biobakery_workflows_halla",dir=output_folder)

    # write the lines for the two feature set files
    feature_set_1 = os.path.join(outfolder,"feature1.tsv")
    feature_set_2 = os.path.join(outfolder,"feature2.tsv")
    for lines, file_name in [[feature_set_1_data,feature_set_1],[feature_set_2_data,feature_set_2]]:
        with open(file_name,"w") as file_handle:
            file_handle.write("\n".join(["\t".join(map(str,l)) for l in lines]))

    # run halla
    halla_command = ["halla","-X", feature_set_1, "-Y", feature_set_2,"--output", outfolder, "--header"]
    if model:
        halla_command += ["-m",model]
    try:
        subprocess.check_call(halla_command)
    except subprocess.CalledProcessError:
        print("Error: Unable to run halla")
    
    if show_table:
        # display the table of associations instead of including the heatmap
        associations_file = os.path.join(outfolder,"associations.txt")
        # check that the file is not empty
        with open(associations_file) as file_handle:
            associations_total = len(file_handle.readlines())-1

        if associations_total > 0:
            columns, row_names, halla_data = document.read_table(associations_file, format_data=str)

            # reduce data to only those columns to display in the table
            cluster_results = [[row[0].replace(";","\n"),row[2].replace(";","\n"),"{0:.7f}".format(float(row[4])),"{0:.3f}".format(float(row[5]))] for row in halla_data]
            # filter any results that are more than the max q_value
            cluster_results = list(filter(lambda row: float(row[-1]) < q_value, cluster_results))
        
            row_labels = list(map(lambda x: str(x)+"   ",range(1,len(cluster_results)+1)))
            column_labels = [axis1_label, axis2_label,"  p-value  ","q-value"]
            document.show_table(cluster_results, row_labels, column_labels, title, font=6)
        else:
            print("No associations found for "+title)
    else:
        # run hallagram
        output_png = os.path.join(outfolder,"hallagram.png")
        hallagram_command = ["hallagram",os.path.join(outfolder,"similarity_table.txt"),
            os.path.join(outfolder,"hypotheses_tree.txt"), os.path.join(outfolder,"associations.txt"),
            "--outfile",output_png,"--strongest",str(strongest),"--axlabels",axis1_label,axis2_label]
        if model:
            hallagram_command += ["--similarity",model]
        try:
            subprocess.check_call(hallagram_command)
        except subprocess.CalledProcessError:
            print("Error: Unable to run hallagram")

        if os.path.isfile(output_png):
            hallagram_png=read_png(output_png)        

            # create a subplot and remove the frame and axis labels
            fig = pyplot.figure()
            subplot = fig.add_subplot(111, frame_on=False)
            subplot.xaxis.set_visible(False)
            subplot.yaxis.set_visible(False)
        
            # set title if provided
            if title:
                fig.suptitle(title, fontsize=16)

            # show but do not interpolate (as this will make the text hard to read)
            try:
                pyplot.imshow(hallagram_png, interpolation="none")
            except TypeError:
                print("Unable to display hallagram plot")
                pass

            # this is needed to increase the image size (to fit in the increased figure)
            pyplot.tight_layout()

    # remove the temp folder
    if not output_folder:
        shutil.rmtree(outfolder) 

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
    columns=[name.replace("Homo_sapiens","hg38") for name in columns]
    columns=[name.replace("human_hg38_refMrna","mRNA") for name in columns]
    
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

def top_average_pathways(document, file, max_sets):
    """ Read the pathways file and get the top average pathways """
    
    # read in the samples and get the data with out the stratification by bug
    samples, pathways, data = document.read_table(file)
    pathway_names = utilities.pathway_names(pathways)
    pathways, data = utilities.remove_stratified_pathways(pathways, 
        data, remove_description=True)
    
    # remove extra identifier from sample name if included in workflow
    samples = [sample.replace("_Abundance","").replace("-RPKs","") for sample in samples]
    
    # get the average abundance for the pathways
    top_pathways, top_data = utilities.top_rows(pathways,
        data, max_sets, function="average")
    
    # get the top names with descriptions
    top_names_and_descriptions = [name+":"+pathway_names[name] for name in top_pathways]
    
    return samples, top_pathways, top_data, top_names_and_descriptions

def show_table_max_rows(document, data, row_labels, column_labels, title, table_file,
    max_rows=20, format_data_comma=None, location="center", font=None, max_columns=7):
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
    document.show_table(data, row_labels, column_labels, 
        title, format_data_comma=format_data_comma, location=location, font=font)
    
    message+="[{file}](data/{file})".format(file=os.path.basename(table_file))
        
    return message

def print_pathways_urls(names, descriptions, total):
    """ List pathways with urls, including descriptions """
    
    print("Detailed functions of the top {} pathways can be found on the following MetaCyc pages:  ".format(total))
    
    print("")
    for pathway, desc in zip(names[:total], descriptions[:total]):
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
        "exploratory figures that summarize HUMAnN2 functional profiling of "+\
        "all samples. HUMAnN2 performs species-specific and species-agnostic "+\
        " quantification of gene families, EC enzyme modules, and pathways, "+\
        "using the UniRef and MetaCyc databases. For more information on "+\
        "functional profiling and the databases used, see websites for "+\
        "[HUMAnN2](http://huttenhower.sph.harvard.edu/humann), "+\
        "[UniRef](http://www.uniprot.org/help/uniref), "+\
        "and [MetaCyc](https://metacyc.org/)."
        
    captions["heatmap_intro"]="Hierarchical clustering of samples "+\
        "and {type}, using top {max_sets} {type} "+\
        "with highest mean relative abundance among samples. "+\
        "The 'average linkage' clustering on the Euclidean "+\
        "distance metric was used to cluster samples. The {type} "+\
        "dendrogram is based on pairwise ( {method} ) correlation between pathways. "+\
        "Samples are columns and pathway are rows. The heatmaps were generated "+\
        "with [Hclust2](https://bitbucket.org/nsegata/hclust2)."
        
    captions["feature_detection"]="Feature detection as a function of sequencing "+\
        "depth. Effect of sample sequencing depth on the ability to detect "+\
        "microbiome functional features in {seq_type} sequence data. HUMAnN2 "+\
        "functional profiling of {seq_short_type} quality filtered reads was performed on "+\
        "individual samples in species-specific mode (blue), i.e. nucleotide "+\
        "alignment against pangenomes of species identified in the sample "+\
        "with MetaPhlAn2, and in combined species-specific and -agnostic "+\
        "(orange) mode, in which reads not matching any pangenome reference "+\
        "sequences were subjected to translated searching against the "+\
        "UniRef90 database. Each profiled sample is represented by a "+\
        "orange and blue point in each plot. Linear regression fit is "+\
        "represented by straight lines in each plot."
        
    captions["pathway_abundance_heatmap"]="Abundances were {norm} transformed "+\
        "prior to clustering. The color bar represents relative abundances on a {norm} scale."  
        
    captions["scatter_reads_aligned"]="Number of aligned reads in species-specific "+\
        "(nucleotide search) and species-agnostic (translated search) HUMAnN2 mode "+\
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
    captions={}

    captions["itsintro"]='Implementing ITS pipeline for resolving sequence variants from ITS region\
         of paired-end sequencing reads, adopting the tutorial from\n\n \
         https://benjjneb.github.io/dada2/ITS_workflow.html \n \
         https://benjjneb.github.io/dada2/tutorial.html and \n \
         https://benjjneb.github.io/dada2/bigdata_paired.html \n\n with minor adjustments. \n \
         "Unlike the 16S rRNA gene, the ITS region is highly variable in length. The commonly amplified ITS1 and ITS2 regions \
         range from 200 - 600 bp in length. This length variation is biological, not technical, and arises from the high rates\
         of insertions and deletions in the evolution of this less conserved gene region. The length variation of the ITS region has \
         significant consequences for the filtering and trimming steps of the standard DADA2 workflow. First, truncation to a fixed \
         length is no longer appropriate, as that approach remove real ITS variants with lengths shorter than the truncation length.\
         Second, primer removal is complicated by the possibility of some, but not all, reads extending into the opposite primer when\
         the amplified ITS region is shorter than the read length." [ITS Tutorial] \n \
         "Critical addition to ITS workflows is the removal of primers on the forward and reverse reads, in a way that accounts\
         for the possibility of read-through into the opposite primer." [ITS Tutorial] \n \
         "Cutadapt" tool is used  for removal of primers from the ITS amplicon sequencing data. After initial step of primers \
         removal, the rest of ITS workflow matches DADA2 workflow.\n \
         Database UNITE is used for taxonomic assignment.\n\n '


    captions["dada2intro"]="Implementing DADA2 pipeline for resolving sequence variants from 16S rRNA \
        gene amplicon paired-end sequencing reads, adopting the tutorial from\n\n \
         https://benjjneb.github.io/dada2/tutorial.html and \n \
         https://benjjneb.github.io/dada2/bigdata_paired.html \n\n with minor adjustments.\
        \n\nThis report captures all the workflow steps necessary to reproduce the analysis. Notes and descriptions\
         of the steps are cited from DADA2 tutorial as well.\
        \n\nMultiple sequence alignment of resolved sequence variants is used to generate a phylogenetic tree,\
        which is required for calculating UniFrac beta-diversity distances between microbiome samples.\n\n\n"
    
    captions["dada2errorintro"]='\n\n "The DADA2 algorithm makes use of a parametric error model (err) and every amplicon dataset has a different set of error rates. \
        \n\nThe learnErrors method learns the error model from the data, by alternating estimation of the error rates and inference of \
        sample composition until they converge on a jointly consistent solution.\n\nAs in many machine-learning problems, the algorithm must \
        begin with an initial guess, for which the maximum possible error rates in this data are used \
        the error rates if only the most abundant sequence is correct and all the rest are errors." [DADA2 Tutorial]\n'
    
    captions["dada2errorinfo"]='"The error rates for each possible transition (eg. A->T,A->G, etc) are shown. \
        Points are the observed error rates for each consensus quality score. \
        \n\nThe black line shows the estimated error rates after convergence. \
         The red line shows the error rates expected under the nominal definition of the Q-value. \
        \n\nIf the black line (the estimated rates) fits the observed rates well, \
        and the error rates drop with increased quality as expected, then everything looks reasonable \
        and can proceed with confidence." [DADA2 Tutorial]\n\n'

    captions["dada2stepsinfo"]='\n\n "Dereplication combines all identical sequencing reads into into unique sequences with a corresponding abundance:\
        the number of reads with that unique sequence .... DADA2 retains a summary of the quality information associated with each unique sequence." [DADA2 Tutorial]\
        \n\nThe consensus quality profile of a unique sequence is the average of the positional qualities from the dereplicated reads.\
        \n\nThese quality profiles inform the error model of the subsequent denoising step, significantly increasing DADA2s accuracy. \
        \n\nThe sample inference step performs the core sequence-variant inference algorithm to the dereplicated data. \
        Spurious sequence variants are further reduced by merging overlapping reads. \n\nThe core function here is mergePairs \
        which depends on the forward and reverse re.samples being in matching order at the time they were dereplicated \
        \n\n "The core dada method removes substitution and indel errors, but chimeras remain.\
        Fortunately, the accuracy of the sequences after denoising makes identifying chimeras simpler than it is when dealing with fuzzy OTUs\
        all sequences which can be exactly reconstructed as a bimera (two-parent chimera) from more abundant sequence." [DADA2 Tutorial]' \
         + '\n\n "Most of reads should remain after chimera removal (it is not uncommon for a majority of sequence variants to be removed though).\
        If most of your reads were removed as chimeric, upstream processing may need to be revisited.\
        In almost all cases this is caused by primer sequences with ambiguous nucleotides that were not removed prior to beginning the DADA2 pipeline." [DADA2 Tutorial]\n'
    
    captions["dada2countsinfo"]='This figure shows the number of reads that made it through each step in the pipeline\
        \n\nThere should no be a step in which a majority of reads are lost, except filtering when it is stringent.\
        \n\n "If a majority of reads failed to merge, you may need to revisit the  truncLen parameter used in the filtering step\
        and make sure that the truncated reads span your amplicon." [DADA2 Tutorial]\
        \n\nIf a majority of reads failed to pass the chimera check, you may need to revisit the removal of primers,\
        as the ambiguous nucleotides in unremoved primers interfere with chimera identification.\n'

    captions["dada2taxinfo"]='"The assignTaxonomy function takes a set of sequences and a training set of taxonomically classified sequences,\
        and outputs the taxonomic assignments with at least minBoot bootstrap confidence." [DADA2 Tutorial]\
        \n\n assignTaxonomy(... ) implements the RDP naive Bayesian classifier method described in Wang et al. 2007.'  \
        + " In short, the kmer profile of the sequences to be classified are compared against the kmer profiles of all sequences in a training set\
        of sequences with assigned taxonomies. The reference sequence with the most similar profile is used to assign taxonomy to the query sequence,\
        and then a bootstrapping approach is used to assess the confidence assignment at each taxonomic level.\n \n"
 
          
    captions["usearchcountsinfo"]="This figure shows counts of reads in three categories: \n \
        \n1) classified: reads that align to OTUs with known taxonomy,\n2) unclassified: reads that align to OTUs of unknown taxonomy,\n3) unmapped: reads that do not align to any OTUs.\n\n The sum of these\
        three read counts for each sample is the total original read count not including filtering prior to OTU clustering.\n"
     
