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

from . import utilities

def feature_counts(document, read_counts_file, feature_counts_file):
    """ Compute feature counts from the humann2 log read counts file and the feature counts file """
    
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
    genefamilies_counts=[utilities.try_log10(row[feature_type.index("humann2_genefamilies_relab_counts")]) for row in sorted_feature_count_data]
    ecs_counts=[utilities.try_log10(row[feature_type.index("humann2_ecs_relab_counts")]) for row in sorted_feature_count_data]
    pathabundance_counts=[utilities.try_log10(row[feature_type.index("humann2_pathabundance_relab_counts")]) for row in sorted_feature_count_data]
    
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
    samples = [sample.replace("_Abundance","") for sample in samples]
    
    # get the average abundance for the pathways
    top_pathways, top_data = utilities.top_rows(pathways,
        data, max_sets, function="average")
    
    # get the top names with descriptions
    top_names_and_descriptions = [name+":"+pathway_names[name] for name in top_pathways]
    
    return samples, top_pathways, top_data, top_names_and_descriptions

def show_table_max_rows(document, data, row_labels, column_labels, title, table_file,
    max_rows=20, format_data_comma=None, location="center", font=None):
    """ For large numbers of samples, only show a reduced table """
    
    table_message="A data file exists of this table: "
    large_table_message="The table is too large to include the full table in this document."+\
        " A partial table is shown which includes only "+str(max_rows)+" rows."+\
        " Please see the data file for the full table: "
        
    if len(row_labels) <= max_rows:
        document.show_table(data, row_labels, column_labels, 
            title, format_data_comma=format_data_comma, location=location, font=font)
    else:
        document.show_table(data[:max_rows], row_labels[:max_rows], 
            column_labels, title+" (partial table)", format_data_comma=format_data_comma, 
            location=location, font=font)
    
    if len(row_labels) <= max_rows:
        message=table_message
    else:
        message=large_table_message
        
    message+="[{file}](data/{file})".format(file=table_file)
        
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
    
    # add captions for the qc templates
    captions["qc_intro"]="This report section contains information about the "+\
        "quality control processing for all {total_samples} {seq_type} fastq input "+\
        "files. These files were run through the "+\
        "[KneadData](http://huttenhower.sph.harvard.edu/kneaddata) QC pipeline."
    
    # add captions for functional data section
    captions["functional_intro"]="This report section contains preliminary "+\
        "exploratory figures that summarize HUMAnN2 functional profiling of "+\
        "all samples. HUMAnN2 performs species-specific and species-agnostic "+\
        " quantification of gene families, EC enzyme modules, and pathways, "+\
        "using the UniRef and MetaCyc databases. For more information on "+\
        "functional profiling and the databases used, see websites for "+\
        "[HUMAnN2](http://huttenhower.sph.harvard.edu/humann2), "+\
        "[UniRef](http://www.uniprot.org/help/uniref), "+\
        "and [MetaCyc](https://metacyc.org/)."
        
    captions["pathway_abundance_intro"]="Hierarchical clustering of samples "+\
        "and pathways, using top {max_sets} pathways "+\
        "with highest mean relative abundance among samples. "+\
        "The 'average linkage' clustering on the Euclidean "+\
        "distance metric was used to cluster samples. The pathway "+\
        "dendrogram is based on pairwise (Spearman) correlation between pathways."+\
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
        
        

