
#+ echo=False
max_sets=50

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

# set the max number of table rows to be shown
# and set the messages to print for each table size
max_table_rows=20
format_table_decimal="{:.3}"
table_message="A data file exists of this table: "
large_table_message="The table is too large to include the full table in this document."+\
    " A partial table is shown which includes only "+str(max_table_rows)+" pathways."+\
    " Please see the data file for the full table: "

#' # Functional Profiling


#' This report section contains preliminary exploratory figures that summarize 
#' HUMAnN2 functional profiling of all samples. HUMAnN2 performs species-specific and species-agnostic 
#' quantification of gene families, EC enzyme modules, and pathways, using the UniRef and MetaCyc databases.
#' For more information on functional profiling and the databases used, see websites for 
#' [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2), 
#' [UniRef](http://www.uniprot.org/help/uniref), 
#' and [MetaCyc](https://metacyc.org/). 

#+ echo=False

# read in the DNA samples and get the data with out the stratification by bug
dna_samples, dna_pathways, dna_data = document.read_table(vars["dna_pathabundance"])
dna_pathway_names = utilities.pathway_names(dna_pathways)
dna_pathways, dna_data = utilities.remove_stratified_pathways(dna_pathways, 
    dna_data, remove_description=True)

# remove extra identifier from sample name if included in workflow
dna_samples = [sample.replace("_Abundance","") for sample in dna_samples]

# get the average abundance for the pathways
dna_top_average_pathways, dna_top_average_data = utilities.top_rows(dna_pathways,
    dna_data, max_sets, function="average")

# get the variance for the pathways
dna_top_variance_pathways, dna_top_variance_data = utilities.top_rows(dna_pathways,
    dna_data, max_sets, function="variance")

#' ## Pathway Abundance

#' Hierarchical clustering of samples and pathways, using top <%= max_sets %> pathways 
#' with highest mean relative abundance among samples. The "average linkage" 
#' clustering on the Euclidean distance metric was used to cluster samples.
#' The pathway dendrogram is based on pairwise (Spearman) correlation between pathways.
#' Samples are columns and pathway are rows. The heatmaps were generated with [Hclust2](https://bitbucket.org/nsegata/hclust2).

#+ echo=False
document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
                      title="Top "+str(max_sets)+" pathways by average abundance")  

#' Abundances were log10 transformed prior to clustering. The color bar 
#' represents relative abundances on a log10 scale.  

#+ echo=False
document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
                      title="Top "+str(max_sets)+" pathways by average abundance",
                      log_scale=False,zscore=True)

#' Abundances were z-score transformed prior to clustering. The color bar 
#' represents relative abundances on a z-score scale.  

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
# get the average abundances and descriptions for the pathways
top_average_pathways_file = os.path.join(document.data_folder,"top_average_pathways_names.tsv")
top_names_and_descriptions = [name+":"+dna_pathway_names[name] for name in dna_top_average_pathways]
# get the average abundances, formatting as a single value per row 
average_abundance_variance=[]
for average, variance in zip(utilities.row_average(dna_top_average_data), utilities.row_variance(dna_top_average_data)):
    average_abundance_variance.append([format_table_decimal.format(average),format_table_decimal.format(variance)])

document.write_table(["# Pathway","Average abundance", "Variance"], top_names_and_descriptions, average_abundance_variance, top_average_pathways_file)

if len(top_names_and_descriptions) <= max_table_rows:
    document.show_table(average_abundance_variance, top_names_and_descriptions, [" Average "," Variance "], 
        "Top "+str(max_sets)+" pathways by average abundance", font=7)
else:
    document.show_table(average_abundance_variance[:max_table_rows], top_names_and_descriptions[:max_table_rows], 
        [" Average "," Variance "], "Top "+str(max_sets)+" pathways by average abundance (partial table)", font=7)
    
#' <% print(large_table_message) if len(dna_top_average_pathways) > max_table_rows else print(table_message) %>
#' [top_average_pathways_names.tsv](data/top_average_pathways_names.tsv)  
        
#' Detailed functions of the top three pathways can be found on the following MetaCyc pages:  

#' <%= "  \r".join([" * ["+full_name+"]("+utilities.metacyc_url(path)+")" for path, full_name in zip(dna_top_average_pathways[:3],top_names_and_descriptions[:3])]) %>  

#' To learn more about other pathways, search for the pathway by name on the [MetaCyc](https://metacyc.org/) website.

#' <% if pdf_format: print("\clearpage") %>

#' ## Features

#' Feature detection as a function of sequencing depth. Effect of sample sequencing depth on the ability to detect 
#' microbiome functional features in metagenomic sequence data. HUMAnN2 functional profiling of 
#' DNA quality filtered reads was performed on individual samples in species-specific mode (blue), i.e. nucleotide 
#' alignment against pangenomes of species identified in the sample with MetaPhlAn2, and in combined species-specific 
#' and -agnostic (orange) mode, in which reads not matching any pangenome reference sequences were subjected to translated 
#' searching against the UniRef90 database. Each profiled sample is represented by a orange and blue point in each plot. 
#' Linear regression fit is represented by straight lines in each plot.

#+ echo=False
import math

# read in the read count and feature count files
read_type, read_samples, read_count_data = document.read_table(vars["read_counts"])
feature_type, feature_samples, feature_count_data = document.read_table(vars["feature_counts"])

# remove any samples for which the prescreen did not find any species so nucleotide search was bypassed
# these samples will have NA read counts but could have non-zero species count (species is the last column)
read_samples, read_count_data = utilities.filter_zero_rows(read_samples, read_count_data, ignore_index=-1)

# get the total reads for samples along with those for nucleotide alignment and translated
# convert values to log10

def try_log10(value):
    """ Try to convert value to log10 """
    
    try:
        new_value = math.log10(value)
    except ValueError:
        new_value = 0
        
    return new_value

total_reads=[try_log10(row[read_type.index("total reads")]) for row in read_count_data]
nucleotide_reads=[try_log10(row[read_type.index("total nucleotide aligned")]) for row in read_count_data]
translated_reads=[try_log10(row[read_type.index("total translated aligned")]) for row in read_count_data]

# sort the feature counts so they are in the same sample order as the read counts
feature_counts={sample:row for sample, row in zip(feature_samples, feature_count_data)}
sorted_feature_count_data=[feature_counts[sample] for sample in read_samples]

# get the counts by each feature type
# convert values to log10
genefamilies_counts=[try_log10(row[feature_type.index("humann2_genefamilies_relab_counts")]) for row in sorted_feature_count_data]
ecs_counts=[try_log10(row[feature_type.index("humann2_ecs_relab_counts")]) for row in sorted_feature_count_data]
pathabundance_counts=[try_log10(row[feature_type.index("humann2_pathabundance_relab_counts")]) for row in sorted_feature_count_data]

# add scatter plots of the data
document.plot_scatter([[total_reads,nucleotide_reads],[total_reads,translated_reads]],title="Read alignment rate",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Input reads)", ylabel="log10(Aligned reads)", trendline=True)

#' Number of aligned reads in species-specific (nucleotide search) and species-agnostic (translated search) HUMAnN2 mode as a function of input reads.

#+ echo=False
document.plot_scatter([[nucleotide_reads,genefamilies_counts],[translated_reads,genefamilies_counts]],title="UniRef90 gene families",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(gene families)", trendline=True)

document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="Enzymes (ECs)",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)

document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="Pathways",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

#' Detection of UniRef90 gene families, enzyme modules, and pathways as a function of aligned reads. 
