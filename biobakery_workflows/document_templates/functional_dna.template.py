
#+ echo=False
max_sets=50

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' # Functional Profiling

#' This report section contains information about the functional profiling run
#' on all samples. These samples were
#' run through [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2).
#' The UniRef90 full database was used for the translated search.

#' The pathways identified in functional profiling are from the MetaCyc database.
#' See the [MetaCyc website](https://metacyc.org/) for detailed information
#' on each pathway.

#+ echo=False

# read in the DNA samples and get the data with out the stratification by bug
dna_samples, dna_pathways, dna_data = document.read_table(vars["dna_pathabundance"])
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

#' ## Heatmaps

#' The top <% print(max_sets) %> pathways based on average relative abundance are
#' shown in the heatmaps. The heatmaps were generated with [Hclust2](https://bitbucket.org/nsegata/hclust2).

#+ echo=False
document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
                      title="Top "+str(max_sets)+" pathways by average abundance")

#+ echo=False
document.show_hclust2(dna_samples,dna_top_variance_pathways,dna_top_variance_data,
                      title="Top "+str(max_sets)+" pathways by variance")

#' <% if pdf_format: print("\clearpage") %>

#' ## Features

#' The total number of reads used for functional profiling along with the total
#' number of reads aligned at the nucleotide and translated search steps are shown.
#' They are plotted against the total number of features identified for each sample.
#' The features include gene families, ECs, and pathways. The feature counts do not
#' include stratification levels.

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
genefamilies_counts=[try_log10(row[feature_type.index("humann2_ecs_relab_counts")]) for row in sorted_feature_count_data]
ecs_counts=[try_log10(row[feature_type.index("humann2_ecs_relab_counts")]) for row in sorted_feature_count_data]
pathabundance_counts=[try_log10(row[feature_type.index("humann2_pathabundance_relab_counts")]) for row in sorted_feature_count_data]

# add scatter plots of the data
document.plot_scatter([[total_reads,nucleotide_reads],[total_reads,translated_reads]],title="Read alignment rate",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Input reads)", ylabel="log10(Aligned reads)", trendline=True)

#+ echo=False
document.plot_scatter([[nucleotide_reads,genefamilies_counts],[translated_reads,genefamilies_counts]],title="UniRef90 gene families",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(gene families)", trendline=True)

document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="Enzymes (ECs)",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)

document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="Pathways",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

