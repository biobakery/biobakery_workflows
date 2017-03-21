
#+ echo=False
max_sets=50

from biobakery_workflows import utilities

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

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

#+ echo=False
document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
                      title="Top "+str(max_sets)+" pathways by average abundance")

#+ echo=False
document.show_hclust2(dna_samples,dna_top_variance_pathways,dna_top_variance_data,
                      title="Top "+str(max_sets)+" pathways by variance")

#' ## Features

#+ echo=False
import math

# read in the read count and feature count files
read_type, read_samples, read_count_data = document.read_table(vars["read_counts"])
feature_type, feature_samples, feature_count_data = document.read_table(vars["feature_counts"])

# get the total reads for samples along with those for nucleotide alignment and translated
# convert values to log10

total_reads=[math.log10(row[read_type.index("total reads")]) for row in read_count_data]
nucleotide_reads=[math.log10(row[read_type.index("total nucleotide aligned")]) for row in read_count_data]
translated_reads=[math.log10(row[read_type.index("total translated aligned")]) for row in read_count_data]

# sort the feature counts so they are in the same sample order as the read counts
feature_counts={sample:row for sample, row in zip(feature_samples, feature_count_data)}
sorted_feature_count_data=[feature_counts[sample] for sample in read_samples]

# get the counts by each feature type
# convert values to log10
genefamilies_counts=[math.log10(row[feature_type.index("humann2_ecs_relab_counts")]) for row in sorted_feature_count_data]
ecs_counts=[math.log10(row[feature_type.index("humann2_ecs_relab_counts")]) for row in sorted_feature_count_data]
pathabundance_counts=[math.log10(row[feature_type.index("humann2_pathabundance_relab_counts")]) for row in sorted_feature_count_data]

# add scatter plots of the data
document.plot_scatter([[total_reads,nucleotide_reads],[total_reads,translated_reads]],title="Read alignment rate",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Input reads)", ylabel="log10(Aligned reads)", trendline=True)

document.plot_scatter([[nucleotide_reads,genefamilies_counts],[translated_reads,genefamilies_counts]],title="UniRef90 gene families",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(gene families)", trendline=True)

document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="Enzymes (ECs)",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)

document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="Pathways",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

