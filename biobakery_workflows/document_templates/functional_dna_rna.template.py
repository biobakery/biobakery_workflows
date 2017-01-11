
#+ echo=False
max_sets=50

from biobakery_workflows import utilities

#' # Functional Profiling

#' This report section contains information about the functional profiling run
#' on all samples. These samples were
#' run through [HUMAnN2](http://huttenhower.sph.harvard.edu/humann2).
#' The taxonomic profiles for the DNA samples were used for both the DNA and RNA
#' samples in this set. The UniRef90 full database was used for the translated search.

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






