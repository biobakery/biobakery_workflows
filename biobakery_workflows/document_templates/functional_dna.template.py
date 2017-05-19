
#+ echo=False
max_sets=50

from biobakery_workflows import utilities
from biobakery_workflows import visualizations

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' # Functional Profiling

#' <%= visualizations.ShotGun.captions["functional_intro"] %>

#+ echo=False

# read in the top average pathways
dna_samples, dna_top_average_pathways, dna_top_average_data, top_names_and_descriptions = visualizations.top_average_pathways(
    document, vars["dna_pathabundance"], max_sets)

#' ## Pathway Abundance

#' <%= visualizations.ShotGun.format_caption("pathway_abundance_intro",max_sets=max_sets) %>

#+ echo=False
document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
                      title="Top "+str(max_sets)+" pathways by average abundance")  

#' <%= visualizations.ShotGun.format_caption("pathway_abundance_heatmap",norm="log10") %> 

#+ echo=False
document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
                      title="Top "+str(max_sets)+" pathways by average abundance",
                      log_scale=False,zscore=True)

#' <%= visualizations.ShotGun.format_caption("pathway_abundance_heatmap",norm="z-score") %> 

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False

# write a table of the pathways average and variance
pathway_file_name="top_average_pathways_names.tsv"
average_abundance_variance=visualizations.write_pathway_average_variance_table(document, pathway_file_name, dna_top_average_data, top_names_and_descriptions)

table_message=visualizations.show_table_max_rows(document, average_abundance_variance, 
    top_names_and_descriptions, [" Average "," Variance "], 
    "Top "+str(max_sets)+" pathways by average abundance", pathway_file_name, font=7)

#' <%= table_message %>

#' <% visualizations.print_pathways_urls(dna_top_average_pathways,top_names_and_descriptions,3) %>

#' <% if pdf_format: print("\clearpage") %>

#' ## Features

#' <%= visualizations.ShotGun.format_caption("feature_detection",seq_type="metagenomic",seq_short_type="DNA") %>

#+ echo=False

# read in the read count and feature count files
total_reads, nucleotide_reads, translated_reads, genefamilies_counts, ecs_counts, pathabundance_counts = visualizations.feature_counts(
    document, vars["read_counts"],vars["feature_counts"])

# add scatter plots of the data
document.plot_scatter([[total_reads,nucleotide_reads],[total_reads,translated_reads]],title="Read alignment rate",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Input reads)", ylabel="log10(Aligned reads)", trendline=True)

#' <%= visualizations.ShotGun.captions["scatter_reads_aligned"] %> 

#+ echo=False
document.plot_scatter([[nucleotide_reads,genefamilies_counts],[translated_reads,genefamilies_counts]],title="UniRef90 gene families",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(gene families)", trendline=True)

document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="Enzymes (ECs)",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)

document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="Pathways",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

#' <%= visualizations.ShotGun.captions["scatter_features"] %>
