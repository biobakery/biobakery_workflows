
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

# set the max number of features for the norm heatmaps
top_norm_pathways=50
top_norm_ecs=100
top_norm_genes=100

#' # Functional Profiling of Metagenomic and Metatranscriptomic Reads

#' <%= visualizations.ShotGun.captions["functional_intro"] %>

#+ echo=False

# read in the top average pathways
dna_samples, dna_top_average_pathways, dna_top_average_data, top_names_and_descriptions = visualizations.top_average_pathways(
    document, vars["dna_pathabundance"], max_sets)

#' ## Pathway Abundance

#' <%= visualizations.ShotGun.format_caption("heatmap_intro",max_sets=max_sets,type="pathways",method="Spearman") %>

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

#+ echo=False
# check for then rna/dna norm files
show_norm_ratio=False
if vars["genefamilies_norm_ratio"] and vars["ecs_norm_ratio"] and vars["paths_norm_ratio"]:
    show_norm_ratio=True 

norm_intro="The heatmaps show RNA normalized to DNA for the three features computed by HUMAnN2: gene families, ecs and pathways."
norm_intro+=" Note the most abundant DNA features are not necessarily those with the highest transcription (RNA) levels."

#' <% if show_norm_ratio and pdf_format: print("\clearpage") %>

#' <% if show_norm_ratio: print("### RNA/DNA Normalized Features") %>

#' <% if show_norm_ratio: print(norm_intro) %>

#+ echo=False
if show_norm_ratio:
    # read in the top average rna normed pathways
    samples, top_pathways, top_pathway_data, top_names_and_descriptions = visualizations.top_average_pathways(
        document, vars["paths_norm_ratio"], top_norm_pathways)
    
    document.show_hclust2(samples,top_pathways,top_pathway_data,
                          title="Top "+str(top_norm_pathways)+" RNA pathways by average abundance")  
 
#' <% if show_norm_ratio and pdf_format: print("\clearpage") %>
 
#+ echo=False
table_message=None
if show_norm_ratio:  
    # write a table of the pathways average and variance
    pathway_file_name="top_average_rna_dna_pathways_names.tsv"
    average_abundance_variance=visualizations.write_pathway_average_variance_table(document, pathway_file_name, top_pathway_data, top_names_and_descriptions)
    
    table_message=visualizations.show_table_max_rows(document, average_abundance_variance, 
        top_names_and_descriptions, [" Average "," Variance "], 
        "Top "+str(top_norm_pathways)+" RNA pathways by average abundance", pathway_file_name, font=7)
    
#' <% if table_message: print(table_message) %>

#' <% if table_message: visualizations.print_pathways_urls(top_pathways,top_names_and_descriptions,3) %>
    
#' <% if show_norm_ratio and pdf_format: print("\clearpage") %>
    
#+ echo=False
if show_norm_ratio:
    # read in the top average rna normed ecs
    samples, top_ecs, top_ec_data, top_names_and_descriptions = visualizations.top_average_pathways(
        document, vars["ecs_norm_ratio"], top_norm_ecs)
    
    document.show_hclust2(samples,top_ecs,top_ec_data,
                          title="Top "+str(top_norm_ecs)+" RNA ECs by average abundance")
#+ echo=False
if show_norm_ratio:
    # read in the top average rna normed genes
    samples, top_genes, top_gene_data, top_names_and_descriptions = visualizations.top_average_pathways(
        document, vars["genefamilies_norm_ratio"], top_norm_genes)
    
    document.show_hclust2(samples,top_genes,top_gene_data,
                          title="Top "+str(top_norm_genes)+" RNA gene families by average abundance")   
    
#+ echo=False
# check if the optional feature files were included
show_dna_features=False
if vars["dna_aligned_read_counts"] and vars["dna_feature_counts"]:
    show_dna_features=True
    
show_rna_features=False
if vars["rna_aligned_read_counts"] and vars["rna_feature_counts"]:
    show_rna_features=True

# determine the wording for the intro based on the input files provided
if show_dna_features and show_rna_features:
    seq_data_type="metagenomic and metatranscriptomic"
    short_type="DNA and RNA"
elif show_dna_features:
    seq_data_type="metagenomic"
    short_type="DNA"
elif show_rna_features:
    seq_data_type="metatranscriptomic"
    short_type="RNA"
else:
    seq_data_type=""
    short_type=""
    
#' <% if show_rna_features and show_norm_ratio and pdf_format: print("\clearpage") %>  

#' <% if show_dna_features or show_rna_features: print("## Features") %>

#' <% if show_dna_features or show_rna_features: print(visualizations.ShotGun.format_caption("feature_detection",seq_type=seq_data_type,seq_short_type=short_type)) %>

#' <% if show_dna_features: print("### DNA Features") %>

#+ echo=False
if show_dna_features:
    total_reads, nucleotide_reads, translated_reads, genefamilies_counts, ecs_counts, pathabundance_counts = visualizations.feature_counts(
        document, vars["dna_aligned_read_counts"],vars["dna_feature_counts"])
        
    # add scatter plots of the data
    document.plot_scatter([[total_reads,nucleotide_reads],[total_reads,translated_reads]],title="DNA Read alignment rate",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Input reads)", ylabel="log10(Aligned reads)", trendline=True)
    
#' <% if show_dna_features: print(visualizations.ShotGun.captions["scatter_reads_aligned"]) %>
#' <% if show_dna_features and pdf_format: print("\clearpage") %>
    
#+ echo=False
if show_dna_features:
    document.plot_scatter([[nucleotide_reads,genefamilies_counts],[translated_reads,genefamilies_counts]],title="DNA UniRef90 gene families",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(gene families)", trendline=True)

#+ echo=False
if show_dna_features:
    document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="DNA Enzymes (ECs)",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)

#+ echo=False
if show_dna_features:
    document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="DNA Pathways",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

#' <% if show_dna_features: print(visualizations.ShotGun.captions["scatter_features"]) %>

#' <% if show_dna_features and pdf_format: print("\clearpage") %>

#' <% if show_rna_features: print("### RNA Features") %>

#+ echo=False
if show_rna_features:
    total_reads, nucleotide_reads, translated_reads, genefamilies_counts, ecs_counts, pathabundance_counts = visualizations.feature_counts(
        document, vars["rna_aligned_read_counts"],vars["rna_feature_counts"])
    
    # add scatter plots of the data
    document.plot_scatter([[total_reads,nucleotide_reads],[total_reads,translated_reads]],title="RNA Read alignment rate",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Input reads)", ylabel="log10(Aligned reads)", trendline=True)
    
#' <% if show_rna_features: print(visualizations.ShotGun.captions["scatter_reads_aligned"]) %>
#' <% if show_rna_features and pdf_format: print("\clearpage") %>
    
#+ echo=False
if show_rna_features:
    document.plot_scatter([[nucleotide_reads,genefamilies_counts],[translated_reads,genefamilies_counts]],title="RNA UniRef90 gene families",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(gene families)", trendline=True)
#+ echo=False
if show_rna_features:
    document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="RNA Enzymes (ECs)",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)
#+ echo=False
if show_rna_features:
    document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="RNA Pathways",
                            row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

#' <% if show_rna_features: print(visualizations.ShotGun.captions["scatter_features"]) %>

#' <% if show_dna_features or show_rna_features: print("\clearpage") %>
