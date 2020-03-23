
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

#' # Functional Profiling of Metagenomic Reads

#' <%= visualizations.ShotGun.captions["functional_intro"] %>

#+ echo=False

# read in the top average pathways
dna_samples, dna_top_average_pathways, dna_top_average_data, top_names_and_descriptions = visualizations.top_average_pathways(
    document, vars["dna_pathabundance"], max_sets)
dna_ecs_samples, dna_top_average_ecs, dna_top_average_ecs_data, top_ecs_names_and_descriptions = visualizations.top_average_pathways(
    document, vars["dna_ecabundance"], max_sets)

#' ## Pathway and ECs Abundance

#' <%= visualizations.ShotGun.format_caption("heatmap_intro",max_sets=max_sets,type="pathways",method="Spearman") %>

#+ echo=False
# update the figure size based on output format for the heatmaps
utilities.change_pweave_figure_size_heatmap(pdf_format)

#+ echo=False
# if there is metadata, add it to the heatmap
def log10_heatmap(dna_samples, dna_top_average_pathways, dna_top_average_data, data_type="pathways"):
    merged_data=[]
    metadata_pathways=[]
    metadata_samples=[]
    if 'metadata' in vars and vars['metadata']:
        merged_data, metadata_samples=utilities.merge_metadata(vars['metadata'], dna_samples, 
            [[dna_top_average_pathways[i]]+dna_top_average_data[i] for i in range(len(dna_top_average_pathways))])
        metadata_pathways=[row.pop(0) for row in merged_data]
        # get the metadata row numbers
        metadata_rows=range(1,len(vars['metadata']))
        document.show_hclust2(metadata_samples, metadata_pathways, merged_data,
            title="Top "+str(max_sets)+" "+data_type+" by average abundance",
            metadata_rows=metadata_rows)
    else:
        document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
            title="Top "+str(max_sets)+" "+data_type+" by average abundance")  

    return merged_data, metadata_pathways, metadata_samples

merged_data, metadata_pathways, metadata_samples=log10_heatmap(dna_samples, dna_top_average_pathways, dna_top_average_data)

#' <%= visualizations.ShotGun.format_caption("pathway_abundance_heatmap",norm="log10") %> 

#+ echo=False
ec_merged_data, ec_pathways, ec_samples=log10_heatmap(dna_ecs_samples, dna_top_average_ecs, dna_top_average_ecs_data, "ecs")

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
def zscore_heatmap(dna_samples, dna_top_average_pathways, dna_top_average_data, merged_data, metadata_pathways, metadata_samples, data_type="pathways"):
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
            title="Top "+str(max_sets)+" "+data_type+" by average abundance",
            log_scale=False,zscore=True,
            metadata_rows=filtered_metadata_rows)
    else:
        document.show_hclust2(dna_samples,dna_top_average_pathways,dna_top_average_data,
            title="Top "+str(max_sets)+" "+data_type+" by average abundance",
            log_scale=False,zscore=True)

zscore_heatmap(dna_samples, dna_top_average_pathways, dna_top_average_data, merged_data, metadata_pathways, metadata_samples)

#' <%= visualizations.ShotGun.format_caption("pathway_abundance_heatmap",norm="z-score") %> 

#+ echo=False
zscore_heatmap(dna_ecs_samples, dna_top_average_ecs, dna_top_average_ecs_data, ec_merged_data, ec_pathways, ec_samples, "ecs")

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
# reset the figure size to the defaults
utilities.reset_pweave_figure_size()

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

#+ echo=False
document.plot_scatter([[nucleotide_reads,ecs_counts],[translated_reads,ecs_counts]],title="Enzymes (ECs)",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(ECs)", trendline=True)

#+ echo=False
document.plot_scatter([[nucleotide_reads,pathabundance_counts],[translated_reads,pathabundance_counts]],title="Pathways",
                       row_labels=["Nucleotide search","Nucleotide + translated search"],xlabel="log10(Aligned reads)", ylabel="log10(Pathways)", trendline=True)

#' <%= visualizations.ShotGun.captions["scatter_features"] %>
