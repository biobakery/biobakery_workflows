
#+ echo=False
import os

from biobakery_workflows import utilities, visualizations, files

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#+ echo=False

# display the heatmap
maaslin_taxonomy_heatmap = vars["maaslin_tasks_info"][0][1]
maaslin_taxonomy_output_folder = os.path.dirname(maaslin_taxonomy_heatmap)

#' # MaAsLin2 Results

#' MaAsLin2 is comprehensive R package for efficiently determining multivariable association between clinical metadata and microbial meta'omic features. MaAsLin2 relies on general linear models to accommodate most modern epidemiological study designs, including cross-sectional and longitudinal, and offers a variety of data exploration, normalization, and transformation methods. More detailed information may be found in the [MaAsLin2 User Manual](https://bitbucket.org/biobakery/maaslin2).

#' See the log file for the exact MaAsLiN2 commands run to generate these outputs. Also please note these are just a subset of the outputs, check out the MaAsLin2 results folders for the complete set of output files.

#' ## Taxonomic Profile

#' This report section contains the results from running the taxonomic profile through MaAsLin2. 

#' ### MaAsLin2 Heatmap
#' <% if os.path.isfile(maaslin_taxonomy_heatmap): print("![Taxonomy heatmap]("+maaslin_taxonomy_heatmap+")\n") %>
#' <% if not os.path.isfile(maaslin_taxonomy_heatmap): print("No significant associations.") %>

#' <% if pdf_format: print("\clearpage") %>

#' ### MaAsLin2 Plots

#' The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.

#+ echo=False

def show_maaslin_metadata_plots(figures_folder, type):
    # show the top plot for each metadata
    images_found = False
    for file_name in os.listdir(figures_folder):
        if file_name.endswith("_1.jpg"):
            images_found = True
            metadata_name=file_name.replace("_1.jpg","")
            print("![Most significant "+metadata_name+" association for "+type+"]("+os.path.join(figures_folder,file_name)+")\n")
    if not images_found:
        print("No significant associations.")

#' <% show_maaslin_metadata_plots(maaslin_taxonomy_output_folder,"taxonomy") %>

#+ echo=False

# check for the pathways results
if len(vars["maaslin_tasks_info"]) > 1:
    maaslin_pathways_heatmap = vars["maaslin_tasks_info"][1][1]
    maaslin_pathways_output_folder = os.path.dirname(vars["maaslin_tasks_info"][1][1])
else:
    maaslin_pathways_heatmap = ""
    maaslin_pathways_output_folder = ""

def display_pathways_heatmap(maaslin_pathways_heatmap):
        if os.path.isfile(maaslin_pathways_heatmap):
            print("![Pathways heatmap]("+maaslin_pathways_heatmap+")\n")
        else:
            print("No significant associations.")

#' <% if maaslin_pathways_output_folder and pdf_format: print("\clearpage") %>

#' <% if maaslin_pathways_output_folder: print("## Pathways\n") %>
#' <% if maaslin_pathways_output_folder: print("This report section contains the results from running the pathways through MaAsLin2.\n") %>

#' <% if maaslin_pathways_output_folder: print("### MaAsLin2 Heatmap\n") %>
#' <% if maaslin_pathways_output_folder: display_pathways_heatmap(maaslin_pathways_heatmap) %>

#' <% if maaslin_pathways_output_folder and pdf_format: print("\clearpage") %>

#' <% if maaslin_pathways_output_folder: print("### MaAsLin2 Plots") %>
#' <% if maaslin_pathways_output_folder: print("The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.") %>

#' <% show_maaslin_metadata_plots(maaslin_pathways_output_folder,"pathways") %>

#+ echo=False

# check for stratified pathways results
filtered_stratified_pathways_plots=[]
for plot_file in vars["stratified_pathways_plots"]:
    if os.path.isfile(plot_file):
        filtered_stratified_pathways_plots.append(plot_file)

def show_stratified_plots(plots):
    # Display each of the plots in the report
    for i, jpg_file in enumerate(plots):
        print("![Pathway #{0} sorted by significance from most to least]({1})\n".format(i+1, jpg_file))

#' <% if filtered_stratified_pathways_plots and pdf_format: print("\clearpage") %>

#' <% if filtered_stratified_pathways_plots: print("# Stratified Pathways Plots") %>

#' <% if maaslin_pathways_output_folder: print("The abundance for each of the "+str(len(filtered_stratified_pathways_plots))+" most significant associations are plotted stratified by species. These plots were generated with the utility script included with HUMAnN2 named humann2_barplot.") %>

#' <% show_stratified_plots(filtered_stratified_pathways_plots) %>

