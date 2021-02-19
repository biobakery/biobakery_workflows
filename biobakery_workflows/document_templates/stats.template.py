
#+ echo=False
import os
import re

from biobakery_workflows import utilities, visualizations, files

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' <% if pdf_format: print("\clearpage") %>

#' <% if vars["mantel_plots"]: print("# All against all: Mantel test") %>

#' <% if vars["mantel_plots"] and os.path.isfile(vars["mantel_plots"][0]): print("![{0}]({1})\n\n\n".format("Mantel test",vars["mantel_plots"][0])) %>
#' <% if pdf_format: print("\clearpage") %>

#' <% if vars["permanova_plots"] or vars["beta_diversity_plots"]["univariate"] or vars["beta_diversity_plots"]["multivariate"]: print("# One against all: Permanova test") %>
#' <% if vars["beta_diversity_plots"]["univariate"]: print("## Univariate") %>

#+ echo=False

def show_all_variate_plots(runtype):
    for filetype in vars["beta_diversity_plots"][runtype]:
        show_univariate_plot(filetype,runtype)

def show_univariate_plot(filetype,runtype):
    if vars["beta_diversity_plots"][runtype] and filetype in vars["beta_diversity_plots"][runtype] and os.path.isfile(vars["beta_diversity_plots"][runtype][filetype]):
        print("![{0} {1}]({2})\n\n\n".format(filetype,runtype,vars["beta_diversity_plots"][runtype][filetype]))
    elif filetype in vars["beta_diversity_plots"][runtype]:
        print("Error generating {0} for filetype {1}".format(runtype,filetype))

def show_all_permanova(permanova_plots):
    for filetype in permanova_plots:
        permanova_file = permanova_plots[filetype]
        if filetype == "all":
            filetype = "Heatmap of all data features"
        if os.path.isfile(permanova_file):
            print("![{0}]({1})\n\n".format(filetype, permanova_file))
        else:
            print("Error generating permanova for filetype {}".format(filetype))

#+ echo=False
#' <% show_all_permanova(vars["permanova_plots"]) %>

#+ echo=False
#' <% show_all_variate_plots("univariate") %>

#' <% if vars["beta_diversity_plots"]["multivariate"] and pdf_format: print("\clearpage") %>

#' <% if vars["beta_diversity_plots"]["multivariate"]: print("## Multivariate") %>

#' <% if vars["beta_diversity_plots"]["multivariate"]: print("For the multivariate model the following covariate equation was provided: 'bray ~ "+vars["covariate_equation"]+"' .") %>

#+ echo=False
#' <% show_all_variate_plots("multivariate") %>
#' <% if vars["beta_diversity_plots"]["multivariate"] and pdf_format: print("\clearpage") %>

#' <% if not vars["bypass_maaslin"]: print("# Each metadata against all data: MaAsLin2") %>

#' <% if not vars["bypass_maaslin"]: print("MaAsLin2 is comprehensive R package for efficiently determining multivariable association between clinical metadata and microbial meta'omic features. MaAsLin2 relies on general linear models to accommodate most modern epidemiological study designs, including cross-sectional and longitudinal, and offers a variety of data exploration, normalization, and transformation methods. More detailed information may be found in the [MaAsLin2 User Manual](https://bitbucket.org/biobakery/maaslin2).") %>

#' <% if not vars["bypass_maaslin"]: print("See the log file for the exact MaAsLiN2 commands run to generate these outputs. Also please note these are just a subset of the outputs, check out the MaAsLin2 results folders for the complete set of output files.") %>

#+ echo=False

def show_heatmaps(heatmap, run_type):
    # display the heatmap if generated

    if os.path.isfile(heatmap):
        print("\n\n!["+run_type+" heatmap]("+heatmap+")\n\n")
    else:
        print("Not enough significant associations for a heatmap.\n\n")

def show_maaslin_tile(figures, type):
    # show the top plots for each metadata

    images_found = False

    for metadata_name in figures:
        if ( os.path.isfile(figures[metadata_name]) and os.path.getsize(figures[metadata_name]) > 0 ):
            images_found = True
            print("#### "+metadata_name+"\n\n")
            print("\n![Top "+metadata_name+" associations for "+type+"]("+figures[metadata_name]+")\n\n\n")
            print("\clearpage \n\n")

    if not images_found:
        print("No significant associations.\n\n")

def show_all_maaslin_run_types(maaslin_tasks_info):
    # search through the maaslin tasks info and show all plots for all data types

    for newtype in maaslin_tasks_info:
        maaslin_heatmap = maaslin_tasks_info[newtype][1]
        maaslin_output_folder = os.path.dirname(maaslin_heatmap)

        # Start tile with upper case
        newtype_title=newtype[0].upper()+newtype[1:]

        print("## {}\n\n".format(newtype_title))
        print("This report section contains the results from running the {} data through MaAsLin2.\n\n".format(newtype))

        show_heatmaps(maaslin_heatmap, newtype)
        print("\clearpage \n\n")

        print("### MaAsLin2 Plots\n\n")
        print("The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.\n\n\n")
        image_files, maaslin_tiles = utilities.get_maaslin_image_files(maaslin_tasks_info)
        
        show_maaslin_tile(maaslin_tiles[newtype], newtype)
        print("\clearpage \n\n")

#' <% if not vars["bypass_maaslin"]: show_all_maaslin_run_types(vars["maaslin_tasks_info"]) %>

#+ echo=False
def show_stratified_plots(plots):
    # Display each of the plots in the report
    no_plots_found = True
    for image_file in sorted(plots, key=lambda x: int(x.replace(".jpg","").split("_")[-1])):
        # get the pathway number and metadata name
        info = image_file.replace(".jpg","").split("_")
        pathway_number = info[-1]
        try:
            metadata_focus = open(image_file.replace(".jpg",".txt")).readline().rstrip()
        except EnvironmentError:
            metadata_focus = "Unknown"

        if os.path.isfile(image_file) and os.path.getsize(image_file) > 0:
            no_plots_found = False
            print("![Pathway #{0} sorted by significance from most to least for metadata focus {1}]({2})\n\n".format(int(pathway_number)+1, metadata_focus, image_file))

    if no_plots_found:
        print("No significant associations for pathways with categorical metadata found.")

#' <% if vars["stratified_pathways_plots"]: print("## MaAsLiN2 stratified pathways plots") %>

#' <% if vars["stratified_pathways_plots"] and "pathways" in vars["maaslin_tasks_info"]: print("The abundance for each of the "+str(len(vars["stratified_pathways_plots"]))+" most significant associations, for categorical features only, are plotted stratified by species. These plots were generated with the utility script included with HUMAnN named humann_barplot.") %>

#' <% if vars["stratified_pathways_plots"]: show_stratified_plots(vars["stratified_pathways_plots"]) %>

#' <% if pdf_format: print("\clearpage") %>

#+ echo=False
#' <% if not vars["bypass_halla"]: print("# Each data type against all other data types: HAllA") %>

#' <% if not vars["bypass_halla"]: print("HAllA (Hierarchical All-against-All Association Testing) discovers densely-associated blocks of features between two high-dimensional 'omics datasets. HAllA was run on each possible set of pairs from the data sets provided. The heatmaps for each run type are shown. For more information from each HAllA run, check out the HAllA results folders for a complete set of output files.") %>

#+ echo=False

def show_halla_results(halla_tasks_info):
    for run_type in halla_tasks_info:
        print("## HAllA "+run_type.replace(" "," vs. ")+"\n\n")
        show_heatmaps(halla_tasks_info[run_type], run_type)
        print("\clearpage \n\n")

#' <% show_halla_results(vars["halla_tasks_info"]) %>

