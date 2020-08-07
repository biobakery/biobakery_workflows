
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
maaslin_taxonomy_heatmap = vars["maaslin_tasks_info"]["taxonomy"][1]
maaslin_taxonomy_output_folder = os.path.dirname(maaslin_taxonomy_heatmap)

#' <% if not vars["bypass_maaslin"]: print("# MaAsLin2 Results") %>

#' <% if not vars["bypass_maaslin"]: print("MaAsLin2 is comprehensive R package for efficiently determining multivariable association between clinical metadata and microbial meta'omic features. MaAsLin2 relies on general linear models to accommodate most modern epidemiological study designs, including cross-sectional and longitudinal, and offers a variety of data exploration, normalization, and transformation methods. More detailed information may be found in the [MaAsLin2 User Manual](https://bitbucket.org/biobakery/maaslin2).") %>

#' <% if not vars["bypass_maaslin"]: print("See the log file for the exact MaAsLiN2 commands run to generate these outputs. Also please note these are just a subset of the outputs, check out the MaAsLin2 results folders for the complete set of output files.") %>

#' <% if not vars["bypass_maaslin"]: print("## Taxonomic Profile") %>

#' <% if not vars["bypass_maaslin"]: print("This report section contains the results from running the taxonomic profile through MaAsLin2.") %>

#' <% if not vars["bypass_maaslin"]: print("### MaAsLin2 Heatmap") %>
#' <% if os.path.isfile(maaslin_taxonomy_heatmap): print("![Taxonomy heatmap]("+maaslin_taxonomy_heatmap+")\n") %>
#' <% if not os.path.isfile(maaslin_taxonomy_heatmap) and not vars["bypass_maaslin"]: print("No significant associations.") %>

#' <% if pdf_format and not vars["bypass_maaslin"]: print("\clearpage") %>

#' <% if not vars["bypass_maaslin"]: print("### MaAsLin2 Plots") %>

#' <% if not vars["bypass_maaslin"]: print("The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.") %>

#+ echo=False

def show_maaslin_metadata_plots(figures_folder, type):
    # show the top plot for each metadata
    images_found = False
    for file_name in os.listdir(figures_folder):
        if file_name.endswith("_1.jpg"):
            images_found = True
            metadata_name=file_name.replace("_1.jpg","")
            print("![Most significant "+metadata_name+" association for "+type+"]("+os.path.join(figures_folder,file_name)+")\n\n")
            if pdf_format:
                print("\clearpage")
    for i in range(2,11):
        for file_name in os.listdir(figures_folder):
            if file_name.endswith("_{}.jpg".format(i)):
                metadata_name=file_name.replace("_{}.jpg".format(i),"")
                print("![Significant #"+str(i)+" "+metadata_name+" association for "+type+"]("+os.path.join(figures_folder,file_name)+")\n\n")
                if pdf_format:
                    print("\clearpage")
    if not images_found:
        print("No significant associations.")

#' <% if not vars["bypass_maaslin"]: show_maaslin_metadata_plots(maaslin_taxonomy_output_folder,"taxonomy") %>

#+ echo=False

# check for the pathways/ecs results
def check_for_maaslin_runs(run_type):
    if run_type in vars["maaslin_tasks_info"] and not vars["bypass_maaslin"]:
        heatmap_file = vars["maaslin_tasks_info"][run_type][1]
        output_folder = os.path.dirname(vars["maaslin_tasks_info"][run_type][1])
    else:
        heatmap_file = ""
        output_folder = ""

    return heatmap_file, output_folder

def display_maaslin_heatmap(maaslin_heatmap, run_type):
        if os.path.isfile(maaslin_heatmap):
            print("!["+run_type+" heatmap]("+maaslin_heatmap+")\n")
        else:
            print("No significant associations.")

maaslin_pathways_heatmap, maaslin_pathways_output_folder = check_for_maaslin_runs("pathways")

#' <% if maaslin_pathways_output_folder and pdf_format: print("\clearpage") %>

#' <% if maaslin_pathways_output_folder: print("## Pathways\n") %>
#' <% if maaslin_pathways_output_folder: print("This report section contains the results from running the pathways through MaAsLin2.\n") %>

#' <% if maaslin_pathways_output_folder: print("### MaAsLin2 Heatmap\n") %>
#' <% if maaslin_pathways_output_folder: display_maaslin_heatmap(maaslin_pathways_heatmap, "Pathways") %>

#' <% if maaslin_pathways_output_folder and pdf_format: print("\clearpage") %>

#' <% if maaslin_pathways_output_folder: print("### MaAsLin2 Plots") %>
#' <% if maaslin_pathways_output_folder: print("The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.") %>

#' <% if maaslin_pathways_output_folder: show_maaslin_metadata_plots(maaslin_pathways_output_folder,"pathways") %>

#+ echo=False

maaslin_ecs_heatmap, maaslin_ecs_output_folder = check_for_maaslin_runs("ecs")

#' <% if maaslin_ecs_output_folder and pdf_format: print("\clearpage") %>

#' <% if maaslin_ecs_output_folder: print("## ECs\n") %>
#' <% if maaslin_ecs_output_folder: print("This report section contains the results from running the ECs through MaAsLin2.\n") %>

#' <% if maaslin_ecs_output_folder: print("### MaAsLin2 Heatmap\n") %>
#' <% if maaslin_ecs_output_folder: display_maaslin_heatmap(maaslin_ecs_heatmap, "ECs") %>

#' <% if maaslin_ecs_output_folder and pdf_format: print("\clearpage") %>

#' <% if maaslin_ecs_output_folder: print("### MaAsLin2 Plots") %>
#' <% if maaslin_ecs_output_folder: print("The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.") %>

#' <% if maaslin_ecs_output_folder: show_maaslin_metadata_plots(maaslin_ecs_output_folder,"ecs") %>

#+ echo=False

# check for stratified pathways results
filtered_stratified_pathways_plots=[]
for plot_file in vars["stratified_pathways_plots"]:
    if os.path.isfile(plot_file):
        filtered_stratified_pathways_plots.append(plot_file)

def show_stratified_plots(plots):
    # Display each of the plots in the report
    for jpg_file in sorted(plots, key=lambda x: int(x.replace(".jpg","").split("_")[-1])):
        # get the pathway number and metadata name
        info = jpg_file.replace(".jpg","").split("_")
        pathway_number = info[-1]
        try:
            metadata_focus = open(jpg_file.replace(".jpg",".txt")).readline().rstrip()
        except EnvironmentError:
            metadata_focus = "Unknown"
        print("![Pathway #{0} sorted by significance from most to least for metadata focus {1}]({2})\n\n".format(int(pathway_number)+1, metadata_focus, jpg_file))

#' <% if filtered_stratified_pathways_plots and pdf_format: print("\clearpage") %>

#' <% if filtered_stratified_pathways_plots: print("# Stratified Pathways Plots") %>

#' <% if maaslin_pathways_output_folder: print("The abundance for each of the "+str(len(filtered_stratified_pathways_plots))+" most significant associations are plotted stratified by species. These plots were generated with the utility script included with HUMAnN named humann_barplot.") %>

#' <% show_stratified_plots(filtered_stratified_pathways_plots) %>

#' <% if pdf_format: print("\clearpage") %>

#' <% if vars["permanova_plots"]: print("# Permanova") %>
#' <% if vars["beta_diversity_plots"]["univariate"]: print("# Univariate") %>

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
        if os.path.isfile(permanova_file):
            print("![{0}]({1})\n\n".format(filetype, permanova_file))
        else:
            print("Error generating permanova for filetype {}".format(filetype))

#+ echo=False
#' <% show_all_permanova(vars["permanova_plots"]) %>

#+ echo=False
#' <% show_all_variate_plots("univariate") %>

#' <% if vars["beta_diversity_plots"]["multivariate"] and pdf_format: print("\clearpage") %>

#' <% if vars["beta_diversity_plots"]["multivariate"]: print("# Multivariate") %>

#' <% if vars["beta_diversity_plots"]["multivariate"]: print("For the multivariate model the following covariate equation was provided: 'bray ~ "+vars["covariate_equation"]+"' .") %>

#+ echo=False
#' <% show_all_variate_plots("multivariate") %>



