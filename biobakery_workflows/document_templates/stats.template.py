
#+ echo=False
import os

from biobakery_workflows import utilities, visualizations, files

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' # Taxonomic Profiling of Metagenomic Reads

#' This report section contains information about the taxonomy
#' for all DNA samples. These samples were
#' run through [MetaPhlAn2](http://huttenhower.sph.harvard.edu/metaphlan2).

#+ echo=False

# read in the taxonomy data
samples, taxonomy, data = document.read_table(vars["taxonomic_profile"])

# display the heatmap
maaslin_taxonomy_output_folder = os.path.join(vars["output_folder"],vars["maaslin_tasks_info"][0][1],"figures")
maaslin_taxonomy_heatmap = os.path.join(maaslin_taxonomy_output_folder,"heatmap.jpg")

#' ## MaAsLin2 Results

#' MaAsLin2 is comprehensive R package for efficiently determining multivariable association between clinical metadata and microbial meta'omic features. MaAsLin2 relies on general linear models to accommodate most modern epidemiological study designs, including cross-sectional and longitudinal, and offers a variety of data exploration, normalization, and transformation methods. More detailed information may be found in the [MaAsLin2 User Manual](https://bitbucket.org/biobakery/maaslin2).

#' See the log file for the exact MaAsLiN2 commands run to generate these outputs. Also please note these are just a subset of the outputs, check out the MaAsLin2 results folders for the complete set of output files.

#' ### MaAsLin2 Heatmap
#' <% if os.path.isfile(maaslin_taxonomy_heatmap): print("![Taxonomy heatmap]("+maaslin_taxonomy_heatmap+")\n") %>
#' <% if not os.path.isfile(maaslin_taxonomy_heatmap): print("No significant associations.") %>

#' <% if pdf_format: print("\clearpage") %>

#' ### MaAsLin2 Plots

#' The most significant association for each metadata are shown. For a complete set of plots, check out the MaAsLin2 results folders.

#+ echo=False

def show_maaslin_metadata_plots(figures_folder):
    # show the top plot for each metadata
    for file_name in os.listdir(figures_folder):
        if file_name.endswith("_1.jpg"):
            metadata_name=file_name.replace("_1.jpg","")
            print("![Most significant "+metadata_name+" association for taxonomy]("+os.path.join(figures_folder,file_name)+")\n")

#' <% show_maaslin_metadata_plots(maaslin_taxonomy_output_folder) %>

