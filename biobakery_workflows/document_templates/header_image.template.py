#' % <% from anadama2 import PweaveDocument; document=PweaveDocument(); vars = document.get_vars(); print(" ")   %>
#'   <% print("  ![{0}]({1})  ".format("",vars["header_image"])) %>
#'   <% print("  "+vars["title"]+" for "+vars["project"]) %>
#' % Author: <% print(vars["author"]) %>
#' % Date: <% import time; print(time.strftime("%m/%d/%Y")+"\n") %>

#+ echo=False
import os
import re
import numpy

from biobakery_workflows import visualizations, utilities, files

# determine the document format
pdf_format=True if vars["format"] == "pdf" else False

# set universal constants
min_abundance=0.01
min_samples=10

max_taxa=15
max_sets_heatmap=25
max_sets_barplot=15

#' # Introduction

#' <% print(vars["introduction_text"]) %>

#' <% if pdf_format: print("\clearpage") %>
