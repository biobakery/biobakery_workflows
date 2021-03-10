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


#' # Introduction

#' <% print(vars["introduction_text"]) %>

#' <% if pdf_format: print("\clearpage") %>
