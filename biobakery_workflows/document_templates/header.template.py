
#+ echo=False
import time

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False
        
#' % <% print(vars["title"]) %>
#' <% if vars["project"]: print("% Project: "+vars["project"]) if pdf_format else print("**Project**: "+vars["project"]) %>  
#' <% print("% Date: "+ time.strftime("%m/%d/%Y")) if pdf_format else print("**Date**: " + time.strftime("%m/%d/%Y")) %>

#' # Introduction

#' <% print(vars["introduction_text"]) %>
