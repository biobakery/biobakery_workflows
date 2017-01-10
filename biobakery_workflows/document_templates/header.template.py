
#+ echo=False
import time

from anadama2 import PweaveDocument

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()
        
#' % <% print(vars["title"]) %>
#' % Project: <% print(vars["project"]) %>
#' % Date: <%= time.strftime("%m/%d/%Y") %>
