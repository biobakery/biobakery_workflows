
#+ echo=False

from anadama2 import PweaveDocument
from anadama2.reporters import LoggerReporter 

document=PweaveDocument()  

# get the variables for this document generation task
vars = document.get_vars()

# determine the document format
pdf_format = True if vars["format"] == "pdf" else False

#' # Data Processing Workflow Information

#' ## Software Versions

#' 
#' <% for version in LoggerReporter.read_log(vars["log"],"versions"): print("  * " + version+"  ") %>
#'

#' ## Tasks Run

#'
#' <% for command in LoggerReporter.read_log(vars["log"],"commands"): print("  * " + command+"  ") %>
#' 
