
#+ echo=False

from anadama2.reporters import LoggerReporter 

#' # Data Processing Workflow Information

#' ## Software Versions

#' 
#' <% for version in LoggerReporter.read_log(vars["log"],"versions"): print("  * " + version+"  \n") %>
#'

#' ## Tasks Run

#'
#' <% for command in LoggerReporter.read_log(vars["log"],"commands"): print("  * " + command+"  \n") %>
#' 
