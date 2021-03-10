
#' # Quality Control

#' ## Forward Reads
#' <% print("![FWD Read](" + vars["readF_qc"] + ")") %>
#' <% if pdf_format: print("\clearpage") %>

#' ## Reverse Reads
#' <% print("![REV Read](" + vars["readR_qc"] + ")") %>
#' <% if pdf_format: print("\clearpage") %>

#' # Error rates
#' <% print(visualizations.Sixteen_S.captions["dada2errorintro"]) %>

#' ## Forward Reads
#' <% print("![FWD Error Rates](" + vars["error_ratesF"] + ")") %>
#' <% if pdf_format: print("\clearpage") %>

#' ## Reverse Reads
#' <% print("![REV Error Rates](" + vars["error_ratesR"] + ")") %>
#' <% if pdf_format: print("\clearpage") %>

#' <% print(visualizations.Sixteen_S.captions["dada2errorinfo"]) %>

