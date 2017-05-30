#' % <% from anadama2 import PweaveDocument; document=PweaveDocument(); vars = document.get_vars(); print(vars["title"]) %>
#' % Project: <% print(vars["project"]) %>
#' % Date: <% import time; print(time.strftime("%m/%d/%Y")) %>

#' # Introduction

#' <% print(vars["introduction_text"]) %>
