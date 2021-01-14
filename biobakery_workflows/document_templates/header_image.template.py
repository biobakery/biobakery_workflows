#' % <% from anadama2 import PweaveDocument; document=PweaveDocument(); vars = document.get_vars(); print(" ")   %>
#'   <% print("  ![{0}]({1})  ".format("",vars["header_image"])) %>
#'   <% print("  "+vars["title"]+" for "+vars["project"]) %>
#' % Author: <% print(vars["author"]) %>
#' % Date: <% import time; print(time.strftime("%m/%d/%Y")+"\n") %>


#' # Introduction

#' <% print(vars["introduction_text"]) %>
