#' % <% from anadama2 import PweaveDocument; document=PweaveDocument(); vars = document.get_vars(); print(vars["title"]) %>
#' % Project: <% print(vars["project"]) %>
#' % Date: <% import time; print(time.strftime("%m/%d/%Y")) %>
#+ echo=False

# get the variable settings from the data processing workflow
from anadama2.reporters import LoggerReporter
from biobakery_workflows import statsvis

try:
    workflow_settings = LoggerReporter.read_log(vars["log"],"variables")
except AttributeError:
    workflow_settings = []

# print a warning if the variables could not be read
if isinstance(workflow_settings, list):
    print("WARNING: Unable to read workflow settings from log file.")
    workflow_settings={}

input_dir= workflow_settings.get("input","UNK")

workflow_data=vars["workflow_data"]

# determine the document format
pdf_format=True if vars["format"] == "pdf" else False


def get_all_images(path):
    import os, fnmatch

    for d in fnmatch.filter(os.listdir(path), '*_images'):
        cur_path=os.path.join(path,d)
        for f in sorted(os.listdir(cur_path)):
            fname=f.split("-")
            ftitle=fname[0].replace("_"," ").capitalize()
            print("## "+ftitle+"\n\n")
            print("![Data Visualization "+ftitle+"]("+os.path.join(cur_path,f)+")\n")
            if pdf_format:
                print("\clearpage")

#' <% if pdf_format: print("\clearpage") %>

#' # Introduction
#+ echo=False
#' <% print(statsvis.Stats.captions["intro"]) %>

#+ echo=False

#' # Permanova

#' <% print(statsvis.Permanova.captions["intro"]) %>
#' <% if pdf_format: print("\clearpage") %>

#' <% get_all_images(vars["output_dir"]+"/permanova") %>

#' # Maaslin

#' <% print(statsvis.Maaslin.captions["intro"]) %>
#' <% if pdf_format: print("\clearpage") %>

#' <% get_all_images(vars["output_dir"]+"/maaslin") %>


#+ echo=False

#' <% if pdf_format: print("\clearpage") %>