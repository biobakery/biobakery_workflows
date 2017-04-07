#!/usr/bin/env python

"""
bioBakery Workflows: 16S workflow for visualization

Copyright (c) 2017 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import os

# import the workflow class from anadama2
from anadama2 import Workflow

# import the document templates from biobakery_workflows
from biobakery_workflows import document_templates

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with multiple input files
workflow = Workflow(version="0.1", remove_options=["input"],
                    description="A workflow for 16S visualization")

# add the custom arguments to the workflow
workflow.add_argument("otu-table",desc="the closed reference taxonomic profile",required=True)
workflow.add_argument("read-count-table",desc="the table of read counts (all, with taxonomy, and unclassified)",required=True)
workflow.add_argument("project-name",desc="the name of the project")
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default=("The samples from this project were run through the standard workflow for 16S sequencing. The workflow "+
             " follows the UPARSE OTU analysis pipeline for OTU calling and taxonomy prediction."+
             " The GreenGenes 16S RNA Gene Database version 13_8 was used for taxonomy prediction."+
             " Reads were filtered for quality control using a MAXEE score of 1. Filtered reads were"+
             " used to generate the OTUs. Reads not passing quality control were kept and used in the step"+
             " assigning reads to OTUs. First these reads were truncated to a max length of 200 bases."))

# get the arguments from the command line
args = workflow.parse_args()

# add the document to the workflow
doc_task=workflow.add_document(
    templates=[document_templates.get_template("header"),
               document_templates.get_template("16S")],
    depends=[args.otu_table, args.read_count_table], 
    targets=workflow.name_output_files("16S_report.pdf"),
    vars={"title":"16S Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "otu_table":os.path.abspath(args.otu_table),
          "read_count_table":os.path.abspath(args.read_count_table)})

# name the archive the same as the output folder
# join with "" so directory always ends with path join so dirname picks up correct directory
archive = os.path.dirname(os.path.join(args.output,""))+".zip"

# add an archive of the document and figures, removing the log file
workflow.add_archive(
    depends=[args.output,doc_task],
    targets=archive,
    archive_software="zip",
    remove_log=True)

# start the workflow
workflow.go()
