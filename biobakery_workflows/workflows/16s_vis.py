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
workflow.add_argument("read-count-table",desc="the table of read counts",required=True)
workflow.add_argument("project-name",desc="the name of the project")
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default="The data was run through the standard workflow for 16S sequencing.")

# get the arguments from the command line
args = workflow.parse_args()

# add the document to the workflow
workflow.add_document(
    templates=[document_templates.get_template("header"),
               document_templates.get_template("16S")],
    depends=[args.otu_table, args.read_count_table], 
    targets=workflow.name_output_files("16S_report.pdf"),
    vars={"title":"16S Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "otu_table":os.path.abspath(args.otu_table),
          "read_count_table":os.path.abspath(args.read_count_table)})

# start the workflow
workflow.go()
