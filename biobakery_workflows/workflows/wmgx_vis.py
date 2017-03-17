#!/usr/bin/env python

"""
bioBakery Workflows: whole metagenome shotgun workflow for visualization

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
                    description="A workflow for whole metagenome shotgun sequence visualization")

# add the custom arguments to the workflow
workflow.add_argument("qc-counts",desc="the kneaddata read counts table",required=True)
workflow.add_argument("taxonomic-profile",desc="the merged taxonomic profiles",required=True)
workflow.add_argument("pathabundance",desc="the pathway abundances",required=True)
workflow.add_argument("read-counts", desc="the read counts from the humann2 log files",required=True)
workflow.add_argument("feature-counts", desc="the counts of features (gene families, ecs, and pathways)", required=True)
workflow.add_argument("project-name",desc="the name of the project")
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default="The data was run through the standard workflow for whole metagenome shotgun sequencing.")

# get the arguments from the command line
args = workflow.parse_args()

# add the document to the workflow
workflow.add_document(
    templates=[document_templates.get_template("header"),
               document_templates.get_template("quality_control_paired_dna"),
               document_templates.get_template("taxonomy"),
               document_templates.get_template("functional_dna")],
    depends=[args.qc_counts, args.taxonomic_profile, args.pathabundance,
             args.read_counts, args.feature_counts], 
    targets=workflow.name_output_files("wmgx_report.pdf"),
    vars={"title":"Metagenome Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "dna_read_counts":os.path.abspath(args.qc_counts),
          "taxonomic_profile":os.path.abspath(args.taxonomic_profile),
          "dna_pathabundance":os.path.abspath(args.pathabundance),
          "read_counts":os.path.abspath(args.read_counts),
          "feature_counts":os.path.abspath(args.feature_counts)})

# start the workflow
workflow.go()
