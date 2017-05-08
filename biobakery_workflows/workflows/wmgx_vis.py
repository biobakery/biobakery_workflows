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

# import the workflow class from anadama2
from anadama2 import Workflow

# import the document templates and utilities from biobakery_workflows
from biobakery_workflows import document_templates
from biobakery_workflows import utilities

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
workflow.add_argument("format",desc="the format for the report, pdf or html", default="pdf")
workflow.add_argument("contaminate-database",desc="the database used for contaminate read filtering for quality control", default="hg38")

# get the arguments from the command line
args = workflow.parse_args()

# select the templates based on the qc data
templates=[document_templates.get_template("header")]

if utilities.is_paired_table(args.qc_counts):
    templates+=[document_templates.get_template("quality_control_paired_dna")]
else:
    templates+=[document_templates.get_template("quality_control_single_dna")]
    
templates+=[document_templates.get_template("taxonomy"),
            document_templates.get_template("functional_dna")]

# add the document to the workflow
doc_task=workflow.add_document(
    templates=templates,
    depends=[args.qc_counts, args.taxonomic_profile, args.pathabundance,
             args.read_counts, args.feature_counts], 
    targets=workflow.name_output_files("wmgx_report."+args.format),
    vars={"title":"Metagenome Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "dna_read_counts":args.qc_counts,
          "taxonomic_profile":args.taxonomic_profile,
          "dna_pathabundance":args.pathabundance,
          "read_counts":args.read_counts,
          "feature_counts":args.feature_counts,
          "format":args.format,
          "contaminate_database":args.contaminate_database})

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
workflow.add_archive(
    depends=[args.output,doc_task],
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
