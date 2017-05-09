#!/usr/bin/env python

"""
bioBakery Workflows: whole metagenome and metatranscriptome shotgun workflow for visualization

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

# import the document templates from biobakery_workflows
from biobakery_workflows import document_templates

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with multiple input files
workflow = Workflow(version="0.1", remove_options=["input"],
                    description="A workflow for whole metagenome and metatranscriptome shotgun sequence visualization")

# add the custom arguments to the workflow
workflow.add_argument("wmgx-qc-counts",desc="the kneaddata read counts table for wmgx samples",required=True)
workflow.add_argument("wmtx-qc-counts",desc="the kneaddata read counts table for wmtx samples",required=True)
workflow.add_argument("taxonomic-profile",desc="the merged taxonomic profiles for the wmgx samples",required=True)
workflow.add_argument("wmgx-pathabundance",desc="the pathway abundances for the wmgx samples",required=True)
workflow.add_argument("wmgx-read-counts", desc="the read counts from the humann2 log files for wmgx samples")
workflow.add_argument("wmtx-read-counts", desc="the read counts from the humann2 log files for wmtx samples")
workflow.add_argument("wmgx-feature-counts", desc="the counts of features (gene families, ecs, and pathways) for wmgx samples")
workflow.add_argument("wmtx-feature-counts", desc="the counts of features (gene families, ecs, and pathways) for wmtx samples")
workflow.add_argument("project-name",desc="the name of the project")
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default="The data was run through the standard workflow for whole metagenome and metatranscriptome shotgun sequencing.")
workflow.add_argument("format",desc="the format for the report, pdf or html", default="pdf")

# get the arguments from the command line
args = workflow.parse_args()

# add the document to the workflow
doc_task=workflow.add_document(
    templates=[document_templates.get_template("header"),
               document_templates.get_template("quality_control_paired_dna_rna"),
               document_templates.get_template("taxonomy"),
               document_templates.get_template("functional_dna_rna")],
    depends=[args.wmgx_qc_counts, args.wmtx_qc_counts,
             args.taxonomic_profile, args.wmgx_pathabundance], 
    targets=workflow.name_output_files("wmgx_wmtx_report."+args.format),
    vars={"title":"Metagenome and Metatranscriptome Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "dna_read_counts":args.wmgx_qc_counts,
          "rna_read_counts":args.wmtx_qc_counts,
          "dna_aligned_read_counts":args.wmgx_read_counts,
          "rna_aligned_read_counts":args.wmtx_read_counts,
          "dna_feature_counts":args.wmgx_feature_counts,
          "rna_feature_counts":args.wmtx_feature_counts,
          "taxonomic_profile":args.taxonomic_profile,
          "dna_pathabundance":args.wmgx_pathabundance,
          "format":args.format})

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
workflow.add_archive(
    depends=[args.output,doc_task],
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
