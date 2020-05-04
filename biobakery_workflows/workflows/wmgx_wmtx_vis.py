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

import os

# import the workflow class from anadama2
from anadama2 import Workflow

# import the document templates from biobakery_workflows
from biobakery_workflows import utilities

# import the files for descriptions and paths
from biobakery_workflows import files

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with multiple input files
workflow = Workflow(version="0.1", remove_options=["input"],
                    description="A workflow for whole metagenome and metatranscriptome shotgun sequence visualization")

# list the required and optional files for the workflow
# these are expected to be included in the input folder
wmgx_input_files={"required":["kneaddata_read_counts","taxonomic_profile","pathabundance_relab"],"optional":["humann_read_counts","feature_counts"]}
wmtx_input_files={"required":["kneaddata_read_counts"],"optional":["humann_read_counts","feature_counts"]}
norm_input_files={"optional":["genefamilies_norm_ratio","ecs_norm_ratio","paths_norm_ratio"]}

# create a custom description for the input argument listing all expected input files
input_desc="A folder containing the final products from the wmgx_mwtx data workflow.\n\nThe input folder should include the following:\n\n"
input_desc+="Whole Metagenome Shotgun\n---------------------------------\n"
input_desc+=files.ShotGun.list_file_path_description(files.ShotGun.wmgx_folder_name,wmgx_input_files)
input_desc+="\n\nWhole Metatranscriptome Shotgun\n---------------------------------\n"
input_desc+=files.ShotGun.list_file_path_description(files.ShotGun.wmtx_folder_name,wmtx_input_files)
input_desc+="\n\nRNA/DNA Norm\n---------------------------------\n"
input_desc+=files.ShotGun.list_file_path_description("",norm_input_files)

# add the custom arguments to the workflow                          
workflow.add_argument("input",desc=input_desc,required=True)
workflow.add_argument("project-name",desc="the name of the project", required=True)
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default="The data was run through the standard workflow for whole metagenome and metatranscriptome shotgun sequencing.")
workflow.add_argument("exclude-workflow-info",desc="do not include data processing task info in report", action="store_true")
workflow.add_argument("format",desc="the format for the report", default="pdf", choices=["pdf","html"])

# get the arguments from the command line
args = workflow.parse_args()

# set the input folders for the wmgx and wmtx files
wmgx_input_folder=os.path.join(args.input,files.ShotGun.wmgx_folder_name)
wmtx_input_folder=os.path.join(args.input,files.ShotGun.wmtx_folder_name)

# get the paths for the required files and check they are found
wmgx_qc_counts=files.ShotGun.path("kneaddata_read_counts",wmgx_input_folder, error_if_not_found=True)
wmtx_qc_counts=files.ShotGun.path("kneaddata_read_counts",wmtx_input_folder, error_if_not_found=True)
taxonomic_profile=files.ShotGun.path("taxonomic_profile",wmgx_input_folder, error_if_not_found=True)
pathabundance=files.ShotGun.path("pathabundance_relab",wmgx_input_folder, error_if_not_found=True)

# get the templates for the report
templates=[utilities.get_package_file("header"),
    utilities.get_package_file("quality_control_paired_dna_rna"),
    utilities.get_package_file("taxonomy"),
    utilities.get_package_file("functional_dna_rna")]

# add the template for the data processing information
log_file=None
if not args.exclude_workflow_info:
    templates+=[utilities.get_package_file("workflow_info")]
    log_file=files.Workflow.path("log", args.input, error_if_not_found=True)

# add the document to the workflow
doc_task=workflow.add_document(
    templates=templates,
    depends=[wmgx_qc_counts, wmtx_qc_counts,
             taxonomic_profile, pathabundance], 
    targets=workflow.name_output_files("wmgx_wmtx_report."+args.format),
    vars={"title":"Metagenome and Metatranscriptome Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "dna_read_counts":wmgx_qc_counts,
          "rna_read_counts":wmtx_qc_counts,
          "dna_aligned_read_counts":files.ShotGun.path("humann_read_counts",wmgx_input_folder, none_if_not_found=True),
          "rna_aligned_read_counts":files.ShotGun.path("humann_read_counts",wmtx_input_folder, none_if_not_found=True),
          "dna_feature_counts":files.ShotGun.path("feature_counts",wmgx_input_folder, none_if_not_found=True),
          "rna_feature_counts":files.ShotGun.path("feature_counts",wmtx_input_folder, none_if_not_found=True),
          "taxonomic_profile":taxonomic_profile,
          "dna_pathabundance":pathabundance,
          "genefamilies_norm_ratio":files.ShotGun.path("genefamilies_norm_ratio",args.input,none_if_not_found=True),
          "ecs_norm_ratio":files.ShotGun.path("ecs_norm_ratio",args.input,none_if_not_found=True),
          "paths_norm_ratio":files.ShotGun.path("paths_norm_ratio",args.input,none_if_not_found=True),
          "format":args.format,
          "log":log_file},
    table_of_contents=True)

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
workflow.add_archive(
    depends=[args.output,doc_task],
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
