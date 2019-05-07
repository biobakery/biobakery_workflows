#!/usr/bin/env python

"""
bioBakery Workflows: statistical analysis workflow

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
#from anadama2.tracked import TrackedDirectory

# import the document templates and utilities from biobakery_workflows
from biobakery_workflows import utilities, statsvis

# import the files for descriptions and paths
from biobakery_workflows import files
import os

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with multiple input files
workflow = Workflow(version="0.1", remove_options=["input"],
                    description="Statistical analysis workflow")

# list the required and optional files for the workflow
# these are expected to be included in the input folder


# create a custom description for the input argument listing all expected input files
input_desc="A folder containing the final products from the wmgx or 16s data workflow.\n\nThe input folder should include the following:\n\n"


# add the custom arguments to the workflow                          
workflow.add_argument("input",desc=input_desc,required=True)

# add the custom arguments to the workflow
workflow.add_argument("project-name",desc="the name of the project", required=True)
workflow.add_argument("input-metadata",desc="the metadata file (samples as columns or rows)")
workflow.add_argument("metadata-categorical",desc="the categorical features", action="append", default=[])
workflow.add_argument("metadata-continuous",desc="the continuous features", action="append", default=[])
workflow.add_argument("metadata-exclude",desc="the features to exclude", action="append", default=[])
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default="The data was run through the standard workflow for whole metagenome shotgun sequencing.")
workflow.add_argument("exclude-workflow-info",desc="do not include data processing task info in report", action="store_true")
workflow.add_argument("format",desc="the format for the report", default="pdf", choices=["pdf","html"])

# get the arguments from the command line
args = workflow.parse_args()


if files.ShotGun.path("taxonomic_profile",args.input, none_if_not_found=True):
    workflow_data = "wmgx"
    # list the required and optional files for the workflow
    # these are expected to be included in the input folder
    input_files={"required":["taxonomic_profile"]}
    input_desc+=files.ShotGun.list_file_path_description("",input_files)
    taxonomic_profile = files.ShotGun.path("taxonomic_profile",args.input, none_if_not_found=True)
    pathabundance_relab = files.ShotGun.path("pathabundance_relab",args.input, none_if_not_found=True)
elif files.SixteenS.path("otu_table_closed_reference",args.input, error_if_not_found=True):
    workflow_data = "16s"
    # list the required and optional files for the workflow
    # these are expected to be included in the input folder
    input_files={"required":["otu_table_closed_reference"]}
    input_desc+=files.SixteenS.list_file_path_description("",input_files)
    taxonomic_profile = files.SixteenS.path("otu_table_closed_reference",args.input, error_if_not_found=True)
else:
    exit("No taxonomic profile found in input folder", args.input)

doc_depends=[taxonomic_profile]
adonis_taxa_dir=os.path.join(args.output,"permanova_taxa")
permanova_taxa_task=statsvis.run_permanova(workflow,taxonomic_profile,args.input,args.output,adonis_taxa_dir,workflow_data)
doc_depends.append(permanova_taxa_task)
split_task_permanova_taxa=statsvis.split_pdfs(workflow,adonis_taxa_dir,permanova_taxa_task)
doc_depends.append(split_task_permanova_taxa)

maaslin_taxa_dir=os.path.join(args.output,"maaslin_taxa")
maaslin_taxa_task=statsvis.run_maaslin(workflow,taxonomic_profile,args.input,maaslin_taxa_dir)
doc_depends.append(maaslin_taxa_task)
split_task_maaslin_taxa=statsvis.split_pdfs(workflow,maaslin_taxa_dir,maaslin_taxa_task)
doc_depends.append(split_task_maaslin_taxa)

if os.path.isfile(pathabundance_relab):
    adonis_path_dir=os.path.join(args.output,"permanova_path")
    permanova_path_task=statsvis.run_permanova(workflow,pathabundance_relab,args.input,args.output,adonis_path_dir,workflow_data)
    doc_depends.append(permanova_path_task)
    split_task_permanova_path=statsvis.split_pdfs(workflow,adonis_path_dir,permanova_path_task)
    doc_depends.append(split_task_permanova_path)

    maaslin_path_dir=os.path.join(args.output,"maaslin_path")
    maaslin_path_task=statsvis.run_maaslin(workflow,pathabundance_relab,args.input,maaslin_path_dir)
    doc_depends.append(maaslin_path_task)
    split_task_maaslin_path=statsvis.split_pdfs(workflow,maaslin_path_dir,maaslin_path_task)
    doc_depends.append(split_task_maaslin_path)


# add the template
templates=[utilities.get_package_file("stats")]


# add the template for the data processing information
log_file=None
if not args.exclude_workflow_info:
    templates+=[utilities.get_package_file("workflow_info")]
    log_file=files.Workflow.path("log", args.input, error_if_not_found=True)

# add the document to the workflow
doc_task=workflow.add_document(
    templates=templates,
    depends=doc_depends,
    targets=workflow.name_output_files("stats_report."+args.format),
    vars={"title":"Statistical Analysis Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "taxonomic_profile":taxonomic_profile,
          "format":args.format,
          "workflow_data": workflow_data,
          "output_dir": args.output,
          "log":log_file
          },
    table_of_contents=True)

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
workflow.add_archive(
    depends=[args.output,doc_task],
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
