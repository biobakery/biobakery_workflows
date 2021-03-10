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

# import the workflow class from anadama2
from anadama2 import Workflow

# import the document templates from biobakery_workflows
from biobakery_workflows import document_templates, utilities, visualizations

# import the files for descriptions and paths
from biobakery_workflows import files
import os

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with multiple input files
workflow = Workflow(version="0.1", remove_options=["input"],
                    description="A workflow for 16S visualization")
                    
# add the custom arguments to the workflow 
# create a custom description for the input argument listing all expected input files
input_desc="A folder containing the final products from the 16s data workflow.\n\nThe input folder must include the following:\n\n"  
                      
workflow.add_argument("input",desc=input_desc,required=True)

# add the custom arguments to the workflow
workflow.add_argument("project-name",desc="the name of the project",required=True)
workflow.add_argument("author-name",desc="the name of the author of the report", required=True)
workflow.add_argument("header-image",desc="the image to add to the report header", default="")
workflow.add_argument("input-metadata",desc="the metadata file (samples as columns or rows)")
workflow.add_argument("input-picard",desc="the folder of picard quality score files")
workflow.add_argument("input-picard-extension",desc="the extensions for the picard quality score files", default="quality_by_cycle_metrics")
workflow.add_argument("metadata-categorical",desc="the categorical features", action="append", default=[])
workflow.add_argument("metadata-continuous",desc="the continuous features", action="append", default=[])
workflow.add_argument("metadata-exclude",desc="the features to exclude", action="append", default=[])
workflow.add_argument("exclude-workflow-info",desc="do not include data processing task info in report", action="store_true")
workflow.add_argument("format",desc="the format for the report", default="pdf", choices=["pdf","html"])
workflow.add_argument("introduction",desc="the introduction to be included in the report [DEFAULT: intro includes information from workflow log]", default="")
workflow.add_argument("print-template",desc="only print the template for the visualization workflow, do not run the workflow", action="store_true")

# get the arguments from the command line
args = workflow.parse_args()

otu_table = files.SixteenS.path("otu_table_closed_reference",args.input, error_if_not_found=True)

# read and label the metadata
metadata=None
metadata_labels=None
if args.input_metadata:
    metadata=utilities.read_metadata(args.input_metadata, otu_table, ignore_features=args.metadata_exclude, otu_table=True)
    metadata_labels, metadata=utilities.label_metadata(metadata, categorical=args.metadata_categorical, continuous=args.metadata_continuous)

# if using a header image then select a different starting template
if args.header_image:
    templates=[utilities.get_package_file("header_image")]
else:
    templates=[utilities.get_package_file("header_author")]

log_file=None
# add the template for the data processing information
if not args.exclude_workflow_info:
    log_file=files.Workflow.path("log", args.input, error_if_not_found=True)

# get the variables, input files, and method depending on the input files provided for the workflow
method_vars, method_depends, input_files, method = utilities.set_variables_for_16s_workflow_based_on_input(args,otu_table,files)

# add additional variables
method_vars["metadata"]=metadata
method_vars["metadata_labels"]=metadata_labels
method_vars["log"]=log_file

# listing all expected input files
input_desc+=files.SixteenS.list_file_path_description("",input_files)

# add the correct QC template based on the method
if method == "usearch":
    templates += [utilities.get_package_file("quality_control_usearch")]
else:
    templates += [utilities.get_package_file("quality_control_dada2")]

# if picard files are present then add to the template
if method_vars["picard"]:
    templates += [utilities.get_package_file("picard")]

# add the correct read count template
if method == "usearch":
    templates += [utilities.get_package_file("read_count_usearch")]
else:
    templates += [utilities.get_package_file("read_count_dada2")]

# add the rest of the 16s template
templates += [utilities.get_package_file("16S")]

if not args.exclude_workflow_info:
    templates += [utilities.get_package_file("workflow_info")]

# get the introduction text if not provided by the user
if not args.introduction:
    method_vars["introduction"]=visualizations.Sixteen_S.compile_default_intro(method_vars)
else:
    method_vars["introduction"]=args.introduction

# add author and image if included
method_vars["author"]=args.author_name
method_vars["header_image"]=args.header_image

if args.print_template:
    # only print the template to stdout
    utilities.print_template(templates)

# add the document to the workflow
doc_task=workflow.add_document(
    templates=templates,
    depends=method_depends, 
    targets=workflow.name_output_files("16S_report."+args.format),
    vars=method_vars,
    table_of_contents=True)

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
workflow.add_archive(
    depends=[args.output,doc_task],
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
