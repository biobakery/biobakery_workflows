#!/usr/bin/env python

"""
bioBakery Workflows: workflow for visualization

Copyright (c) 2021 Harvard School of Public Health

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
from biobakery_workflows import files, data
import os
import sys

# create a workflow instance, providing the version number and description
# remove the input folder option as it will be replaced with multiple input files
workflow = Workflow(version="3.1", remove_options=["input"],
                    description="A workflow for visualization")
                    
# add the custom arguments to the workflow 
input_desc=utilities.get_vis_input_description(files)
workflow.add_argument("input",desc=input_desc,required=True)

# add the custom arguments to the workflow
workflow.add_argument("project-name",desc="the name of the project", default="")
workflow.add_argument("author-name",desc="the name of the author of the report", default="")
workflow.add_argument("header-image",desc="the image to add to the report header", default="")
workflow.add_argument("input-metadata",desc="the metadata file (samples as columns or rows)",default="")
workflow.add_argument("input-file-type",desc="the file type for an input file formatted as 'filename,filetype'", action="append", default=[])
workflow.add_argument("input-picard",desc="the folder of picard quality score files")
workflow.add_argument("input-picard-extension",desc="the extensions for the picard quality score files", default="quality_by_cycle_metrics")
workflow.add_argument("metadata-categorical",desc="the categorical features", action="append", default=[])
workflow.add_argument("metadata-continuous",desc="the continuous features", action="append", default=[])
workflow.add_argument("metadata-exclude",desc="the features to exclude", action="append", default=[])
workflow.add_argument("min-abundance",desc="the min abundance to use for filtering", default=0.01)
workflow.add_argument("min-samples",desc="the min samples to use for filtering (as a percentage)", default=3)
workflow.add_argument("max-sets-heatmap",desc="the max sets to show for a heatmap", default=25)
workflow.add_argument("max-missing",desc="the max percentage of missing values for a metadata to not be filtered", default="20.0")
workflow.add_argument("max-sets-barplot",desc="the max sets to show for a barplot", default=15)
workflow.add_argument("max-groups-barplot",desc="the max number of grouped barplots to show for a single metadata variable", default=5)
workflow.add_argument("correlation-threshold",desc="the min value to use with the spearman correlation for filtering features on heatmap", default=0.7)
workflow.add_argument("format",desc="the format for the report", default="pdf", choices=["pdf","html"])
workflow.add_argument("introduction-text",desc="the introduction to be included in the report [DEFAULT: intro includes information from workflow log]", default="")
workflow.add_argument("print-template",desc="only print the template for the visualization workflow, do not run the workflow", action="store_true")
workflow.add_argument("use-template",desc="provide a report template to use instead of using that which is automatically generated", default="")

# get the arguments from the command line
args = workflow.parse_args()

templates=[utilities.get_package_file("universal_vis")]

# check for the log file
log_file=files.Workflow.path("log", args.input, none_if_not_found=True, search_for_file=True)

# get the paths for the required files from the set of all input files
data_files=utilities.identify_data_files(files,args.input,args.input_file_type,args.input_metadata)

# error if no input files are found
if len(data_files.keys()) < 1:
    sys.exit("ERROR: No data files found in the input folder.")

# check the required taxonomic profile is provided
if data_files.get("wmgx_taxonomy",""):
    workflow_type = "WGX"
elif "16s_taxonomy_asv" in data_files or "16s_taxonomy_otu" in data_files: 
    workflow_type = "16S"
else:
    print(input_desc)
    sys.exit(1)

if workflow_type == "16S" :
    # get the variables, input files, and method depending on the input files provided for the workflow
    template_variables, template_depends, method, otu_table = utilities.set_variables_for_16s_workflow_based_on_input(args,data_files)

    # read and label the metadata
    metadata=None
    metadata_labels=None
    if args.input_metadata:
        metadata, samples_missing_metadata=utilities.read_metadata(args.input_metadata, otu_table, ignore_features=args.metadata_exclude, otu_table=True)
        metadata_labels, metadata=utilities.label_metadata(metadata, categorical=args.metadata_categorical, continuous=args.metadata_continuous)

    # get the introduction text if not provided by the user
    if not args.introduction_text:
        template_variables["log"]=log_file
        if not log_file:
            sys.exit("When running the workflow without a log file, please provide the introduction text with the option '--introduction-text <txt>'")
        template_variables["introduction_text"]=visualizations.Sixteen_S.compile_default_intro(template_variables)
    else:
        template_variables["introduction_text"]=args.introduction_text

    workflow_targets=workflow.name_output_files("16S_report."+args.format)

    # if metadata are provided then generate alpha diversity plots
    template_variables["alpha_diversity_plots"],alpha_task=utilities.generate_alpha_diversity_plots(workflow,"16S",args.output,args.input_metadata,otu_table,args.max_missing)
else:
    # set default introduction text
    if not args.introduction_text:
        with open(data.get_file("wmgx_methods.txt")) as file_handle:
            intro_methods=file_handle.readlines()[0]
        args.introduction_text = "![]("+utilities.get_package_file("wms_workflow","image")+"){#id .class width=675px height=505px}\n\n"+intro_methods

    # get the paths for the required files and check they are found
    qc_counts=utilities.find_data_file(data_files,"wmgx_qc_readcounts")
    taxonomic_profile=utilities.find_data_file(data_files,"wmgx_taxonomy")
    pathabundance=utilities.find_data_file(data_files,"both_function_pathway")
    ecsabundance=utilities.find_data_file(data_files,"wmgx_function_ec")
    read_counts=utilities.find_data_file(data_files,"wmgx_humann_counts")
    feature_counts=utilities.find_data_file(data_files,"wmgx_feature_counts")

    # read and label the metadata
    metadata=None
    metadata_labels=None
    if args.input_metadata:
        metadata, samples_missing_metadata=utilities.read_metadata(args.input_metadata, taxonomic_profile,
            name_addition="_taxonomic_profile", ignore_features=args.metadata_exclude)
        metadata_labels, metadata=utilities.label_metadata(metadata, categorical=args.metadata_categorical, continuous=args.metadata_continuous)

    template_depends=[taxonomic_profile]
    workflow_targets=workflow.name_output_files("wmgx_report."+args.format)

    # if the ecs do not include the names then add them
    ecsabundance=visualizations.add_ec_names(workflow,ecsabundance,args.output,template_depends)

    template_variables={"title":"Metagenome Report",
          "project":args.project_name,
          "introduction_text":args.introduction_text,
          "dna_read_counts":qc_counts,
          "is_paired":utilities.is_paired_table(qc_counts) if qc_counts else False,
          "taxonomic_profile":taxonomic_profile,
          "dna_pathabundance":pathabundance,
          "dna_ecabundance": ecsabundance,
          "read_counts":read_counts,
          "feature_counts":feature_counts,
          "log":log_file,
          "metadata":metadata,
          "metadata_labels":metadata_labels}

    # if metadata are provided then generate alpha diversity plots
    template_variables["alpha_diversity_plots"],alpha_task=utilities.generate_alpha_diversity_plots(workflow,"wmgx",args.output,args.input_metadata,taxonomic_profile,args.max_missing)

# add author and image if included
template_variables["author"]=args.author_name
template_variables["header_image"]=args.header_image
template_variables["study_type"]=workflow_type

# add formatting and plotting settings
template_variables["pdf_format"]=True if args.format == "pdf" else False
template_variables["min_abundance"] = float(args.min_abundance)
template_variables["min_samples"] = int(args.min_samples)
template_variables["max_sets_heatmap"] = int(args.max_sets_heatmap)
template_variables["max_sets_barplot"] = int(args.max_sets_barplot)
template_variables["max_groups_barplot"] = int(args.max_groups_barplot)

# add additional variables
template_variables["metadata"]=metadata
template_variables["metadata_labels"]=metadata_labels
template_variables["log"]=log_file

template_variables["correlation_threshold"]=float(args.correlation_threshold)

if args.print_template:
    # only print the template to stdout
    utilities.print_template(templates)

# use the template from the user if provided
if args.use_template:
    templates=[args.use_template]

# add the alpha task if needed
if alpha_task:
    template_depends+=[alpha_task]

# add the document to the workflow
doc_task=workflow.add_document(
    templates=templates,
    depends=template_depends, 
    targets=workflow_targets,
    vars=template_variables,
    table_of_contents=True)

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
# add the input metadata file if provided
optional_files=[]
if args.input_metadata:
    optional_files=[args.input_metadata]

workflow.add_archive(
    depends=[args.output,doc_task]+optional_files,
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
