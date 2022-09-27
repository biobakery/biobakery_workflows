#!/usr/bin/env python

"""
bioBakery Workflows: stats visualization workflow

Copyright (c) 2019 Harvard School of Public Health

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
import sys
import re

# import the workflow class from anadama2
from anadama2 import Workflow

# import the document templates and utilities from biobakery_workflows
from biobakery_workflows import utilities

# import the task to convert from biom to tsv
from biobakery_workflows.tasks.sixteen_s import convert_from_biom_to_tsv_list

# import the files for descriptions and paths
from biobakery_workflows import files, visualizations

# create a workflow instance, providing the version number and description
workflow = Workflow(version="3.1", remove_options=["input"],
                    description="A workflow for stats on wmgx and 16s data sets")

# add the custom arguments to the workflow                          
workflow.add_argument("input",desc="the folder containing taxonomy and functional data files",required=True)

# add the custom arguments to the workflow
workflow_vis = visualizations.Stats()
workflow.add_argument("project-name",desc="the name of the project", default="")
workflow.add_argument("author-name",desc="the name of the author of the report", default="")
workflow.add_argument("header-image",desc="the image to add to the report header", default="")
workflow.add_argument("input-metadata",desc="the metadata file (samples as columns or rows)", required=True)
workflow.add_argument("transform",desc="the transform to apply to the data with MaAsLin2 (default is the MaAsLin2 default transform)", default="")
workflow.add_argument("adonis-method",desc="the method to apply for the adonis, default is based on file type (bray for relab)", default="bray", choices=["manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", "chao", "cao", "mahalanobis"])
workflow.add_argument("fixed-effects",desc="the fixed effects to use in the models", default="")
workflow.add_argument("multivariable-fixed-effects",desc="the fixed effects that are multivariable (ordered first in covariate equation)", default="")
workflow.add_argument("random-effects",desc="the random effects to use in the models", default="")
workflow.add_argument("bypass-maaslin",desc="bypass running MaAsLiN", action="store_true")
workflow.add_argument("bypass-halla",desc="bypass running HAllA", action="store_true")
workflow.add_argument("maaslin-options",desc="additional MaAsLiN options", default="")
workflow.add_argument("halla-options",desc="additional HAllA options", default="")
workflow.add_argument("permutations",desc="the total number of permutations to apply to the permanova or mantel tests", default="4999")
workflow.add_argument("static-covariates",desc="the covariates, comma-delimited, that do not change per individual (to permutate within in permanova)", default="")
workflow.add_argument("scale",desc="the scale to apply with the permanova", default="100")
workflow.add_argument("min-abundance",desc="the min abundance to apply for filtering", default="0.0001")
workflow.add_argument("min-prevalence",desc="the min prevalence to apply for filtering", default="0.1")
workflow.add_argument("max-missing",desc="the max percentage of missing values for a metadata to not be filtered", default="20.0")
workflow.add_argument("format",desc="the format for the report", default="pdf", choices=["pdf","html"])
workflow.add_argument("top-pathways",desc="the top N significant pathways/metadata variables to plot stratified abundance", default=3)
workflow.add_argument("metadata-categorical",desc="the categorical features (for the plot stratified pathways)", action="append", default=[])
workflow.add_argument("metadata-continuous",desc="the continuous features (for the plot stratified pathways)", action="append", default=[])
workflow.add_argument("metadata-exclude",desc="the features to exclude (for the plot stratified pathways)", action="append", default=[])
workflow.add_argument("input-file-type",desc="the file type for an input file formatted as 'filename,filetype'", action="append", default=[])
workflow.add_argument("introduction-text",desc="the text to include in the intro of the report",
    default=workflow_vis.captions["intro"])
workflow.add_argument("print-template",desc="only print the template for the visualization workflow, do not run the workflow", action="store_true")
workflow.add_argument("use-template",desc="provide a report template to use instead of using that which is automatically generated", default="")

# get the arguments from the command line
args = workflow.parse_args()

# get the paths for the required files from the set of all input files
data_files=utilities.identify_data_files(files,args.input,args.input_file_type,args.input_metadata)

# error if no input files are found
if len(data_files.keys()) < 1:
    sys.exit("ERROR: No data files found in the input folder.")

# get the study type
study_type=utilities.get_study_type(data_files)

# get inputs based on study type
taxonomic_profile,pathabundance,other_data_files,study_type=utilities.get_input_files_for_study_type(data_files,study_type,workflow="stats")

# check for any biom files that need to be converted to txt
taxonomic_profile,pathabundance=convert_from_biom_to_tsv_list(workflow,[taxonomic_profile,pathabundance],args.output)
# add tasks to convert biom to tsv if needed
other_data_files=convert_from_biom_to_tsv_list(workflow,other_data_files,args.output)

# get metadata variables and check sample names
metadata_variables, samples_as_columns=utilities.get_metadata_variables(args.input_metadata,taxonomic_profile)

# check the options include valid variables
utilities.check_effects_are_included_in_metadata(args.fixed_effects, args.random_effects, metadata_variables)

## Add tasks to the workflow 

## 0. See if names need to be added to ecs
template_depends=[]
try:
    other_data_files_ecs=list(filter(lambda x: x[1] == "ec", other_data_files.items()))[0]
except IndexError:
    other_data_files_ecs=[]

if other_data_files_ecs:
    new_ecs_file=visualizations.add_ec_names(workflow,other_data_files_ecs[0],args.output,template_depends)
    del other_data_files[other_data_files_ecs[0]]
    other_data_files[new_ecs_file]="ec"

## 1. Add tasks to create feature tables from all input files to be used. Feature tables will be used in downstream tasks                   .

feature_tasks_info=utilities.create_feature_table_inputs(workflow,study_type,args.output,taxonomic_profile,pathabundance,other_data_files)

## 2. Add tasks to run mantel tests on all data file pairs

mantel_tasks,mantel_plots=utilities.run_mantel_tests(workflow,feature_tasks_info,args.input_metadata,args.output,args.permutations)

## 3. Add tasks to run MaAsLiN2 on all input data files (feature tables)

all_maaslin_tasks=[]
if not args.bypass_maaslin:
    maaslin_tasks=utilities.run_maaslin_on_input_file_set(workflow,feature_tasks_info,args.input_metadata,args.transform,args.fixed_effects,args.random_effects,args.maaslin_options)
    maaslin_tiles_task=workflow.add_task(
        utilities.partial_function(utilities.generate_tiles_of_maaslin_figures, feature_tasks_info=feature_tasks_info),
        depends=maaslin_tasks)
    all_maaslin_tasks=maaslin_tasks+[maaslin_tiles_task]

## 4. Add tasks to run halla on all sets of data files            

halla_tasks_info=[]
halla_tasks=[]
if not args.bypass_halla:
    halla_tasks,halla_tasks_info=utilities.run_halla_on_input_file_set(workflow,feature_tasks_info,args.input_metadata,args.output,args.halla_options, samples_as_columns)

## 5. Add tasks to run stratified pathways plots, if pathways provided  

stratified_plots_tasks=[]
stratified_pathways_plots=[]
if not args.bypass_maaslin:
    stratified_pathways_plots,stratified_plots_tasks=utilities.create_stratified_pathways_plots(workflow,study_type,pathabundance,args.input_metadata,args.metadata_exclude,args.metadata_categorical,args.metadata_continuous,args.top_pathways,feature_tasks_info,args.output)

## 6. Add tasks to run permanova if longitudinal (if random effects are set) and if not run beta diversity script

additional_stats_tasks=[]
permanova_plots=[]
beta_diversity_plots={"univariate": {}, "multivariate": {}, "pairwise": {}}
covariate_equation=""

if args.random_effects:
    additional_stats_tasks,permanova_plots=utilities.run_permanova(workflow,args.static_covariates,feature_tasks_info,args.input_metadata,args.scale,args.min_abundance,args.min_prevalence,args.permutations,args.output,additional_stats_tasks)
else:
    additional_stats_tasks,beta_diversity_plots,covariate_equation=utilities.run_beta_diversity(workflow,feature_tasks_info,args.input_metadata,args.min_abundance,args.min_prevalence,args.max_missing,[args.multivariable_fixed_effects,args.fixed_effects],args.output,additional_stats_tasks,args.random_effects,metadata_variables,args.adonis_method)

templates=[utilities.get_package_file("stats")]

if args.print_template:
    # only print the template to stdout
    utilities.print_template(templates)

# use the template from the user if provided
if args.use_template:
    templates=[args.use_template]

# update the intro text if halla is not run and the default intro is used
if args.bypass_halla and args.introduction_text == workflow_vis.captions["intro"]:
    args.introduction_text = workflow_vis.captions["intro_bypass_halla"]

# add the document to the workflow
doc_task=workflow.add_document(
    templates=templates,
    depends=all_maaslin_tasks+stratified_plots_tasks+[taxonomic_profile]+additional_stats_tasks+halla_tasks+mantel_tasks+template_depends, 
    targets=workflow.name_output_files("stats_report."+args.format),
    vars={"title":"Statistics report",
          "project":args.project_name,
          "author":args.author_name,
          "header_image":args.header_image,
          "introduction_text":args.introduction_text,
          "taxonomic_profile":taxonomic_profile,
          "mantel_plots":mantel_plots,
          "feature_tasks_info":feature_tasks_info,
          "halla_tasks_info":halla_tasks_info,
          "bypass_maaslin":args.bypass_maaslin,
          "bypass_halla":args.bypass_halla,
          "stratified_pathways_plots":stratified_pathways_plots,
          "top_pathways":args.top_pathways,
          "permanova_plots":permanova_plots,
          "beta_diversity_plots":beta_diversity_plots,
          "covariate_equation":covariate_equation,
          "pdf_format":True if args.format == "pdf" else False},
    table_of_contents=True)

# add an archive of the document and figures, removing the log file
# the archive will have the same name and location as the output folder
workflow.add_archive(
    depends=[args.output,doc_task,args.input_metadata],
    targets=args.output+".zip",
    remove_log=True)

# start the workflow
workflow.go()
