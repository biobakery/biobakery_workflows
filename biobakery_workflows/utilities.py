"""
bioBakery Workflows: utilities module
Utility functions for workflows and tasks

Copyright (c) 2016 Harvard School of Public Health

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
import math
import functools
import time
import collections
import re

from anadama2.tracked import TrackedDirectory

# try to import urllib.request.urlretrieve for python3
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

MIN_SAMPLES_DATA_FILE = 3
TAXONOMY_DELIMITER = "|"
MAX_METADATA_CATEGORIES = 10

def run_mantel_tests(workflow,maaslin_tasks_info,output,nperm):

    # name the output file and folder
    mantel_plots=[name_files("mantel_plot.png",output,subfolder="mantel_test",create_folder=True)]

    input_files=[data[1][0] for data in maaslin_tasks_info.items()]
    workflow.add_task(
        "[args[0]] '[args[1]]' [targets[0]] --permutations [args[2]]",
        depends=input_files,
        targets=mantel_plots,
        name="mantel_test",
        args=[get_package_file("mantel_test", "Rscript"),",".join(input_files),nperm])

    return mantel_plots

def get_metadata_variables(input_metadata, taxonomic_profile):
    # get the metadata variables (might be columns or rows)

    with open(taxonomic_profile) as file_handle:
        samples = file_handle.readline().rstrip().split("\t")[1:-1]

    row_names = []
    col_names = []
    with open(input_metadata) as file_handle:
        for line in file_handle:
            if not col_names:
                col_names=line.rstrip().split("\t")[1:]
            else:
                row_names.append(line.rstrip().split("\t")[0])

    if set(samples).intersection(set(col_names)):
        metadata_variables = row_names
    else:
        metadata_variables = col_names

    # check if samples start with numbers
    start_number = [x for x in samples if x[0].isdigit()]
    if start_number:
        sys.exit("ERROR: Sample names start with a number: "+",".join(start_number)+".")

    # check for duplicate samples
    if len(list(set(samples))) != len(samples):
        duplicate=[x for x, count in collections.Counter(samples).items() if count > 1]
        sys.exit("ERROR: Duplicate samples in the taxonomic profile: "+",".join(duplicate)+".")

    return metadata_variables

def run_permanova(workflow,static_covariates,maaslin_tasks_info,input_metadata,scale,min_abundance,min_prevalence,permutations,output,additional_stats_tasks):
    # if longitudinal run the permanova

    permanova_plots = {}
    if static_covariates:
        optional_args=" --static_covariates "+static_covariates
    else:
        sys.exit("ERROR: Please provide the individual covariates when running with longitudinal metadata (ie --individual-covariates='age,gender')")

    input_files=[]
    permanova_plots["all"]=name_files("permanova.png",output,subfolder="permanova",create_folder=True)

    for filetype in maaslin_tasks_info.keys():
        input_files.append(maaslin_tasks_info[filetype][0])
        permanova_plots[filetype]=name_files("permanova_{}.png".format(filetype),output,subfolder="permanova",create_folder=True)

    permanova_script_path = get_package_file("permanova_hmp2", "Rscript")

    additional_stats_tasks.append(
        workflow.add_task(
            "[args[0]] '[args[1]]' [depends[0]] [targets[0]] --scale [args[2]] --min_abundance [args[3]] --min_prevalence [args[4]] --permutations [args[5]] [args[6]]",
            depends=[input_metadata]+input_files,
            targets=list(permanova_plots.values()),
            args=[permanova_script_path,",".join(input_files),scale,min_abundance,min_prevalence,permutations,optional_args],
            name="hmp2_permanova"))

    return additional_stats_tasks,permanova_plots


def run_beta_diversity(workflow,maaslin_tasks_info,input_metadata,min_abundance,min_prevalence,max_missing,fixed_effects,output,additional_stats_tasks,random_effects,metadata_variables,adonis_method):
    # if not longitudinal then run univariate plus multivariate if set

    # if set, add the adnois method
    optional_args=""
    if adonis_method:
        optional_args=" --adonis_method "+adonis_method

    # construct the equation for the model based on the fixed effects provided
    ordered_fixed_effects=list(collections.OrderedDict.fromkeys(",".join(fixed_effects).split(",")).keys())
    covariate_equation=""
    if len(ordered_fixed_effects) > 1:
        covariate_equation=" + ".join(ordered_fixed_effects)

    # determine the covariate equation from the metadata if not provided
    if not covariate_equation and len(metadata_variables) > 1:
        # discard subject if found
        metadata_variables_set = set(metadata_variables)
        metadata_variables_set.discard("subject")
        fixed_effects = list(metadata_variables_set.difference(set(random_effects.split(","))))
        # only include a multivariate equation if more than one variable
        if len(fixed_effects) > 1:
            covariate_equation=" + ".join(fixed_effects)

    beta_diversity_plots = {"univariate": {}, "multivariate": {}}
    univariate_script_path = get_package_file("beta_diversity", "Rscript")
    for filetype in maaslin_tasks_info.keys():
        univariate=name_files(filetype+"_univariate.png",output,subfolder="beta_diversity",create_folder=True)

        additional_stats_tasks.append(
            workflow.add_task(
                "[args[0]] [depends[0]] [depends[1]] [targets[0]] --min_abundance [args[1]] --min_prevalence [args[2]] --max_missing [args[3]]"+optional_args,
                depends=[maaslin_tasks_info[filetype][0],input_metadata],
                targets=univariate,
                args=[univariate_script_path,min_abundance,min_prevalence,max_missing],
                name="beta_diversity_univarite_"+filetype))
        beta_diversity_plots["univariate"][filetype]=univariate

    if covariate_equation:
        for filetype in maaslin_tasks_info.keys():
            multivariate=name_files(filetype+"_multivariate.png",output,subfolder="beta_diversity",create_folder=True)

            additional_stats_tasks.append(
                workflow.add_task(
                    "[args[0]] [depends[0]] [depends[1]] [targets[0]] --min_abundance [args[1]] --min_prevalence [args[2]] --max_missing [args[3]] --covariate_equation='[args[4]]'"+optional_args,
                    depends=[maaslin_tasks_info[filetype][0],input_metadata],
                    targets=multivariate,
                    args=[univariate_script_path,min_abundance,min_prevalence,max_missing,covariate_equation],
                    name="beta_diversity_multivariate_"+filetype))
            beta_diversity_plots["multivariate"][filetype]=multivariate

    return additional_stats_tasks,beta_diversity_plots,covariate_equation


def create_stratified_pathways_plots(workflow,study_type,pathabundance,input_metadata,metadata_exclude,metadata_categorical,metadata_continuous,top_pathways,maaslin_tasks_info,output):
    # if pathways are provided then generate stratified plots

    stratified_pathways_plots = []
    stratified_plots_tasks = []

    if pathabundance and study_type=="wmgx":
        # read in the metadata to merge with the data for the barplot script
        metadata=read_metadata(input_metadata, pathabundance,
            name_addition="_Abundance", ignore_features=metadata_exclude)

        metadata_labels, metadata=label_metadata(metadata, categorical=metadata_categorical, continuous=metadata_continuous)
        # get all continuous or samples ids and remove (as they are not to be used for the plots)
        metadata_exclude=metadata_exclude+[x for x,y in filter(lambda x: x[1] == "con", metadata_labels.items())]
        for metadata_row in metadata[1:]:
            if len(list(set(metadata_row[1:]))) > MAX_METADATA_CATEGORIES:
                metadata_exclude+=[metadata_row[0]]
        metadata_exclude=list(set(metadata_exclude))
        metadata=read_metadata(input_metadata, pathabundance,
            name_addition="_Abundance", ignore_features=metadata_exclude)
        metadata_labels, metadata=label_metadata(metadata, categorical=metadata_categorical, continuous=metadata_continuous)

        humann_barplot_input = name_files("merged_data_metadata_input.tsv", output, subfolder="stratified_pathways", create_folder=True)
        workflow.add_task(
            partial_function(create_merged_data_file, metadata=metadata, name_addition="_Abundance"),
            depends=pathabundance,
            targets=humann_barplot_input)

        # only include the categorical metadata
        metadata_row_names=[row[0] for row in metadata[1:] if row[0] in metadata_labels.keys()]
        metadata_end=metadata_row_names[-1]
        for i in range(top_pathways):

            new_pathways_plot=name_files("stratified_pathways_{0}.jpg".format(i), output, subfolder="stratified_pathways")
            stratified_plots_tasks.append(workflow.add_task(
                partial_function(run_humann_barplot, number=i, metadata_end=metadata_end, categorical=list(metadata_labels.keys())),
                depends=[maaslin_tasks_info["pathways"][2],humann_barplot_input],
                targets=new_pathways_plot,
                name="run_humann_barplot_pathway_{0}".format(i)))
            stratified_pathways_plots.append(new_pathways_plot)

    return stratified_pathways_plots,stratified_plots_tasks


def generate_tile_of_images(input_files, output_file):
    # for the set of images, generate a single tile with a table of images
    import matplotlib.pyplot as pyplot

    rows=3
    columns=3

    figure = pyplot.figure(figsize=(8,8),dpi=300)
    for index in range(1, columns*rows+1):
        try:
            new_file=input_files[index]
        except IndexError:
            break

        image = pyplot.imread(new_file)
        subplot=figure.add_subplot(rows, columns, index, frame_on=False)
        subplot.xaxis.set_visible(False)
        subplot.yaxis.set_visible(False)

        pyplot.imshow(image, interpolation="none")
    pyplot.tight_layout()
    pyplot.draw()

    pyplot.savefig(output_file, dpi=300)


def get_maaslin_image_files(maaslin_tasks_info):
    # Get a list of all the generaged maaslin2 image files
    maaslin_tiles={}
    metadata_images={}
    for datatype in maaslin_tasks_info:
        # get all the available images and sort by rank
        figures_folder = os.path.dirname(maaslin_tasks_info[datatype][1])
        ranked_files = [ (filename, re.findall(r'\d+', filename)[0]) for filename in os.listdir(figures_folder) if re.search(r'_\d+.png$',filename)]
        ordered_files = sorted( ranked_files, key=lambda x: int(x[1]))

        # group images by metadata type
        metadata_images[datatype]={}
        for file_name, rank in ordered_files:
            if file_name.endswith("_{}.png".format(rank)):
                images_found = True
                metadata_name=file_name.replace("_{}.png".format(rank),"")
                if not metadata_name in metadata_images[datatype]:
                    metadata_images[datatype][metadata_name]=[]
                metadata_images[datatype][metadata_name].append(os.path.join(figures_folder,file_name))

        maaslin_tiles[datatype]={}
        for metadata_name in metadata_images[datatype]:
            new_image="{}_tiled.png".format(os.path.join(figures_folder,metadata_name))
            maaslin_tiles[datatype][metadata_name]=new_image

    return metadata_images, maaslin_tiles

def generate_tiles_of_maaslin_figures(task, maaslin_tasks_info):
    # Generate a tile of the top plots for each metadata for each maaslin run to be displayed in the report

    metadata_images, maaslin_tiles = get_maaslin_image_files(maaslin_tasks_info)
    for datatype in maaslin_tasks_info:
        for metadata_name in metadata_images[datatype]:
            run_task(
                "create_image_tile.py --input '[args[0]]' --output [targets[0]]",
                depends=metadata_images[datatype][metadata_name],
                targets=maaslin_tiles[datatype][metadata_name],
                args=",".join(metadata_images[datatype][metadata_name]))

def run_halla_on_input_file_set(workflow,maaslin_tasks_info,output,halla_options=""):
    # Run maaslin on all files in input set
    
    halla_tasks=[]
    halla_tasks_info={}
    for run_type, infiles in maaslin_tasks_info.items():
        for run_type2, infiles2 in maaslin_tasks_info.items():
            if run_type == run_type2 or run_type2+" "+run_type in halla_tasks_info or "gene" in run_type or "gene" in run_type2:
                continue

            current_target=os.path.join(os.path.abspath(output),"halla_"+run_type+"_"+run_type2,"hallagram.png")
            halla_tasks.append(
                workflow.add_task(
                    "halla -x [depends[0]] -y [depends[1]] -o [args[0]]",
                    depends=[infiles[0], infiles2[0]],
                    targets=current_target,
                    args=os.path.dirname(current_target),
                    name="HAllA_{0}_{1}".format(run_type,run_type2)))

            halla_tasks_info[run_type+" "+run_type2]=current_target

    return halla_tasks, halla_tasks_info

def run_maaslin_on_input_file_set(workflow,maaslin_tasks_info,input_metadata,transform,fixed_effects,random_effects,maaslin_options=""):
    # Run maaslin on all files in input set
    
    maaslin_tasks=[]
    maaslin_optional_args=maaslin_options
    # add comma if not included
    if maaslin_optional_args and not maaslin_optional_args.startswith(","):
        maaslin_optional_args=","+maaslin_optional_args

    if transform:
        maaslin_optional_args+=",transform='"+transform+"'"
    if fixed_effects:
        maaslin_optional_args+=",fixed_effects='"+fixed_effects+"'"
    if random_effects:
        maaslin_optional_args+=",random_effects='"+random_effects+"'"

    for run_type, (maaslin_input_file, maaslin_heatmap, maaslin_results_table) in maaslin_tasks_info.items():
        maaslin_tasks.append(
            workflow.add_task(
                "R -e \"library('Maaslin2'); results <- Maaslin2('[depends[0]]','[depends[1]]','[args[0]]'"+maaslin_optional_args+")\"",
                depends=[maaslin_input_file, input_metadata],
                targets=maaslin_results_table,
                args=os.path.dirname(maaslin_results_table),
                name="R_Maaslin2_{}".format(run_type)))

    return maaslin_tasks

def create_maaslin_feature_table_inputs(workflow,study_type,output,taxonomic_profile,pathabundance,other_data_files):
    # For all input files based on type create feature tables for input to maaslin

    taxon_feature=name_files("taxonomy_features.txt",output,subfolder="features",create_folder=True)
    create_feature_table_tasks_info=[]
    if study_type == "wmgx":
        create_feature_table_tasks_info=[(taxonomic_profile,taxon_feature,"--sample-tag-column '_taxonomic_profile' --reduce-stratified-species-only")]
    else:
        # reformat this table to move the taxonomic column and sum for species if 16s data
        workflow.add_task(
            "trim_taxonomy.py --input [depends[0]] --output [targets[0]] --end-taxonomy-column 0",
             depends=taxonomic_profile,
             targets=taxon_feature)

    maaslin_tasks_info={"taxonomy":(taxon_feature,name_files("heatmap.png", output, subfolder=os.path.join("maaslin2_taxa","figures")),
        name_files("significant_results.tsv", output, subfolder="maaslin2_taxa"))}

    if pathabundance:
        pathabundance_feature=name_files("pathways_features.txt",output,subfolder="features",create_folder=True)
        create_feature_table_tasks_info.append((pathabundance,pathabundance_feature,"--sample-tag-column '_Abundance' --remove-stratified"))
        maaslin_tasks_info["pathways"]=(pathabundance_feature,name_files("heatmap.png", output, subfolder=os.path.join("maaslin2_pathways","figures")),
            name_files("significant_results.tsv", output, subfolder="maaslin2_pathways"))

    for newfile in other_data_files:
        newfile_type = other_data_files[newfile]
        new_feature=name_files(newfile_type+"_features.txt",output,subfolder="features",create_folder=True)
        new_subfolder="maaslin2_"+newfile_type
        create_feature_table_tasks_info.append((newfile,new_feature,"--sample-tag-column '_Abundance' --remove-stratified"))
        maaslin_tasks_info[newfile_type]=(new_feature,name_files("heatmap.png", output, subfolder=os.path.join(new_subfolder,"figures")),
            name_files("significant_results.tsv", output, subfolder=new_subfolder))

    for input_file, output_file, options in create_feature_table_tasks_info:
        workflow.add_task(
            "create_feature_table.py --input [depends[0]] --output [targets[0]] [args[0]]",
            depends=input_file,
            targets=output_file,
            args=[options])

    return maaslin_tasks_info

def get_input_files_for_study_type(data_files, study_type):
    # based on the type of study, find the input files in the input folder

    # if study is of type "both" then first look for wmgx taxonomy file
    taxonomic_profile=find_data_file(data_files,"wmgx_taxonomy",required=False)

    if study_type=="wmgx" or (taxonomic_profile and study_type=="both"):
        study_type="wmgx"

        # get the paths for the optional files from the set of input files
        pathabundance=find_data_file(data_files, "function_pathway", required=False)
        other_data_files=dict([(filename[0], type.split("_")[-1]) for type, filename in data_files.items() if not filename[0] in [taxonomic_profile,pathabundance]])

    else:
        taxonomic_profile=find_data_file(data_files,"16s_taxonomy",required=True)

        # get the paths for the optional files from the set of input files
        pathabundance=find_data_file(data_files,"function_pathway", required=False)
        other_data_files=dict([(filename[0], type.split("_")[-1]) for type, filename in data_files.items() if not filename[0] in [taxonomic_profile,pathabundance]])

    return taxonomic_profile,pathabundance,other_data_files,study_type

# create a merged metadata table to be used as input for humann_barplot
def create_merged_data_file(task, metadata, name_addition):
    # read in the pathabundance file
    data = []
    with open(task.depends[0].name) as file_handle:
        samples = [i.split(name_addition)[0] for i in file_handle.readline().rstrip().split("\t")[1:]]
        for line in file_handle:
            line=line.rstrip().split("\t")
            data.append(line)    

    merged_data, metadata_samples=merge_metadata(metadata, samples, data)

    with open(task.targets[0].name, "w") as file_handle:
        header = "\t".join(["Feature"]+metadata_samples)+"\n"
        file_handle.write(header)
        for line in merged_data:
            file_handle.write("\t".join(line)+"\n")


# gather the top pathways to plot from maaslin2 outputs
# only gather categorical metadata features
def gather_top_N_associations_maaslin2_results(filename, N, categorical):
    associations=[]
    with open(filename) as file_handle:
        for line in file_handle:
            data = line.rstrip().split("\t")
            associations.append([data[0], data[1]])

    # use N+1 to allow for the header value
    try:
        selected_pathway, metadata_focus = associations[N+1]
        while not metadata_focus in categorical:
            N=N+1
            selected_pathway, metadata_focus = associations[N+1]
    except IndexError:
        selected_pathway, metadata_focus = "", ""

    return selected_pathway, metadata_focus

def run_humann_barplot(task, number, metadata_end, categorical):
    # determine the pathway name
    try:
        original_selected_pathway, metadata_focus = gather_top_N_associations_maaslin2_results(task.depends[0].name, number, categorical)
    except IndexError:
        original_selected_pathway, metadata_focus = None

    if original_selected_pathway:
        # only use pathway name and replace periods if present with dash
        selected_pathway = original_selected_pathway
        if not "-" in selected_pathway:
            selected_pathway = "-".join(original_selected_pathway.split(".")[0:2])

        if not "pwy" in selected_pathway.lower():
            selected_pathway = "-".join(original_selected_pathway.split(".")[0:3])

        run_task(
            "humann_barplot --input [depends[1]] --focal-feature [args[0]] --output [targets[0]] --last-metadatum [args[1]] --focal-metadatum [args[2]] --sort [args[3]] && echo '[args[2]]' > [targets[1]]",
            depends=task.depends,
            targets=task.targets+[task.targets[0].name.replace(".jpg",".txt")],
            args=[selected_pathway, metadata_end, metadata_focus, "metadata"])
    else:
        run_task("touch [targets[0]] && touch [targets[1]]",
            depends=task.depends,
            targets=task.targets+[task.targets[0].name.replace(".jpg",".txt")])

def find_data_file(data_files, type, required=False):
    """ Return an error if the file of that type has not been found """

    not_found = False
    try:
        file_name = data_files[type][0]
    except KeyError:
        not_found = True

    if not_found:
        # look for a partial key
        for data_type in data_files.keys():
            if data_type.startswith(type) or data_type.endswith(type):
                file_name = data_files[data_type][0]
                not_found = False

    if not_found:
        if required:
            sys.exit("ERROR: A required file of type "+type+" can not be found in the input folder")
        else:
            file_name=""

    return file_name

def get_study_type(data_files):
    """ Determine the type of study based on the data files """

    types = [file_type.split("_")[0] for file_type in data_files.keys()]

    type = set(types)
    type.discard("both")

    if "16s" in type and "wmgx" in type:
        files = [key+"\t"+value[0] for key,value in data_files.items()]
        sys.exit("ERROR: Input files found of multiple study types:" + ",".join(types)+"\n"+"\n".join(files))
    elif "16s" in type:
        return "16s"
    elif "wmgx" in type:
        return "wmgx"
    elif "both" in types:
        return "both"

    return types[0]
    

def identify_data_files(folder,input_file_type,metadata_input):
    """ For all files in the folder and subfolders, return all tab delimited files with their data type 

        Args:
            folder (string): The path to the main folder
            input_file_type (array): The user provided file types
            metadata_input (string): Full path to the user provided metadata file

        Returns:
            dict : A dictionary of data files organised by type 
    """

    class openfile_txt_and_biom_gz():
        def __init__(self, filename):
            self.filename = filename

        def readline(self):
            try:
                line = self.lines.pop(0)
            except IndexError:
                line = ""
            return line

        def __enter__(self):
            if self.filename.endswith(".biom"):
                import biom
                self.lines = biom.load_table(self.filename).to_tsv().split("\n")[:100]
            elif self.filename.endswith(".gz"):
                import gzip
                self.lines = gzip.open(self.filename,"rt").readlines()[:100]
            else:
                self.lines = open(self.filename).readlines()[:100]
            return self

        def __exit__(self, exc_type, exc_value, traceback):
            self.lines = []

    # remove the path from the metadata input file
    metadata_input = os.path.basename(metadata_input)

    # set the user provided data file types
    known_filetypes = dict([set.split(",") for set in input_file_type])

    # get all of the files of the tab delimited type
    data_files = []
    for root, dir, files in os.walk(folder):
        for file_name in files:
            if file_name.endswith(".tsv") or file_name.endswith(".txt") or file_name.endswith(".biom") or file_name.endswith(".gz"):
                data_files.append(os.path.join(root,file_name))

    # determine the type of each file
    data_files_types = {}
    for file in data_files:

        # ignore the metadata input file
        if os.path.basename(file) == metadata_input:
            continue

        with openfile_txt_and_biom_gz(file) as file_handle:
            file_type = None

            # ignore all comment lines and header
            header = file_handle.readline()
            next_line = file_handle.readline()
            while next_line.startswith("#"):
                next_line = file_handle.readline()
            header = next_line

            # check this is a data file for at least a few samples
            first_data_line = file_handle.readline()
            data_info = first_data_line.split("\t")
            if len(data_info) > MIN_SAMPLES_DATA_FILE + 1:
                # remove stratification if present
                if TAXONOMY_DELIMITER in data_info[0]:
                    data_info[0]=data_info[0].split(TAXONOMY_DELIMITER)[0]                    
                
                # check the first column for the data type
                if "." in data_info[0] and data_info[0].split(":")[0].replace(".","").isdigit():
                    file_type = "wmgx_function_ec"
                elif "pwy" in data_info[0].lower():         
                    file_type = "both_function_pathway"
                elif "k__" in data_info[0] or "s__" in data_info[0] or "g__" in data_info[0]:         
                    file_type = "wmgx_taxonomy"
                elif data_info[0].startswith("K0"):
                    file_type = "16s_function_gene"
                elif data_info[0].startswith("EC"):
                    file_type = "16s_function_ec"
                elif data_info[0].startswith("ko"):
                    file_type = "16s_function_pathway"
                elif data_info[0].startswith("M0"):
                    file_type = "16s_function_module"
                elif data_info[0].startswith("ASV") and data_info[-1].lower().startswith("k__"):
                    file_type = "16s_taxonomy_asv"
                elif data_info[-1].lower().startswith("k__"):
                    file_type = "16s_taxonomy_otu"

                # replace with user provided type, if set
                if file in known_filetypes:
                    # check for valid file type
                    file_type = known_filetypes[file]       
                    file_type_split_info = file_type.split("_")
                    if not (file_type_split_info[0] in ["wmgx","16s"] and file_type_split_info[1] in ["function","taxonomy"] and len(file_type_split_info)==3):
                        sys.exit("Please provide a valid file type of the format [wmgx|16s]_[function|taxonomy]_[type*] (valid functions types [ec|pathway|gene|module] and valid taxonomy types for 16s [otu|asv]) replacing the input provided of '"+file_type+"'.")

                if not file_type:
                    # set the file type to the name of the file if not identified (remove underscore as those are used for file type info)
                    file_type = os.path.basename(file).split(".")[0].replace("_","")

                if not file_type in data_files_types:
                    data_files_types[file_type]=[]
                data_files_types[file_type].append(file)
  
    return data_files_types

def get_package_file(basename, type="template"):
    """ Get the full path to a file included in the installed python package.

        Args:
            basename (string) : The basename of the file
            type (string) : The type of file to find (template or Rscript)

        Returns: 
            string : The full path to the file

    """

    if type == "template":
        subfolder = "document_templates"
        extension = ".template.py"
    else:
        subfolder = "Rscripts"
        extension = ".R"

    # get all of the templates in this folder
    package_install_folder=os.path.join(os.path.dirname(os.path.realpath(__file__)), subfolder)
    found_files=list(filter(lambda file: file.endswith(extension),os.listdir(package_install_folder)))

    # return the template with the name
    matching_file=list(filter(lambda file: file.startswith(basename+extension), found_files))

    if matching_file:
        matching_file=os.path.join(package_install_folder,matching_file[0])
    else:
        matching_file=""

    return matching_file
    


def change_pweave_figure_size_heatmap(pdf_format):
    """ Change the figure size for heatmaps based on the output format"""
    fig_size = (4,4) if pdf_format else (2.5,2.5)
    change_pweave_figure_size(fig_size)
    
def reset_pweave_figure_size():
    """ Set the pweave figure size back to the default """
    change_pweave_figure_size((8,6))
    
def change_pweave_figure_size(fig_size):
    """ Change the pweave default figure size """
    import pweave
    pweave.rcParams["chunk"]["defaultoptions"].update({'f_size': fig_size})

def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """
    
    return byte / (1024.0**2)

class ReportHook():
    def __init__(self):
        self.start_time=time.time()
        
    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        
        if blocknum == 0:
            self.start_time=time.time()
            if total_size > 0:
                print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
        else:
            total_downloaded=blocknum*block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))
                    
            if total_size > 0:
                percent_downloaded=total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stdout to overwrite stdout
                try:
                    download_rate=total_downloaded/(time.time()-self.start_time)
                    estimated_time=(total_size-total_downloaded)/download_rate
                except ZeroDivisionError:
                    download_rate=0
                    estimated_time=0
                estimated_minutes=int(estimated_time/60.0)
                estimated_seconds=estimated_time-estimated_minutes*60.0
                status +="{:3.2f}".format(percent_downloaded) + " %  " + \
                    "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
                    "{:2.0f}".format(estimated_minutes) + " min " + \
                    "{:2.0f}".format(estimated_seconds) + " sec "
            status+="        \r"
            sys.stdout.write(status)

def read_file_n_lines(file,n):
    """ Read a file n lines at a time """

    line_set=[]
    with open(file) as file_handle:
        for line in file_handle:
            if len(line_set) == n:
                yield line_set
                line_set=[]
            line_set.append(line)

    # yield the last set
    if len(line_set) == n:
        yield line_set
            
def read_file_catch(file, delimiter="\t"):
    """ Try to read the file, catch on error. Split data by delimiter. """
    
    try:
        handle=open(file)
        lines=handle.readlines()
        handle.close()
    except EnvironmentError:
        sys.exit("Error: Unable to read file: "+file)
    
    # split rows by delimiter
    new_lines=[line.rstrip().split(delimiter) for line in lines]    
    
    return new_lines
            
def read_metadata(metadata_file, taxonomy_file, name_addition="", ignore_features=[], otu_table=False):
    """ Read in the metadata file. Samples can be rows or columns.
    Ignore features if set. 
        
    Args:
        metadata_file (string): The path to the metadata file.
        taxonomy_file (string): The path to a taxonomy file (or any file with
            the sample names as the column names).
        name_addition (string): Any strings that should be removed from the
            names in the taxonomy files to get the sample names.
        ignore_features (list): A list of strings of features to ignore
            when reading the metadata file.
        otu_table (bool): Set if otu table is used for taxonomy to allow for different
            header format.
        
    Returns:
        (list): A list of lists of the metadata (samples as columns)
    """
    
    # read in a taxonomy file to get the sample names from the columns
    if otu_table:
        samples=set([name.replace(name_addition,"") for name in read_file_catch(taxonomy_file)[0][1:-1]])
    else:
        samples=set([name.replace(name_addition,"") for name in read_file_catch(taxonomy_file)[0][1:]])
        
    # read in the metadata file
    data=read_file_catch(metadata_file)
    # check if the columns or rows are samples
    possible_samples=data[0][1:]
    overlap=samples.intersection(possible_samples)
    if len(list(overlap)) == 0:
        # the samples must be the rows so invert the data
        data=[list(a) for a in zip(*data)]
        possible_samples=data[0][1:]
        overlap=samples.intersection(possible_samples)
    
    # check for samples not included in metadata
    if len(list(overlap)) < len(list(samples)):
        sys.exit("ERROR: Not all of the samples in the data set have"+
            " metadata. Please review the metadata file. The following samples"+
            " were not found: "+",".join(list(samples.difference(possible_samples))))

    # remove any features that should be ignored
    new_data=[]
    for row in data:
        if not row[0] in ignore_features:
            new_data.append(row)
        else:
            ignore_features.remove(row[0])

    # check for any features that do not vary
    check_diff = {}
    for row in new_data[1:]:
        for value in row[1:]:
            if not row[0] in check_diff:
                check_diff[row[0]]=[]
            check_diff[row[0]]=check_diff[row[0]]+[value]

    for feature_name, feature_types in check_diff.items():
        if len(list(set(feature_types))) == 1:
            sys.exit("ERROR: Please remove feature '"+feature_name+"' as it"+
                " only includes a single type of '"+feature_types[0]+"'.")

    # check for any features that were not found
    if ignore_features:
        sys.exit("ERROR: Unable to find features that should be ignored: "+
            ",".join(ignore_features))
        
    return new_data

def label_metadata(data, categorical=[], continuous=[]):
    """ Label the metadata type. All numerical is continous.
    
    Args:
        data (lists of lists): A list of metadata 
        categorical (list): A list of categorical features.
        continuous (list): A list of continuous features.
        
    Returns:
        (dict): A dictionary of metadata labels
        (list): A list of lists of the metadata (converted to floats if continuous)
    """
    
    # add labels to the metadata, ignore sample names, convert to nan misisng values if continuous
    missing = ["Unknown", "unknown", "NA", "na", "nan", "NaN", "NAN", " "]
    labeled_data=[data[0]]
    labels={}
    for row in data[1:]:
        # check if there are missing values in row, convert to nan
        for index, item in enumerate(row):
            if item in missing or item == "":
                row[index] = "nan"
        # apply specific labels if set
        label=None
        if row[0] in continuous:
            try:
                row[1:] = map(float, row[1:])
                label="con"
            except ValueError:
                label="cat"
            continuous.remove(row[0])
        if row[0] in categorical:
            label="cat"
            categorical.remove(row[0])
        if not label:
            try:
                row[1:] = map(float, row[1:])
                label="con"
            except ValueError:
                label="cat"
        if label == "cat":
            # if label is categorical convert misisng values to Unknown
            for index, item in enumerate(row):
                if item == "nan":
                    row[index] = "NA"
        labeled_data.append(row)
        labels[row[0]]=label
        
    # check for remaining labels
    if categorical:
        sys.exit("ERROR: Unable to find and label categorical feature in metadata: "+
            ",".join(categorical))
    if continuous:
        sys.exit("ERROR: Unable to find and label continuous feature in metadata: "+
            ",".join(continuous))
        
    return labels, labeled_data

def filter_metadata_categorical(metadata, metadata_labels):
    """ Return only the metadata that is categorical 
    
    Args:
        metadata (lists of lists): A list of metadata.
        metadata_labels (dict): The labels for the metadata.
        
    Returns:
        metadata (lists of lists): Only the categorical metadata.
    """
    
    new_metadata=[]
    for row in metadata:
        if metadata_labels.get(row[0],"con") == "cat":
            new_metadata.append(row)
            
    return new_metadata

def group_samples_by_metadata(metadata, data, samples):
    """ Return the samples grouped by the metadata. The data and metadata
    should have the same ordering of sample columns.
    
    Args:
        metadata (list): A single metadata list.
        data (list): A single data set.
        samples (list): The samples (corresponding to the columns in the 
            data/metadata).
    
    Return:
        data (list of lists): The data organized by the metadata groups.
    """
    
    # get the samples for each metadata type
    sorted_samples_grouped={}
    for name, type in zip(samples,metadata[1:]):
        sorted_samples_grouped[type]=sorted_samples_grouped.get(type,[])+[name]
        
    sorted_data_grouped={}
    for row in data:
        sorted_temp={}
        for data_point, type in zip(row, metadata[1:]):
            sorted_temp[type]=sorted_temp.get(type,[])+[data_point]
            
        for key, value in sorted_temp.items():
            if not key in sorted_data_grouped:
                sorted_data_grouped[key]=[]
            sorted_data_grouped[key].append(value)
    
    return sorted_data_grouped, sorted_samples_grouped
    
def merge_metadata(metadata, samples, values, values_without_names=None):
    """ Merge the metadata and values into a single set. Samples are columns. 

    Args:
        metadata (lists of lists): A list of metadata. 
        samples (list): A list of samples that correspond with value columns.
        values (lists of lists): A list of values to merge with the metadata.
        values_without_names (bool): Set if the values do not have row names.
    Returns:
        (list): A list of lists of the merged data.
        (list): A list of the samples (might be a subset based on metadata available).
    """
    
    # get the indexes for the samples in the data that match with the metadata
    sample_index=[]
    metadata_index=[]
    samples_found=[]
    for index, name in enumerate(metadata[0][1:]):
        if name in samples:
            samples_found.append(name)
            sample_index.append(samples.index(name))
            metadata_index.append(index)
            
    # warn if no metadata samples match the values samples
    if len(sample_index) == 0:
        print("Warning: Metadata does not match samples.")
        return values, samples
    
    if len(samples) > len(metadata[0][1:]):
        print("Warning: Metadata only provided for a subset of samples.")
    
    # add metadata to the new data
    new_data=[]
    for row in metadata[1:]:
        # add only matching samples
        new_data.append([row[0]]+[row[i+1] for i in metadata_index])
        
    # add abundance values to the new data
    for row in values:
        if values_without_names:
            new_data.append([row[i] for i in sample_index])
        else:
            new_data.append([row[0]]+[row[i+1] for i in sample_index])
        
    return new_data, samples_found

def download_file(url, download_file):
    """
    Download a file from a url
    Create folder for downloaded file if it does not exist
    """

    create_folders(os.path.dirname(download_file))

    try:
        print("Downloading "+url)
        file, headers = urlretrieve(url,download_file,reporthook=ReportHook().report)
        # print final return to start new line of stdout
        print("\n")
    except EnvironmentError:
        print("WARNING: Unable to download "+url)

def try_log10(value):
    """ Try to convert value to log10 """
    
    try:
        new_value = math.log10(value)
    except ValueError:
        new_value = 0
        
    return new_value

def name_task(sample,software):
    """ Name the task based on the sample name and software """
    
    return software+"____"+os.path.basename(sample)

def add_to_list(items,new_item):
    """ Add the value to the list/tuple. If the item is not a list, create a new
        list from the item and the value 
        
    Args:
        items (list, string or tuple): Single or multiple items
        new_item (string): The new value
        
    Returns:
        (list): A list of all values
    """
    
    if isinstance(items,tuple):
        items=[i for i in items]
    
    if not isinstance(items,list):
        items=[items]
        
    return items+[new_item]

def metacyc_url(pathway):
    """ Return the url for the pathway on the MetaCyc website 
    
    Args:
        pathway (string): The MetaCyc pathway
        
    Returns
       (string): The url to the website for the pathway
       
    """
    
    return "http://metacyc.org/META/NEW-IMAGE?type=NIL&object="+pathway

def run_task(command, **keywords):
    """ Run the task command, formatting command with keywords. The command stdout
        and stderr are written to the workflow log.
    
    Args:
        command (string): A string to execute on the command line. It can be
            formatted the same as a task command.
       
    Returns:
        (int): Return code from command.     
    """

    from anadama2.helpers import format_command
    from anadama2.helpers import sh
    
    # format the command to include the items for this task
    command=format_command(command, **keywords)
    
    # run the command
    return_code = sh(command)()
    
    return return_code

def partial_function(function, **keywords):
    """ Return a partial function, setting function name attribute
    
    Args:
        function (function): A function
        keywords: One or more keywords to be applied to the function
        
    Returns:
        (function): A partial function
        
    """
    
    partial = functools.partial(function, **keywords)
    partial.__name__ = function.__name__
    
    return partial
        

def paired_files(files, extension, pair_identifier=None):
    """ Find sets of paired-end reads
    
    This function will find sets of paired end reads from a list of files.
    
    Args:
        files (list): A list of files (with or without the full paths)
        extension (string): The extension for all files.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
            
    Requires:
        None
        
    Returns:
        list: A list of paired files.
        
    Example:
        paired_set = paired_reads(["1.R1.fq", "1.R2.fq"],".fq")
        
    """
    
    # add period to extension if not included
    if not extension.startswith("."):
        extension="."+extension
    
    if pair_identifier is None:
        pair_identifier=".R1"
    
    # check for the one in the pair identifier
    if not "1" in pair_identifier:
        sys.exit("Please provide the identifier for the first pair set (ie R1).")
    
    pair_identifier2=pair_identifier.replace("1","2",1)

    input_pair1 = list(filter(lambda file: os.path.basename(file).replace(extension,"").endswith(pair_identifier), files))
    input_pair2 = list(filter(lambda file: os.path.basename(file).replace(extension,"").endswith(pair_identifier2), files))
    
    # only return matching pairs of files in the same order
    paired_file_set = [[],[]]
    for file1 in sorted(input_pair1):
        # find the matching file in the second set
        name1=sample_names(file1, extension, pair_identifier)
        for file2 in input_pair2:
            name2=sample_names(file2, extension, pair_identifier2)
            if name1 and name1 == name2:
                paired_file_set[0].append(file1)
                paired_file_set[1].append(file2)
                input_pair2.remove(file2)
                break
    
    return paired_file_set

def sample_names(files,extension,pair_identifier=None):
    """ Return the basenames of the files, without any extensions, as the sample names
    
    Args:
        files (list): A list of files (with or without the full paths)
        extension (string): The extension for all files.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        
    Requires:
        None
        
    Returns:
        list: A list of sample names (file basenames)

    Example:
        names = sample_names(["1.R1.fq", "1.R2.fq"],".fq")
        
    """
    
    # add period to extension if not included
    if not extension.startswith("."):
        extension="."+extension
    
    # if files is a string, convert to a list
    convert=False
    if isinstance(files,str):
        files=[files]
        convert=True
        
    samples=[os.path.basename(file).replace(extension,"") for file in files]
    
    # remove the pair_idenifier from the sample name, if provided
    if pair_identifier:
        # only remove the last instance of the pair identifier
        samples=[pair_identifier.join(sample.split(pair_identifier)[:-1]) if pair_identifier in sample else sample for sample in samples]
    
    if convert:
        samples=samples[0]
    
    return samples

def find_files(folder, extension=None, exit_if_not_found=None):
    """ Return the files in the given folder with the extension if provided
    
    Args:
        folder (string): A path to a folder
        extension (string): The file extension to search for (optional)
        exit_if_not_found (bool): Indicator to check if files exist (optional) 
        
    Requires:
        None
        
    Returns:
        list: A list of files in the folder

    Example:
        files = find_files("examples","fastq")
    """
    
    # get all of the files in the folder
    files=[os.path.join(folder,file) for file in os.listdir(folder)]
    files=list(filter(lambda file: os.path.isfile(file),files))
    
    # filter to only files with extension
    if extension:
        files=list(filter(lambda file: file.endswith(extension), files))
    
    if exit_if_not_found:
        if not files:
            message="ERROR: No files were found in the folder " + folder
            if extension:
                message+=" with extension "+extension
            sys.exit(message+" .\n")
            
    return files

def name_files(names, folder, subfolder=None, tag=None, extension=None, create_folder=None):
    """ Return a list of file names based on the names and folders provided
    
    Args:
        names (list or string): A list of basenames or files.
        folder (string): The path to the folder.
        subfolder (string): The subfolder to use with the files (optional).
        tag (string): The tag to add to the file basenames (optional).
        extension (string): The extension to use for the files (optional).
        create_folder (bool): Create the folder and subfolder if they do not exist (optional).

    Requires:
        None
        
    Returns:
        list: A list of file names (or string if input is string).
        
    Example:
        files = name_files(["file1","file2"], "output")
    """
    
    # if names is a list, convert to string
    was_string=False
    if isinstance(names, str):
        was_string=True
        names=[names]
    
    # get the basenames from the files
    names=[os.path.basename(name) for name in names]
    
    # use the full path to the folder
    folder=os.path.abspath(folder)
    
    # get the name of the full folder plus subfolder if provided
    if subfolder:
        folder=os.path.join(folder,subfolder)

    # add the extension if provided, and replace existing
    if extension:
        names=[os.path.splitext(name)[0]+"."+extension for name in names]

    # add the tag to the names, if provided
    if tag:
        names=[os.path.splitext(name)[0]+"_"+tag+os.path.splitext(name)[1] for name in names]
        
    files=[os.path.join(folder,name) for name in names]
    
    if create_folder:
        create_folders(os.path.dirname(files[0]))
        
    # if the input was originally a string, convert from list
    if was_string:
        files=files[0]
        
    return files

def create_folders(folder):
    """ Create folder if it does not exist
    
    Args:
        folder (string): The full path to the folder.
        
    Requires:
        None
        
    Returns:
        None
        
    Example:
        create_folders("new_folder")
    """
    
    try:
        if not os.path.exists(folder):
            os.makedirs(folder)
    except EnvironmentError:
        print("Warning: Unable to create folder: "+ folder)
        
def match_files(files1,files2,mapping):
    """ Match files from two sets using the mapping provided
    
    Args:
        files1 (list): A list of files for set1
        files2 (list): A list of files for set2
        mapping (string): The file with the mapping information. This file
            should be tab delimited. It can have headers starting with "#". 
            It should have the basenames for the files for set1 and set2 with
            each line as "fileA\tfileB" with fileA in set files1 and fileB in
            set files2.
            
    Requires:
        None
    
    Returns:
        (list): An ordered list of the first set of files
        (list): An ordered list of the second set of files
        The two lists will be the same length with pairs having the same index
            in each list.
        
    Example:
        match_files(["wts_1.fastq","wts_2.fastq"],["wms_1.tsv","wms_2.tsv"],"mapping.tsv")
        
        mapping.tsv contains:
        # wts   wms
        wts_1   wms_1
        wts_2   wms_2

    """
    
    # read in the mapping file
    set_mappings={}
    try:
        file_handle=open(mapping,"r")
        lines=file_handle.readlines()
        file_handle.close()                
    except EnvironmentError:
        sys.exit("ERROR: Unable to read mapping file: " + mapping)
    
    for line in lines:
        if not line.startswith("#"):
            data=line.rstrip().split("\t")
            if len(data) > 1:
                item1=data[0]
                item2=data[1]
                if "." in item1 or "." in item2:
                    sys.exit("ERROR: Sample names should not contain file extensions ('.fastq'),"+
                    "pair identifiers ('.R1.'), or periods ('.') as part of the name.")
                # check for duplicate mappings
                if item1 in set_mappings:
                    print("Warning: Duplicate mapping in file: " + item1)
                set_mappings[item1]=item2
    
    pair1=[]
    pair2=[]
    for item1,item2 in set_mappings.items():
        file1=list(filter(lambda file: os.path.basename(file).startswith(item1),files1))
        file2=list(filter(lambda file: os.path.basename(file).startswith(item2), files2))
        if len(file1) == 1 and len(file2) == 1:
            # check for the pair
            pair1.append(file1[0])
            pair2.append(file2[0])
        elif len(file1) == 0:
            print("Warning: Unable to find file with key, " + item1 + " in folder " + os.path.dirname(files1))
        elif len(file2) == 0:
            print("Warning: Unable to find file with key, " + item2 + " in folder " + os.path.dirname(files2))
        else:
            print("Warning: Duplicate files found for mapping keys: " + item1 + " " + item2)

    if len(pair1) != len(files1):
        print("Warning: Unable to find matches for all of the files in the set.")
    
    return pair1, pair2

def row_average(data):
    """ Compute the average of each row in a data set
        
        Args:
            data (list of lists): Each list in data represents a row of data. 
                
        Requires:
            None
        
        Returns:
            (list): A list of averages, one for each row in the original data.
            
        Example:
            row_average([[1,2,3],[4,5,6]])
    """    
    
    return [sum(row)/(len(row)*1.0) for row in data]

def row_variance(data):
    """ Compute the variance of each row in a data set
        
        Args:
            data (list of lists): Each list in data represents a row of data. 
                
        Requires:
            None
        
        Returns:
            (list): A list of variances, one for each row in the original data.
            
        Example:
            row_variance([[1,2,3],[4,5,6]])
    """
      
    data_averages=row_average(data)
    data_variances=[]
    for average, row in zip(data_averages,data):
        data_variances.append(sum((i-average)**2 for i in row)/(len(row)*1.0))
        
    return data_variances

def relative_abundance(data, percent=False):
    """ Compute the relative abundance values for a set of data 
        
        Args:
            data (list of lists): Each list in data represents a row of data. 
            percent (bool): Abundance is a percent of 100 (30 for 30% instead of 0.3)
                
        Requires:
            None
        
        Returns:
            (list of lists): Each list in data represents a row of data with relative abundance values.
            
        Example:
            relative_abundance([[1,2,3],[4,5,6]])   
    """ 

    # compute the sum for each column
    sums=[0.0]*len(data[0])
    for i in range(len(data[0])):
        for row in data:
            sums[i]+=float(row[i])
            
    relab=[]
    for row in data:
        new_row=[]
        for i, value in enumerate(row):
            try:
                new_value=value/sums[i]
            except ZeroDivisionError:
                new_value=0
            new_row.append(new_value)
        if percent:
            new_row = list(map(lambda x: x * 100.0, new_row))
        relab.append(new_row)
        
    return relab

def filter_zero_rows(taxa, data, ignore_index=None):
    """ Remove any taxa and data rows from the lists if the data sum for a row is zero.
        
        Args:
            taxa (list): The list of taxa.
            data (list of lists): Each list in data represents a row of data.
            ignore_index (int): An index to ignore in each row in computing the sum.
                
        Requires:
            None
        
        Returns:
            (list): A list of labels for the non-zero rows.
            (list of lists): Each list in data represents a row of data that is non-zero.  
    """ 
    new_taxa=[]
    new_data=[]
    for taxon, row in zip(taxa, data):
        if ignore_index is not None:
            temp_row = row
            del temp_row[ignore_index]
            row_sum = sum(temp_row)
        else:
            row_sum = sum(row)
        if row_sum != 0:
            new_taxa.append(taxon)
            new_data.append(row)
            
    return new_taxa, new_data

def taxa_shorten_name(taxa, level, remove_identifier=None):
    """ Shorten the taxa name by removing the levels indicated (useful for plotting)
        
        Args:
            taxa (list): The list of taxa.
            level (int): The level to filter.
            remove_identifier (bool): If set remove the [k|p|c|r|f|g|s|t__]) from the name. 
                
        Requires:
            None
        
        Returns:
            (list): The list of taxa after removing the unclassified names.  
    """ 

    new_names=[]
    for taxon in taxa:
        name=taxon.split(";")[level]
        if remove_identifier:
            name=name.split("__")[-1]
        new_names.append(name)
        
    return new_names

def top_rows(row_labels, data, max_sets, function):
    """ Get the top rows in the data based on the metric provided 
        
        Args:
            row_labels (list): A list of labels for each row.
            data (list of lists): Each list in data represents a row of data. 
            max_sets (int): Total number of top rows to return.
            function (string): The function to run to get the top values (average or variance)
                
        Requires:
            None
        
        Returns:
            (list): A list of labels for the top rows.
            (list of lists): Each list in data represents a row of data for the top data.
            
        Example:
            top_rows(["row1","row2"],[[1,2,3],[4,5,6]],1)
    """
    
    # get the data after applying the metric function
    if function == "variance":
        stats_data=row_variance(data)
    else:
        stats_data=row_average(data)
    
    # sort the numbers by decreasing order
    sorted_indexes=sorted(range(len(stats_data)),key=lambda i: stats_data[i],reverse=True)
    
    # reduce max sets if the number of rows is less than max sets
    if len(row_labels) < max_sets:
        max_sets = len(row_labels)
    
    top_labels=[]
    top_data=[]
    for i in range(max_sets):
        top_labels.append(row_labels[sorted_indexes[i]])
        top_data.append(data[sorted_indexes[i]])
        
    return top_labels, top_data

def remove_stratified_pathways(pathways, data, remove_description=None):
    """ Remove the stratified pathways from the data set.
        Also remove the unintegrated and unmapped values. Remove the descriptions
        from the pathway names if set.
    
        Args:
            pathways (list): A list of pathway names for each row.
            data (list of lists): Each list in data represents a row of data. 
            remove_description (bool): If set, remove the pathway description
                from the names returned.
                
        Requires:
            None
        
        Returns:
            (list): A list of pathway names.
            (list): A list of lists of the data.
            
        Example:
            remove_stratified_pathways(["pwy1","pwy1|bug1"],[[1,2,3],[4,5,6]])
    """
    
    new_pathways=[]
    new_data=[]
    
    for path, row in zip(pathways, data):
        if not "|" in path and not "UNINTEGRATED" in path and not "UNMAPPED" in path:
            if remove_description:
                path=path.split(":")[0]
            new_pathways.append(path)
            new_data.append(row)
            
    return new_pathways, new_data 

def pathway_names(pathways):
    """ Split the pathway names and descriptions
    
        Args:
            pathways (list): A list of pathway names and descriptions 
                (from a pathway abundance or coverage file)  
        
        Requires:
            None
            
        Returns:
            (dict): A dictionary of pathway names to descriptions
            
    """
    
    path_names = {}
    for path in pathways:
        # ignore stratified pathways
        if not "|" in path:
            try:
                description = path.split(":")
                name = description.pop(0)
                description=":".join(description)
            except ValueError:
                continue
            
            path_names[name]=description
        
    return path_names
        

def filter_taxa_abundance(taxonomy, data, min_abundance, min_samples):
    """ Remove the taxons by min abundance and min samples.
    
        Args:
            taxonomy (list): A list of taxonomy strings for each row.
            data (list of lists): Each list in data represents a row of data. 
            min_abundance (float): Remove data without min abundance. 
            min_samples (float): Remove data not in min samples.
                
        Requires:
            None
        
        Returns:
            (list): A list of species names.
            (list): A list of lists of the data.
            
        Example:
            filter_taxa_abundance(["g__ABC","s__DEF"],[[1,2,3],[4,5,6]],10,2)
    """ 

    filtered_data=[]
    filtered_taxonomy=[]
    # compute the min samples required for this data set
    min_samples_required=math.ceil(len(data[0])*(min_samples/100.0))
    for taxon, data_row in zip(taxonomy, data):
        # filter the species to only include those with min abundance in min of samples
        total_samples_pass_filter=len(list(filter(lambda x: x>min_abundance, data_row)))
        if total_samples_pass_filter >= min_samples_required: 
            filtered_taxonomy.append(taxon)
            filtered_data.append(data_row)
    
    return filtered_taxonomy, filtered_data

def filter_taxa_level_metaphlan_format(taxonomy, data, min_abundance=None, min_samples=None, level=6):
    """ Remove the taxons that are not a species level (or set a different level with keyword) from the data set.
        Also filter the species if filters are provided. Metaphlan2 format with "|" delimiters and tiered
        abundances (so genus level is split and repeated stratified by species).
    
        Args:
            taxonomy (list): A list of taxonomy strings for each row.
            data (list of lists): Each list in data represents a row of data. 
            min_abundance (float): If set, remove data without min abundance. To
                be used with min_samples.
            min_samples (float): If set, remove data not in min samples.
            level (int): Taxonomic level (default set to species)
                
        Requires:
            None
        
        Returns:
            (list): A list of species names.
            (list): A list of lists of the data.
            
        Example:
            filter_taxa_level_metaphlan_format(["g__ABC","s__DEF"],[[1,2,3],[4,5,6]])
    """

    taxonomic_levels=["|k__","|p__","|c__","|o__","|f__","|g__","|s__","|t__"]

    # get the level indicator to search for and the next level to filter
    search_taxa = taxonomic_levels[level]
    remove_taxa = taxonomic_levels[level+1] if level+1 < len(taxonomic_levels) else "\n"

    species_data=[]
    species_taxonomy=[]
    # identify the species data in the data set
    # filter out those with species and strain information
    for taxon, data_row in zip(taxonomy, data):
        if search_taxa in taxon and not remove_taxa in taxon:
            species_taxonomy.append(taxon.split("|")[-1].replace(search_taxa[1:],"").replace("_"," "))
            species_data.append(data_row)

    # if filters are provided, then filter the data by both min abundance
    # and min samples
    if min_abundance is not None and min_samples is not None:
        species_taxonomy, species_data = filter_taxa_abundance(species_taxonomy, species_data, min_abundance, min_samples)

    return species_taxonomy, species_data

def read_otu_table(file):
    """ Read in an otu table. Remove extra brackets from taxonomy names if present.
    
        Args:
            file (string): A file containing the otu table (tsv format).
                
        Requires:
            None
        
        Returns:
            (list): A list of samples.
            (list): A list of otu ids.
            (list): A list of taxons.
            (list): A list of lists of data.
            
        Example:
            samples, ids, taxonomy, data = read_otu_table("otu_table.tsv")
    """
        
    data=[]
    samples=[]
    taxonomy=[]
    ids=[]
    with open(file) as file_handle:
        samples = file_handle.readline().rstrip().split("\t")[1:-1]
        for line in file_handle:
            data_points=line.rstrip().split("\t")
            ids.append(data_points.pop(0))
            taxonomy.append(data_points.pop().replace("[","").replace("]",""))
            data.append([float(i) for i in data_points])
            
    return samples, ids, taxonomy, data
    

def sort_data(data, samples, sort_by_name=False, sort_by_name_inverse=False):
    """ Sort the data with those with the largest values first or by sample name

        Args:
            data (list): The data points for each sample.
            samples (list): The sample names that correspond to each data point. 
            sort_by_name (bool): If true, sort by sample name
            sort_by_name_inverse (bool): If true, sort by the inverse of the name (so the reverse of the string)
                this is useful for samples with sample name plus features

        Requires:
            None
        
        Returns:
            (list): The data points for each sample sorted.
            (list): The sample names that correspond to each data point sorted.
            
    """
    import numpy
    
    # if the data is a list of lists of single values, then convert to a list of values
    if isinstance(data[0], list) and max([len(row) for row in data]) == 1:
        data=[row[0] for row in data]
       
    # if set, sort by sample name (samples as columns in the data)
    if sort_by_name or sort_by_name_inverse:
        data_by_sample={sample:data_point for sample,data_point in zip(samples,numpy.transpose(data))}
        sorted_samples=sorted(samples)
        if sort_by_name_inverse:
            # use the reverse of the sample names to sort
            sorted_samples=sorted(samples, key=lambda x: x[::-1])
        sorted_data_transpose=[data_by_sample[sample] for sample in sorted_samples]
        sorted_data = numpy.transpose(sorted_data_transpose)
    else:
        data_by_sample={sample:data_point for sample,data_point in zip(samples,data)}
        sorted_samples=sorted(data_by_sample,key=data_by_sample.get, reverse=True)
        sorted_data=[data_by_sample[sample] for sample in sorted_samples]
   
    return sorted_samples, sorted_data

def is_paired_table(file):
    """ Check if a file contains paired read counts using the header information.
    
        Args:
            file (string): A file of read counts.
                
        Requires:
            None
        
        Returns:
            (bool): True if the file contains paired read counts
            
    """
    
    # read in the first line in the file
    try:
        with open(file) as file_handle:
            header=file_handle.readline()
    except EnvironmentError:
        sys.exit("Unable to read file: " + file)
        
    paired = True if "pair" in header.lower() else False
    
    return paired

def microbial_read_proportion_multiple_databases(data, columns, orphan_data=None, rna=None):
    """ Compute microbial read proportions from the KneadData read counts for 
        multiple databases 
        
        Args:
            data (list of lists): The single or paired data read counts for each sample.
            columns (list): The names of the columns corresponding to the paired data. These
                columns include the reference database names.
            orphan_data (list of lists): The orphan data (if paired end reads)
            rna (bool): If set, this data set is RNA so compute additional ratio.
        
        Requires:
            None
        
        Returns:
            (list of lists): A list of ratios for each sample.
            (list): A list of strings with labels for the ratios.
    """
    
    # compute ratios for each database used for qc
    dna_microbial_reads=[]
    dna_microbial_labels=[]
    for index, qc_database in enumerate(columns[2:]):
        # get a subset of the data for this ratio
        data_subset=[row[:2]+[row[index+2]] for row in data]

        # create subset of orphan data if provided
        orphan_subset=None
        if orphan_data and len(orphan_data[0]) > 2:
            orphan_subset=[row[:2]+[row[index+2]] for row in orphan_data]
        else:
            orphan_subset=orphan_data
        
        reads_ratio, ratio_labels = microbial_read_proportion(data_subset, 
            orphan_data=orphan_subset, database_name=qc_database, rna=rna)
        dna_microbial_labels+=ratio_labels
        if not dna_microbial_reads:
            dna_microbial_reads=reads_ratio
        else:
            dna_microbial_reads=[row1+row2 for row1, row2 in zip(dna_microbial_reads,reads_ratio)]
            
    return dna_microbial_reads, dna_microbial_labels

def microbial_read_proportion(paired_data, orphan_data=None, rna=None, database_name=None):
    """ Compute microbial read proporations from the KneadData read counts.
    
        Args:
            paired_data (list of lists): The paired data read counts for each sample.
            orphan_data (list of lists): The orphan data read counts for each sample. 
            rna (bool): If set, this data set is RNA so compute additional ratio.
            database_name (string): The name of the contaminate database.
                
        Requires:
            None
        
        Returns:
            (list of lists): A list of ratios for each sample.
            (list): A list of strings with labels for the ratios.
            
    """
    
    # if the database name is not set, use the default
    if database_name is None:
        database_name="hg38"
    
    # if the orphan reads are not provided, create an empty set of data
    if orphan_data is None:
        orphan_data=[]
        for i in range(len(paired_data)):
            orphan_data.append([0,0,0,0])
    
    proportion_decontaminated = []
    for paired_row, orphan_row in zip(paired_data, orphan_data):
        decontaminated_sum = 2.0 * paired_row[-1] + orphan_row[-1] + orphan_row[-2]
        decon_trim = decontaminated_sum / (2.0 * paired_row[1] + orphan_row[0] + orphan_row[1])
        decon_raw = decontaminated_sum / (2.0 * paired_row[0])
        if rna:
            decon_ratio = decontaminated_sum / (2.0 * paired_row[-2] + orphan_row[-3] + orphan_row[-4])
            proportion_decontaminated.append(["{0:.5f}".format(i) for i in [decon_trim, decon_ratio, decon_raw]])
        else:
            proportion_decontaminated.append(["{0:.5f}".format(i) for i in [decon_trim, decon_raw]])

    if rna:
        labels=[database_name+" mRNA / Trim",database_name+" mRNA / "+database_name,database_name+" mRNA / Raw"]
    else:
        labels=[database_name+" / Trim",database_name+" / Raw"]
        
    return proportion_decontaminated, labels


def taxa_remove_unclassified(taxa, delimiter=";"):
    """ Rename the taxa to remove the unclassified levels
    
        Args:
            taxa (list): The list of taxa.
            delimiter (str): The string delimiter (usually pipe or semi-colon).
                
        Requires:
            None
        
        Returns:
            (list): The list of taxa after removing the unclassified names.
    """
    
    # remove any levels where the name is unknown (ie empty)
    for taxon in taxa:
        new_name=[]
        for level in taxon.replace(" ","").split(delimiter):
            try:
                rank, name = level.split("__")
            except ValueError:
                # ignore identities like "unclassified" if present
                continue
            if name:
                new_name.append(level)
            else:
                break
        yield delimiter.join(new_name)
        
def taxonomy_trim(taxa):
    """ Trim the taxonomy name to include the most specific known name followed by unclassified
    
        Args:
            taxa (list): The list of taxa.
                
        Requires:
            None
        
        Returns:
            (list): The list of taxa after trimming.
    """
    
    # remove any spaces from the taxonomy
    taxa = [taxon.replace(" ","") for taxon in taxa]

    # determine the delimiter (pipe or semi-colon)
    delimiter="|" if "|" in taxa[0] else ";"
    
    # get the taxa with unclassified levels removed
    taxa_unclassified_removed = taxa_remove_unclassified(taxa,delimiter)
    
    trimmed_taxa=[]
    for taxon_full, taxon_reduced in zip(taxa, taxa_unclassified_removed):
        # if the taxon is specific to species level, then 
        # return the genus and species level
        if taxon_full == taxon_reduced:
            data = taxon_full.split(delimiter)
            trimmed_taxa.append(data[-2]+"."+data[-1])
            
        else:
            most_specific_clade = taxon_reduced.split(delimiter)[-1]
            if not most_specific_clade:
                trimmed_taxa.append(taxon_reduced.replace(delimiter,".")) 
            else:
                data = taxon_full.split(most_specific_clade)
                trimmed_taxa.append(most_specific_clade+data[-1].replace(delimiter,"."))
            
    return trimmed_taxa
        
def terminal_taxa(taxa, data):
    """ Reduce the list of taxa to just those that represent the terminal nodes. If there
        are duplicate terminal nodes, then sum the duplicates.
    
        Args:
            taxa (list): The list of taxa (in the same order as the data).
            data (list of lists): The data points for all samples for each taxa.
                
        Requires:
            None
        
        Returns:
            (list): The list of taxa (terminal node only)
            (list of lists): The data after reducing to terminal node taxa.
    """    
    
    terminal_node_taxa=[]
    # check the taxa by level, starting with the most specific level of strain
    # use a full match with strain instead of just startswith to allow for unclassified
    # strains to not match with classified strains
    # if strains are not present, then run at a species level instead
    
    # check for the most specific taxonomy level (ie strain or species)
    max_taxonomy_level=max([len(taxon.split(";")) for taxon in taxa])
    
    taxa_for_level, data_level=taxa_by_level(taxa, data, level=max_taxonomy_level-1, keep_unclassified=True)
    for taxon in taxa_for_level:
        matching_taxa=list(filter(lambda x: x.replace(" ","") == taxon.replace(" ",""), terminal_node_taxa))
        if len(matching_taxa) == 0:
            terminal_node_taxa.append(taxon)
    
    for level in reversed(range(max_taxonomy_level-1)):
        taxa_for_level, data_level=taxa_by_level(taxa, data, level, keep_unclassified=True)
        for taxon in taxa_for_level:
            # check if part of this taxon is already included
            matching_taxa=list(filter(lambda x: x.replace(" ","").startswith(taxon.replace(" ","")), terminal_node_taxa))
            if len(matching_taxa) == 0:
                terminal_node_taxa.append(taxon)
                
    # create a set of terminal node taxa and data
    new_taxa={}
    for taxon, row in zip(taxa, data):
        if taxon in terminal_node_taxa:
            if taxon in new_taxa:
                new_taxa[taxon]=[a+b for a,b in zip(new_taxa[taxon],row)]
            else:
                new_taxa[taxon]=row
               
    new_taxa_list=sorted(new_taxa.keys())
    new_data_list=[new_taxa[i] for i in new_taxa_list] 
               
    return new_taxa_list, new_data_list
                    
def taxa_by_level(taxa, data, level, keep_unclassified=None):
    """ Combine the data to represent the taxa by a specific level
    
        Args:
            taxa (list): The list of taxa (in the same order as the data).
            data (list of lists): The data points for all samples for each taxa.
            level (int): The level to sum the taxa (zero is kingdom level).
            keep_unclassified (bool): If set, keep unclassified taxa.
                
        Requires:
            None
        
        Returns:
            (list): The list of taxa (all to the level specified)
            (list of lists): The data after summing to the taxa level specified.
    """    

    # first remove any unclassified levels
    if not keep_unclassified:
        taxa=taxa_remove_unclassified(taxa)
    
    # sum the taxa by the level provided
    data_sum={}
    for taxon, taxon_data in zip(taxa, data):
        split_taxon=taxon.split(";")
        if len(split_taxon) < (level+1):
            # do not include those taxa that are not specified to the level requested
            continue
        new_taxon_level=";".join(split_taxon[:(level+1)])
        if new_taxon_level in data_sum:
            data_sum[new_taxon_level]=[a+b for a,b in zip(data_sum[new_taxon_level],taxon_data)]
        else:
            data_sum[new_taxon_level]=taxon_data
        
    new_taxa=[]
    new_data=[]    
    for taxon, taxon_data in data_sum.items():
        new_taxa.append(taxon)
        new_data.append(taxon_data)
        
    return new_taxa, new_data

def format_data_comma(data):
    """ Format the numbers in the string to include commas.
    
    Args:
        data (string or list): A text string.
        
    Requires:
        None
        
    Returns:
        (string): A text string.
        
    """
    
    if not isinstance(data,list):
        data=data.split()

    new_string=[]
    for token in data:
        try:
            new_token="{:,}".format(int(token))
        except ValueError:
            new_token=token
        new_string.append(new_token)
        
    return " ".join(new_string)

def read_eestats2(file):
    """ Read the eestats2 file which is an ascii table.
    
    Args:
        file (string): The path to the eestats file.
        
    Requires:
        None
        
    Returns:
        (list): The table rows.
        (list): The table columns.
        (list): The table data.
        (string): The summary.
        
    """
    
    with open(file) as file_handle:
        eestats_lines = file_handle.readlines()
        
    # read in the overall stats line
    overall_stats = format_data_comma(eestats_lines[1].rstrip())
    # read in the maxee values from the columns
    columns = list(filter(lambda x: x.strip() and not x in ["Length","MaxEE"],eestats_lines[3].rstrip().split()))
    columns = [column+" maxee" for column in columns]
    rows = []
    data = []
    # read through the data table
    for line in eestats_lines[5:]:
        stats = list(filter(lambda x: x.strip(),line.strip().split("   ")))
        # move spaces in data values and percents
        stats = [stat.replace("(  ","(").replace("( ","(").replace("("," (") for stat in stats]
        rows.append(stats.pop(0)+" nt")
        data.append([format_data_comma(stat) for stat in stats])
        
    return rows, columns, data, overall_stats

def get_files(folder, extension):
    """ Return paths to all files in a folder with a given extension """
    
    for file in os.listdir(folder):
        file=os.path.join(folder, file)
        if os.path.isfile(file) and file.endswith(extension):
            yield file

def read_picard(file, threshold=20):
    """ Read the picard file which is an ascii table of quality scores per base.
    
    Args:
        file (string): The path to the picard file.
        
    Requires:
        None
        
    Returns:
        (list): A list of tuples of base / quality score. 
        (bool): True if any quality score in the set does not meet the threshold.
    """
    
    with open(file) as file_handle:
        picard_lines = file_handle.readlines()
        
    # read through file to get quality scores ignoring all lines with comments/headers
    data=[]
    below_threshold=False
    for line in picard_lines:
        if not (line.startswith("#") or line.startswith("CYCLE")):
            new_data=line.rstrip().split("\t")
            # check if the quality score passes the threshold
            try:
                new_data=(int(new_data[0]),float(new_data[1]))
                if new_data[1] < threshold:
                    below_threshold=True
                data.append(new_data)
            except ValueError:
                pass
        
    return data, below_threshold

def rank_species_average_abundance(file, id_index=-1, only_species=True):
    """ Read in a taxonomy file, and sort species by average abundance 
    
    Args:
        file (string): The path to the merged taxonomy file generated by MetaPhlAn
        
    Requires:
        None
        
    Returns:
        (list): A list of species ordered by decreasing average abundance
    """
    
    def try_format_data(value):
        """ Try to format the data in a file """
        try:
            value = float(value)
        except ValueError:
            value= 0
            
        return value
    
    species=collections.OrderedDict()
    with open(file) as file_handle:
        try:
            column_names = file_handle.readline().rstrip().split("\t")[1:]
        except IndexError:
            column_names = []
        for line in file_handle:
            line=line.rstrip().split("\t")
            taxonomy=line.pop(0).split("|")[id_index]
            data=[try_format_data(i) for i in line]
            try:
                average=sum(data)/(len(data)*1.0)
            except ZeroDivisionError:
                average=0
            # only store values for species
            if only_species:
                if taxonomy.startswith("s__"):
                    species[taxonomy]=average
            else:
                species[taxonomy]=average
                
    # sort the species from highest to lowest average abundance
    sorted_species = sorted(species, key=species.get, reverse=True)
    # if abundances are not provided then use original ordering
    if sum(species.values()) == 0:
        sorted_species = list(species.keys())
    
    return sorted_species

def order_clade_list(task,clade_list,abundance_file,output_file):
    """ Using the average abundances order the clade list from strainphlan 
    
    Args:
        task (anadama2.task): An instance of the task class.
        clade_list (string): The path to the file containing the strainphlan list of clades.
        abundance_file (string): The path to the merged abundance file from metaphlan.
        output_file (string): The file to write the ordered clade list.
        
    Requires:
        none
        
    Returns:
        none
    """
    
    # get the species listed by average abundance
    species_ranked = rank_species_average_abundance(abundance_file)
    
    # read in the clade list
    clades=set()
    with open(clade_list) as file_handle:
        for line in file_handle:
            if "s__" in line:
                clades.add(line.strip().split("\t")[1].split(": in ")[0])
        
    # write out ordered species also included in clade list
    with open(output_file,"w") as file_handle:
        for taxon in species_ranked:
            if taxon in clades:
                file_handle.write(taxon+"\n")
    

def sort_fastq_file(task):
    """Sorts a FASTQ file by name (sequence identifier contents).

    Args:
        task (anadama2.task): An instance of the task class.

    Requires:
        None

    Returns:
        None
    """
    sample_name = os.path.basename(task.depends[0].name)
    output_dir = os.path.dirname(task.targets[0].name)
    temp_dir = os.path.join(output_dir, "%s.tmp" % sample_name)

    if task.depends[0].name.endswith('.gz'):
        sort_command = "zcat "
    else:
        sort_command = "cat "        
    
    sort_command += ("[depends[0]] | paste - - - - | sort -T [depends[1]] -k1,1 | "
                    "tr '\t' '\n' > [targets[0]]")

    run_task('mkdir -p [targets[0]]',
             depends=[TrackedDirectory(output_dir)],
             targets=[TrackedDirectory(temp_dir)])

    run_task(sort_command,
             depends=task.depends + [TrackedDirectory(temp_dir), TrackedDirectory(output_dir)],
             targets=task.targets)

    run_task("rm -rf [depends[0]]",
             depends=[TrackedDirectory(temp_dir)] + task.targets)


def extract_orphan_reads(task):
    """Extracts orphan reads from the provided input files. Orphan reads are saved into a separate file 
    for further downstream analysis.

    Args:
        task (anadama2.task): An instance of the task class.

    Requires:
        seqtk v1.2+: A fast and lightweight tool for processing sequences in the FASTA
             or FASTQ format

    Returns:
        None
    """
    sample_name = os.path.basename(os.path.dirname(task.depends[0].name))
    orphans_dir = os.path.dirname(task.targets[0].name)

    run_task("seqtk dropse [depends[0]] > [targets[0]]",
             depends=task.depends,
             targets=task.targets[0])

    run_task("extract_orphan_reads.py -r [depends[0]] -b [depends[1]] -o [depends[2]]",
             depends=[task.depends[0], task.targets[0], TrackedDirectory(orphans_dir)],
             targets=[task.targets[1]],
             args=[".fastq"])


def is_paired_end(input_files, extension, pair_identifier):
    """Returns true if the provided set of input files are paired-end.

    Args:
        input_files (list): A list of files to verify if paired-end data.
        extension (string): The extension for all files.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set.

    Requires:
        None

    Returns:
        bool: True if paired-end datasets; False otherwise.
    """
    input_pair1, input_pair2 = paired_files(input_files, extension, pair_identifier)
    
    return True if input_pair1 else False
