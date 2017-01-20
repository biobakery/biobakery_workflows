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

def paired_files(files, pair_identifier=None):
    """ Find sets of paired-end reads
    
    This function will find sets of paired end reads from a list of files.
    
    Args:
        files (list): A list of files (with or without the full paths)
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
            
    Requires:
        None
        
    Returns:
        list: A list of paired files.
        
    Example:
        paired_set = paired_reads(["1.R1.fq", "1.R2.fq"])
        
    """
    
    if pair_identifier is None:
        pair_identifier=".R1."
    
    pair_identifier2=pair_identifier.replace("1","2",1)

    input_pair1 = list(filter(lambda file: pair_identifier in os.path.basename(file), files))
    input_pair2 = list(filter(lambda file: pair_identifier2 in os.path.basename(file), files))
    
    # only return matching pairs of files in the same order
    paired_file_set = [[],[]]
    for file1, file2 in zip(sorted(input_pair1), sorted(input_pair2)):
        if sample_names(file1,pair_identifier) == sample_names(file2,pair_identifier2):
            paired_file_set[0].append(file1)
            paired_file_set[1].append(file2)
    
    return paired_file_set

def sample_names(files,pair_identifier=None):
    """ Return the basenames of the files, without any extensions, as the sample names
    
    Args:
        files (list): A list of files (with or without the full paths)
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set (optional).
        
    Requires:
        None
        
    Returns:
        list: A list of sample names (file basenames)

    Example:
        names = sample_names(["1.R1.fq", "1.R2.fq"])
        
    """
    
    # if files is a string, convert to a list
    convert=False
    if isinstance(files,basestring):
        files=[files]
        convert=True
        
    samples=[os.path.basename(file).split(".")[0] for file in files]
    
    # remove the pair_idenifier from the sample name, if provided
    if pair_identifier:
        samples=[sample.replace(pair_identifier,"") for sample in samples]
    
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
        list: A list of file names.
        
    Example:
        files = name_files(["file1","file2"], "output")
    """
    
    # if names is a list, convert to string
    if isinstance(names, basestring):
        names=[names]
    
    # get the basenames from the files
    names=[os.path.basename(name) for name in names]
    
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
            (list): A list of variances, one for each row in the original data.
            
        Example:
            row_variance([[1,2,3],[4,5,6]])
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

def filter_species(taxonomy, data, min_abundance=None, min_samples=None):
    """ Remove the taxons that are not a species level from the data set.
        Also filter the species if filters are provided.
    
        Args:
            taxonomy (list): A list of taxonomy strings for each row.
            data (list of lists): Each list in data represents a row of data. 
            min_abundance (float): If set, remove data without min abundance. To
                be used with min_samples.
            min_samples (float): If set, remove data not in min samples.
                
        Requires:
            None
        
        Returns:
            (list): A list of species names.
            (list): A list of lists of the data.
            
        Example:
            filter_species(["g__ABC","s__DEF"],[[1,2,3],[4,5,6]])
    """
    
    species_data=[]
    species_taxonomy=[]
    # identify the species data in the data set
    # filter out those with species and strain information
    for taxon, data_row in zip(taxonomy, data):
        if "|s__" in taxon and not "|t__" in taxon:
            species_taxonomy.append(taxon.split("|")[-1].replace("s__","").replace("_"," "))
            species_data.append(data_row)

    # if filters are provided, then filter the data by both min abundance
    # and min samples
    if min_abundance is not None and min_samples is not None:
        filtered_data=[]
        filtered_taxonomy=[]
        # compute the min samples required for this data set
        min_samples_required=math.ceil(len(species_data[0])*(min_samples/100.0))
        for taxon, data_row in zip(species_taxonomy, species_data):
            # filter the species to only include those with min abundance in min of samples
            total_samples_pass_filter=len(list(filter(lambda x: x>min_abundance, data_row)))
            if total_samples_pass_filter >= min_samples_required: 
                filtered_taxonomy.append(taxon)
                filtered_data.append(data_row)
        species_taxonomy=filtered_taxonomy
        species_data=filtered_data

    return species_taxonomy, species_data



def microbial_read_proportion(paired_data, orphan_data, rna=None):
    """ Compute microbial read proporations from the KneadData read counts.
    
        Args:
            paired_data (list of lists): The paired data read counts for each sample.
            orphan_data (list of lists): The orphan data read counts for each sample. 
            rna (bool): If set, this data set is RNA so compute additional ratio.
                
        Requires:
            None
        
        Returns:
            (list of lists): A list of ratios for each sample.
            (list): A list of strings with labels for the ratios.
            
    """
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
        labels=["hg38 mRNA / Trim","hg38 mRNA / hg38","hg38 mRNA / Raw"]
    else:
        labels=["hg38 / Trim","hg38 / Raw"]
        
    return proportion_decontaminated, labels


