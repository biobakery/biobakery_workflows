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

def paired_files(files, pair_identifier=None):
    """ Find sets of paired-end reads
    
    This function will find sets of paired end reads from a list of files.
    
    Args:
        files (list): A list of files (with or without the full paths)
        pair_identifer (string): The string in the file basename to identify
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
    
    input_pair1 = list(filter(lambda file: pair_identifier in os.path.basename(file), files))
    input_pair2 = list(filter(lambda file: pair_identifier.replace("1","2") in os.path.basename(file), files))
    
    # only return matching pairs of files in the same order
    paired_files = [[],[]]
    for file1, file2 in zip(sorted(input_pair1), sorted(input_pair2)):
        if sample_names(file1) == sample_names(file2):
            paired_files[0].append(file1)
            paired_files[1].append(file2)
    
    return [input_pair1, input_pair2]

def sample_names(files):
    """ Return the basenames of the files, without any extensions, as the sample names
    
    Args:
        files (list): A list of files (with or without the full paths)
        
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
