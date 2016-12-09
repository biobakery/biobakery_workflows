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
