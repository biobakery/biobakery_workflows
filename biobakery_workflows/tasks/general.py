"""
bioBakery Workflows: general tasks module
A collection of common tasks for workflows with shotgun sequences

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
import sys

from anadama2.tracked import TrackedExecutable, TrackedDirectory

from biobakery_workflows import utilities


def demultiplex(workflow, input_files, extension, output_folder, barcode_file, index_files, min_phred, pair_identifier):
    """Demultiplex the files (single end or paired)
    
    Args:
        workflow (anadama2.workflow): An instance of the workflow class.
        input_files (list): A list of paths to fastq files for input to ea-utils.
        extension (string): The extension for all files.
        output_folder (string): The path of the output folder.
        barcode_file (string): A file of barcodes.
        index_files (string): A list of paths to the index files.
        min_phred (int): The min phred quality score to use in the demultiplex command.
        pair_identifier (string): The string in the file basename to identify
            the first pair in the set.
        
    Requires:
        ea-utils fastq-multx: A tool to demultiplex fastq files.
        
    Returns:
        list: A list of the demultiplexed files
        string: output folder of demultiplexed files
        
    """
    
    # error if there is more than one index file
    if len(index_files) > 1:
        sys.exit("ERROR: Only one index file expected for demultiplexing step.")
    
    # read the barcode file to get the expected output files 
    try:
        file_handle=open(barcode_file)
        lines=file_handle.readlines()
        file_handle.close()
    except EnvironmentError:
        sys.exit("ERROR: Unable to read barcode file: " + barcode_file)
        
    samples=set()
    for line in lines:
        # ignore headers or comment lines
        if not line.startswith("#"):
            sample_name=line.rstrip().split("\t")[0]
            if sample_name:
                samples.add(sample_name)
            
    # get the names of the expected output files
    demultiplex_fastq_files = utilities.name_files(samples,output_folder,subfolder="demultiplex",extension="fastq")
    
    # name the barcode file with the reverse complement barcodes added
    expanded_barcode_file = utilities.name_files("expanded_barcode_file.txt",output_folder,subfolder="demultiplex",create_folder=True)
    
    # create a file that includes the reverse complements of the barcodes
    workflow.add_task(
        "reverse_compliment_barcodes.py --input [depends[0]] --output [targets[0]]",
        depends=barcode_file,
        targets=expanded_barcode_file)
    
    # check for paired input files
    input_pair1, input_pair2 = utilities.paired_files(input_files, extension, pair_identifier)
    
    # capture the demultiplex stats in output files, one for each set of input files
    if input_pair1:
        demultiplex_log = utilities.name_files(input_pair1[0],output_folder,subfolder="demultiplex",extension="log")
    else:
        demultiplex_log = utilities.name_files(input_files[0],output_folder,subfolder="demultiplex",extension="log")
        
    # get the output folder for all files
    demultiplex_output_folder = os.path.dirname(demultiplex_log)
    
    # get the basenames of the output files, one for each sample
    demultiplex_output_basenames = utilities.name_files(samples,output_folder,subfolder="demultiplex")
    
    # create a tracked executable
    fastq_multx_tracked = TrackedExecutable("fastq-multx",version_command="echo 'fastq-multx' `fastq-multx 2>&1 | grep Version`")
    
    if input_pair1 and input_pair2:
        # this run has paired input files
        # get the second pair identifier
        pair_identifier2=pair_identifier.replace("1","2",1)
        # get the names of the expected output files
        demultiplex_fastq_files_R1 = [file+pair_identifier+".fastq" for file in demultiplex_output_basenames]
        demultiplex_fastq_files_R2 = [file+pair_identifier2+".fastq" for file in demultiplex_output_basenames]
        demultiplex_fastq_files = demultiplex_fastq_files_R1+demultiplex_fastq_files_R2
        
        if index_files:
            # this run has index files
            workflow.add_task(
                "fastq-multx -l [depends[0]] [depends[1]] [depends[2]] [depends[3]] -o [args[1]]/%_I1_001.fastq [args[1]]/%[args[2]].fastq [args[1]]/%[args[3]].fastq -q [args[0]] > [targets[0]]",
                depends=[expanded_barcode_file, index_files[0], input_pair1[0], input_pair2[0], fastq_multx_tracked],
                args=[min_phred, demultiplex_output_folder, pair_identifier, pair_identifier2],
                targets=demultiplex_log,
                name="demultiplex")
            
        else:
            workflow.add_task(
                "fastq-multx -l [depends[0]] [depends[1]] [depends[2]] -o [args[1]]/%[args[2]].fastq [args[1]]/%[args[3]].fastq -q [args[0]] > [targets[0]]",
                depends=[expanded_barcode_file, input_pair1[0], input_pair2[0], fastq_multx_tracked],
                args=[min_phred, demultiplex_output_folder, pair_identifier, pair_identifier2],
                targets=demultiplex_log,
                name="demultiplex")
        
    else:
        # this run has single end input files
        # get the names of the expected output files
        demultiplex_fastq_files = [file+pair_identifier+".fastq" for file in demultiplex_output_basenames]
        
        if index_files:
            # this run has index files
            workflow.add_task(
                "fastq-multx -l [depends[0]] [depends[1]] [depends[2]] -o [args[1]]/%_I1_001.fastq [args[1]]/%[args[2]].fastq -q [args[0]] > [targets[0]]",
                depends=[expanded_barcode_file, index_files[0], input_files[0], fastq_multx_tracked],
                args=[min_phred, demultiplex_output_folder, pair_identifier],
                targets=demultiplex_log,
                name="demultiplex")
            
        else:
            workflow.add_task(
                "fastq-multx -l [depends[0]] [depends[1]] -o [args[1]]/%[args[2]].fastq -q [args[0]] > [targets[0]]",
                depends=[expanded_barcode_file, input_files[0]],
                args=[min_phred, demultiplex_output_folder, pair_identifier, fastq_multx_tracked],
                targets=demultiplex_log,
                name="demultiplex")

    demultiplex_fastq_files = demultiplex_check(workflow, demultiplex_log, demultiplex_fastq_files)


    return demultiplex_fastq_files, demultiplex_output_folder


def demultiplex_dual(workflow, output_folder, input_files, extension,
            barcode_files, dual_barcode_path, min_phred, pair_identifier):

    """Demultiplex the files (dual indexed paired)

        Args:
            workflow (anadama2.workflow): An instance of the workflow class.
            input_files (list): A list of paths to fastq(gz) files for input to ea-utils.
            extension (string): The extension for all files.
            output_folder (string): The path of the output folder.
            barcode_files (list): A list of barcode files.
            dual_index_path (string): A paths to the dual index file.
            min_phred (int): The min phred quality score to use in the demultiplex command.
            pair_identifier (string): The string in the file basename to identify
                the first pair in the set.

        Requires:
            ea-utils fastq-multx: A tool to demultiplex fastq files.

        Returns:
            list: A list of the demultiplexed files
            string: output folder of demultiplexed files

        """

    # capture the demultiplex stats in log file, one for each set of input files
    demultiplex_log = utilities.name_files(input_files[0],output_folder,subfolder="demultiplex",extension="log",create_folder=True)
    demultiplex_output_folder = os.path.dirname(demultiplex_log)

    # create a tracked executable
    fastq_multx_tracked = TrackedExecutable("fastq-multx",
                                            version_command="echo 'fastq-multx' `fastq-multx 2>&1 | grep Version`")

    # check for paired input files
    input_pair1, input_pair2 = utilities.paired_files(input_files, extension, pair_identifier)

    # get barcode files
    barcode1, barcode2 = utilities.paired_files(barcode_files, extension, pair_identifier)

    # get the second pair identifier
    pair_identifier2 = pair_identifier.replace("1", "2", 1)

    try:
        file_handle = open(dual_barcode_path)
        lines = file_handle.readlines()
        file_handle.close()
    except EnvironmentError:
        sys.exit("ERROR: Unable to read dual barcode file: " + dual_barcode_path)

    run_name = os.path.basename(input_pair1[0]).replace(pair_identifier, "").replace("." + extension, "")
    demultiplex_files = set()
    for line in lines:
        # ignore headers or comment lines
        if not line.startswith("#"):
            sample_name = line.split("\t")[0]

            if sample_name:
                nm1 = demultiplex_output_folder + "/" + run_name + "_" + sample_name + pair_identifier + "." + extension
                nm2 = demultiplex_output_folder + "/" + run_name + "_" + sample_name + pair_identifier2 + "." + extension
                demultiplex_files.add(nm1)
                demultiplex_files.add(nm2)

    # get the names of the expected output files
    # demultiplex_files = utilities.name_files(samples, demultiplex_output_folder, extension=extension)

    workflow.add_task(
        "fastq-multx -B [depends[0]] [depends[1]] [depends[2]] [depends[3]] [depends[4]]\
         -o n/a -o n/a -o [args[0]]/[args[5]]_%[args[3]].[args[1]] -o [args[0]]/[args[5]]_%[args[4]].[args[1]]\
         -q [args[2]] > [targets[0]]",
        depends=[dual_barcode_path, barcode1[0], barcode2[0], input_pair1[0], input_pair2[0]],
        args=[demultiplex_output_folder, extension, min_phred, pair_identifier, pair_identifier2, run_name, fastq_multx_tracked],
        targets=[demultiplex_log, TrackedDirectory(demultiplex_output_folder)],
        name="demultiplex_dual")

    demultiplex_files = demultiplex_check(workflow, demultiplex_log, demultiplex_files)


    return demultiplex_files, demultiplex_output_folder


def demultiplex_check(workflow, demultiplex_log, demultiplex_files):

    # fastq-multx only creates files for those samples with reads mapping to barcodes
    # if a sample in the barcode file does not have any reads, the expected files for
    # that sample will not be created. This task group will create empty files for
    # any samples that do not have reads so that all expected files exist.
    workflow.add_task_group(
        "bash -c \"[ -e [targets[0]] ] || touch [targets[0]]\"",
        depends=[demultiplex_log] * len(demultiplex_files),
        targets=demultiplex_files,
        name="check_demultiplex")

    return demultiplex_files

def generate_dual_barcode(barcode_files, dual_barcode_file):

    """Generate dual barcode file for demultiplexing

        Args:
            workflow (anadama2.workflow): An instance of the workflow class.
            barcode_files (list): A list of barcode files.
            dual_barcode_file (string): A paths to the dual barcode file.

        Requires:
            None

        Returns:
            Creates dual index file (path dual_index_file)
        """

    import itertools

    allbarcodes = set()
    for barcode_file in barcode_files:
        try:
            file_handle = open(barcode_file)
            lines = file_handle.readlines()
            file_handle.close()
        except EnvironmentError:
            sys.exit("ERROR: Unable to read barcode: " + barcode_file)
        allbarcodes.update(lines[1::4])

    dual_indexes_all = list(itertools.combinations(allbarcodes, 2))
    dual_indexes = set(dual_indexes_all)

    fh = open(dual_barcode_file, "w")
    i = 0
    for ind in dual_indexes:
        i+= 1
        ind1 = ind[0].replace("\n", "")
        ind2 = ind[1].replace("\n", "")
        print(ind1,ind2)
        write_string = "Sample_" + str(i) + "\t" + ind1 + "-" + ind2 + "\tNextra\n"
        fh.write(write_string)
    fh.close()

    print("Dual barcode file " + dual_barcode_file + " has been generated")