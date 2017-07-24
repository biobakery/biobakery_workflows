#!/usr/bin/env python

"""
biobakery_workflows: databases module
Download databases for workflows

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

import sys
    
from . import utilities
from . import config
from . import data

import argparse
import os
import subprocess
import shutil

def check_dependencies(depends):
    """ Check the required software is installed """
    
    for tool in depends:
        try:
            subprocess.check_output(["which",tool])
        except subprocess.CalledProcessError:
            sys.exit("ERROR: Unable to find tool required for database install: "+tool)

def default_install_location():
    """ Select the default install location based on the permissions """
    
    for folder in config.install_locations():
        # check the base folder is writeable
        if os.access(os.path.dirname(folder), os.W_OK):
            return folder
        
    return ""

def try_create_folder(folder):
    """ Try to create a folder, exit on error """

    try:
        if not os.path.exists(folder):
            os.makedirs(folder)
    except EnvironmentError:
        sys.exit("Unable to create folder: "+ folder)

def run_command(command,shell=False):
    """ Run the command, exit if error """
    
    try:
        id=subprocess.check_call(command, shell=shell)
    except subprocess.CalledProcessError:
        sys.exit("Unable to install database. Error running command: "+" ".join(command))
        
def create_strainphlan_db(location):
    """ Create the strainphlan fasta file from the bowtie2 indexes """
    
    # get the default fasta install folder
    install_folder=os.path.join(location,config.ShotGun.vars["strainphlan_db_markers"].default_folder)
    install_folder_reference=os.path.join(location,config.ShotGun.vars["strainphlan_db_reference"].default_folder)
    
    # create the install folder if it does not exist
    try_create_folder(install_folder)
    try_create_folder(install_folder_reference)
    
    # find the strainphlan db folder and index files
    try:
        strainphlan_db=os.path.join(os.path.dirname(subprocess.check_output(["which","strainphlan.py"])),"db_v20","mpa_v20_m200")
    except subprocess.CalledProcessError:
        sys.exit("Unable to find strainphlan install.")
        
    # generate the fasta marker files
    print("Generating strainphlan fasta database")
    run_command(" ".join(["bowtie2-inspect",strainphlan_db,">",os.path.join(install_folder,"all_markers.fasta")]),shell=True)

def parse_arguments(args):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "bioBakery workflows databases\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--available", 
        action="store_true",
        help="print the available databases\n")
    parser.add_argument(
        "--install",
        metavar="<workflow>",
        choices=["wmgx","wmgx_demo","wmgx_wmtx","16s"],
        help="install the databases for the selected workflow\n")
    parser.add_argument(
        "--location", 
        default=default_install_location(),
        help="location to install databases [DEFAULT: "+default_install_location()+")]\n")
    
    return parser.parse_args()

def main():
    # Parse arguments from the command line
    args=parse_arguments(sys.argv)
    
    # if no option is set, show available databases
    if not args.install and not args.available:
        args.available = True
        
    if args.available:
        print("There are four available database sets each corresponding to a data "+
            "processing workflow.\n"+
            "wmgx: The full databases for the whole metagenome workflow\n"+
            "wmgx_demo: The demo databases for the whole metagenome workflow\n"+
            "wmgx_wmtx: The full databases for the whole metagenome and metatranscriptome workflow\n"+
            "16s: The full databases for the 16s workflow")
        sys.exit(0)
    
    # Check the install location
    if not args.location:
        sys.exit("Please select the install location as the default options are not writable")
        
    if not args.install:
        sys.exit("Please select a database set to install")
        
    # try to create the base install folder
    try_create_folder(args.location)
    
    # check required dependencies are installed for wmgx installs
    if "wmgx" in args.install:
        dependencies=["humann2_databases","kneaddata_database","strainphlan.py","bowtie2-inspect"]
    else:
        dependencies=["usearch"]    
    check_dependencies(dependencies)
        
    # install humann2 utility dbs for all shotgun workflows
    humann2_install_folder=os.path.join(args.location,"humann2")
    if "wmgx" in args.install:
        print("Installing humann2 utility mapping database")
        run_command(["humann2_databases","--download","utility_mapping","full",humann2_install_folder])
        
        # create the strainphlan fasta database of markers
        create_strainphlan_db(args.location)
        
    # install the databases based on the workflow selected
    if args.install in ["wmgx","wmgx_wmtx"]:
        # install the full chocophlan and uniref90
        print("Installing humann2 nucleotide and protein databases")
        run_command(["humann2_databases","--download","chocophlan","full",humann2_install_folder])
        run_command(["humann2_databases","--download","uniref","uniref90_diamond",humann2_install_folder])
        
        # install the two kneaddata databases
        print("Installing hg and rRNA kneaddata databases")
        run_command(["kneaddata_database","--download","human_genome","bowtie2",
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_human_genome"].default_folder)])
        run_command(["kneaddata_database","--download","ribosomal_RNA","bowtie2",
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_rrna"].default_folder)])
        
    elif args.install == "wmgx_demo":
        # install the demo chocophlan and demo uniref90
        print("Installing humann2 DEMO nucleotide and protein databases")
        run_command(["humann2_databases","--download","chocophlan","DEMO",humann2_install_folder])
        run_command(["humann2_databases","--download","uniref","DEMO_diamond",humann2_install_folder])
        
        # Install demo kneaddata databases from examples folder to install folder
        print("Installing DEMO hg and rRNA kneaddata databases")
        shutil.copytree(data.get_kneaddata_hg_demo_folder(),
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_human_genome"].default_folder))
        shutil.copytree(data.get_kneaddata_silva_demo_folder(),
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_rrna"].default_folder))
        
    elif args.install == "16s":
        # download the green genes fasta and taxonomy files
        print("Downloading green genes database files")
        usearch_fasta_install_path=os.path.join(args.location,config.SixteenS.vars["greengenes_fasta"].default_path)
        utilities.download_file(config.SixteenS.vars["greengenes_fasta"].url,
            usearch_fasta_install_path)
        utilities.download_file(config.SixteenS.vars["greengenes_taxonomy"].url,
            os.path.join(args.location,config.SixteenS.vars["greengenes_taxonomy"].default_path))
        
        # create the usearch database from the fasta file
        print("Creating usearch green genes database")
        usearch_db_folder=os.path.join(args.location,config.SixteenS.vars["greengenes_usearch"].default_folder)
        try_create_folder(usearch_db_folder)
        run_command(["usearch","-makeudb_usearch",usearch_fasta_install_path,
            "-output",os.path.join(args.location,config.SixteenS.vars["greengenes_usearch"].default_path)])

    # if metatranscriptome workflow, install the additional kneaddata database
    if args.install == "wmgx_wmtx":
        print("Installing mRNA kneaddata database")
        run_command(["kneaddata_database","--download","human_transcriptome","bowtie2",
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_human_metatranscriptome"].default_folder)])
        
    # Check for a custom install location
    if args.location != default_install_location():
        print("\n\nA custom install location was selected. Please set the "+
            "environment variable $"+config.Workflow.base_environment_variable+" to the install location.")
    
