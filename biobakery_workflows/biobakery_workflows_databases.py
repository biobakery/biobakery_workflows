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
import zipfile

def check_dependencies(depends):
    """ Check the required software is installed """
    
    for tool, package in depends:
        try:
            subprocess.check_output(["which",tool],stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            sys.exit("ERROR: Unable to find tool required for database install: "+tool+
                ". Please install "+package+".")

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
        if isinstance(command, list):
            command = " ".join(command)
        sys.exit("Unable to install database. Error running command: "+command)
        
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
        import metaphlan
        strainphlan_db=os.path.join(os.path.dirname(metaphlan.__file__),"metaphlan_databases","mpa_v30_CHOCOPhlAn_201901")
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
        choices=["wmgx","wmgx_demo","wmgx_wmtx","16s_usearch","16s_dada2","16s_its","isolate_assembly"],
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
        print("There are five available database sets each corresponding to a data "+
            "processing workflow.\n"+
            "wmgx: The full databases for the whole metagenome workflow\n"+
            "wmgx_demo: The demo databases for the whole metagenome workflow\n"+
            "wmgx_wmtx: The full databases for the whole metagenome and metatranscriptome workflow\n"+
            "16s_usearch: The full databases for the 16s workflow\n" +
            "16s_dada2: The full databases for the dada2 workflow\n"+
            "16s_its: The unite database for the its workflow\n"+
            "isolate_assembly: The eggnog-mapper databases for the assembly workflow")
        sys.exit(0)
    
    # Check the install location
    if not args.location:
        sys.exit("Please select the install location as the default options are not writable")
        
    if not args.install:
        sys.exit("Please select a database set to install")
        
    # try to create the base install folder
    try_create_folder(args.location)
    
    # check required dependencies are installed for wmgx installs
    dependencies=[]
    if "wmgx" in args.install:
        dependencies=[("humann_databases","HUMAnN2"),("kneaddata_database","Kneaddata"),
            ("strainphlan","StrainPhlAn"),("bowtie2-inspect","Bowtie2")]
    elif "usearch" in args.install:
        dependencies=[("usearch","usearch")]    
    if dependencies:
        check_dependencies(dependencies)
        
    # install humann utility dbs for all shotgun workflows
    humann_install_folder=os.path.join(args.location,"humann")
    if "wmgx" in args.install:
        print("Installing humann utility mapping database")
        run_command(["humann_databases","--download","utility_mapping","full",humann_install_folder])
        
        # create the strainphlan fasta database of markers
        create_strainphlan_db(args.location)
        
    # install the databases based on the workflow selected
    if args.install in ["wmgx","wmgx_wmtx"]:
        # install the full chocophlan and uniref90
        print("Installing humann nucleotide and protein databases")
        run_command(["humann_databases","--download","chocophlan","full",humann_install_folder])
        run_command(["humann_databases","--download","uniref","uniref90_diamond",humann_install_folder])
        
        # install the two kneaddata databases
        print("Installing hg kneaddata database")
        run_command(["kneaddata_database","--download","human_genome","bowtie2",
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_human_genome"].default_folder)])
        
    elif args.install == "wmgx_demo":
        # install the demo chocophlan and demo uniref90
        print("Installing humann DEMO nucleotide and protein databases")
        run_command(["humann_databases","--download","chocophlan","DEMO",humann_install_folder])
        run_command(["humann_databases","--download","uniref","DEMO_diamond",humann_install_folder])
        
        # Install demo kneaddata databases from examples folder to install folder
        print("Installing DEMO hg kneaddata database")
        shutil.copytree(data.get_kneaddata_hg_demo_folder(),
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_human_genome"].default_folder))
        
    elif args.install == "16s_usearch":
        # download the green genes fasta and taxonomy files
        print("Downloading green genes database files")
        usearch_fasta_install_path=os.path.join(args.location,config.SixteenS.vars["greengenes_fasta"].default_path)
        usearch_database_install_path=os.path.join(args.location,config.SixteenS.vars["greengenes_usearch"].default_path)
        utilities.download_file(config.SixteenS.vars["greengenes_fasta"].url,
            usearch_fasta_install_path)
        # use the fasta file also as the database so install works for both usearch and vsearch
        try_create_folder(os.path.dirname(usearch_database_install_path))
        shutil.copy(usearch_fasta_install_path,usearch_database_install_path)
        utilities.download_file(config.SixteenS.vars["greengenes_taxonomy"].url,
            os.path.join(args.location,config.SixteenS.vars["greengenes_taxonomy"].default_path))
        
    elif args.install == "16s_dada2":
        # download the green genes fasta and taxonomy files
        print("Downloading dada2 green genes database files")
        dada2_install_path=os.path.join(args.location,config.SixteenS.vars["greengenes_dada2"].default_path)
        utilities.download_file(config.SixteenS.vars["greengenes_dada2"].url,
            dada2_install_path)
        utilities.download_file(config.SixteenS.vars["rdp_dada2"].url,
            os.path.join(args.location,config.SixteenS.vars["rdp_dada2"].default_path))
        utilities.download_file(config.SixteenS.vars["silva_dada2"].url,
            os.path.join(args.location,config.SixteenS.vars["silva_dada2"].default_path))
        utilities.download_file(config.SixteenS.vars["rdp_species_dada2"].url,
            os.path.join(args.location,config.SixteenS.vars["rdp_species_dada2"].default_path))
        utilities.download_file(config.SixteenS.vars["silva_species_dada2"].url,
            os.path.join(args.location,config.SixteenS.vars["silva_species_dada2"].default_path))

    elif args.install == "16s_its":
        # download unite database for its workflow
        print("Downloading UNITE database files")
        its_install_path = os.path.join(args.location, config.SixteenS.vars["unite_zip"].default_path)
        utilities.download_file(config.SixteenS.vars["unite_zip"].url,
                                its_install_path)
        zip_ref = zipfile.ZipFile(os.path.join(args.location,config.SixteenS.vars["unite_zip"].default_path), 'r')
        zip_ref.extractall(os.path.join(args.location,config.SixteenS.vars["unite"].default_folder))
        zip_ref.close()

    elif args.install == "isolate_assembly":
        print("Downloading eggnog mapper databases")
        eggnog_install_path = os.path.join(args.location, "eggnog_mapper/")
        try_create_folder(eggnog_install_path)
        run_command(["download_eggnog_data.py","--data_dir",eggnog_install_path,"-y"])
        run_command(["ln","-s",eggnog_install_path+"/eggnog.db","/opt/conda/bin/data/"])
        run_command(["ln","-s",eggnog_install_path+"/eggnog_proteins.dmnd","/opt/conda/bin/data/"])

    if "16s" in args.install:
        # install the picrust databases
        run_command(["download_picrust_files.py"])
       
    # if metatranscriptome workflow, install the additional kneaddata database
    if args.install == "wmgx_wmtx":
        print("Installing rRNA and mRNA kneaddata database")
        run_command(["kneaddata_database","--download","ribosomal_RNA","bowtie2",
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_rrna"].default_folder)])
        run_command(["kneaddata_database","--download","human_transcriptome","bowtie2",
            os.path.join(args.location,config.ShotGun.vars["kneaddata_db_human_metatranscriptome"].default_folder)])
        
    # Check for a custom install location
    if args.location != default_install_location():
        print("\n\nA custom install location was selected. Please set the "+
            "environment variable $"+config.Workflow.base_environment_variable+" to the install location.")
    
