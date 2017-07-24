"""
bioBakery Workflows: config module
Configuration settings for workflows and tasks

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

def get_home_directory():
    """ Try to get the users home directory """
    
    try:
        directory=os.path.expanduser("~")
    except EnvironmentError:
        directory=None
        
    return directory
    
def install_locations():
    """ Get the possible default install locations """
    
    main_folder=Workflow.base_environment_variable.lower()
    options=[os.path.join(os.path.sep,"opt",main_folder)]
    home_dir=get_home_directory()
    if home_dir:
        options+=[os.path.join(home_dir,main_folder)]
        
    return options

class DBInfo(object):
    def __init__(self, name, description, url=None, file_name=None, default_path=None):
        self.name = name
        self.description = description
        self.default_folder = name.lower()
        self.url = url
        
        # if the file name is not set, get it from the url
        if not file_name and url:
            file_name = os.path.basename(url)
        self.file_name = file_name
        
        # set default path if not provided to default folder/file names
        if file_name and not default_path:
            default_path = os.path.join(self.default_folder, file_name)
        elif not default_path:
            default_path = self.default_folder
            
        self.default_path = default_path
        
def get_environment_variable(name):
    """ Try to get the environment variable """

    variable = None
    try:
        variable = os.environ[name]
    except KeyError:
        pass
    
    return variable
    

class Workflow(object):
    base_environment_variable = "BIOBAKERY_WORKFLOWS_DATABASES"
    
    def __getattr__(self, name):
        """ Try to get the database location """
        
        # first check for the environment variable for the database
        variable = get_environment_variable(self.vars[name].name)
        
        if not variable:
            # next check for the overall location variable
            base_folder = get_environment_variable(self.base_environment_variable)
            if base_folder:
                variable = os.path.join(base_folder,self.vars[name].default_path)
                
        if variable and not (os.path.isfile(variable) or os.path.isdir(variable)):
            # if this is not a valid folder/file, then resume search
            variable = None
        
        # check in the default folders, if not found with env variables
        if not variable:
            for folder in install_locations():
                if os.path.isdir(folder):
                    variable = os.path.join(folder,self.vars[name].default_path)
                    if os.path.isfile(variable) or os.path.isdir(variable):
                        return variable
        
        if not variable:
            sys.exit("ERROR: Unable to find database "+self.vars[name].name+
                ". "+self.vars[name].description+" Unable to find in default"+
                " install folders or with environment variables.")
            
        return variable

class ShotGun(Workflow):
    vars={}
    vars["kneaddata_db_human_genome"]=DBInfo("KNEADDATA_DB_HUMAN_GENOME",
        description="This is the KneadData bowtie2 database of the human genome. "+\
            "This database can be downloaded with KneadData.")
    vars["kneaddata_db_human_metatranscriptome"] = DBInfo("KNEADDATA_DB_HUMAN_TRANSCRIPTOME",
        description="This is the folder containing the KneadData bowtie2 database for the human transcriptome."+\
            "This database can be downloaded with KneadData.")
    vars["kneaddata_db_rrna"] = DBInfo("KNEADDATA_DB_RIBOSOMAL_RNA",
        description="This is the folder containing the KneadData bowtie2 database of SILVA rRNA."+\
            "This database can be downloaded with KneadData.")
    vars["strainphlan_db_reference"] = DBInfo("STRAINPHLAN_DB_REFERENCE",
        description="This is the folder containing the reference genomes used "+\
            "when running StrainPhlAn.")
    vars["strainphlan_db_markers"] = DBInfo("STRAINPHLAN_DB_MARKERS",
        description="This is the folder containing the StrainPhlAN marker files.")

class SixteenS(Workflow):
    vars={}
    vars["greengenes_fasta"] = DBInfo("GREEN_GENES_FASTA_DB",
        description="This is the GreenGenes fasta file.",
        url="ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/rep_set/97_otus.fasta")
    vars["greengenes_usearch"] = DBInfo("GREEN_GENES_USEARCH_DB",
        description="This is the GreenGenes fasta file formatted for Usearch. "+\
            "A fasta formatted file can also be used.",
        file_name=vars["greengenes_fasta"].file_name.replace(".fasta",".udb"))
    vars["greengenes_taxonomy"] = DBInfo("GREEN_GENES_TAXONOMY_DB",
        description="This is the GreenGenes taxonomy file matching the fasta file.",
        url="ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt")



