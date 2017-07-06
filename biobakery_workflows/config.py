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

class DBInfo(object):
    def __init__(self, name, description):
        self.name = name
        self.description = description

class Workflow(object):
    def __getattr__(self, name):
        """ Try to get the environment variable """
        try:
            variable = os.environ[self.vars[name].name]
        except KeyError:
            message="ERROR: Please set the environment variable $" + self.vars[name].name
            message+="\n"+self.vars[name].description
            sys.exit(message)
            
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
    vars["greengenes_usearch"] = DBInfo("GREEN_GENES_USEARCH_DB",
        description="This is the GreenGenes fasta file formatted for Usearch. "+\
            "A fasta formatted file can also be used.")
    vars["greengenes_fasta"] = DBInfo("GREEN_GENES_FASTA_DB",
        description="This is the GreenGenes fasta file.")
    vars["greengenes_taxonomy"] = DBInfo("GREEN_GENES_TAXONOMY_DB",
        description="This is the GreenGenes taxonomy file matching the fasta file.")



