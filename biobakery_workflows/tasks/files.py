"""
bioBakery Workflows: tasks.files module
A collection of file names used by tasks

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
import copy

from biobakery_workflows import utilities

def settings(name=None, subfolder=None, tag=None, extension=None):
    keywords={"names":name, "subfolder":subfolder, "tag":tag, "extension":extension}
    return {key:value for key, value in keywords.items() if value}

class Workflow(object):
    @classmethod
    def file(cls, name, main_folder, **keywords):
        merged_keywords = copy.copy(keywords)
        merged_keywords.update(cls.file_settings[name])
        return utilities.name_files(folder=main_folder, **merged_keywords)

class ShotGun(Workflow):
    file_settings={}
    
    # set the kneaddata file name
    file_settings["kneaddata_read_counts"]=settings("kneaddata_read_count_table.tsv",subfolder="counts")

    # set the taxonomy file names
    file_settings["taxonomic_profile"]=settings("taxonomic_profiles.tsv")
    file_settings["species_counts"]=settings("metaphlan2_species_counts_table.tsv",subfolder="counts")
    
    # set the merged feature file names
    file_settings["genefamilies"]=settings("genefamilies.tsv")
    file_settings["ecs"]=settings("ecs.tsv")
    file_settings["pathabundance"]=settings("pathabundance.tsv")
    
    # set the normed feature file names
    file_settings["genefamilies_relab"]=settings("genefamilies_relab.tsv")
    file_settings["ecs_relab"]=settings("ecs_relab.tsv")
    file_settings["pathabundance_relab"]=settings("pathabundance_relab.tsv")
    
    # set the feature count file names
    file_settings["genefamilies_relab_counts"]=settings("humann2_genefamilies_relab_counts.tsv", subfolder="counts")
    file_settings["ecs_relab_counts"]=settings("humann2_ecs_relab_counts.tsv", subfolder="counts")
    file_settings["pathabundance_relab_counts"]=settings("humann2_pathabundance_relab_counts.tsv", subfolder="counts")
    
    # set the all feature counts file names
    file_settings["feature_counts"]=settings("humann2_feature_counts.tsv", subfolder="counts")
    file_settings["humann2_read_counts"]=settings("humann2_read_and_species_count_table.tsv",subfolder="counts")
    
    # set the names for the rna/dna normed files
    file_settings["genefamilies_norm_ratio"]=settings("rna_dna_relative_expression_unstratified.tsv",subfolder=os.path.join("norm","genes"))
    file_settings["ecs_norm_ratio"]=settings("rna_dna_relative_expression_unstratified.tsv",subfolder=os.path.join("norm","ecs"))
    file_settings["paths_norm_ratio"]=settings("rna_dna_relative_expression_unstratified.tsv",subfolder=os.path.join("norm","paths"))

    
    