"""
bioBakery Workflows: statistical analysis visualizations module
A collection of utilities for the statistical analysis visualizations

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


import subprocess
import utilities
from anadama2.tracked import TrackedDirectory
import os


def run_R(script,abandance,input_dir):
    subprocess.call(['Rscript', script, abandance,input_dir])

def split_pdfs_do(task):
    import fnmatch
    from anadama2.util import get_name

    pdfs_path=get_name(task.depends[1])
    for f in fnmatch.filter(os.listdir(pdfs_path), '*.pdf'):
        pdf_name=f.replace(".pdf", "")
        pdf_images=str(pdfs_path)+"/"+pdf_name+"_images"
        command="convert "+ str(pdfs_path)+"/"+str(f)+" "+str(pdf_images)+"/"+str(pdf_name)+"-%04d.png"
        utilities.run_task("mkdir "+pdf_images,depends=task.depends[1], targets=task.targets)
        utilities.run_task(command, depends=task.depends[1], targets=task.targets)

def split_pdfs(workflow,pdfs_path,perv_task):

    split_task=workflow.add_task(
        split_pdfs_do,
        depends=[perv_task,TrackedDirectory(pdfs_path)],
        targets=[TrackedDirectory(pdfs_path)],
        name="split_pdfs"
    )

    return split_task

def run_permanova(workflow,abundance,input_dir,output_dir,adonis_dir,workflow_data):

    adonis_dirname=os.path.basename(adonis_dir)
    metadata_path = input_dir + "/metadata.tsv"
    options_file=input_dir+"/"+adonis_dirname+"_config.txt"
    options=""
    if os.path.isfile(options_file):
        with open(options_file) as f:
            options_string = f.read().splitlines()
        options=" ".join(options_string)

    script_path = utilities.get_package_file("permanova", "Rscript")

    permanova_task=workflow.add_task(
        "[vars[0]] --input_dir=[args[0]] --abundance=[depends[0]] --metadata=[depends[1]]\
         --output_dir=[args[1]] --adonis_dir=[args[2]] --workflow=[args[3]] [vars[1]]",
        depends=[abundance,metadata_path],
        targets=[TrackedDirectory(adonis_dir)],
        vars=[script_path,options],
        args=[input_dir,output_dir,adonis_dir,workflow_data],
        name=adonis_dirname
    )

    return permanova_task

def run_maaslin(workflow,abundance,input_dir,maaslin_dir):

    maaslin_dirname = os.path.basename(maaslin_dir)
    metadata_path=input_dir+"/metadata.tsv"
    options_file=input_dir+"/"+maaslin_dirname+"_config.txt"
    options=""
    if os.path.isfile(options_file):
        with open(options_file) as f:
            options_string = f.read().splitlines()
        options=" ".join(options_string)

    script_path = utilities.get_package_file("maaslin_call", "Rscript")

    maaslin_task=workflow.add_task(
        "[vars[0]]  --input_dir=[args[0]] --abundance=[depends[0]] --metadata=[depends[1]] --maaslin_output=[args[1]]",
        depends=[abundance,metadata_path],
        targets=[TrackedDirectory(maaslin_dir)],
        vars=[script_path,options],
        args=[input_dir,maaslin_dir],
        name=maaslin_dirname
    )

    return maaslin_task

class Workflow(object):

    @classmethod
    def format_caption(cls, name, **keywords):
        return cls.captions[name].format(**keywords)

class Stats(Workflow):
    captions={}
    
    # add captions for functional data section
    captions["intro"]="This is Statistical Analysis Workflow introduction"

class Permanova(Workflow):
    captions={}

    captions["intro"]="This is Permanova introduction"

class Maaslin(Workflow):
    captions={}

    captions["intro"]="This is Maaslin introduction"

