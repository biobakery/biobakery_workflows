#!/usr/bin/env python

"""
bioBakery Workflows scripts : A collection of AnADAMA2 workflows

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


import sys
import argparse
import os
import subprocess

from . import utilities

VERSION = "3.0.0-alpha.4"
SCRIPTS_FOLDER="Rscripts"
VIS_SCRIPT_EXTENSION="_vis.R"
SCRIPT_EXTENSION=".R"

def find_vis_scripts():
    """ Search for installed vis scripts """
    
    scripts={}
    scripts_folder=os.path.join(os.path.dirname(os.path.abspath(__file__)),SCRIPTS_FOLDER)
    for file in os.listdir(scripts_folder):
        if file.endswith(VIS_SCRIPT_EXTENSION):
            scripts[file.replace(SCRIPT_EXTENSION,"")]=file.replace(SCRIPT_EXTENSION,"")
    
    return scripts

def parse_arguments(args,scripts):
    """ 
    Parse the arguments from the user
    """
    parser = argparse.ArgumentParser(
        description= "bioBakery workflows visualization scripts\n",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog="biobakery_workflows_scripts")
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s v"+VERSION)
    parser.add_argument(
        "script",
        choices=scripts,
        help="visualization script to run")
    
    return parser.parse_args(args)

def run_script(args, script):
    """ Run the script with the arguments provided """

    # check for visualization script
    visualization_script = utilities.get_package_file(script, "Rscript")

    if not visualization_script:
        sys.exit("ERROR: Unable to find script "+script)
    
    try:
        command = [visualization_script]+args[2:]
        subprocess.call(command)
    except ( subprocess.CalledProcessError, EnvironmentError):
        sys.exit("Error: Unable to run: " +" ".join(command))


def main():
    scripts=find_vis_scripts()
    
    # parse the arguments (only the first two as the rest are for the script)
    args=parse_arguments(sys.argv[1:2],scripts.keys())
    
    # run the script (providing all of the arguments)
    run_script(sys.argv,scripts[args.script])
    
    
if __name__ == "__main__":
    main()
    
