#!/usr/bin/env python

""" This script will remove a file if it exists.
    It will report error messages from trying to remove the file. """
    
import sys
import shutil

try:
    import argparse
except ImportError:
    sys.exit("Please upgrade to at least python v2.7")
    
import os

def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Remove file if it exists.")
    parser.add_argument('file',help="File (or folder) to remove")
    parser.add_argument('--is-folder',help="Indicate if input is a folder", action="store_true")

    return parser.parse_args()    

    
def main():
    # parse arguments
    args = parse_arguments(sys.argv)
    
    # check if the file exists
    if os.path.isfile(args.file):
        # if it exists, then try to remove
        # Report error messages if unable to remove file that exists
        try:
            os.unlink(args.file)
        except EnvironmentError:
            sys.exit("ERROR: Unable to remove file: " + args.file)
    if args.is_folder:
        try:
            shutil.rmtree(args.file)    
        except EnvironmentError:
            sys.exit("ERROR: Unable to remove folder: " + args.file)

    
if __name__ == "__main__":
    main()
