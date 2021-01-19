#!/usr/bin/env python

""" This script will take a set of images and create a single image tile. """
   
import sys
import argparse

from biobakery_workflows import utilities

def parse_arguments(args):
    """ Parse the arguments from the user"""
    
    parser=argparse.ArgumentParser(description="Create a tile of images.")
    parser.add_argument('--input',help="Comma delimited list of input image files", required=True)
    parser.add_argument('--output',help="Output file name", required=True)

    return parser.parse_args()    
    
def main():
    # parse arguments
    args = parse_arguments(sys.argv)
    
    utilities.generate_tile_of_images(args.input.split(","),args.output)    

if __name__ == "__main__":
    main()
