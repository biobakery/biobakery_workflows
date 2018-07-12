#!/usr/bin/env python

import sys
import re
import argparse

from biobakery_workflows import utilities

# This script will check a fastq file is of the expected format.
# To run: $ check_fastq_format.py --input input.fastq

def parse_arguments(args):
    # Parse the arguments from the user
    
    parser = argparse.ArgumentParser(
        description= "Check fastq format\n",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-i", "--input",
        help="the fastq file\n[REQUIRED]",
        metavar="<input.fastq>",
        required=True)
    parser.add_argument(
        "-s","--stop-first-error",
        help="stop after the first formatting error\n[DEFAULT: print all formatting errors]",
        action='store_true')

    return parser.parse_args() 

def write_error(message,id,sequence,qual_id,quality_scores,stop,errors_found):
    # write the error message including the full sequence 
    print("FORMAT ERROR: " + message)
    print("".join([id,sequence,qual_id,quality_scores]))
    if stop:
        sys.exit()

    return errors_found+1

def main():
    # parse command line arguments
    args=parse_arguments(sys)

    # this is a regex of the allowable quality score ascii
    regex = re.compile(r"[a-zA-Z0-9#;/>@?<=\-\\&,\'\$%\*\+\(\)]")

    # check for possible fastq formatting errors
    errors_found=0
    for lines in utilities.read_file_n_lines(args.input,4):
        seq_id, sequence, qual_id, quality_scores = lines
        remainder = regex.sub("",re.escape(quality_scores.rstrip()))
        if len(remainder) > 0:
            errors_found=write_error("Unexpected character in quality scores (found: {})".format(remainder),
                seq_id,sequence, qual_id, quality_scores, args.stop_first_error, errors_found)
        if len(sequence) != len(quality_scores):
            errors_found=write_error("Sequence and quality score lengths differ",
                seq_id,sequence, qual_id, quality_scores, args.stop_first_error, errors_found)
        if len(quality_scores.rstrip()) == 0 or len(sequence.rstrip()) == 0:
            errors_found=write_error("Empty sequence", seq_id, sequence, 
                qual_id, quality_scores, args.stop_first_error, errors_found)
        if seq_id[0] != "@":
            errors_found=write_error("Unexpected format in sequence id (first character is not @)",
                seq_id,sequence, qual_id, quality_scores, args.stop_first_error, errors_found)
        if qual_id[0] != "+":
            errors_found=write_error("Unexpected format in quality id (first character is not +)",
                seq_id,sequence, qual_id, quality_scores, args.stop_first_error, errors_found)

    if errors_found == 0:
        print("No errors identified in fastq format")
    else:
        print("{} errors identified in fastq format".format(errors_found))

if __name__ == "__main__":
    main()

