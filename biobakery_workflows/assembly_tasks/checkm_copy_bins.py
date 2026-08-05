import argparse
from pathlib import Path
from datetime import date
import os
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-a", help="Copy bins from this directory", type=str, required=True)
parser.add_argument("-b", help="Copy bins to this directory", type=str, required=True)
parser.add_argument("-n", help="Number of bins per sub-directory", type=int, default=1000)
args = parser.parse_args()

####################################
# list bins in the input directory #
####################################

all_bins = args.a

all_paths = Path(all_bins).rglob('*.bin.*')

all_files = []

for path in all_paths:
  file = path.as_posix()
  if not "unbinned" in file and not "tooShort" in file and not "lowDepth" in file and not "final.contigs" in file:
    all_files.append(file)

##################################################
# list pre-existing bins in the output directory #
##################################################

old_bins = args.b

old_paths = Path(old_bins).rglob('*.bin.*')

old_files = []

for path in old_paths:
  file = path.as_posix().split("/")[-1]
  if not "unbinned" in file:
    old_files.append(file)

##################################################
# list bins in input that are absent from output #
##################################################

new_files = []

for file in all_files:
  base = file.split("/")[-1]
  if not base in old_files:
    new_files.append(file)

###########################
# copy the bins to output #
###########################

out_dir = old_bins

# number of bins

n = len(new_files)

# number of bins to add to each sub-directory

size = args.n

# subset the list of bins

subsets = [new_files[i:i+size] for i in range(0, n, size)]

# create the sub-directories and copy bins

for i in range(len(subsets)):
  sub_dir = out_dir + "_" + i + "/"
  if not os.path.isdir(sub_dir):
    os.makedirs(sub_dir)
  sub_files = subsets[i]
  for file in sub_files:
    command = "cp " + file + " " + sub_dir
    subprocess.run(command, shell=True)
