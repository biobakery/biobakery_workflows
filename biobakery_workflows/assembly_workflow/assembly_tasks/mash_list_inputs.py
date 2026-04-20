import argparse
import pandas as pd
from pathlib import Path
import os
import multiprocess as mp
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--checkm", help="CheckM QA and n50 merged table", type=str, required=True)
parser.add_argument("--phylophlan", help="Phylophlan metagenomic relabeled", type=str, required=True)
parser.add_argument("--bins", help="Directory containing bins", type=str, required=True)
parser.add_argument("--mash", help="Write output to this directory", type=str, required=True)
parser.add_argument("--threads", help="Number of processes to run in parallel", type=int, default=4)
args = parser.parse_args()

checkm = args.checkm
bins = "/" + args.bins.strip("/") + "/"
mash = "/" + args.mash.strip("/") + "/"
threads = args.threads

# function to recursively list files

def list_files(directory, string):
	paths = Path(directory).rglob(string)
	files = []
	for path in paths:
		file = path.as_posix()
		files.append(file)
	return files

##################
# list QC'd MAGs #
##################

data = pd.read_csv(checkm, sep="\t")
phylophlan = pd.read_csv(args.phylophlan, sep="\t")

# Keep quality mags that are also unplaced by phylophlan
mags = set(data.loc[data["keep"] == "keep"]["bin_id"]).intersection(set(phylophlan[phylophlan['taxon'] == 'UNKNOWN']['mag']))

##############################
# copy QC'd MAGs to qc_bins/ #
##############################

qc_bins = bins + "qc_bins/"

if not os.path.isdir(qc_bins):
	os.makedirs(qc_bins)

files_to_copy = list_files(bins, "*.fa")

def copy(files):
	for file in files:
		name = file.split("/")[-1].split(".fa")[0]
		if name in mags and not os.path.exists(qc_bins + name + ".fa"):
			command = "cp " + file + " " + qc_bins
			subprocess.run(command, shell=True)

if __name__ == "__main__":
	pool = mp.Pool(threads)
	n = len(files_to_copy)
	for i in range(n):
		subset = [files_to_copy[i]]
		pool.apply_async(copy, (subset,))
	pool.close()
	pool.join()

###############################
# generate a map to QC'd MAGs #
###############################

if not os.path.isdir(mash):
	os.makedirs(mash)

output = mash + "mags_filepaths.txt"

files = list_files(qc_bins, "*.fa")

def print_name(file):
	return file

if __name__ == "__main__":
	p = mp.Pool(threads)
	with open(output, "w") as f:
		for file in p.imap(print_name, files):
			f.write(file + "\n")
