import argparse
from pathlib import Path
import os
import multiprocess as mp
import pandas as pd
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Directory containing bins", type=str, required=True)
parser.add_argument("-o", help="Write output to this directory", type=str, required=True)
parser.add_argument("-t", help="Number of processes to run in parallel", type=int, default=4)
args = parser.parse_args()

# input
# abspath rather than forcing a leading "/": a relative -i such as "." became
# "/./" and walked the whole filesystem
in_dir = os.path.abspath(args.i) + os.sep

paths = Path(in_dir).rglob('*.fa')

files = []

for path in paths:
	files.append(path.as_posix())

# output
out_dir = os.path.abspath(args.o) + os.sep

if not os.path.exists(out_dir + "tmp/"):
	os.makedirs(out_dir + "tmp/")

# threads
threads = args.t

##############################
# calculate N50 for each MAG #
##############################

def n50(file):
	name = file.replace(".fa", "").split("/")[-1:][0]
	command = "assembly-stats -t -u " + file + " | cut -f1,9 > " + out_dir + "tmp/" + name + ".tsv"
	subprocess.run(command, shell=True)

if __name__ == "__main__":
	pool = mp.Pool(threads)
	for file in files:
		pool.map(n50, [file])

#######################
# combine the outputs #
#######################

tsvs = Path(out_dir + "tmp/").rglob('*.tsv')

dfs = []

for tsv in tsvs:
	df = pd.read_csv(tsv.as_posix(), sep="\t", header=None)
	dfs.append(df)

if dfs:
	df = pd.concat(dfs, axis=0, ignore_index=True)
	df.rename(columns={ df.columns[0]: "MAG", df.columns[1]: "N50" }, inplace = True)
else:
	# no bins to measure, e.g. every sample assembled no contigs above the
	# minimum length; still write the table so downstream tasks have their input
	df = pd.DataFrame(columns=["MAG", "N50"])

df.to_csv(out_dir + "mags_n50.tsv", sep="\t", index=False)

#
