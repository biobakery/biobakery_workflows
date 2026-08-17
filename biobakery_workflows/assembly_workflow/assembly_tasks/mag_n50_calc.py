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
in_dir = "/" + args.i.strip("/") + "/"

paths = Path(in_dir).rglob('*.fa')

files = []

for path in paths:
	files.append(path.as_posix())

# output
out_dir = "/" + args.o.strip("/") + "/"

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

df = pd.concat(dfs, axis=0, ignore_index=True)

df.rename(columns={ df.columns[0]: "MAG", df.columns[1]: "N50" }, inplace = True)

df.to_csv(out_dir + "mags_n50.tsv", sep="\t", index=False)

#
