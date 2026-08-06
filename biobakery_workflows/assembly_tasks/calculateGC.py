import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import os
import glob
import multiprocessing as mp
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--in-dir", help="input directory")
parser.add_argument("--out-file", help="output directory")
parser.add_argument("--threads", help="threads to use", type=int, default = 1)
args = parser.parse_args()

in_dir = args.in_dir
out_file = args.out_file
threads = args.threads

exclude = [".unbinned.", ".tooShort.", ".lowDepth.", "final.contigs"]

# abspath, and os.walk(followlinks=True) rather than rglob, which will not
# descend into symlinked directories such as staged workflow inputs
paths = []
for _root, _dirs, _names in os.walk(os.path.abspath(in_dir), followlinks=True):
	for _n in _names:
		if _n.endswith(".fa"):
			paths.append(Path(os.path.join(_root, _n)))

files = []
for path in paths:
    if not any(x in str(path) for x in exclude):
        files.append(str(path))

def calculate_gc(fasta):
    total_GC = 0
    total_length = 0
    for record in SeqIO.parse(fasta,'fasta'):
        total_GC += GC(record.seq) * len(record)
        total_length += len(record)

    if total_length > 0:
        gc_content = total_GC/total_length
    else:
        gc_content = 0
    gc_content=round(gc_content,2)
    filename = fasta.split("/")[-1].split(".")[0]

    return (filename, gc_content, total_length)

with mp.Pool(processes = threads) as p:
    gc_out = p.map(calculate_gc, files)

if not os.path.isdir(out_file.rsplit("/", 1)[0]):
    os.makedirs(out_file.rsplit("/", 1)[0])

pd.DataFrame(gc_out, columns=['filename', 'gc_contents', 'length']).to_csv(out_file, sep='\t', index=False)
