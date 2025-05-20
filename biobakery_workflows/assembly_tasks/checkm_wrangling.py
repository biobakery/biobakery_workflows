import argparse
import glob
import os
import subprocess
import shutil
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(prog="checkm_wrangling.py")
parser.add_argument("--checkm-qa", help="checkm qa table", type=str, required=True)
parser.add_argument("--n50", help="n50 table", type=str, required=True)
parser.add_argument("--out_file", help="output file", type=str, required=True)
parser.add_argument("--completeness", help="completeness", type=float, default=0.5)
parser.add_argument("--contamination", help="contamination", type=float, default=0.1)
args = parser.parse_args()

checkm_qa = args.checkm_qa
n50 = args.n50
out_file = args.out_file

data = pd.read_csv(checkm_qa, sep='\t', header=0)
data = data.iloc[:, [0, 1, 2]]

data = data.rename(columns={'Name': "bin_id", "Completeness": "completeness", "Contamination": "contamination"})
data = data[~data["bin_id"].str.contains("unbinned|tooShort|lowDepth")]
data['quality'] = np.where(
    (data['completeness'] >= 90) & (data['contamination'] <= 5), "high_quality", np.where(
    ((data['completeness'] >= 90) & (data['contamination'] >= 5) & (data['contamination'] <= 10)) | ((data['completeness'] <= 90) & (data['completeness'] >= 50) & (data['contamination'] <= 10)), "medium_quality", "low_quality"))
data['keep'] = np.where((data['completeness'] >= args.completeness) & (data['contamination'] <= args.contamination), "keep", "reject")

n50 = pd.read_csv(n50, sep='\t', header=0)
n50 = n50.rename(columns={'MAG': "bin_id", "N50": "n50"})
n50['bin_id'] = [name.split("/")[-1].rsplit(".",1)[0] for name in n50['bin_id']]
output = data.merge(n50, left_on='bin_id', right_on='bin_id')

output.to_csv(out_file, sep='\t', index=False)
