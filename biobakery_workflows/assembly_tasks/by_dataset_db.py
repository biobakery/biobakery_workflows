import argparse
import glob
import os
import subprocess
import shutil

parser = argparse.ArgumentParser(prog="by_dataset_db.py")
parser.add_argument("--mag_dir", help="Directory containing MAGs", type=str, required=True)
parser.add_argument("--out_dir", help="Output directory", type=str, required=True)
args = parser.parse_args()

mag_dir = "/" + args.mag_dir.strip("/") + "/"
out_dir = "/" + args.out_dir.strip("/") + "/"

# create tmp directory
tmp = out_dir + "tmp/"

if not os.path.exists(tmp):
	os.makedirs(tmp)

files = glob.glob(mag_dir + "*/bins/*.fa")

for file in files:
	copy = "cp " + file + " " + tmp
	subprocess.run(copy, shell=True)
	name = file.split("/")[-1:][0]
	tmp_file = tmp + name
	part_a = "sed -i 's/>" + name + "_/>/g' " + tmp_file
	part_b = "sed -i 's/>/>" + name + "_/g' " + tmp_file
	header = part_a + "; " + part_b
	subprocess.run(header, shell=True)

merged = out_dir + "mag_db.fa"

if not os.path.exists(merged):
	subprocess.run("cat " + tmp + "*.fa > " + merged, shell=True)

# create the bowtie2 indices
index = out_dir + "mag_db.1.bt2"

if not os.path.exists(index):
	subprocess.run("bowtie2-build " + merged + " " + out_dir + "mag_db", shell=True)
