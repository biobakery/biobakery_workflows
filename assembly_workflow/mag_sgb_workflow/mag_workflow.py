#!/usr/bin/env python

from anadama2 import Workflow
import glob
import os
import math
import shutil
from anadama2.tracked import TrackedDirectory
from pathlib import Path

# Parse arguments
workflow = Workflow(version="0.4", description="MAG and SGB workflow")
workflow.add_argument("cores", desc="The number of CPU cores allocated to the job", type=int, default=4)
workflow.add_argument("mem", desc="The maximum memory in megabytes allocated to run any individual command", type=int, default=180000)
workflow.add_argument("input-extension", desc="the input file extension", default="fastq.gz")
workflow.add_argument("time", desc="The maximum time in minutes allocated to run any individual command", type=int, default=10000)
workflow.add_argument("paired", desc="Whether the inputs are \"unpaired\", \"paired\", or \"concatenated\"", default="paired")
workflow.add_argument("phylophlan-database", desc="Database name", default="SGB.Jul20")
workflow.add_argument("phylophlan-database-folder", desc="Folder with the phylophlan database files such as SGB.Jul20/, SGB.Jul20.txt.bz2, etc.")
workflow.add_argument("checkm-data-path", desc='CHECKM_DATA_PATH with selected_marker_sets.tsv, taxon_marker_sets.tsv, etc. if not already set')
workflow.add_argument("checkm-path", desc='path to checkm2 directory with bin/')
workflow.add_argument("checkm-predict-options", desc='checkm2 predict options as a text string with quotes', default="")
workflow.add_argument("skip-contigs", desc="Whether to skip MEGAHIT, contigs should be in $OUTPUT_DIRECTORY/assembly/main/$SAMPLE_NAME/$SAMPLE_NAME.final.contigs.fa")
workflow.add_argument("min-contig-length", desc='MEGAHIT --min-contig-length and MetaBAT -m parameter', default=1500)
workflow.add_argument("megahit-options", desc='MEGAHIT options as a text string with quotes', default="")
workflow.add_argument("metabat-options", desc='MetaBAT options as a text string with quotes', default="")
workflow.add_argument("checkm-coverage-options", desc='checkm coverage options as a text string with quotes', default="")
workflow.add_argument("checkm-predict-options", desc='checkm2 predict options as a text string with quotes', default="")
workflow.add_argument("abundance-type", desc='by_sample or by_dataset', default="by_sample")
workflow.add_argument("completeness", desc='completeness threshold for retaining bins', default=50)
workflow.add_argument("contamination", desc='contamination threshold for retaining bins', default=10)
workflow.add_argument("phylophlan-metagenomic-options", desc='PhyloPhlAn metagenomic options as a text string with quotes', default="")
workflow.add_argument("gc-length-stats", desc='calculate GC and length stats for each bin')
workflow.add_argument("mash-sketch-options", desc='Mash sketch options as a text string with quotes', default="")
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

# Set CHECKM_DATA_PATH
if args.checkm_data_path:
	os.environ["CHECKM_DATA_PATH"] = args.checkm_data_path

try:
	os.environ["CHECKM_DATA_PATH"]
except:
	raise ValueError("CHECKM_DATA_PATH not provided or set")

if not args.checkm_path:
	raise ValueError("No checkm_path provided")

database = args.phylophlan_database
database_folder = args.phylophlan_database_folder

# Check valid input extension
input_extension = args.input_extension
try:
	assert input_extension in ['fastq', 'fastq.gz', 'fq', 'fq.gz']
except:
	raise ValueError("--input-extension must be fastq, fastq.gz, fq, or fq.gz")

paired = args.paired
try:
	assert paired in ['paired', 'unpaired', 'concatenated']
except:
	raise ValueError("--paired must be paired, unpaired, or concatenated")

abundance_type = args.abundance_type
try:
	assert abundance_type in ['by_sample', 'by_dataset']
except:
	raise ValueError("--abundance must be by_sample or by_dataset")

# output
output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
	os.makedirs(output)

# make necessary directories
def make_directory(path):
	if not os.path.isdir(path):
		os.makedirs(path)

# scratch directory
scratch = "/" + args.grid_scratch.strip("/") + "/"
make_directory(scratch)

scratch_searched = scratch + "searched/"
make_directory(scratch_searched)

scratch_deconcatenated = scratch + "deconcatenated/"
make_directory(scratch_deconcatenated)

deconcatenated_dir = output + "deconcatenated/"
make_directory(deconcatenated_dir)

assembly_dir = output + "assembly/"
make_directory(assembly_dir)

megahit_scratch = scratch + "megahit/"
make_directory(megahit_scratch)

contigs_dir = assembly_dir + "main/"
make_directory(contigs_dir)

mags_scratch = scratch + "bins/"
make_directory(mags_scratch)

bowtie2_tmp = scratch + "bowtie2_tmp/"
make_directory(bowtie2_tmp)

bowtie2_global_dir = scratch + "bowtie2_global_dir/"
make_directory(bowtie2_global_dir)

depths_dir = assembly_dir + "contig_depths/"
make_directory(depths_dir)

bins_dir = output + "bins/"
make_directory(bins_dir)

abundance_dir = output + "abundance_" + abundance_type + "/"
make_directory(abundance_dir)

checkm_dir = output + "checkm/"
make_directory(checkm_dir)

checkm_n50_dir = checkm_dir + "n50/"
make_directory(checkm_n50_dir)

qa_dir = checkm_dir + "qa/"
make_directory(qa_dir)

phylophlan_dir = output + "phylophlan/"
make_directory(phylophlan_dir)

#phylophlan_dir_scratch = scratch + "phylophlan/"

sgb_dir = output + "sgbs/"
make_directory(sgb_dir)

mash_dir = sgb_dir + "mash/"
make_directory(mash_dir)

sgbs_scratch = scratch + "sgbs/"
make_directory(sgbs_scratch)

sgbs = sgb_dir + "sgbs/"

map_out = sgbs_scratch + "abundances/"
if not os.path.isdir(map_out):
	os.makedirs(map_out)

# grid
memory = args.mem
cores = args.cores
partition = args.grid_partition
max_time = args.time
local_jobs = args.jobs

# list the input fastq files
in_dir = args.input

# Get filepath without paired/unpaired ending
if paired == "paired":
	paths = glob.glob("/" + in_dir.strip("/") + "/" + '*.' + input_extension)
	names = set(file.split('_paired_1')[0].split('_paired_2')[0].split('_unmatched_1')[0].split('_unmatched_2')[0] for file in paths)
else:
	paths = glob.glob("/" + in_dir.strip("/") + "/" + '*.' + input_extension)
	names = set(file.split("." + input_extension)[0] for file in paths)

#######################################
# function to calculate tool runtimes #
#######################################

def calculate_time(name, step, paired):
	time = 0
	if step == "deconcatenate":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		time = 40 * n_gigabytes
	elif paired == "paired":
		n_gigabytes = math.ceil(os.path.getsize(name + "_paired_1." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "megahit":
			time = 210 * n_gigabytes + 3
		elif step == "align":
			time = 60 * n_gigabytes + 3
		elif step == "metabat":
			time = 30 * n_gigabytes + 3
		elif step == "abundance":
			time = 20 * n_gigabytes + 3
	elif paired == "unpaired":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "megahit":
			time = 105 * n_gigabytes + 3
		elif step == "align":
			time = 30 * n_gigabytes + 3
		if step == "metabat":
			time = 15 * n_gigabytes + 3
		if step == "abundance":
			time = 10 * n_gigabytes + 3
	elif paired == "concatenated":
		n_gigabytes = math.ceil(os.path.getsize(name + "." + input_extension) / (1024 * 1024 * 1024.0))
		if step == "megahit":
			time = 210 * n_gigabytes + 3
		elif step == "align":
			time = 60 * n_gigabytes + 3
		if step == "metabat":
			time = 30 * n_gigabytes + 3
		if step == "abundance":
			time = 20 * n_gigabytes + 3
	if 'gz' in input_extension:
		time = time * 6
	if time > max_time:
		time = max_time
	return int(time)

#################################
# function to list dependencies #
#################################

def list_depends(name, step, paired):
	if step == "deconcatenate":
		return [name + "." + input_extension].extend([scratch_searched + name.split("/")[-1] + "_searched.log" for name in names])
	elif step == "megahit":
		if paired == "paired":
			return [str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension), str(name + "_unmatched_1." + input_extension), str(name + "_unmatched_2." + input_extension)]
		elif paired == "concatenated":
			return [str(deconcatenated_dir + name.split("/")[-1] + ".done"), list_paired]
		else:
			return [str(name + "." + input_extension)]
	elif step == "align":
		return [str(contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa")]
	elif step == "metabat":
		return [str(depths_dir + name.split("/")[-1] + ".contig_depths.txt")]
	elif step == "abundance_sample":
		if paired == "paired":
			return [str(bins_dir + name.split("/")[-1] + ".done"), str(name + "_paired_1." + input_extension), str(name + "_paired_2." + input_extension), str(name + "_unmatched_1." + input_extension), str(name + "_unmatched_2." + input_extension)]
		else:
			return [str(bins_dir + name.split("/")[-1] + ".done"), str(name + "." + input_extension)]

############################
# function to list targets #
############################

def list_targets(name, step, paired):
	if step == "deconcatenate":
		return [deconcatenated_dir + name.split("/")[-1] + ".done"]
	elif step == "megahit":
		return [str(contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa")]
	elif step == "align":
		return [str(depths_dir + name.split("/")[-1] + ".contig_depths.txt")]
	elif step == "metabat":
		return [str(bins_dir + name.split("/")[-1] + ".done")]
	elif step == "abundance":
		targets0 = abundance_dir + name.split("/")[-1] + ".coverage.tsv"
		targets1 = abundance_dir + name.split("/")[-1] + ".abundance.tsv"
		targets2 = abundance_dir + name.split("/")[-1] + ".mapped_read_num.txt"
		return [str(targets0), str(targets1), str(targets2)]

#######################################################################
# function to detect paired concatenated files and deconcatenate them #
#######################################################################

def deconcatenate(name):
	if input_extension in ["fastq.gz", "fq.gz"]:
		unzipped_name = scratch_deconcatenated + name.split("/")[-1] + ".fastq"
		command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {k}'''.format(
			a = "if grep -q -m 1 " + name + " " + list_paired + "; then  gunzip -c " + name + "." + input_extension + " > " + unzipped_name,
			b = "python " + this_folder + "deconcatenate.py " + unzipped_name + " " + scratch_deconcatenated + name.split("/")[-1],
			c = "rm " + unzipped_name,
			d = "gzip -f " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_1.fastq",
			e = "gzip -f " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_2.fastq",
			f = "cat " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip -f > " + scratch_deconcatenated + name.split("/")[-1] + "_1_sorted.fastq.gz",
			g = "cat " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip -f > " + scratch_deconcatenated + name.split("/")[-1] + "_2_sorted.fastq.gz",
			h = "rm " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq",
			i = "rm " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq",
			j = "mv " + scratch_deconcatenated + name.split("/")[-1] + "_1_sorted.fastq.gz " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq.gz",
			k = "mv " + scratch_deconcatenated + name.split("/")[-1] + "_2_sorted.fastq.gz " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq.gz; fi && touch [targets[0]]"
			)
	else:
		command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i}'''.format(
			a = "if grep -q -m 1 " + name + " " + list_paired + "; then  python " + this_folder + "deconcatenate.py " + name + "." + input_extension + " " + scratch_deconcatenated + name.split("/")[-1],
			b = "gzip -f " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_1.fastq",
			c = "gzip -f " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_2.fastq",
			d = "cat " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip -f > " + scratch_deconcatenated + name.split("/")[-1] + "_1_sorted.fastq.gz",
			e = "cat " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' | gzip -f > " + scratch_deconcatenated + name.split("/")[-1] + "_2_sorted.fastq.gz",
			f = "rm " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq",
			g = "rm " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq",
			h = "mv " + scratch_deconcatenated + name.split("/")[-1] + "_1_sorted.fastq.gz " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq.gz",
			i = "mv " + scratch_deconcatenated + name.split("/")[-1] + "_2_sorted.fastq.gz " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq.gz; fi && touch [targets[0]]"
			)
	return command

if paired == "concatenated":
	list_paired = scratch + "paired_list.txt"
	for name in names:
		if input_extension in ["fastq.gz", "fq.gz"]:
			command = "if zgrep -q -m 1 /2$ " + name + "." + input_extension + "; then echo " + name + "." + input_extension + " >> " + list_paired + "; fi && touch " + scratch_searched + name.split("/")[-1] + "_searched.log"
			workflow.add_task(command, depends=[name + "." + input_extension], targets = [scratch_searched + name.split("/")[-1] + "_searched.log", list_paired])
		else:
			command = "if grep -q -m 1 /2$ " + name + "." + input_extension + "; then echo " + name + "." + input_extension + " >> " + list_paired + "; fi && touch " + scratch_searched + name.split("/")[-1] + "_searched.log"
			workflow.add_task(command, depends=[name + "." + input_extension], targets = [scratch_searched + name.split("/")[-1] + "_searched.log", list_paired])

	for name in names:
		workflow.add_task_gridable(deconcatenate(name),
			depends=list_depends(name=name, step="deconcatenate", paired="concatenated"),
			targets=list_targets(name=name, step="deconcatenate", paired="concatenated"),
			time=calculate_time(name=name, step="deconcatenate", paired="concatenated"),
			mem=memory,
			cores=cores,
			partition=partition
			)

###################################################
# function to call MEGAHIT and calculate coverage #
###################################################

def megahit(name, paired):
	if paired == "paired":
		command = '''{a} && {b} && {c}'''.format(
			a = "megahit -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -r " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " -o " + megahit_scratch + name.split("/")[-1] + " --out-prefix " + name.split("/")[-1] + " -m 0.99 -t " + str(cores) + " --continue --min-contig-len " + str(args.min_contig_length) + " " + args.megahit_options,
			b = "mkdir -p " + contigs_dir + name.split("/")[-1],
			c = "mv " + megahit_scratch + name.split("/")[-1] + "/" + name.split("/")[-1] + ".contigs.fa " + contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
			)
	if paired == "unpaired":
		command = '''{a} && {b} && {c}'''.format(
			a = "megahit -r " + name + "." + input_extension + " -o " + megahit_scratch + name.split("/")[-1] + " --out-prefix " + name.split("/")[-1] + " -m 0.99 -t " + str(cores) + " --continue --min-contig-len " + str(args.min_contig_length) + " " + args.megahit_options,
			b = "mkdir -p " + contigs_dir + name.split("/")[-1],
			c = "mv " + megahit_scratch + name.split("/")[-1] + "/" + name.split("/")[-1] + ".contigs.fa " + contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
			)
	if paired == "concatenated":
		command = '''{a} && {b} && {c}'''.format(
			a = "if grep -q -m 1 " + name + " " + list_paired + "; then megahit -1 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq.gz" + " -2 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq.gz" +
				" -r " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_1.fastq.gz" + "," + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_2.fastq.gz" + " -o " + megahit_scratch + name.split("/")[-1] + " --out-prefix " + name.split("/")[-1] + " -m 0.99 -t " + str(cores) +
				" --continue --min-contig-len " + str(args.min_contig_length) + " " + args.megahit_options + "; else megahit -r " + name + "." + input_extension + " -o " + megahit_scratch + name.split("/")[-1] + " --out-prefix " + name.split("/")[-1] + " -m 0.99 -t " + str(cores) + " --continue --min-contig-len " + str(args.min_contig_length) + " " + args.megahit_options + "; fi",
			b = "mkdir -p " + contigs_dir + name.split("/")[-1],
			c = "mv " + megahit_scratch + name.split("/")[-1] + "/" + name.split("/")[-1] + ".contigs.fa " + contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
			)
	return str(command)

if not args.skip_contigs:
	for name in names:
		if not os.path.isfile(list_targets(name=name, step="megahit", paired=paired)[0]):
			workflow.add_task_gridable(actions=megahit(name, paired),
				depends=list_depends(name=name, step="megahit", paired=paired),
				targets=list_targets(name=name, step="megahit", paired=paired),
				time=calculate_time(name=name, step="megahit", paired=paired),
				mem=memory,
				cores=cores,
				partition=partition
				)

###############################################
# function to align reads and calculate depth #
###############################################

def align(name, paired):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = mags_scratch + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_dir + name.split("/")[-1]
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"

	if paired == "paired":
		command = '''{a} && {b} && {c} && {d} && {e} && {f}'''.format(
			a = "mkdir -p " + bowtie2_dir,
			b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
			c = "bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -U " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
			d = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
			e = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
			f = "jgi_summarize_bam_contig_depths --outputDepth " + depths_dir + name.split("/")[-1] + ".contig_depths.txt " + bam_sorted,
			)
	elif paired == "unpaired":
		command = '''{a} && {b} && {c} && {d} && {e} && {f}'''.format(
			a = "mkdir -p " + bowtie2_dir,
			b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
			c = "bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
			d = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
			e = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
			f = "jgi_summarize_bam_contig_depths --outputDepth " + depths_dir + name.split("/")[-1] + ".contig_depths.txt " + bam_sorted,
			)
	elif paired == "concatenated":
		command = '''{a} && {b} && {c} && {d} && {e} && {f}'''.format(
			a = "mkdir -p " + bowtie2_dir,
			b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2-build " + contigs + " " + index,
			c = "if grep -q -m 1 " + name + " " + list_paired + "; then bowtie2 -x " + index + " -1 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq.gz" + " -2 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq.gz" + " -U " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_1.fastq.gz" +
				"," + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_2.fastq" + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal; else bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " +
				sam + " -p " + str(cores) + " --very-sensitive-local --no-unal; fi",
			d = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
			e = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
			f = "jgi_summarize_bam_contig_depths --outputDepth " + depths_dir + name.split("/")[-1] + ".contig_depths.txt " + bam_sorted,
			)

	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="align", paired=paired)[0]):
		workflow.add_task_gridable(actions=align(name, paired),
			depends=list_depends(name=name, step="align", paired=paired),
			targets=list_targets(name=name, step="align", paired=paired),
			time=calculate_time(name=name, step="align", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

#############################
# function to run MetaBAT 2 #
#############################

def metabat(name):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	depth = depths_dir + name.split("/")[-1] + ".contig_depths.txt"
	metabat_tmp = mags_scratch + name.split("/")[-1] + "/bins/" + name.split("/")[-1]
	metabat_out = bins_dir + name.split("/")[-1] + "/bins/"
	command = '''{a} && {b} && {c} && {d}'''.format(
		a = "if [ ! -s " + contigs + " ]; then mkdir -p " + metabat_tmp + " && touch " + metabat_tmp + ".bin.lowDepth.fa && touch " + metabat_tmp + ".bin.tooShort.fa && touch " + metabat_tmp + ".bin.unbinned.fa; else metabat2 -i " + contigs + " -a " + depth + " -o " + metabat_tmp + ".bin --unbinned -m 1500 -t " + str(cores) + " " + args.metabat_options + "; fi",
		b = "mkdir -p " + metabat_out,
		c = "cp " + metabat_tmp + "*.fa " + metabat_out,
		d = "touch [targets[0]]"
		)
	return str(command)

for name in names:
	if not os.path.isfile(list_targets(name=name, step="metabat", paired=paired)[0]):
		workflow.add_task_gridable(actions=metabat(name),
			depends=list_depends(name=name, step="metabat", paired=paired),
			targets=list_targets(name=name, step="metabat", paired=paired),
			time=calculate_time(name=name, step="metabat", paired=paired),
			mem=memory,
			cores=cores,
			partition=partition
			)

###################################
# function to calculate abundance #
###################################

def abundance_sample(name, paired):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = mags_scratch + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_dir + name.split("/")[-1]
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"
	bam_index = bowtie2_dir + name.split("/")[-1] + ".sorted.bam.bai"
	bin = bins_dir + name.split("/")[-1] + "/bins"

	if paired == "paired":
		if input_extension in ["fastq.gz", "fq.gz"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(zcat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(zcat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				g = "unpaired1=$(echo $(zcat " + name + "_unmatched_1." + input_extension + "|wc -l)/4|bc)",
				h = "unpaired2=$(echo $(zcat " + name + "_unmatched_2." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)
		elif input_extension in ["fastq", "fq"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "paired1=$(echo $(cat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				f = "paired2=$(echo $(cat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				g = "unpaired1=$(echo $(cat " + name + "_unmatched_1." + input_extension + "|wc -l)/4|bc)",
				h = "unpaired2=$(echo $(cat " + name + "_unmatched_2." + input_extension + "|wc -l)/4|bc)",
				i = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				j = "rm -r " + bowtie2_dir
				)
	elif paired in ["unpaired", "concatenated"]:
		if input_extension in ["fastq.gz", "fq.gz"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} '''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "echo $(zcat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				f = "rm -r " + bowtie2_dir
				)
		elif input_extension in ["fastq", "fq"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} '''.format(
				a = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				b = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				c = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				d = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				e = "echo $(cat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				f = "rm -r " + bowtie2_dir
				)
	return str(command)

def rebuild_bowtie2_db():
	command = '''{a} && {b} && {c} && {d}'''.format(
		a = "rm -r " + bowtie2_global_dir,
		b = "mkdir -p " + bowtie2_global_dir,
		c = "python " + this_folder + "by_dataset_db.py --mag_dir " + bins_dir + " --out_dir " + bowtie2_global_dir,
		d = "touch " + abundance_dir + "built.done"
		)
	return str(command)

def abundance_dataset(name, paired):
	contigs = contigs_dir + name.split("/")[-1] + "/" + name.split("/")[-1] + ".final.contigs.fa"
	bowtie2_dir = mags_scratch + name.split("/")[-1] + "/bowtie2/"
	index = bowtie2_global_dir + "mag_db"
	sam = bowtie2_dir + name.split("/")[-1] + ".sam"
	bam_unsorted = bowtie2_dir + name.split("/")[-1] + ".unsorted.bam"
	bam_sorted = bowtie2_dir + name.split("/")[-1] + ".sorted.bam"
	bam_index = bowtie2_dir + name.split("/")[-1] + ".sorted.bam.bai"
	bin = bowtie2_global_dir + "tmp/"

	if paired == "paired":
		if input_extension in ["fastq.gz", "fq.gz"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {k} && {l} && {m}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -U " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "paired1=$(echo $(zcat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				j = "paired2=$(echo $(zcat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				k = "unpaired1=$(echo $(zcat " + name + "_unmatched_1." + input_extension + "|wc -l)/4|bc)",
				l = "unpaired2=$(echo $(zcat " + name + "_unmatched_2." + input_extension + "|wc -l)/4|bc)",
				m = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				)
		elif input_extension in ["fastq", "fq"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j} && {k} && {l} && {m}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -1 " + name + "_paired_1." + input_extension + " -2 " + name + "_paired_2." + input_extension + " -U " + name + "_unmatched_1." + input_extension + "," + name + "_unmatched_2." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "paired1=$(echo $(cat " + name + "_paired_1." + input_extension + "|wc -l)/4|bc)",
				j = "paired2=$(echo $(cat " + name + "_paired_2." + input_extension + "|wc -l)/4|bc)",
				k = "unpaired1=$(echo $(cat " + name + "_unmatched_1." + input_extension + "|wc -l)/4|bc)",
				l = "unpaired2=$(echo $(cat " + name + "_unmatched_2." + input_extension + "|wc -l)/4|bc)",
				m = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				)
	elif paired == "unpaired":
		if input_extension in ["fastq.gz", "fq.gz"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi " + args.checkm_coverage_options,
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "echo $(zcat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				j = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				)
		elif input_extension in ["fastq", "fq"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r " + args.checkm_coverage_options,
				g = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "echo $(cat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				j = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				)

	elif paired == "concatenated":
		if input_extension in ["fastq.gz", "fq.gz"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else if grep -q -m 1 " + name + " " + list_paired + "; then bowtie2 -x " + index + " -1 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq.gz" + " -2 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq.gz" + " -U " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_1.fastq.gz" +
					"," + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_2.fastq" + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal; else bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " +
					sam + " -p " + str(cores) + " --very-sensitive-local --no-unal; fi",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r",
				g = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "echo $(zcat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				j = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				)
		elif input_extension in ["fastq", "fq"]:
			command = '''{a} && {b} && {c} && {d} && {e} && {f} && {g} && {h} && {i} && {j}'''.format(
				a = "mkdir -p " + bowtie2_dir,
				b = "if [ ! -s " + contigs + " ]; then touch " + bam_sorted + "; else if grep -q -m 1 " + name + " " + list_paired + "; then bowtie2 -x " + index + " -1 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_1.fastq.gz" + " -2 " + scratch_deconcatenated + name.split("/")[-1] + "_paired_2.fastq.gz" + " -U " + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_1.fastq.gz" +
					"," + scratch_deconcatenated + name.split("/")[-1] + "_unmatched_2.fastq" + " -S " + sam + " -p " + str(cores) + " --very-sensitive-local --no-unal; else bowtie2 -x " + index + " -U " + name + "." + input_extension + " -S " +
					sam + " -p " + str(cores) + " --very-sensitive-local --no-unal; fi",
				c = "samtools view -bS -F 4 " + sam + " > " + bam_unsorted,
				d = "samtools sort " + bam_unsorted + " -o " + bam_sorted + " --threads " + str(cores) + "; fi",
				e = "if [ ! -s " + contigs + " ]; then echo -e \"Sequence Id\tBin Id\tSequence length (bp)\tBam Id\tCoverage\tMapped reads\" > [targets[0]] && echo -e \"Bin Id\tBin size (Mbp)\t" + name.split("/")[-1] + ".sorted: mapped reads\t" + name.split("/")[-1] + ".sorted: % mapped reads\t" + name.split("/")[-1] + ".sorted: % binned populations\t" + name.split("/")[-1] + ".sorted: % community\" > [targets[1]] && echo 0 > [targets[2]]; else samtools index " + bam_sorted + " -@ " + str(cores) + " " + bam_index,
				f = "checkm coverage " + bin + " [targets[0]] " + bam_sorted + " -x fa -t " + str(cores) + " -r",
				g = "checkm profile [targets[0]] --tab_table -f [targets[1]]",
				h = "samtools view -c -F 260 " + bam_sorted + " -o [targets[2]]; fi",
				i = "echo $(cat " + name + "." + input_extension + "|wc -l)/4|bc &>> [targets[2]]",
				j = "echo $((paired1+paired2+unpaired1+unpaired2)) &>> [targets[2]]",
				)
	return str(command)

if abundance_type == "by_sample":
	for name in names:
		if not os.path.isfile(list_targets(name=name, step="abundance", paired=paired)[0]):
			workflow.add_task_gridable(actions=abundance_sample(name, paired),
				depends=list_depends(name=name, step="abundance_sample", paired=paired),
				targets=list_targets(name=name, step="abundance", paired=paired),
				time=calculate_time(name=name, step="abundance", paired=paired),
				mem=memory,
				cores=cores,
				partition=partition
				)
else:
	workflow.add_task_gridable(actions=rebuild_bowtie2_db(),
		depends=[str(bins_dir + name.split("/")[-1] + ".done") for name in names],
		targets=abundance_dir + "built.done",
		time=30 + 10 * len(names),
		mem=memory,
		cores=1,
		partition=partition
		)
	for name in names:
		if not os.path.isfile(list_targets(name=name, step="abundance", paired=paired)[0]):
			workflow.add_task_gridable(actions=abundance_dataset(name, paired),
				depends=abundance_dir + "built.done",
				targets=list_targets(name=name, step="abundance", paired=paired),
				time=calculate_time(name=name, step="abundance", paired=paired),
				mem=memory,
				cores=cores,
				partition=partition
				)

################
# Calculate GC #
################

if args.gc_length_stats and not os.path.isdir(qa_dir + "GC_content.tsv"):
	command = "python " + this_folder + "calculateGC.py --in-dir " + bins_dir + " --out-file " + qa_dir + "GC_content.tsv" + " --threads " + str(local_jobs)
	workflow.add_task(command, targets=qa_dir + "GC_content.tsv", depends=[item for name in names for item in list_targets(name=name, step="abundance", paired=paired)])

######################################
# run checkm and phylophlan workflow #
######################################

if not os.path.isfile(qa_dir + "checkm_qa_and_n50.tsv") or not os.path.isfile(phylophlan_dir + "phylophlan_out.tsv"):
	add_string = ""
	if args.checkm_predict_options != "":
		add_string = add_string + " --checkm-predict-options \"" + args.checkm_predict_options + "\""
	if args.phylophlan_metagenomic_options != "":
		add_string = add_string + " --phylophlan-metagenomic-options \"" + args.phylophlan_metagenomic_options + "\""
	if args.checkm_data_path:
		add_string = add_string + " --checkm-data-path " + args.checkm_data_path
	if args.grid_options:
		add_string = add_string + " --grid-options=\"" + " ".join(args.grid_options) + "\""
	command = "python " + this_folder + "checkm_phylophlan_workflow.py -o " + output + "checkm_phylophlan_steps/" + " --n 300 --grid-scratch " + scratch + "checkm_phylophlan_steps/" + " --grid-jobs " + str(args.grid_jobs) + " --grid-partition " + partition + " --cores " + str(cores) + " --mem " + str(memory) + " --completeness " + str(args.completeness) + " --contamination " + str(args.contamination) + " --phylophlan-database " + database + " --phylophlan-database-folder " + database_folder + " --checkm-path " + args.checkm_path + add_string
	workflow.add_task(actions=command,
		depends=[str(bins_dir + name.split("/")[-1] + ".done") for name in names],
		targets=[qa_dir + "checkm_qa_and_n50.tsv", phylophlan_dir + "phylophlan_out.tsv"])

if not os.path.isfile(phylophlan_dir + "phylophlan_relab.tsv"):
	workflow.add_task("python " + this_folder + "phylophlan_add_tax_assignment.py --table " + phylophlan_dir + "phylophlan_out.tsv" + " --output " + phylophlan_dir + "phylophlan_relab.tsv",
		depends=phylophlan_dir + "phylophlan_out.tsv",
		targets=phylophlan_dir + "phylophlan_relab.tsv")

########
# SGBs #
########

############################
# list MAGs to run Mash on #
############################

if not os.path.isfile(mash_dir + "mags_filepaths.txt"):
	list_inputs = "python " + this_folder + "mash_list_inputs.py --checkm " + qa_dir + "checkm_qa_and_n50.tsv --phylophlan " + phylophlan_dir + "phylophlan_relab.tsv" + " --bins " + bins_dir + " --mash " + mash_dir + " --threads " + str(local_jobs)
	workflow.add_task(list_inputs, depends=[phylophlan_dir + "phylophlan_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv"], targets=mash_dir + "mags_filepaths.txt")

###############
# mash sketch #
###############

if not os.path.isfile(mash_dir + "sketches.msh"):
	sketch = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "sketches.msh" + "; else mash sketch -p " + str(local_jobs) + " -l " + mash_dir + "mags_filepaths.txt" + " -o " + mash_dir + "sketches " + args.mash_sketch_options + "; fi"
	workflow.add_task(sketch, depends=mash_dir + "mags_filepaths.txt", targets=mash_dir + "sketches.msh")

##############
# mash paste #
##############

if not os.path.isfile(mash_dir + "references.msh"):
	paste = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "references.msh" + "; else mash paste " + mash_dir + "references " + mash_dir + "sketches.msh; fi"
	workflow.add_task(paste, depends=mash_dir + "sketches.msh", targets=mash_dir + "references.msh")

#############
# mash dist #
#############

if not os.path.isfile(mash_dir + "mash_dist_out.tsv"):
	dist = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then touch " + mash_dir + "mash_dist_out.tsv" + "; else mash dist -p " + str(local_jobs) + " -t " + mash_dir + "references.msh " + mash_dir + "sketches.msh > " + mash_dir + "mash_dist_out.tsv; fi"
	workflow.add_task(dist, depends=mash_dir + "references.msh", targets=mash_dir + "mash_dist_out.tsv")

#################
# identify SGBS #
#################

depends_list = [item for name in names for item in list_targets(name=name, step="abundance", paired=paired)]
depends_list.extend([phylophlan_dir + "phylophlan_relab.tsv", qa_dir + "checkm_qa_and_n50.tsv", sgb_dir + "sgbs/SGB_info.tsv"])

# groups MAGs into SGBs using (1) Mash and (2) fastANI
if not os.path.isfile(sgb_dir + "fastANI/SGB_list.txt") or not os.path.isfile(sgb_dir + "sgbs/SGB_info.tsv"):
	cluster = "if [ ! -s " + mash_dir + "mags_filepaths.txt" + " ]; then mkdir -p " + sgb_dir + "fastANI/" + " && touch " + sgb_dir + "fastANI/SGB_list.txt && mkdir -p " + sgb_dir + "sgbs/" + " && echo -e \"cluster\tgenome\tcluster_members\tn_genomes\tcompleteness\tcontamination\tstrain_heterogeneity\tn50\tquality\tkeep\tcluster_name\tsgb\" > " + sgb_dir + "sgbs/SGB_info.tsv; else Rscript " + this_folder + "mash_clusters.R --mash " + mash_dir + "mash_dist_out.tsv --checkm " + qa_dir + "checkm_qa_and_n50.tsv" + " --phylo " + phylophlan_dir + "phylophlan_relab.tsv --out_dir " + sgb_dir + " --threads " + str(local_jobs) + " --mag_dir " + bins_dir + "qc_bins/; fi"
	workflow.add_task(cluster, depends=[mash_dir + "mash_dist_out.tsv", phylophlan_dir + "phylophlan_relab.tsv"], targets=[sgb_dir + "sgbs/SGB_info.tsv"])

if not os.path.isfile(output + "final_profile.tsv"):
	merge = "Rscript " + this_folder + "merge_tax_and_abundance.R" + " -i " + abundance_dir + " --tax " + phylophlan_dir + "phylophlan_relab.tsv" + " --qa " + qa_dir + "checkm_qa_and_n50.tsv --sgbs " + sgb_dir + "sgbs/SGB_info.tsv" + " -o " + output + "final_profile_" + abundance_type + ".tsv"
	workflow.add_task(merge, depends=depends_list, targets = output + "final_profile_" + abundance_type + ".tsv")

####################
# run the workflow #
####################

workflow.go()

#
