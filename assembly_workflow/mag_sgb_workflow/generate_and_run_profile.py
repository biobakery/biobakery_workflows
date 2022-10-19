from anadama2 import Workflow
import os
import pandas as pd
from datetime import date

workflow = Workflow(version = "0.4")

workflow.add_argument("start-date", desc="YYYY/MM/DD")
workflow.add_argument("end-date", desc="YYYY/MM/DD")
workflow.add_argument("tmp", desc="TMP_DIRECTORY that will be deleted")
workflow.add_argument("strain-n", desc="number of strains per species")
workflow.add_argument("n", desc="number of species", default=100)
workflow.add_argument("samples", desc="number of samples")
workflow.add_argument("camisim-path", desc="folder with metagenomesimulation.py")
workflow.add_argument("sample-size", desc="gigabases")
workflow.add_argument("email", desc="email for Entrez")
workflow.add_argument("api-key", desc="api key for Entrez")
workflow.add_argument("gc", desc="GC content from 0 to 1", default="unrestricted")
workflow.add_argument("genome-size", desc="genome-size k", default="unrestricted")
workflow.add_argument("sigma", desc="sigma", default=0.5)
args = workflow.parse_args()

this_folder = os.path.realpath(__file__).rsplit("/", 1)[0] + "/"

current_date = date.today().strftime("%Y%m%d")

output = "/" + args.output.strip("/") + "/"
if not os.path.isdir(output):
    os.makedirs(output)

profiles = output + "profiles/"
if not os.path.isdir(profiles):
    os.makedirs(profiles)

config_dir = output + "config/"
if not os.path.isdir(config_dir):
    os.makedirs(config_dir)

inputs_dir = output + "inputs/"
if not os.path.isdir(output + "inputs/"):
    os.makedirs(output + "inputs/")

tmp = "/" + args.tmp.strip("/") + "/"
if not os.path.isdir(tmp):
    os.makedirs(tmp)

camisim_path = "/" + args.camisim_path.strip("/") + "/"

gc = args.gc
genome_size_k = args.genome_size

if (gc != "unrestricted") == (genome_size_k != "unrestricted"):
    raise ValueError("GC content and genome size cannot both be specified or unspecified")

profile = pd.read_csv(args.input, sep='\t')
if not profile.columns.tolist() == ["genome_ID", "NCBI_ID", "filepath", "abundance", "novelty_category"]:
    raise ValueError("Wrong column names")

if not abs(sum(profile['abundance']) - 1) < 0.001:
    raise ValueError("Abundances don't sum to 1")

if any([item == "" or item is None or pd.isnull(item) for item in profile.loc[["BIN_" in item for item in profile['genome_ID']], 'filepath'].tolist()]):
    raise ValueError("Assembled bins need file paths")

if any([item == "" or item is None or pd.isnull(item) for item in profile.loc[["UNKNOWN" not in item for item in profile['genome_ID']], 'novelty_category'].tolist()]):
    raise ValueError("Non-UNKNOWN entries need a novelty category")

if any([item == "" or item is None or pd.isnull(item) for item in profile['abundance'].tolist()]):
    raise ValueError("All entries need an abundance")

if any([item == "" or item is None or pd.isnull(item) for item in profile.loc[["UNKNOWN" not in item for item in profile['genome_ID']], 'NCBI_ID'].tolist()]):
    raise ValueError("Non-UNKNOWN entries need an NCBI ID")

workflow.add_task("python " + this_folder + "download_new_genomes.py --out-folder " + inputs_dir + " --tax-dump-folder " + inputs_dir + "NCBI_tax/ --start-date " + args.start_date + " --end-date " + args.end_date + " --current-date " + current_date + " --n all --email " + args.email + " --api-key " + args.api_key,
    targets=[inputs_dir + "NCBI_tax/ncbi-taxonomy_" + current_date + ".tar.gz", inputs_dir + "new_genomes_metadata.tsv"])

profile[["UNKNOWN" not in item and "BIN_" not in item for item in profile['genome_ID']]].to_csv(inputs_dir + "reference_profile.tsv", index = False, sep = "\t")

workflow.add_task("python " + this_folder + "download_known_genomes.py --reference-profile " + inputs_dir + "reference_profile.tsv --out-folder " + inputs_dir + " --email " + args.email + " --api-key " + args.api_key,
    targets=inputs_dir + "reference_metadata.tsv", depends=inputs_dir + "reference_profile.tsv")

files_to_move = profile[["BIN_" in item for item in profile['genome_ID']]]['filepath'].tolist()

if not os.path.isdir(inputs_dir + "genomes/"):
    os.makedirs(inputs_dir + "genomes/")

for file in files_to_move:
    workflow.add_task("cp " + file + " " + inputs_dir + "genomes/" + file.split("/")[-1],
        targets=inputs_dir + "genomes/" + file.split("/")[-1], depends=file)

if not os.path.isfile(inputs_dir + "GC_and_lengths.tsv"):
    command = "python " + this_folder + "gunzip_files.py -i " + inputs_dir + "genomes/ --local-jobs " + str(args.jobs) + " -o " + inputs_dir + "genomes/ && python " + this_folder + "calculateGC.py --in-dir " + inputs_dir + "genomes/" + " --out-file " + inputs_dir + "GC_and_lengths.tsv" + " --threads " + str(args.jobs)
    depends = [inputs_dir + "genomes/" + file.split("/")[-1] for file in files_to_move]
    depends.extend([inputs_dir + "reference_metadata.tsv", inputs_dir + "new_genomes_metadata.tsv"])
    workflow.add_task(command, targets=inputs_dir + "GC_and_lengths.tsv", depends=depends)

# Create tsvs for CAMISIM
depends = [inputs_dir + "genomes/" + file.split("/")[-1] for file in files_to_move]
depends.extend([inputs_dir + "reference_profile.tsv", inputs_dir + "new_genomes_metadata.tsv", inputs_dir + "GC_and_lengths.tsv"])

workflow.add_task("python " + this_folder + "make_config_files.py --profile " + args.input + " --inputs-folder " + inputs_dir + " --configs-folder " + config_dir + " --strain-n " + str(args.strain_n) + " --local-jobs " + str(args.jobs) +
    " --samples " + str(args.samples) + " --tmp-dir " + tmp + " --profiles-dir " + profiles + " --camisim-path " + camisim_path + " --sample-size " + str(args.sample_size) + " --current-date " + current_date + " --gc-len-file " + inputs_dir + "GC_and_lengths.tsv",
    targets=[config_dir + "metadata.tsv", config_dir + "genome_to_id.tsv", config_dir + "config.ini"], depends=depends)

if genome_size_k == "unrestricted":
    workflow.add_task("Rscript get_GC_abundances.R --GC_target " + str(gc) + " --num_species " + str(args.n) + " --genome_to_id " + config_dir + "genome_to_id.tsv" + " --profile " + args.input + " --GC_and_lengths " + inputs_dir + "GC_and_lengths.tsv" + " -o " + config_dir + "abundances.tsv",
        targets=[config_dir + "abundances.tsv"], depends=[config_dir + "metadata.tsv", config_dir + "genome_to_id.tsv", config_dir + "config.ini"])

if gc == "unrestricted":
    workflow.add_task("Rscript get_genome_size_abundances.R --k " + str(genome_size_k) + " --sigma " + str(args.sigma) + " --num_species " + str(args.n) + " --genome_to_id " + config_dir + "genome_to_id.tsv" + " --profile " + args.input + " --GC_and_lengths " + inputs_dir + "GC_and_lengths.tsv" + " -o " + config_dir + "abundances.tsv",
        targets=[config_dir + "abundances.tsv"], depends=[config_dir + "metadata.tsv", config_dir + "genome_to_id.tsv", config_dir + "config.ini"])

# Run CAMISIM
"""
workflow.add_task("python " + camisim_path + "metagenomesimulation.py " + config_dir + "config.ini",
    targets=[profiles + "taxonomic_profile_" + str(i) + ".txt" for i in range(int(args.samples))],
    depends=[config_dir + "metadata.tsv", config_dir + "genome_to_id.tsv", config_dir + "abundances.tsv", config_dir + "config.ini"])
"""

workflow.go()
